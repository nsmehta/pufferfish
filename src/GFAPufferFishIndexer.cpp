#include "FastxParser.hpp"
#include <cmath>
#include <iostream>
#include <iterator>
#include <type_traits>
#include <vector>
#include <bitset>
#include <sstream>

#include "cereal/archives/json.hpp"
#include "CanonicalKmer.hpp"
#include "OurGFAReader.hpp"
#include "PufferFS.hpp"
#include "ScopedTimer.hpp"
#include "Util.hpp"
#include "jellyfish/mer_dna.hpp"
#include "sdsl/int_vector.hpp"
#include "sdsl/rank_support.hpp"
#include "sdsl/select_support.hpp"
#include "spdlog/spdlog.h"
#include "PufferfishIndex.hpp"
#include "sparsepp/spp.h"
//#include "gfakluge.hpp"

uint64_t swap_uint64(uint64_t val) {
  val = ((val << 8) & 0xFF00FF00FF00FF00ULL) |
        ((val >> 8) & 0x00FF00FF00FF00FFULL);
  val = ((val << 16) & 0xFFFF0000FFFF0000ULL) |
        ((val >> 16) & 0x0000FFFF0000FFFFULL);
  return (val << 32) | (val >> 32);
}

// adapted from :
// http://stackoverflow.com/questions/34875315/implementation-my-own-list-and-iterator-stl-c
class ContigKmerIterator {
public:
  typedef ContigKmerIterator self_type;
  typedef uint64_t value_type;
  typedef value_type& reference;
  typedef value_type* pointer;
  typedef std::forward_iterator_tag iterator_category;
  typedef int64_t difference_type;

  ContigKmerIterator(sdsl::int_vector<>* storage, sdsl::bit_vector* rank,
                     uint8_t k, uint64_t startAt)
      : storage_(storage), rank_(rank), k_(k), curr_(startAt) {
    if (curr_ + k_ <= rank_->size()) {
      mer_.fromNum(storage_->get_int(2 * curr_, 2 * k_));
      // mer_.word__(0) = storage_->get_int(2 * curr_, 2 * k_);
    }
    // rcMer_ = mer_.get_reverse_complement();
  }

  ContigKmerIterator&
  operator=(ContigKmerIterator& other) { //}= default;//(sdsl::int_vector<>*
                                         //storage, sdsl::bit_vector* rank,
                                         //uint8_t k, uint64_t startAt) :
    storage_ = other.storage_;
    rank_ = other.rank_;
    k_ = other.k_;
    curr_ = other.curr_;
    mer_ = other.mer_;
    // rcMer_ = other.rcMer_;
    word_ = other.word_;
    return *this;
  }

  ContigKmerIterator operator++() {
    ContigKmerIterator i = *this;
    advance_();
    return i;
  }

  ContigKmerIterator operator++(int) {
    advance_();
    return *this;
  }

  reference operator*() {
    // word_ = (mer_.word(0) < rcMer_.word(0)) ? mer_.word(0) : rcMer_.word(0);
    word_ = mer_.getCanonicalWord();
    return word_;
  }

  difference_type pos() { return curr_; }

  bool isEndKmer() {
	  size_t endPos = curr_ + k_ - 1;
	  if((*rank_)[endPos] == 1)
		  return true;
	  else
		  return false;
  }


  bool isCanonical(){
	  return mer_.fwWord() == mer_.getCanonicalWord() ;
  }

  pointer operator->() {
    word_ = mer_.getCanonicalWord(); //(mer_.word(0) < rcMer_.word(0)) ?
                                     //mer_.word(0) : rcMer_.word(0);
    return &word_;
  }
  bool operator==(const self_type& rhs) { return curr_ == rhs.curr_; }

  bool operator!=(const self_type& rhs) { return curr_ != rhs.curr_; }

  bool operator<(const self_type& rhs) { return curr_ < rhs.curr_; }

  bool operator<=(const self_type& rhs) { return curr_ <= rhs.curr_; }

private:
  void advance_() {
    size_t endPos = curr_ + k_ - 1;
    if (endPos + 1 < rank_->size() and (*rank_)[endPos] == 1) {
      curr_ += k_;
      mer_.fromNum(storage_->get_int(2 * curr_, 2 * k_));
    } else {
      if (curr_ + k_ < rank_->size()) {
        int c = (*storage_)[curr_ + k_];
        mer_.shiftFw(c);
      } else {
        mer_.fromNum(storage_->get_int(2 * (rank_->size() - k_), 2 * k_));
      }
      ++curr_;
    }
  }
  sdsl::int_vector<>* storage_{nullptr};
  sdsl::bit_vector* rank_{nullptr};
  uint8_t k_{0};
  uint64_t curr_{0};
  CanonicalKmer mer_;
  uint64_t word_{0};
};


int pufferfishTest(util::TestOptions& testOpts) {
  (void)testOpts;
  std::cerr << "this command is not yet implemented\n";
  return 1;
}

void sampledPositions(size_t tlen, int sampleSize, std::vector<size_t>& sampledInds){
  sampledInds.clear() ;
  auto numOfKmers = tlen - 31 ;
  size_t lastCovered = 0 ;


  for(size_t j = 0 ; j <= numOfKmers; j++){
    if(j > lastCovered){
      auto next_samp = std::min(j + sampleSize/2 - 1 ,numOfKmers) ;
      sampledInds.push_back(next_samp) ;
      lastCovered = next_samp + sampleSize/2 + 1;
    }
  }
  if(lastCovered == numOfKmers)
    sampledInds.push_back(lastCovered) ;
}

std::string packedToString(sdsl::int_vector<>& seqVec, uint64_t offset, uint32_t len) {
  std::stringstream s;
  for (size_t i = offset; i < offset + len; ++i) {
    auto c = seqVec[i];
    s << my_mer::rev_code(c);
  }
  auto st = s.str();
  return st;
}


enum class NextSampleDirection : uint8_t { FORWARD = 0, REVERSE=1 };

uint32_t getEncodedExtension(sdsl::int_vector<>& seqVec, uint64_t firstSampPos, uint64_t distToSamplePos,
                             uint32_t maxExt, NextSampleDirection dir, bool print=false) {
  uint32_t encodedNucs{0};
  uint32_t bitsPerCode{3};
  std::vector<uint32_t> charsToExtend;
  size_t i = 0;
  for (; i < distToSamplePos; ++i) {
    if ( firstSampPos + i >= seqVec.size() ) {
      std::cerr << "seqVec.size() " << seqVec.size() << ", looking for index " << firstSampPos + i << "\n";
      std::cerr << "dist to sample is " << distToSamplePos << ", i = " << i << "\n";
    }
    auto c = seqVec[firstSampPos + i];
    charsToExtend.push_back(c);
  }
  if (dir == NextSampleDirection::REVERSE) {
    std::reverse(charsToExtend.begin(), charsToExtend.end());
  }
  if (i < maxExt) {
    charsToExtend.push_back(0x4);
  }

  if (print) {
  std::cerr << "ext = [";
  }
  for (size_t j = 0; j < charsToExtend.size(); ++j) {
    auto c = charsToExtend[j];
    if (print) { std::cerr << c << ", "; }
    encodedNucs |= (c << (bitsPerCode  * (maxExt - j - 1)));
    std::bitset<32> pIt(encodedNucs) ;
    if(print){std::cerr << " j:  " << j << " bits: " << pIt << "\n" ; }
  }
  if (print) {std::cerr << "]\n";}
  return encodedNucs;
}


int pufferfishIndex(util::IndexOptions& indexOpts) {
  uint32_t k = indexOpts.k;
  std::string gfa_file = indexOpts.gfa_file;
  std::string rfile = indexOpts.rfile;
  std::string outdir = indexOpts.outdir;

  // If the user included the '/' in the output directory path, remove
  // it here
  if (outdir.back() == '/') { outdir.pop_back(); }
  std::vector<std::string> read_file = {rfile};

  auto console = spdlog::stderr_color_mt("console");

  size_t tlen{0};
  size_t numKmers{0};
  size_t nread{0};
  CanonicalKmer::k(k);

  puffer::fs::MakeDir(outdir.c_str());

  PosFinder pf(gfa_file.c_str(), k - 1);
  pf.parseFile();
  pf.mapContig2Pos();
  pf.serializeContigTable(outdir);

  {
    auto& cnmap = pf.getContigNameMap();
    for (auto& kv : cnmap) {
      auto& r1 = kv.second;
      tlen += r1.length();
      numKmers += r1.length() - k + 1;
      ++nread;
    }
    console->info("# segments = {}", nread);
    console->info("total length = {}", tlen);
  }


  // now we know the size we need --- create our bitvectors and pack!
  size_t gpos{0};
  size_t w = std::log2(tlen) + 1;
  console->info("positional integer width = {}", w);
  sdsl::int_vector<> seqVec(tlen, 0, 2);

  sdsl::bit_vector rankVec(tlen);
  tlen = 0;

  size_t numContigs{0};
  size_t nkeys{0};
  {
    auto& cnmap = pf.getContigNameMap();
    for (auto& kv : cnmap) {
      ++numContigs;
      size_t len{0};
      const auto& r1 = kv.second;
      CanonicalKmer mer;
      for (size_t i = 0; i < r1.length(); ++i) {
        auto offset = i; // r1.length() - i - 1;
        // NOTE: Having to add things in the reverse order here is strange
        // we should make sure that this doesn't mess with the positions we
        // end up storing!
        auto c = my_mer::code(r1[i]);
        seqVec[gpos + offset] = c;
        mer.shiftFw(c);
        ++len;
        if (len >= k) {
#ifdef PUFFER_DEBUG
          my_mer mm;
          uint64_t num = seqVec.get_int(2 * (gpos + len - k), 2 * k);
          mm.word__(0) = num; // mm.canonicalize();
          if (!(mm == mer.fwMer() or
                mm == mer.rcMer())) { //}mm != mer.get_canonical()) {
            std::cerr << "num & 0x3 = " << (num & 0x3) << "\n";
            std::cerr << "i = " << i << "\n";
            std::cerr << mer.to_str() << ", " << mm.to_str() << "\n";
          }
#endif // PUFFER_DEBUG
          ++nkeys;
        }
        //++gpos;
      }
      gpos += r1.length();
      tlen += r1.length();
      rankVec[tlen - 1] = 1;
    }
  }
  std::cerr << "seqSize = " << sdsl::size_in_mega_bytes(seqVec) << "\n";
  std::cerr << "rankSize = " << sdsl::size_in_mega_bytes(rankVec) << "\n";
  //std::cerr << "posSize = " << sdsl::size_in_mega_bytes(posVec) << "\n";
  std::cerr << "num keys = " << nkeys << "\n";
  typedef boomphf::SingleHashFunctor<uint64_t> hasher_t;
  typedef boomphf::mphf<uint64_t, hasher_t> boophf_t;

  ContigKmerIterator kb(&seqVec, &rankVec, k, 0);
  ContigKmerIterator ke(&seqVec, &rankVec, k, seqVec.size() - k + 1);

#ifdef PUFFER_DEBUG
  auto ks = kb;
  size_t nkeyIt{0};
  for (; ks < ke; ++ks) {
    nkeyIt++;
  }
  std::cerr << "num keys (iterator)= " << nkeyIt << "\n";
#endif //PUFFER_DEBUG

  auto keyIt = boomphf::range(kb, ke);
  boophf_t* bphf = new boophf_t(nkeys, keyIt, 16, 2.5); // keys.size(), keys, 16);
  std::cerr << "mphf size = " << (bphf->totalBitSize() / 8) / std::pow(2, 20)
            << "\n";

  sdsl::store_to_file(seqVec, outdir + "/seq.bin");
  sdsl::store_to_file(rankVec, outdir + "/rank.bin");

  //sdsl::int_vector<> posVec(tlen, 0, w);

  

   int sampleSize = 5;
   int extensionSize = 2 ;
   size_t sampledKmers{0};


   std::vector<uint32_t> startSamplePositions ;
   std::vector<std::vector<size_t> > allSampledPositions ;
   std::vector<size_t> contigLengths ;
   //track start positions
   /*{
 	  auto& cnmap = pf.getContigNameMap() ;
 	  for(auto& kv : cnmap){
 		  auto& r1 = kv.second ;
 		  sampledKmers += ((int)(r1.length()-k)/sampleSize + 1) ;
 		  startSamplePositions.push_back((r1.length()-k)%sampleSize) ;
 	  }
 	  console->info("# sampled kmers ={}", sampledKmers) ;
 	  console->info("# skipped kmers ={}", numKmers - sampledKmers) ;
    }*/

  //fill up optimal positions
   {
     auto& cnmap = pf.getContigNameMap() ;
     for(auto& kv : cnmap){
       auto& r1 = kv.second ;
       std::vector<size_t> sampledInds ;
       sampledPositions(r1.length(), sampleSize, sampledInds) ;
       sampledKmers += sampledInds.size() ;
       allSampledPositions.push_back(sampledInds) ;
       contigLengths.push_back(r1.length()) ;
       //startSamplePositions.push_back((r1.length()-k)%sampleSize) ;
     }
     console->info("# sampled kmers = {}", sampledKmers) ;
     console->info("# skipped kmers = {}", numKmers - sampledKmers) ;
   }



  //fill up sampledPosVec
  //store data for sampled information
  //sdsl::bit_vector samplePosVec(sampledKmers*w+3*(numKmers-sampledKmers)) ;
  sdsl::bit_vector presenceVec(nkeys);
  sdsl::bit_vector canonicalNess(nkeys);
  //check if we find the kmers
  //in the way we thought we stored them
  //need these two to fetch

  //feel up the canonicalNess vector
  sdsl::int_vector<> auxInfo((numKmers-sampledKmers),0,3*extensionSize) ;
  sdsl::bit_vector direction(numKmers - sampledKmers) ;
  sdsl::int_vector<> samplePosVec(sampledKmers, 0, w);
  
  {
    ContigKmerIterator kb1(&seqVec, &rankVec, k, 0);
    ContigKmerIterator ke1(&seqVec, &rankVec, k, seqVec.size() - k + 1);
    for(;kb1!=ke1;++kb1){
    	auto idx = bphf->lookup(*kb1) ;
    	canonicalNess[idx] = kb1.isCanonical() ;
    }
  }

  //old presenceVec presenceVec 
  /*
  {
    size_t i = 0 ;
    ContigKmerIterator kb1(&seqVec, &rankVec, k, 0);
    ContigKmerIterator ke1(&seqVec, &rankVec, k, seqVec.size() - k + 1);
    size_t contigId{0};
    int sampleCounter = 0 ;
    while(kb1 != ke1){
      my_mer r;
      r.word__(0) = *kb1;
      if(r.to_str() == "CAAGGACTCTTAGTCTCTCTGGGTCTTTTTA" or r.to_str() == "ATTTTTCTGGGTCTCTCTGATTCTCAGGAAC" or r.to_str() == "TAAAAAGACCCAGAGAGACTAAGAGTCCTTG" or r.to_str() == "GTTCCTGAGAATCAGAGAGACCCAGAAAAAT"){
        std::cerr << "should be set to 0" <<"\n" ;
        std::cerr << "local pos " << kb1.pos() <<"\n" ;
        std::exit(1) ;
      }
      kb1++ ;

    }

	std::cerr << "i = " << i
			 << " sampledKmers = " << sampledKmers << "\n" ;
	if(i != sampledKmers)
		std::exit(1) ;

    }*/


  // new presence Vec
  {
    std::cerr << "\nFilling presence Vector \n" ;

    size_t i = 0 ;
    ContigKmerIterator kb1(&seqVec, &rankVec, k, 0);
    ContigKmerIterator ke1(&seqVec, &rankVec, k, seqVec.size() - k + 1);
    size_t contigId{0};
    int sampleCounter = 0 ;

    //debug flags
    int loopCounter = 0;
    size_t ourKeys = 0 ;
    spp::sparse_hash_set<size_t> idxset;
    while(kb1 != ke1){
    	auto sampledPositions = allSampledPositions[contigId] ;
    	contigId++;
      size_t ite = 0;
      size_t skip = 0 ;


      //std::cerr<<loopCounter << ": "<<sampledPositions.size() << " contig length: "<< contigLengths[contigId-1]
      //       << " #kmers stored: "<< i << " #sampled kmers: "
      //       << sampledKmers << " #nkeys: " << nkeys << " #ourkeys: " << ourKeys <<"\n" ;
        loopCounter++ ;


      
        if(contigLengths[contigId-1] == 61){
          //for(auto sp : sampledPositions){ std::cout << sp << " " ; }
        }
        //std::cout<<"\n" ;

        my_mer r;

        auto zeroPos = kb1.pos();
        auto clen = contigLengths[contigId-1];
        auto nextSampIter = sampledPositions.begin();
        auto prevSamp = *nextSampIter;
        auto skipLen = kb1.pos() - zeroPos;
        bool didSample = false;
        bool done = false;
        
        for (size_t j = 0; j < clen - k + 1; ++kb1, ++j) {
          skipLen = kb1.pos() - zeroPos;
          if (!done and skipLen == *nextSampIter) {
            auto idx = bphf->lookup(*kb1);
            if (idxset.contains(idx)) {
              std::cerr << "that's some BS right there\n";
              r.word__(0) = *kb1;
              std::string theKmer = r.to_str();
              std::reverse(theKmer.begin(), theKmer.end());
              std::cerr <<  theKmer << "\n";
            }
            presenceVec[idx] = 1 ;
            idxset.insert(idx);
            i++ ;
            didSample = true;
            prevSamp = *nextSampIter;
            ++nextSampIter;
            if (nextSampIter == sampledPositions.end()) {
              done = true;
            }
          }


          r.word__(0) = *kb1;
          if(r.to_str() == "CAAGGACTCTTAGTCTCTCTGGGTCTTTTTA" or r.to_str() == "ATTTTTCTGGGTCTCTCTGATTCTCAGGAAC" or r.to_str() == "TAAAAAGACCCAGAGAGACTAAGAGTCCTTG" or r.to_str() == "GTTCCTGAGAATCAGAGAGACCCAGAAAAAT"){
            if (didSample) {
              std::cerr << "I skipped , I should not skip " << skipLen << " " << prevSamp <<"\n" ;
              std::cerr << "contig id and length " << contigId -1 << " " << contigLengths[contigId-1] << "\n" ;
              auto cseq = packedToString(seqVec,zeroPos,contigLengths[contigId-1]);
              std::cerr << "cotig: " << cseq <<"\n" ;
              //std::exit(1) ;
            }
          }

          didSample = false;
        }
        if (nextSampIter != sampledPositions.end()) {
          std::cerr << "I didn't sample " << std::distance(nextSampIter, sampledPositions.end()) << " samples for contig " << contigId - 1 << "\n";
          std::cerr << "last sample is " << sampledPositions.back() << "\n";
          std::cerr << "contig length is " << contigLengths[contigId-1] << "\n";
        }
        //++kb1;


        /*
      for(auto sp : sampledPositions){
        //std::cout<<"skip: "<<skip<<"\n" ;
      auto skipLen = kb1.pos() - zeroPos;
        while(skipLen < sp) {
          //skip++ ;
          r.word__(0) = *kb1;
          if(r.to_str() == "CAAGGACTCTTAGTCTCTCTGGGTCTTTTTA" or r.to_str() == "ATTTTTCTGGGTCTCTCTGATTCTCAGGAAC" or r.to_str() == "TAAAAAGACCCAGAGAGACTAAGAGTCCTTG" or r.to_str() == "GTTCCTGAGAATCAGAGAGACCCAGAAAAAT"){
            std::cerr << "I skipped , I should not skip " << skipLen << " " << sp <<"\n" ;
            std::cerr << "contig id and length " << contigId -1 << " " << contigLengths[contigId-1] << "\n" ;
            std::cerr << "cotig: " << packedToString(seqVec,zeroPos,contigLengths[contigId-1])<<"\n" ;
            std::exit(1) ;
          }
          if(!kb1.isEndKmer())
            kb1++ ;
          ourKeys++ ;
          skipLen = kb1.pos() - zeroPos;
        }
        auto idx = bphf->lookup(*kb1) ;

        if(idx >= presenceVec.size()){
          std::cerr << "i: " << i
                    << "idx: "<< idx
                    << "nkeys: " << presenceVec.size() << "\n" ;
        }
        r.word__(0) = *kb1;
        if(r.to_str() == "CAAGGACTCTTAGTCTCTCTGGGTCTTTTTA" or r.to_str() == "ATTTTTCTGGGTCTCTCTGATTCTCAGGAAC" or r.to_str() == "TAAAAAGACCCAGAGAGACTAAGAGTCCTTG" or r.to_str() == "GTTCCTGAGAATCAGAGAGACCCAGAAAAAT"){
          std::cerr << "should be set to 0" <<"\n" ;
          std::exit(1) ;
        }
        presenceVec[idx] = 1 ;
        i++ ;
        }

      while(!kb1.isEndKmer()){
        kb1++ ;
        ourKeys++ ;
      }
      kb1++ ;
      ourKeys++ ; */
    }

    std::cerr << "i = " << i
              << " sampledKmers = " << sampledKmers
              << " Loops = "<< loopCounter
              << " Contig array = "<<contigLengths.size()  
              << "\n" ;
    //if(i != sampledKmers)
    //  std::exit(1) ;

  }

  sdsl::bit_vector::rank_1_type realPresenceRank(&presenceVec) ;
  sdsl::bit_vector::select_1_type realPresenceSelect(&presenceVec) ;

  std::cerr << " num ones in presenceVec = " << realPresenceRank(presenceVec.size()-1) << "\n" ;
  //bidirectional version 1
  {

    ContigKmerIterator kb1(&seqVec, &rankVec, k, 0);
    ContigKmerIterator ke1(&seqVec, &rankVec, k, seqVec.size() - k + 1);

    size_t contigId{0} ;
    size_t coveredKeys{0} ;
    size_t totalKmersIshouldSee{0} ;

    // For every valid k-mer (i.e. every contig)
    while(kb1 != ke1){
      auto sampledPositions = allSampledPositions[contigId] ;
      auto thisContigLength = contigLengths[contigId] ;

      totalKmersIshouldSee += (thisContigLength - 31 + 1);

      contigId++ ;

      size_t skip = 0 ;

         my_mer r;

        auto zeroPos = kb1.pos();
        auto clen = contigLengths[contigId-1];
        auto nextSampIter = sampledPositions.begin();
        auto prevSampIter = sampledPositions.end();//*nextSampIter;
        auto skipLen = kb1.pos() - zeroPos;
        NextSampleDirection sampDir = NextSampleDirection::FORWARD;
        bool done = false;
        for (size_t j = 0; j < clen - k + 1; ++kb1, ++j) {
          int64_t nextSampPos = (nextSampIter != sampledPositions.end()) ? *nextSampIter : -1;
          int64_t prevSampPos = (prevSampIter != sampledPositions.end()) ? *prevSampIter : -1;
          uint64_t distToNext = (nextSampPos >= 0) ? nextSampPos - j : std::numeric_limits<uint64_t>::max();
          uint64_t distToPrev = (prevSampPos >= 0) ? j - prevSampPos : std::numeric_limits<uint64_t>::max();

          if (distToNext == std::numeric_limits<uint64_t>::max() and
              distToPrev == std::numeric_limits<uint64_t>::max()) {
            std::cerr << "We have fucked up royally\n";
          }

          sampDir = (distToNext < distToPrev) ? NextSampleDirection::FORWARD : NextSampleDirection::REVERSE;
          skipLen = kb1.pos() - zeroPos;
          // If this is a sampled position
          if (!done and skipLen == *nextSampIter) {
            prevSampIter = nextSampIter;
            ++nextSampIter;
            if (nextSampIter == sampledPositions.end()) {
              done = true;
            }
            auto idx = bphf->lookup(*kb1);
            auto rank = (idx == 0) ? 0 : realPresenceRank(idx);
            samplePosVec[rank] = kb1.pos();
          } else { // not a sampled position
            uint32_t ext = 0;
            size_t firstSampPos = 0;
            if (sampDir == NextSampleDirection::FORWARD) {
              firstSampPos = zeroPos + j + k;

              bool doPrint = false;
              my_mer r1 ;
              r1.word__(0) = *kb1 ;
              if(r1.to_str() == "ATGGTGACTGAATTCATTTTTCTGGGTCTCT" or r1.to_str() == "AGAGACCCAGAAAAATGAATTCAGTCACCAT" or r1.to_str() == "TCTCTGGGTCTTTTTACTTAAGTCAGTGGTA" or r1.to_str() == "TACCACTGACTTAAGTAAAAAGACCCAGAGA"){
                doPrint = true;
              }
              ext = getEncodedExtension(seqVec, firstSampPos, distToNext, extensionSize, sampDir, doPrint);
            } else if (sampDir == NextSampleDirection::REVERSE) {
              firstSampPos = zeroPos + prevSampPos;
              ext = getEncodedExtension(seqVec, firstSampPos, distToPrev, extensionSize, sampDir);
            } else {
              std::cerr << "go home, you're drunk!\n";
              std::exit(1);
            }
            auto idx = bphf->lookup(*kb1);
            auto rank = (idx == 0) ? 0 : realPresenceRank(idx);

            //debug print 
            my_mer r1 ;
            r1.word__(0) = *kb1 ;
            if(r1.to_str() == "ATGGTGACTGAATTCATTTTTCTGGGTCTCT" or r1.to_str() == "AGAGACCCAGAAAAATGAATTCAGTCACCAT" or r1.to_str() == "TCTCTGGGTCTTTTTACTTAAGTCAGTGGTA" or r1.to_str() == "TACCACTGACTTAAGTAAAAAGACCCAGAGA"){
              std::bitset<6> ext1(ext) ;
              std::string contigSeq = packedToString(seqVec, zeroPos, contigLengths[contigId-1]);
              std::cerr << "contig = " << contigSeq << "\n";
              std::cerr << "zeroPos = " << zeroPos << ", j = " << j << ", firstSampPos = " << firstSampPos << "\n";
              std::cerr << "backward nucl  " << ext1 << "\n" ;
              //std::exit(1) ;
            }

            auxInfo[idx - rank] = ext;
            direction[idx - rank] = (sampDir == NextSampleDirection::FORWARD) ? 1 : 0;
          }
        }
    }
    /////// ======== end of new code ========= //////

    // For every valid k-mer (i.e. every contig)
    /*
    while(kb1 != ke1){
      auto sampledPositions = allSampledPositions[contigId] ;
      auto thisContigLength = contigLengths[contigId] ;

      //std::cerr << "\nContig to cover: "<<contigId <<" last coveredKeys: "  <<coveredKeys
      //      << " ContigLength " << contigLengths[contigId]
      //      << " numSamples " << sampledPositions.size()
      //      << " numKmers: "<<nkeys << "\n" ;

      totalKmersIshouldSee += (thisContigLength - 31 + 1);

      contigId++ ;

      auto kbStampP = kb1 ; //previous kmer that is present
      auto kbStampN = kb1 ; //next kmer to look at unless contig is exhausted

      size_t skip = 0 ;

      if(sampledPositions.size() == 1){
        auto sampledKmerIndex = sampledPositions[0] ;
        if(sampledKmerIndex > 0){

          //fill up forward
          uint32_t extendNucl = seqVec.get_int(2*kb1.pos() + 2*k, 2*sampledKmerIndex) ;

          uint32_t appendNucl = extendNucl & 0x3 ;

          extendNucl = extendNucl >> 2 ;

          for(int j = 1 ; j < sampledKmerIndex ; ++j){
            appendNucl = appendNucl << 3 ;
            appendNucl = appendNucl | (extendNucl & 0x3) ;
            extendNucl = extendNucl >> 2 ;
          }

          auto rawAppendNucl = appendNucl ;

          if(sampledKmerIndex < extensionSize){
            appendNucl = appendNucl << 3 ;
            appendNucl = appendNucl | 0x4 ;
          }
          //pad it more in case it is notFound

          int e2 = sampledKmerIndex + 1 ;
          while(e2 < extensionSize){
            appendNucl = appendNucl << 3 ;
            e2++;
          }

          //add for the first extension 
          {
            auto tidx = bphf->lookup(*kb1);
            auto trank = realPresenceRank(tidx);
            auxInfo[tidx - trank] = appendNucl ;
            direction[tidx - trank] = 1 ;
            kb1++ ;
            coveredKeys++ ;

          }
          //append the rest on the way
          int gap = sampledKmerIndex ;
          gap = gap -1 ;
          

          while(gap > 0){
            auto tidx = bphf->lookup(*kb1) ;
            auto trank = realPresenceRank(tidx) ;
            appendNucl = appendNucl << 3 ;
            appendNucl = appendNucl | 0x4 ;
            gap = gap - 1 ;
            auxInfo[tidx - trank] = appendNucl ;
            direction[tidx - trank] = 1 ;
            kb1++ ;
            coveredKeys++ ;
          }

          auto idx = bphf->lookup(*kb1) ;
          auto rank = realPresenceRank(idx) ;
          if(presenceVec[idx] == 1){
            samplePosVec[rank] = kb1.pos() ;
          }else{
            std::cerr << "\nThis should not happen\n" ;
            std::exit(1) ;
          }

          if(!kb1.isEndKmer()){
            //fill up backward
            gap = 0 ;
            auto kbStart = kb1 ;
            while(!kb1.isEndKmer()){
              gap++ ;
              kb1++ ;
              coveredKeys++ ;
            }
            //I am standing at end kmer this nucliotide is not required

            uint32_t backwardNucl = seqVec.get_int(2*kbStart.pos(),2*gap) ;
            //make a 3 bit encoded format from it
            uint32_t appendNucl = backwardNucl & 0x3 ;
            backwardNucl = backwardNucl >> 2 ;
            auto skip = 1 ;
            while(skip < gap){
              appendNucl = appendNucl << 3 ;
              appendNucl = appendNucl | (backwardNucl & 0x3) ;
              backwardNucl = backwardNucl >> 2 ;
              skip++ ;
            }
            //let us move
            kbStart++ ;

            size_t innerInd = 0 ;

            while(innerInd < gap){
              size_t pushed = 0 ;
              auto appendNuclCopy = appendNucl ;

              auto shifts = gap - innerInd  - 1;
              while(shifts>0){
                appendNuclCopy = appendNuclCopy >> 3 ;
                shifts-- ;
              }

              auto numOfL = innerInd + 1 ;
              //rotate
              uint32_t appendNuclFinal{0} ;
              size_t e1 = 0 ;
              while(e1 < numOfL){
                appendNuclFinal = appendNuclFinal << 3 ;
                appendNuclFinal = appendNuclFinal | (appendNuclCopy & 0x7) ;
                appendNuclCopy = appendNuclCopy >> 3 ;
                e1++ ;
              }

              size_t e2 =  innerInd + 1;
              if(e2 < extensionSize){
                appendNuclFinal = appendNuclFinal << 3 ;
                appendNuclFinal = appendNuclFinal | 0x4 ;
                e2++ ;
              }
              while(e2 < extensionSize){
                appendNuclFinal = appendNuclFinal << 3 ;
                e2++ ;
              }

              auto tidx = bphf->lookup(*kbStart) ;
              auto trank = realPresenceRank(tidx) ;
              auxInfo[tidx - trank] = appendNuclFinal ;
              direction[tidx - trank] = 0 ;
              innerInd++;
              kbStart++ ;

            }

            if(kb1.isEndKmer()){ kb1++; coveredKeys++ ;}

          }else{
            kb1++ ;
            coveredKeys++ ;
          }
        }else{
          //The first kmer is stored

          auto idx = bphf->lookup(*kb1) ;
          auto rank = realPresenceRank(idx) ;

          if(presenceVec[idx] == 1){
            samplePosVec[rank] = kb1.pos() ;
          }else{
            std::cerr << "\nWTF\n" ;
            //std::cerr << "\nI should be on the first kmer\n" ;
          }
            //fill up backward
          size_t gap = 0 ;
          auto kbStart = kb1 ;
          while(!kb1.isEndKmer()){
            gap++ ;
            kb1++ ;
            coveredKeys++ ;
          }
          //I am standing at end kmer this nucliotide is not required
          //store backward

          uint32_t backwardNucl = seqVec.get_int(2*kbStart.pos(),2*gap) ;
          //make a 3 bit encoded format from it
          uint32_t appendNucl = backwardNucl & 0x3 ;
          backwardNucl = backwardNucl >> 2 ;
          auto skip = 1 ;
          while(skip < gap){
            appendNucl = appendNucl << 3 ;
            appendNucl = appendNucl | (backwardNucl & 0x3) ;
            backwardNucl = backwardNucl >> 2 ;
            skip++ ;
          }
          //let us move
          kbStart++ ;

          size_t innerInd = 0 ;

          while(innerInd < gap){
            size_t pushed = 0 ;
            auto appendNuclCopy = appendNucl ;

            auto shifts = gap - innerInd  - 1;
            while(shifts>0){
              appendNuclCopy = appendNuclCopy >> 3 ;
              shifts-- ;
            }

            auto numOfL = innerInd + 1 ;
            //rotate
            uint32_t appendNuclFinal{0} ;
            size_t e1 = 0 ;
            while(e1 < numOfL){

              appendNuclFinal = appendNuclFinal << 3 ;
              appendNuclFinal = appendNuclFinal | (appendNuclCopy & 0x7) ;
              appendNuclCopy = appendNuclCopy >> 3 ;
              e1++ ;
            }

            size_t e2 =  innerInd + 1;
            if(e2 < extensionSize){
              appendNuclFinal = appendNuclFinal << 3 ;
              appendNuclFinal = appendNuclFinal | 0x4 ;
              e2++ ;
            }
            while(e2 < extensionSize){
              appendNuclFinal = appendNuclFinal << 3 ;
              e2++ ;
            }

            auto tidx = bphf->lookup(*kbStart) ;
            auto trank = realPresenceRank(tidx) ;
            auxInfo[tidx - trank] = appendNuclFinal ;
            direction[tidx - trank] = 0 ;
            innerInd++;
            kbStart++ ;
           
          }

          if(kb1.isEndKmer()) { kb1++ ; coveredKeys++ ;}
        }

        continue ;

      }

      //otherwise we are in between two
      //present kmers

      auto sampledKmerIndex = sampledPositions[0] ;

      if(sampledKmerIndex > 0){

        //std::cerr<<"\nI entered extension before "<<sampledKmerIndex<<"\n" ;

        uint32_t extendNucl = seqVec.get_int(2*kb1.pos() + 2*k, 2*sampledKmerIndex) ;

        uint32_t appendNucl = extendNucl & 0x3 ;
        extendNucl = extendNucl >> 2 ;

        for(int j = 1 ; j < sampledKmerIndex ; ++j){
          appendNucl = appendNucl << 3 ;
          appendNucl = appendNucl | (extendNucl & 0x3) ;
          extendNucl = extendNucl >> 2 ;
        }

        auto rawAppendNucl = appendNucl ;

        if(sampledKmerIndex < extensionSize){
          appendNucl = appendNucl << 3 ;
          appendNucl = appendNucl | 0x4 ;
        }
        //pad it more in case it is notFound

        int e2 = sampledKmerIndex + 1 ;
        while(e2 < extensionSize){
          appendNucl = appendNucl << 3 ;
          e2++;
        }

        //add for the first extension 
        {
          auto tidx = bphf->lookup(*kb1);
          auto trank = realPresenceRank(tidx);
          auxInfo[tidx - trank] = appendNucl ;
          direction[tidx - trank] = 1 ;
          kb1++ ;
          coveredKeys++ ;

        }
        //append the rest on the way
        int gap = sampledKmerIndex ;
        gap = gap -1 ;


        while(gap > 0){
          auto tidx = bphf->lookup(*kb1) ;
          auto trank = realPresenceRank(tidx) ;
          appendNucl = appendNucl << 3 ;
          //auto thisAppendNucl = rawAppendNucl ;
          appendNucl = appendNucl | 0x4 ;
          gap = gap - 1 ;
          auxInfo[tidx - trank] = appendNucl ;
          direction[tidx - trank] = 1 ;
          kb1++ ;
          coveredKeys++;

        }
        //std::cerr<<"\nAfter  "<<sampledKmerIndex<<" kmers covered " << coveredKeys <<"\n" ;

      }

      //Now in between two present vector 
      
      for(size_t i = 0; i < sampledPositions.size() - 1; i++){
        auto idx = bphf->lookup(*kb1) ;
        auto rank = realPresenceRank(idx) ;

        auto prevIndex = sampledPositions[i] ;
        auto nextIndex = sampledPositions[i+1] ;

        auto jumpSize = nextIndex - prevIndex - 1 ;

        if(presenceVec[idx] == 1){
          samplePosVec[rank] = kb1.pos() ;
        }else{
          std::cerr<<"\nWTF\n" ;
          //std::cerr << "I should be standing at the sampled position in for loop \n" ;
          //std::exit(1) ;
        }

        //we have to fill backward
        {
          size_t gap = extensionSize ;
          auto kbStart = kb1 ;

          uint32_t backwardNucl = seqVec.get_int(2*kbStart.pos(),2*gap) ;

          //parent kmer = AGTTCCTGAGAATCAGAGAGACCCAGAAAAA
          //make a 3 bit encoded format from it
          my_mer r1 ;
          r1.word__(0) = *kbStart ;
          if(r1.to_str() == "AGTTCCTGAGAATCAGAGAGACCCAGAAAAA" or r1.to_str() == "AAAAAGACCCAGAGAGACTAAGAGTCCTTGA" or r1.to_str() == "TTTTTCTGGGTCTCTCTGATTCTCAGGAACT" or r1.to_str() == "TCAAGGACTCTTAGTCTCTCTGGGTCTTTTT"){
            std::bitset<6> ext1(backwardNucl) ;
            std::cerr << "backward nucl  " << ext1 << "\n" ;
            //std::exit(1) ;
          }

          uint32_t appendNucl = backwardNucl & 0x3 ;
          backwardNucl = backwardNucl >> 2 ;
          auto skip = 1 ;
          while(skip < gap){
            appendNucl = appendNucl << 3 ;
            appendNucl = appendNucl | (backwardNucl & 0x3) ;
            backwardNucl = backwardNucl >> 2 ;
            skip++ ;
          }
          //let us move
          kbStart++ ;


          size_t innerInd = 0 ;

          while(innerInd < gap){
            size_t pushed = 0 ;
            auto appendNuclCopy = appendNucl ;

            auto shifts = gap - innerInd  - 1;
            while(shifts>0){
              appendNuclCopy = appendNuclCopy >> 3 ;
              shifts-- ;
            }
            auto numOfL = innerInd + 1 ;
            //rotate
            uint32_t appendNuclFinal{0} ;
            size_t e1 = 0 ;
            while(e1 < numOfL){
              appendNuclFinal = appendNuclFinal << 3 ;
              appendNuclFinal = appendNuclFinal | (appendNuclCopy & 0x7) ;
              appendNuclCopy = appendNuclCopy >> 3 ;
              e1++ ;
            }

            size_t e2 =  innerInd + 1;
            if(e2 < extensionSize){
              appendNuclFinal = appendNuclFinal << 3 ;
              appendNuclFinal = appendNuclFinal | 0x4 ;
              e2++ ;
            }
            while(e2 < extensionSize){
              appendNuclFinal = appendNuclFinal << 3 ;
              e2++ ;
            }

            my_mer r ;
            r.word__(0) = *kbStart ;
            if(r.to_str() == "GTTCCTGAGAATCAGAGAGACCCAGAAAAAT" or r.to_str() == "TAAAAAGACCCAGAGAGACTAAGAGTCCTTG" or r.to_str() == "ATTTTTCTGGGTCTCTCTGATTCTCAGGAAC" or r.to_str() == "CAAGGACTCTTAGTCTCTCTGGGTCTTTTTA"){
              std::bitset<6> ext1(appendNuclFinal) ;
              std::cerr << "ext in for: " << ext1 << "\n" ;
              std::cerr << "gap: " << gap << " InnerInd " << innerInd << "\n" ;
              //std::exit(1) ;
            }



            auto tidx = bphf->lookup(*kbStart) ;
            auto trank = realPresenceRank(tidx) ;
            auxInfo[tidx - trank] = appendNuclFinal ;
            direction[tidx - trank] = 0 ;
            innerInd++;
            kbStart++ ;
           
          }


          //forward pointer by extensionSize
          size_t j = 0 ;
          while(j < extensionSize) {kb1++ ; j++ ; coveredKeys++ ;}
          //one more to stand on the neucl of interest
          kb1++ ;
          coveredKeys++ ;

        }
        //we have to fill forward
        if(jumpSize > extensionSize)
        {
          //fill up forward
          size_t skip = jumpSize - extensionSize ;
          uint32_t extendNucl = seqVec.get_int(2*kb1.pos() + 2*k, 2*skip) ;

          my_mer r;
          r.word__(0) = *kb1;
          if(r.to_str() == "TGGCGCAGGCTGGGTGGAGCCGTCCCCCCAT" or r.to_str() == "TACCCCCCTGCCGAGGTGGGTCGGACGCGGT" or r.to_str() == "ATGGGGGGACGGCTCCACCCAGCCTGCGCCA" or r.to_str() == "ACCGCGTCCGACCCACCTCGGCAGGGGGGTA"){
            std::bitset<4> ext1(extendNucl) ;
            std::cerr << "ext: " << ext1 << "\n" ;
            //std::exit(1) ;
          }


          uint32_t appendNucl = extendNucl & 0x3 ;
          extendNucl = extendNucl >> 2 ;

          for(int j = 1 ; j < skip ; ++j){
            appendNucl = appendNucl << 3 ;
            appendNucl = appendNucl | (extendNucl & 0x3) ;
            extendNucl = extendNucl >> 2 ;
          }

          if(r.to_str() == "TGGCGCAGGCTGGGTGGAGCCGTCCCCCCAT" or r.to_str() == "TACCCCCCTGCCGAGGTGGGTCGGACGCGGT" or r.to_str() == "ATGGGGGGACGGCTCCACCCAGCCTGCGCCA" or r.to_str() == "ACCGCGTCCGACCCACCTCGGCAGGGGGGTA"){
            //if(r.to_str() == "AGAGAGACCCAGAAAAATGAATTCAGTCACC" or r.to_str() == "CCACTGACTTAAGTAAAAAGACCCAGAGAGA" or r.to_str() == "GGTGACTGAATTCATTTTTCTGGGTCTCTCT" or r.to_str() == "TCTCTCTGGGTCTTTTTACTTAAGTCAGTGG"){
            std::bitset<6> ext1(appendNucl) ;
            std::cerr << "ext: " << ext1 << "\n" ;
            //std::exit(1) ;
          }

          auto rawAppendNucl = appendNucl ;

          //add pad for the first extension
          if(skip < extensionSize){
            appendNucl = appendNucl << 3 ;
            appendNucl = appendNucl | 0x4 ;
          }
          //pad it more in case it is notFound

          size_t e2 = skip + 1 ;
          while(e2 < extensionSize){
            appendNucl = appendNucl << 3 ;
            e2++;
          }

          {
            auto tidx = bphf->lookup(*kb1);
            auto trank = realPresenceRank(tidx);
            auxInfo[tidx - trank] = appendNucl ;
            direction[tidx - trank] = 1 ;
            kb1++ ;
            coveredKeys++ ;

          }

          if(r.to_str() == "TGGCGCAGGCTGGGTGGAGCCGTCCCCCCAT" or r.to_str() == "TACCCCCCTGCCGAGGTGGGTCGGACGCGGT" or r.to_str() == "ATGGGGGGACGGCTCCACCCAGCCTGCGCCA" or r.to_str() == "ACCGCGTCCGACCCACCTCGGCAGGGGGGTA"){
            //if(r.to_str() == "AGAGAGACCCAGAAAAATGAATTCAGTCACC" or r.to_str() == "CCACTGACTTAAGTAAAAAGACCCAGAGAGA" or r.to_str() == "GGTGACTGAATTCATTTTTCTGGGTCTCTCT" or r.to_str() == "TCTCTCTGGGTCTTTTTACTTAAGTCAGTGG"){
            std::bitset<6> ext1(appendNucl) ;
            std::cerr << "ext: " << ext1 << "\n" ;
            //std::exit(1) ;
          }
          //append the rest on the way
          size_t gap = skip ;
          gap = gap -1 ;
          

          while(gap > 0){
            auto tidx = bphf->lookup(*kb1) ;
            auto trank = realPresenceRank(tidx) ;

            r.word__(0) = *kb1 ;

            appendNucl = appendNucl << 3 ;
            appendNucl = appendNucl | 0x4 ;
            gap = gap - 1 ;
            auxInfo[tidx - trank] = appendNucl ;

            if(r.to_str() == "TGGCGCAGGCTGGGTGGAGCCGTCCCCCCAT" or r.to_str() == "TACCCCCCTGCCGAGGTGGGTCGGACGCGGT" or r.to_str() == "ATGGGGGGACGGCTCCACCCAGCCTGCGCCA" or r.to_str() == "ACCGCGTCCGACCCACCTCGGCAGGGGGGTA"){
              //            if(r.to_str() == "GAGAGACCCAGAAAAATGAATTCAGTCACCA" or r.to_str() == "ACCACTGACTTAAGTAAAAAGACCCAGAGAG" or r.to_str() == "TGGTGACTGAATTCATTTTTCTGGGTCTCTC" or r.to_str() == "CTCTCTGGGTCTTTTTACTTAAGTCAGTGGT"){
              std::bitset<6> ext1(appendNucl) ;
              std::cerr << "ext: " << ext1 << "\n" ;
              //std::exit(1) ;
            }
            direction[tidx - trank] = 1 ;
            kb1++ ;
            coveredKeys++ ;
          }

        }

       }


      //std::cerr<<"\nAfter  "<<sampledPositions[sampledPositions.size()-2]<<" kmers covered " << coveredKeys <<"\n" ;
      //this is the last sampled kmer
      //The last remaining battle 
      if(!kb1.isEndKmer()){
          //store the current kmer

          //have to store backward

          auto idx = bphf->lookup(*kb1) ;
          auto rank = realPresenceRank(idx) ;

          if(presenceVec[idx] == 1){
            samplePosVec[rank] = kb1.pos() ;
          }else{
            std::cerr << "\nWTF\n" ;
            //std::cerr << "sampled kmer is not end kmer \n" ;

            //std::exit(1) ;
          }

          //fill up backward
          size_t gap = 0 ;
          auto kbStart = kb1 ;
          while(!kb1.isEndKmer()){
            gap++ ;
            kb1++ ;
            coveredKeys++ ;
          }
          //I am standing at end kmer this nucliotide is not required

          uint32_t backwardNucl = seqVec.get_int(2*kbStart.pos(),2*gap) ;
          //make a 3 bit encoded format from it
          uint32_t appendNucl = backwardNucl & 0x3 ;
          backwardNucl = backwardNucl >> 2 ;
          auto skip = 1 ;
          while(skip < gap){
            appendNucl = appendNucl << 3 ;
            appendNucl = appendNucl | (backwardNucl & 0x3) ;
            backwardNucl = backwardNucl >> 2 ;
            skip++ ;
          }
          //let us move
          kbStart++ ;

          size_t innerInd = 0 ;

          while(innerInd < gap){
            size_t pushed = 0 ;
            auto appendNuclCopy = appendNucl ;

            auto shifts = gap - innerInd  - 1;
            while(shifts>0){
              appendNuclCopy = appendNuclCopy >> 3 ;
              shifts-- ;
            }

            auto numOfL = innerInd + 1 ;
            //rotate
            uint32_t appendNuclFinal{0} ;
            size_t e1 = 0 ;
            while(e1 < numOfL){
              appendNuclFinal = appendNuclFinal << 3 ;
              appendNuclFinal = appendNuclFinal | (appendNuclCopy & 0x7) ;
              appendNuclCopy = appendNuclCopy >> 3 ;
              e1++ ;
            }

            size_t e2 =  innerInd + 1;
            if(e2 < extensionSize){
              appendNuclFinal = appendNuclFinal << 3 ;
              appendNuclFinal = appendNuclFinal | 0x4 ;
              e2++ ;
            }
            while(e2 < extensionSize){
              appendNuclFinal = appendNuclFinal << 3 ;
              e2++ ;
            }

            auto tidx = bphf->lookup(*kbStart) ;
            auto trank = realPresenceRank(tidx) ;
            auxInfo[tidx - trank] = appendNuclFinal ;
            direction[tidx - trank] = 0 ;
            innerInd++;
            kbStart++ ;
           
          }

          if(kb1.isEndKmer()) {kb1++ ; coveredKeys++ ;}
      }else{
        auto idx = bphf->lookup(*kb1) ;
        auto rank = realPresenceRank(idx) ;

        if(presenceVec[idx] == 1){
          samplePosVec[rank] = kb1.pos() ;
        }else{
          std::cerr << "\nWTF\n" ;
          //std::cerr << "The end kmer \n" ;
          //std::exit(1) ;
        }
        kb1++ ;
        coveredKeys++ ;
        
      }
      //std::cerr<<"\nAfter  "<<sampledPositions[sampledPositions.size()-1]<<" kmers covered " << coveredKeys <<"\n" ;
      

    }
    */



  }

  //bidirectional sampling
  /*
  {
    ContigKmerIterator kb1(&seqVec, &rankVec, k, 0);
    ContigKmerIterator ke1(&seqVec, &rankVec, k, seqVec.size() - k + 1);


    while(kb1 != ke1){
      //There are two cases
      //I will make sure
      //that I will traverse the entire
      //contig before I end the iteration
       
      auto kbStampP = kb1 ; //previous kmer that is present
      auto kbStampN = kb1 ; //next kmer to look at unless contig is exhausted

      while(!kb1.isEndKmer()){
        //case 1: kmer is present
        auto idx = bphf->lookup(*kb1) ;
        auto rank = realPresenceRank(idx) ;
        if(presenceVec[idx] == 1){
          kbStampP = kb1 ;
          samplePosVec[rank] = kb1.pos() ;
          kb1++ ;
        }else{
          //case 2: kmer not present
          //We are going to fill the
          //extension before we see another present kmer
          int gap = 0 ;
          auto kbStart = kb1 ; //start of the absent kmer

          // measure the gap to the next present kmer

          while(!kb1.isEndKmer()){
            auto tidx = bphf->lookup(*kb1) ;
            if(presenceVec[tidx] == 1){
              kbStampN = kb1 ;
              break ;
            }
            gap++ ;
            kb1++ ;
          }

          if(gap > extensionSize){
            //there has to be a previous and next present kmer
            int gap2 = gap - extensionSize ;

            //half the way will be filled up backward
            {
              uint32_t extendNucl = seqVec.get_int(2*kbStampP.pos(), 2*extensionSize) ;
              auto kbIt = kbStart ;
              for(int j = 0 ; j < extensionSize ; j++){
                int i = j ;
                auto extendNuclCopy = extendNucl  ;
                uint32_t appendNucl{0} ;
                int pushed = 0;
                while(i >= 0){
                  appendNucl = appendNucl | (extendNuclCopy & 0x3) ;
                  appendNucl = appendNucl << 3 ;
                  extendNuclCopy = extendNuclCopy >> 2 ;
                  i-- ;
                  pushed++ ;
                }
                
                if(pushed < extensionSize){
                  appendNucl = appendNucl << 3 ;
                  appendNucl = appendNucl | 0x4 ;
                  pushed++ ;
                }
                while(pushed < extensionSize){
                  appendNucl = appendNucl << 3 ;
                  pushed++ ;
                }

                auto tidx = bphf->lookup(*kbIt) ;
                auto trank = realPresenceRank(tidx) ;
                auxInfo[tidx - trank] = appendNucl ;
                direction[tidx - trank] = 0 ;
                kbIt++ ;

              }
              kbStart = kbIt ;
              
            }
            //halfWay would be filled up forward
            {

              uint32_t extendNucl = seqVec.get_int(2*kbStart.pos() + 2*k, 2*gap2) ;
              uint32_t appendNucl = extendNucl & 0x3 ;
              
              extendNucl = extendNucl >> 2 ;
              for(int j = 1 ; j < gap2 ; ++j){
                appendNucl = appendNucl << 3 ;
                appendNucl = appendNucl | (extendNucl & 0x3) ;
                extendNucl = extendNucl >> 2 ;
              }
              bool delm = false ;
              auto rawAppendNucl = appendNucl ;
              int skip = gap2 ;

              if(skip < extensionSize){
                delm = true ;
                appendNucl = appendNucl << 3 ;
                appendNucl = appendNucl | 0x4 ;
              }
              //pad it more in case it is notFound

              int e2 = gap + 1 ;
              while(e2 < extensionSize){
                appendNucl = appendNucl << 3 ;
                e2++;
              }

              //add for the first extension 
              {
                auto tidx = bphf->lookup(*kbStart);
                auto trank = realPresenceRank(tidx);
                auxInfo[tidx - trank] = appendNucl ;
                direction[tidx - trank] = 1 ;
                kbStart++ ;

              }
              //append the rest on the way
              while(kbStart != kb1){
                auto tidx = bphf->lookup(*kbStart) ;
                auto trank = realPresenceRank(tidx) ;
                skip = skip - 1;
                rawAppendNucl = rawAppendNucl >> 3 ;
                auto thisAppendNucl = rawAppendNucl ;
                thisAppendNucl = thisAppendNucl << 3 ;
                thisAppendNucl = thisAppendNucl | 0x4 ;
                int e3 = skip + 1 ;
                while(e3 < extensionSize){
                  thisAppendNucl = thisAppendNucl << 3 ;
                  e3++ ;
                }
                auxInfo[tidx - trank] = thisAppendNucl ;
                direction[tidx - trank] = 1 ;
                kbStart++ ;
              }


            }

          }else if(gap > 0){
            //this is the boundary case
            //if I am standing on the last kmer
            //then the present kmer is before me
            if(!kb1.isEndKmer()){

              //I am standing at the start of
              //the contig

              uint32_t extendNucl = seqVec.get_int(2*kbStart.pos() + 2*k, 2*gap) ;
              uint32_t appendNucl = extendNucl & 0x3 ;
              extendNucl = extendNucl >> 2 ;
              for(int j = 1 ; j < gap ; ++j){
                appendNucl = appendNucl << 3 ;
                appendNucl = appendNucl | (extendNucl & 0x3) ;
                extendNucl = extendNucl >> 2 ;
              }
              bool delm = false ;
              auto rawAppendNucl = appendNucl ;

              if(gap  < extensionSize){
                delm = true ;
                appendNucl = appendNucl << 3 ;
                appendNucl = appendNucl | 0x4 ;
              }
              //pad it more in case it is notFound

              int e2 = gap + 1 ;
              while(e2 < extensionSize){
                appendNucl = appendNucl << 3 ;
                e2++;
              }

              //add for the first extension 
              {
                auto tidx = bphf->lookup(*kbStart);
                auto trank = realPresenceRank(tidx);
                auxInfo[tidx - trank] = appendNucl ;
                direction[tidx - trank] = 1 ;
                kbStart++ ;

              }
              //append the rest on the way
              int skip = gap ;
              while(kbStart != kb1){
                auto tidx = bphf->lookup(*kbStart) ;
                auto trank = realPresenceRank(tidx) ;
                skip = skip - 1;
                rawAppendNucl = rawAppendNucl >> 3 ;
                auto thisAppendNucl = rawAppendNucl ;
                thisAppendNucl = thisAppendNucl << 3 ;
                thisAppendNucl = thisAppendNucl | 0x4 ;
                int e3 = skip + 1 ;
                while(e3 < extensionSize){
                  thisAppendNucl = thisAppendNucl << 3 ;
                  e3++ ;
                }
                auxInfo[tidx - trank] = thisAppendNucl ;
                direction[tidx - trank] = 1 ;
                kbStart++ ;
              }

              
            }else{
              //I am standing at the end kmer
              //Here we have to add previous if
              //presenceVec is not 1 on this bit
              //this case is handled at the end
              //two case
              //case 1 : last kmer is stored
              auto tidx = bphf->lookup(*kb1) ;
              auto trank = realPresenceRank(tidx) ;
              if(presenceVec[tidx] == 1){
                //store things forward
                samplePosVec[trank] = kb1.pos() ;
                uint32_t extendNucl = seqVec.get_int(2*kbStart.pos() + 2*k, 2*gap) ;
                uint32_t appendNucl = extendNucl & 0x3 ;
                extendNucl = extendNucl >> 2 ;
                for(int j = 1 ; j < gap ; ++j){
                  appendNucl = appendNucl << 3 ;
                  appendNucl = appendNucl | (extendNucl & 0x3) ;
                  extendNucl = extendNucl >> 2 ;
                }
                bool delm = false ;
                auto rawAppendNucl = appendNucl ;

                if(gap < extensionSize){
                  delm = true ;
                  appendNucl = appendNucl << 3 ;
                  appendNucl = appendNucl | 0x4 ;
                }
                //pad it more in case it is notFound

                int e2 = gap + 1 ;
                while(e2 < extensionSize){
                  appendNucl = appendNucl << 3 ;
                  e2++;
                }

                //add for the first extension 
                {
                  auto tidx = bphf->lookup(*kbStart);
                  auto trank = realPresenceRank(tidx);
                  auxInfo[tidx - trank] = appendNucl ;
                  direction[tidx - trank] = 1 ;
                  kbStart++ ;

                }
                //append the rest on the way
                //int skip = gap ;
                while(kbStart != kb1){
                  auto tidx = bphf->lookup(*kbStart) ;
                  auto trank = realPresenceRank(tidx) ;
                  gap = gap - 1;
                  rawAppendNucl = rawAppendNucl >> 3 ;
                  auto thisAppendNucl = rawAppendNucl ;
                  thisAppendNucl = thisAppendNucl << 3 ;
                  thisAppendNucl = thisAppendNucl | 0x4 ;
                  int e3 = gap + 1 ;
                  while(e3 < extensionSize){
                    thisAppendNucl = thisAppendNucl << 3 ;
                    e3++ ;
                  }
                  auxInfo[tidx - trank] = thisAppendNucl ;
                  direction[tidx - trank] = 1 ;
                  kbStart++ ;
                }
                kb1++ ;
                break ;

              }else{
                //store things backward
                uint32_t extendNucl = seqVec.get_int(2*kbStampP.pos(), 2*extensionSize) ;
                auto kbIt = kbStart ;
                for(int j = 0 ; j < extensionSize ; j++){
                  int i = j ;
                  auto extendNuclCopy = extendNucl  ;
                  uint32_t appendNucl{0} ;
                  int pushed = 0;
                  while(i >= 0){
                    appendNucl = appendNucl | (extendNuclCopy & 0x3) ;
                    appendNucl = appendNucl << 3 ;
                    extendNuclCopy = extendNuclCopy >> 2 ;
                    i-- ;
                    pushed++ ;
                  }

                  if(pushed < extensionSize){
                    appendNucl = appendNucl << 3 ;
                    appendNucl = appendNucl | 0x4 ;
                    pushed++ ;
                  }
                  while(pushed < extensionSize){
                    appendNucl = appendNucl << 3 ;
                    pushed++ ;
                  }

                  auto tidx = bphf->lookup(*kbIt) ;
                  auto trank = realPresenceRank(tidx) ;
                  auxInfo[tidx - trank] = appendNucl ;
                  direction[tidx - trank] = 0 ;
                  kbIt++ ;
                }
                kb1++ ;
                break ;
             }
          }

        }
      } 

      }//inner while loop traversing over one contig 

    }//end of outer while each loop = a contig 

    }*/// end of scope

  /*
  {
    size_t i = 0;
    ContigKmerIterator kb1(&seqVec, &rankVec, k, 0);
    ContigKmerIterator ke1(&seqVec, &rankVec, k, seqVec.size() - k + 1);
    int sampleCounter = 0;

    while(kb1 != ke1){
        auto kbStamp = kb1 ;

    	auto idx = bphf->lookup(*kb1) ;
    	auto rank = realPresenceRank(idx) ;

        if(presenceVec[idx] == 1){
			  if (idx >= presenceVec.size()) {
				std::cerr << "i =  " << i << ", size = " << seqVec.size()
						  << ", idx = " << idx << ", size = " << presenceVec.size() << "\n";
			  }

#ifdef PUFFER_DEBUG
				std::cerr << "i =  " << i << ", size = " << seqVec.size()
						  << ", rank = " << rank << ", size = " << samplePosVec.size() << "\n";
#endif
				samplePosVec[rank] = kb1.pos() ;

    		if(idx < rank){
    			std::cerr << "idx = " << idx
    					<< "rank = " << rank << ", size = " << presenceVec.size() << "\n" ;
    		}
    		kb1++ ;
    		i++ ;
    	}else{
            bool debugExt = false ;

			auto kbIt = kb1 ;
			int extendLength = 0;
			while(extendLength < extensionSize){
				auto nIdx = bphf->lookup(*kb1) ;
                if(kb1.isEndKmer()){
                    //std::cerr<<"breaking b/c end kmer "<<"\n" ;
					break ;
                }
                if(presenceVec[nIdx] == 1){
                    //std::cerr<<"breaking b/c end presenceVec "<<"\n" ;
					break ;
                }
				extendLength++;
				kb1++;
			}
			uint32_t extendNucl ;
			if(extendLength > 0){
				extendNucl = seqVec.get_int(2 * kbIt.pos() + 2* k, 2*extendLength) ;

                
                //std::bitset<32> ext1(extendNucl);
                //std::cout<<"extension bits "<<ext1<<"\n" ;
                //std::cout<<"extension length "<<extendLength<<"\n" ;
                
                auto extendNucl1 = extendNucl ;

				uint32_t appendNucl = extendNucl & 0x3 ;
				extendNucl = extendNucl >> 2 ;
				//do something like store these with proper encoding
				for(int j = 1; j < extendLength ; ++j){
					appendNucl = appendNucl << 3 ;
					appendNucl = appendNucl | (extendNucl & 0x3) ;
					extendNucl = extendNucl >> 2 ;
				}
				//add delimeter
				if(extendLength < extensionSize){
					appendNucl = appendNucl << 3 ;
					appendNucl = appendNucl | 0x4 ;
				}
                int e2 = extendLength + 1;
                while(e2 < extensionSize){
                    appendNucl = appendNucl << 3;
                    e2++;
                }

                
                //std::bitset<16> ext2(appendNucl);
                //std::cout<<"appended bits "<<ext2<<"\n" ;
                //std::exit(1) ;
                
                if(extendNucl1 != 0x0 and appendNucl == 0x0){
                    std::cerr<<"this should not happen\n";
                }

				if(idx < rank){
					std::cerr << "idx = " << idx
							<< "rank = " << rank << ", size = " << presenceVec.size() << "\n" ;
				}
#ifdef PUFFER_DEBUG
				std::cerr << "i =  " << i << ", size = " << seqVec.size()
						  << ", idx - rank = " << (idx - rank) << ", size = " << auxInfo.size() << "\n";
#endif
				auxInfo[idx - rank] = appendNucl ;
			}
            kbStamp++ ;
            kb1 = kbStamp ;
    	}

  // validate
#ifdef PUFFER_DEBUG
      uint64_t kn = seqVec.get_int(2 * kb1.pos(), 2 * k);
      CanonicalKmer sk;
      sk.fromNum(kn);
      if (sk.isEquivalent(*kb1) == KmerMatchType::NO_MATCH) {
        my_mer r;
        r.word__(0) = *kb1;
        std::cerr << "I thought I saw " << sk.to_str() << ", but I saw " << r.to_str() << "\n";
      }
#endif
    }

	std::cerr << "i = " << i
			 << " sampledKmers = " << sampledKmers << "\n" ;
	if(i != sampledKmers)
		std::exit(1) ;
  }
  */






  /** Write the index **/
  std::ofstream descStream(outdir + "/info.json");
  {
    cereal::JSONOutputArchive indexDesc(descStream);
    std::string sampStr = "Sparse";
    indexDesc(cereal::make_nvp("sampling_type", sampStr));
    indexDesc(cereal::make_nvp("sample_size", sampleSize));
    indexDesc(cereal::make_nvp("extension_size", extensionSize));
    indexDesc(cereal::make_nvp("k", k));
    indexDesc(cereal::make_nvp("num_kmers", nkeys));
    indexDesc(cereal::make_nvp("num_sampled_kmers",sampledKmers));
    indexDesc(cereal::make_nvp("num_contigs", numContigs));
    indexDesc(cereal::make_nvp("seq_length", tlen));
  }
  descStream.close();

  //sdsl::store_to_file(posVec, outdir + "/pos.bin");
  std::ofstream hstream(outdir + "/mphf.bin");
  sdsl::store_to_file(presenceVec, outdir + "/presence.bin");
  sdsl::store_to_file(samplePosVec, outdir + "/sample_pos.bin");
  sdsl::store_to_file(auxInfo, outdir + "/extension.bin");
  sdsl::store_to_file(canonicalNess, outdir + "/canonical.bin");
  sdsl::store_to_file(direction, outdir + "/direction.bin");
  bphf->save(hstream);
  hstream.close();

   return 0;

}
