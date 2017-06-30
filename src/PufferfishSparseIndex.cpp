#include <fstream>
#include <iostream>
#include <bitset>

#include "cereal/archives/json.hpp"
#include "cereal/archives/binary.hpp"
#include "CLI/Timer.hpp"

#include "PufferfishSparseIndex.hpp"
#include "PufferFS.hpp"

PufferfishSparseIndex::PufferfishSparseIndex() {}

PufferfishSparseIndex::PufferfishSparseIndex(const std::string& indexDir) {
  if (!puffer::fs::DirExists(indexDir.c_str())) {
    std::cerr << "The index directory " << indexDir << " does not exist!\n";
    std::exit(1);
  }

  {
    std::ifstream infoStream(indexDir + "/info.json");
    cereal::JSONInputArchive infoArchive(infoStream);
    infoArchive(cereal::make_nvp("k", k_));
    infoArchive(cereal::make_nvp("num_kmers", numKmers_));
    std::cerr << "k = " << k_ << '\n';
    std::cerr << "num kmers = " << numKmers_ << '\n';
    infoStream.close();
  }

  //std::cerr << "loading contig table ... ";
  {
    CLI::AutoTimer timer {"Loading contig table", CLI::Timer::Big};
    std::ifstream contigTableStream(indexDir + "/ctable.bin");
    cereal::BinaryInputArchive contigTableArchive(contigTableStream);
    contigTableArchive(refNames_);
    contigTableArchive(contigTable_);
    contigTableStream.close();
  }
  numContigs_ = contigTable_.size();

  {
    CLI::AutoTimer timer {"Loading eq table", CLI::Timer::Big};
    std::ifstream eqTableStream(indexDir + "/eqtable.bin");
    cereal::BinaryInputArchive eqTableArchive(eqTableStream);
    eqTableArchive(eqClassIDs_);
    eqTableArchive(eqLabels_);
    eqTableStream.close();
  }
  //std::cerr << "done\n";

  {
    CLI::AutoTimer timer {"Loading mphf table", CLI::Timer::Big};
    std::string hfile = indexDir + "/mphf.bin";
    std::ifstream hstream(hfile);
    hash_.reset(new boophf_t);
    hash_->load(hstream);
    hstream.close();
  }

  {
    CLI::AutoTimer timer {"Loading contig boundaries", CLI::Timer::Big};
    std::string bfile = indexDir + "/rank.bin";
    sdsl::load_from_file(contigBoundary_, bfile);
    contigRank_ = decltype(contigBoundary_)::rank_1_type(&contigBoundary_);
    contigSelect_ = decltype(contigBoundary_)::select_1_type(&contigBoundary_);
  }

  {
    CLI::AutoTimer timer {"Loading sequence", CLI::Timer::Big};
    std::string sfile = indexDir + "/seq.bin";
    sdsl::load_from_file(seq_, sfile);
  }
  /*
  {
    CLI::AutoTimer timer {"Loading positions", CLI::Timer::Big};
    std::string pfile = indexDir + "/pos.bin";
    sdsl::load_from_file(pos_, pfile);
  }*/

  {
    CLI::AutoTimer timer {"Loading presence vector", CLI::Timer::Big};
    std::string bfile = indexDir + "/presence.bin";
    sdsl::load_from_file(presenceVec_, bfile);
    presenceRank_ = decltype(presenceVec_)::rank_1_type(&presenceVec_);
    presenceSelect_ = decltype(presenceVec_)::select_1_type(&presenceVec_);
  }
  {
    CLI::AutoTimer timer {"Loading canonical vector", CLI::Timer::Big};
    std::string pfile = indexDir + "/canonical.bin";
    sdsl::load_from_file(canonicalNess_, pfile);
  }
  {
    CLI::AutoTimer timer {"Loading sampled positions", CLI::Timer::Big};
    std::string pfile = indexDir + "/sample_pos.bin";
    sdsl::load_from_file(sampledPos_, pfile);
  }

  {
    CLI::AutoTimer timer {"Loading extension vector", CLI::Timer::Big};
    std::string pfile = indexDir + "/extension.bin";
    sdsl::load_from_file(auxInfo_, pfile);
  }

  {
    CLI::AutoTimer timer {"Loading direction vector", CLI::Timer::Big};
    std::string pfile = indexDir + "/direction.bin";
    sdsl::load_from_file(directionVec_, pfile);
  }

}

PufferfishSparseIndex::EqClassID PufferfishSparseIndex::getEqClassID(uint32_t contigID) {
  return eqClassIDs_[contigID];
}

const PufferfishSparseIndex::EqClassLabel& PufferfishSparseIndex::getEqClassLabel(uint32_t contigID) {
  return eqLabels_[getEqClassID(contigID)];
}

//auto endContigMap() -> decltype(contigTable_.begin()) { return contigTable_.end(); }
uint64_t PufferfishSparseIndex::getRawPos(CanonicalKmer& mer)  {
  auto km = mer.getCanonicalWord();
  size_t res = hash_->lookup(km);
  uint64_t pos =
    (res < numKmers_) ? pos_[res] : std::numeric_limits<uint64_t>::max();
  if (pos <= seq_.size() - k_) {
    uint64_t fk = seq_.get_int(2 * pos, 2 * k_);
    my_mer fkm;
    fkm.word__(0) = fk;
    auto keq = mer.isEquivalent(fkm);
    if (keq != KmerMatchType::NO_MATCH) {
      //}mer.fwWord() == fkm.word(0) or mer.rcWord() == fkm.word(0)) {
      return pos;
    }
  }
  pos = std::numeric_limits<uint64_t>::max();
  return pos;
}

bool PufferfishSparseIndex::contains(CanonicalKmer& mer) {
  return isValidPos(getRawPos(mer));
}

bool PufferfishSparseIndex::isValidPos(uint64_t pos) {
  return pos != std::numeric_limits<uint64_t>::max();
}

uint32_t PufferfishSparseIndex::contigID(CanonicalKmer& mer) {
    auto km = mer.getCanonicalWord();
    size_t res = hash_->lookup(km);
    uint64_t pos =
      (res < numKmers_) ? pos_[res] : std::numeric_limits<uint64_t>::max();
    if (pos <= seq_.size() - k_) {
      uint64_t fk = seq_.get_int(2 * pos, 2 * k_);
      my_mer fkm;
      fkm.word__(0) = fk;
      auto keq = mer.isEquivalent(fkm);
      if (keq != KmerMatchType::NO_MATCH) {
        auto rank = contigRank_(pos);
        return rank;
      }
    }
    return std::numeric_limits<uint32_t>::max();
  }

auto PufferfishSparseIndex::getRefPos(CanonicalKmer mern) -> util::ProjectedHits {

  using IterT = std::vector<util::Position>::iterator;//decltype(contigTable_.begin());

  auto km = mern.getCanonicalWord() ;
  CanonicalKmer mer;
  mer.fromNum(mern.getCanonicalWord()) ;
  //bool searchCanon = (mer.fwWord() == mer.getCanonicalWord()) ;
  size_t idx = hash_->lookup(km);
  int extensionSize = 2;
  //here the logic for searching the sparse
  //index comes
        //two options, either found or not
	uint64_t pos{0} ;
    //hack
	//uint64_t posT =
	//	   (idx < numKmers_) ? pos_[idx] : std::numeric_limits<uint64_t>::max();

    /*
    if(mern.to_str() == "CTTAGTCTCTCTGGGTCTTTTTACTTAAGTC")
	    std::cerr << " real position for kmer = " << mer.to_str() << " Canonicalness = "<< canonicalNess_[idx]<< " pos " << posT << "\n" ;
    */
	auto currRank = (idx == 0) ? 0 : presenceRank_(idx) ;

	//std::cerr << "idx = " << idx << ", size = " << presenceVec_.size() << "\n";

  //if(mern.to_str() == "GTCTCTCTGGGTCTTTTTACTTAAGTCAGTG"){
  std::cerr << "the kmer we are searching "<<mer.to_str()<< " original " << mern.to_str() << "\n" ;
  std::cerr << "index of k-mer " << idx << ", rank of kmer " << currRank << "\n";
            std::cerr<<"presence vec: " << presenceVec_[idx] << "\n" ;
            if (presenceVec_[idx] != 1) {std::cerr<<"direction vec: " << directionVec_[idx-currRank]<<"\n" ;}
            std::cerr<<"Canonicalness: " << canonicalNess_[idx] << "\n" ;
            if (presenceVec_[idx] != 1) {
              uint32_t extensionWord = auxInfo_[idx-currRank] ;
              std::bitset<6> ext1(extensionWord) ;
              std::cerr << " extension bits " <<ext1<<"\n" ;
            }
            //std::exit(1);
            //    }


	if(presenceVec_[idx] == 1){
		std::cerr << "currRank = " << currRank << ", size = " << sampledPos_.size() << "\n";
		pos = sampledPos_[currRank];
	}else{
		size_t shiftp{0};
    size_t shiftm{0};
		int inLoop = 0 ;
		do{
			auto extensionPos = idx - currRank ;
		    //std::cerr << "idx - currRank = " <<(idx - currRank) << ", size = " << auxInfo_.size() << "\n";
			uint32_t extensionWord = auxInfo_[extensionPos] ;
            std::bitset<6> ext1(extensionWord) ;
                
            if(mern.to_str() == "GTCTCTCTGGGTCTTTTTACTTAAGTCAGTG")
                    std::cerr << " extension bits " <<ext1<<"\n" ;
                
            if(!canonicalNess_[idx])
                mer.swap();
            /*
            CanonicalKmer mer2_ ;
            mer2_.fromNum(extensionWordKmer);
            std::cerr << " extension bits " <<mer2_.to_str()<<"\n" ;
            */

            auto merStamp = mer ;

			uint32_t mask{0};
			mask = mask | 0x7 ;
			int i = 0;
			while(i < extensionSize - 1){
				mask = mask <<  3 ;
				i++ ;
			}
			i = extensionSize ;
			int movements = 0;

			while(i > 0){
				auto currCode = extensionWord & mask ;
				int j = 0;

				while(j < i-1){
					currCode = currCode >> 3;
					j++ ;
				}
				if(currCode >= 4){
					break ;
				}
				else{
          /*
					if(canonicalNess_[idx]){
						mer.shiftFw((int)(currCode & 0x3)) ;
                        shiftp++;
                    }
					else{
                        //mer.swap();
                        //std::cerr << "I am here"<< "\n";
						mer.shiftFw((int)(currCode & 0x3)) ;
                        shiftm++;
                        }*/

          if(directionVec_[extensionPos] == 1){
            mer.shiftFw((int)(currCode & 0x3)) ;
            shiftp++ ;
          }else{
            mer.shiftBw((int)(currCode & 0x3)) ;
            shiftm++ ;
          }
				}
                
        if(mern.to_str() == "GTCTCTCTGGGTCTTTTTACTTAAGTCAGTG"){
                    std::cerr << "after shift the kmer " << mer.to_str() << "\n" ;
        }
                

				i--;
				mask = mask >> 3 ;
				movements++;
			}

            /*
			if(canonicalNess_[idx]){
				std::cerr << "moving forward " << movements << "\n" ;
			}else{
				std::cerr << "moving backward " << movements << "\n" ;
			}*/
			if(movements > 10)
				std::exit(1) ;

			km = mer.getCanonicalWord() ;
			idx = hash_->lookup(km) ;
			//std::cerr << "idx = " << idx << "\n" ;

      if(mern.to_str() == "GTCTCTCTGGGTCTTTTTACTTAAGTCAGTG"){
        std::cerr << "after shift the kmer " << mer.to_str() << "\n" ;
        std::cerr << "presenceVec " << presenceVec_[idx] << "\n" ;
      }


      if(idx >= numKmers_)
				std::cerr << "this kmer is not found " << mer.to_str() << "\n" ;
			currRank = presenceRank_(idx) ;

			if(inLoop > 5){
				std::cerr<<"I will break from the while loop now" << "\n" ;
				std::cerr<<" mer = "<<mern.to_str() << "\n" ;
				std::exit(1) ;
			}else{
				//std::cerr<<" mer = "<<mer.to_str() << "\n" ;
			}

			inLoop++ ;
		}while(presenceVec_[idx] != 1) ;

		if(presenceVec_[idx]){
			//std::cerr << "idx = " << idx << "\n" ;
			auto sampledPos = sampledPos_[presenceRank_(idx)] ;
      if(shiftp > 0){
        pos = sampledPos - shiftp ;
      }else if(shiftm > 0){

        pos = sampledPos + shiftm;
      }
            //pos = sampledPos - shiftm ;
               // if(mern.to_str() == "CTTAGTCTCTCTGGGTCTTTTTACTTAAGTC")
	             //   std::cerr << " sampled pos " << sampledPos << " pos_ = " << pos << " shiftm "<< shiftm << " shiftp "<< shiftp << "\n" ;
              //  if(pos != posT){
                //    std::cerr << posT << " " << pos<<"\n" ;
               // }
		}
	}
	//end of sampling based pos detection


  //uint64_t pos =
    //(res < numKmers_) ? pos_[res] : std::numeric_limits<uint64_t>::max();
  if (pos <= seq_.size() - k_) {
    uint64_t fk = seq_.get_int(2 * pos, 2 * k_);
    my_mer fkm;
    fkm.word__(0) = fk;
    auto keq = mern.isEquivalent(fkm);
    if (keq != KmerMatchType::NO_MATCH) {
      auto rank = contigRank_(pos);
      auto& pvec = contigTable_[rank];
      // start position of this contig
      uint64_t sp = (rank == 0) ? 0 : static_cast<uint64_t>(contigSelect_(rank)) + 1;
      uint32_t relPos = static_cast<uint32_t>(pos - sp);
      // start position of the next contig - start position of this one
      auto clen = static_cast<uint64_t>(contigSelect_(rank + 1) + 1 - sp);
      bool hitFW = keq == KmerMatchType::IDENTITY_MATCH;
      return {relPos, hitFW, clen, k_, core::range<IterT>{pvec.begin(), pvec.end()}};
    } else {
      //std::cerr << " not found " << mern.to_str() << "\n";
      //std::exit(1) ;
      return {std::numeric_limits<uint32_t>::max(), true, 0, k_, core::range<IterT>{}};
    }
  }

  return {std::numeric_limits<uint32_t>::max(), true, 0, k_, core::range<IterT>{}};
}

uint32_t PufferfishSparseIndex::k() { return k_; }

/**
 * Return the position list (ref_id, pos) corresponding to a contig.
 */
const std::vector<util::Position>& PufferfishSparseIndex::refList(uint64_t contigRank) {
  return contigTable_[contigRank];
}

const std::string& PufferfishSparseIndex::refName(uint64_t refRank) {
  return refNames_[refRank];
}
