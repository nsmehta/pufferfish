#include "OurGFAReader.hpp"
#include "CanonicalKmer.hpp"
#include "cereal/archives/binary.hpp"
#include "xxhash.h"
#include <algorithm>
#include <string>
#include "Util.hpp"

PosFinder::PosFinder(const char* gfaFileName, size_t input_k) {
  filename_ = std::string(gfaFileName);
  std::cerr << "Reading GFA file " << gfaFileName << "\n";
  file.reset(new zstr::ifstream(gfaFileName));
  k = input_k;
}

size_t PosFinder::fillContigInfoMap_() {
  std::string ln;
  std::string tag, id, value;
  size_t contig_ctr{0};
  //size_t contig_len{0};
  while (std::getline(*file, ln)) {
    char firstC = ln[0];
    if (firstC != 'S')
      continue;
    stx::string_view lnview(ln);
    std::vector<stx::string_view> splited = util::split(lnview, '\t');
    try {
      auto nid = std::stoull(splited[1].to_string());
      (void)nid;
      auto clen = splited[2].length();
      contigid2seq[nid] = {contig_ctr, 0, static_cast<uint32_t>(clen), 0};
      ++contig_ctr;
      // contigid2seq[id] = value;
      // contig_cnt++;
    } catch (std::exception& e) {
      // not numeric contig
    }
  }
  size_t total_len{0};
  for (auto& kv : contigid2seq) {
    kv.second.offset = total_len;
    total_len += kv.second.length;
  }
  return total_len;
}

void PosFinder::encodeSeq(sdsl::int_vector<2>& seqVec, size_t offset,
                          stx::string_view str) {
  for (size_t i = 0; i < str.length(); ++i) {
    auto c = my_mer::code(str[i]);
    seqVec[offset + i] = c;
  }
}

sdsl::int_vector<2>& PosFinder::getContigSeqVec() { return seqVec_; }

void PosFinder::parseFile() {
  size_t total_len = fillContigInfoMap_();
  file.reset(new zstr::ifstream(filename_));
  std::cerr << "total contig length = " << total_len << "\n";
  std::cerr << "packing contigs into contig vector\n";
  seqVec_ = sdsl::int_vector<2>(total_len, 0);

  std::string ln;
  std::string tag, id, value;
  size_t contig_cnt{0};
  size_t ref_cnt{0};
  while (std::getline(*file, ln)) {
    char firstC = ln[0];
    if (firstC != 'S' and firstC != 'P')
      continue;
    stx::string_view lnview(ln);
    std::vector<stx::string_view> splited = util::split(lnview, '\t');
    tag = splited[0].to_string();
    id = splited[1].to_string();
    // value = splited[2].to_string();
    // A segment line
    if (tag == "S") {
      try {
        uint64_t nid = std::stoll(id);
        encodeSeq(seqVec_, contigid2seq[nid].offset, splited[2]);
        // contigid2seq[nid] = value;
      } catch (std::exception& e) {
        // not a numeric contig id
      }
      contig_cnt++;
    }
    // A path line
    if (tag == "P") {
      auto pvalue = splited[2];
      std::vector<std::pair<uint64_t, bool>> contigVec = util::explode(pvalue, ',');
      // parse value and add all conitgs to contigVec

      path[ref_cnt] = contigVec;
      refMap.push_back(id);
      ref_cnt++;

	  maxTxpKmerCnt_ = std::max(contigVec.size()-k+1, maxTxpKmerCnt_);
	  for (auto& c : contigVec) {
		  contigid2seq[c.first].txpCnt++;	
		  totalnTxp_++;
	  }
    }
  }
  std::cerr << " Total # of Contigs : " << contig_cnt
            << " Total # of numerical Contigs : " << contigid2seq.size()
            << "\n\n";
}

// spp::sparse_hash_map<uint64_t, std::string>& PosFinder::getContigNameMap() {
spp::sparse_hash_map<uint64_t, util::PackedContigInfo>&
PosFinder::getContigNameMap() {
  return contigid2seq;
}

std::vector<std::string>& PosFinder::getRefIDs() { return refMap; }

void PosFinder::mapContig2Pos() {
  uint64_t pos = 0;
  uint64_t accumPos;
  uint64_t currContigLength = 0;
  uint64_t total_output_lines = 0;
  
  size_t contigOffsetlen = log(contigid2seq[contigid2seq.size()-1].offset);
  size_t txpIdWlen = totalnTxp_*log(path.size());
  size_t txpPosWlen = totalnTxp_*log(maxTxpKmerCnt_);

  contigOffset = sdsl::int_vector<>(contigid2seq.size(), 0, contigOffsetlen);
  txpID = sdsl::int_vector<>(totalnTxp_, 0, txpIdWlen);
  txpPos = sdsl::int_vector<>(totalnTxp_, 0, txpPosWlen);
  contigOri = sdsl::bit_vector(totalnTxp_, 0);

  uint64_t offset{0};
  for (size_t i = 0; i < contigid2seq.size(); i++) {
		  contigOffset[i] = offset;
		  offset += contigid2seq[i].txpCnt;
  }

  for (auto const& ent : path) {
    const uint64_t& tr = ent.first;
    const std::vector<std::pair<uint64_t, bool>>& contigs = ent.second;
    accumPos = 0;
    for (size_t i = 0; i < contigs.size(); i++) {
      if (contig2pos.find(contigs[i].first) == contig2pos.end()) {
        contig2pos[contigs[i].first] = {};
        total_output_lines += 1;
      }
      if (contigid2seq.find(contigs[i].first) == contigid2seq.end()) {
        std::cerr << "ERROR: Couldn't find the contig in the path : " << contigs[i].first << "\n";
      }
      pos = accumPos;
      currContigLength = contigid2seq[contigs[i].first].length;
      accumPos += currContigLength - k;
//      (contig2pos[contigs[i].first])
//          .push_back(util::Position(tr, pos, contigs[i].second));
	  size_t idx = contigOffset[contigs[i].first]+contigid2seq[i].txpCnt-1;
	  txpID[idx] = tr;
	  txpPos[idx] = pos;
	  contigOri[idx] = 1;
	  contigid2seq[i].txpCnt--;
    }
  }
  std::cerr << "\nTotal # of segments we have position for : "
            << total_output_lines << "\n";
}

void PosFinder::clearContigTable() {
  refMap.clear();
  contig2pos.clear();
}

// Note : We assume that odir is the name of a valid (i.e., existing) directory.
void PosFinder::serializeContigTable(const std::string& odir) {
  std::string ofile = odir + "/ctable.bin";
  std::string eqfile = odir + "/eqtable.bin";
  std::ofstream ct(ofile);
  std::ofstream et(eqfile);
  cereal::BinaryOutputArchive ar(ct);
  cereal::BinaryOutputArchive eqAr(et);
  {
    // We want to iterate over the contigs in precisely the
    // order they appear in the contig array (i.e., the iterator
    // order of contigid2seq).
    std::vector<std::string> refNames;
    refNames.reserve(refMap.size());
    for (size_t i = 0; i < refMap.size(); ++i) {
      refNames.push_back(refMap[i]);
    }
    ar(refNames);

    class VecHasher {
    public:
      size_t operator()(const std::vector<uint32_t>& vec) const {
        return XXH64(const_cast<std::vector<uint32_t>&>(vec).data(),
                     vec.size() * sizeof(decltype(vec.front())), 0);
      }
    };

    spp::sparse_hash_map<std::vector<uint32_t>, uint32_t, VecHasher> eqMap;
    std::vector<uint32_t> eqIDs;
    std::vector<std::vector<util::Position>> cpos;

    for (auto& kv : contigid2seq) {
      cpos.push_back(contig2pos[kv.first]);
      std::vector<uint32_t> tlist;
      for (auto& p : contig2pos[kv.first]) {
        tlist.push_back(p.transcript_id());
      }
      std::sort(tlist.begin(), tlist.end());
      tlist.erase(std::unique(tlist.begin(), tlist.end()), tlist.end());
      size_t eqID = eqMap.size();
      if (eqMap.contains(tlist)) {
        eqID = eqMap[tlist];
      } else {
        eqMap[tlist] = eqID;
      }
      eqIDs.push_back(eqID);
    }
    std::cerr << "there were " << eqMap.size() << " equivalence classes\n";
    eqAr(eqIDs);
    eqIDs.clear();
    eqIDs.shrink_to_fit();
    std::vector<std::vector<uint32_t>> eqLabels;
    eqLabels.reserve(eqMap.size());
    for (auto& kv : eqMap) {
      eqLabels.push_back(kv.first);
    }
    std::sort(eqLabels.begin(), eqLabels.end(),
              [&](const std::vector<uint32_t>& l1,
                  const std::vector<uint32_t>& l2) -> bool {
                return eqMap[l1] < eqMap[l2];
              });
    eqAr(eqLabels);
    ar(cpos);
  }
}


