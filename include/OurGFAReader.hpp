#ifndef OUR_GFA_READER_H
#define OUR_GFA_READER_H

#include "FatPufferGraph.hpp"
#include "Util.hpp"
#include "cereal/types/string.hpp"
#include "cereal/types/vector.hpp"
#include "sparsepp/spp.h"
#include "string_view.hpp"
#include "zstr/zstr.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <vector>

class PosFinder {
private:
  std::string filename_;
  std::unique_ptr<zstr::ifstream> file;
  size_t k;
  struct Contig {
    std::string seq;
    std::string id;
  };

  spp::sparse_hash_map<uint64_t, util::PackedContigInfo> contigid2seq; 
  spp::sparse_hash_map<uint64_t, std::vector<std::pair<uint64_t, bool>>> path;
  spp::sparse_hash_map<uint64_t, uint32_t> refIDs;
  std::vector<std::string> refMap;
  sdsl::int_vector<2> seqVec_;

  size_t maxTxpKmerCnt_{0};
  size_t totalnTxp_{0};
  util::Positions positions;

  size_t fillContigInfoMap_();
  void encodeSeq(sdsl::int_vector<2>& seqVec, size_t offset,
                 stx::string_view str);

public:
  spp::sparse_hash_map<uint64_t, std::vector<util::Position>> contig2pos;
  PosFinder(const char* gfaFileName, size_t input_k);
  spp::sparse_hash_map<uint64_t, util::PackedContigInfo>& getContigNameMap();

  spp::sparse_hash_map<std::string, std::string>& getContigIDMap();
  std::vector<std::string>& getRefIDs();
  sdsl::int_vector<2>& getContigSeqVec();

  void parseFile();
  void mapContig2Pos();
  void clearContigTable();
  void serializeContigTable(const std::string& odir);
};

#endif
