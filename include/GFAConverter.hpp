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

enum class BoundaryType: uint8_t {
	STARTER = 0,
	TERMINATOR = 1,
	BOTH = 2,
	NONE = 3
};


struct SegInfo {
	SegInfo(std::string seqIn) : seq(seqIn) {
		boundaryType = BoundaryType::NONE; //is neither start nor end
	}
	void add_oldId(uint64_t oldId, bool ori) {
		oldIds.emplace_back(oldId, ori);
	}
	void set_start() {
			if (boundaryType == BoundaryType::NONE) 
				boundaryType = BoundaryType::STARTER;
			else if (boundaryType == BoundaryType::TERMINATOR)
				boundaryType = BoundaryType::BOTH;			
	}
	void set_end() {
			if (boundaryType == BoundaryType::NONE)
					boundaryType = BoundaryType::TERMINATOR;
			else if (boundaryType == BoundaryType::STARTER)
					boundaryType = BoundaryType::BOTH;
	}
	bool is_start() {return boundaryType == BoundaryType::STARTER or boundaryType == BoundaryType::BOTH;}
	bool is_end() {return boundaryType == BoundaryType::TERMINATOR or boundaryType == BoundaryType::BOTH;}

	void disregard() {
		seq = "";
		oldIds.clear();
	}
	bool is_valid() {return seq != "";}
	std::vector<std::pair<uint64_t, bool> >& get_oldIds() {return oldIds;}
	std::string& get_seq() {return seq;}
	std::string seq;
	BoundaryType boundaryType;
	std::vector<std::pair<uint64_t, bool> > oldIds;
};

class GFAConverter {
private:
  std::string filename_;
  std::unique_ptr<zstr::ifstream> file;
  size_t k;
  struct Contig {
    std::string seq;
    std::string id;
  };

  std::vector<SegInfo> newSegs;
/*  spp::sparse_hash_map<
      uint64_t,
      std::pair<std::string, std::vector<std::pair<uint64_t, bool>>>>
      new2seqAoldids;*/
  std::vector<uint64_t> ksizeContig;
  spp::sparse_hash_map<uint64_t, std::vector<std::pair<uint64_t, bool>>>
      old2newids;
  // path maps each transcript_id to a pair of <contig_id, orientation>
  // orientation : +/true main, -/false reverse
  spp::sparse_hash_map<std::string, std::vector<std::pair<uint64_t, bool>>>
      path;

//  std::map<std::pair<uint64_t, bool>, bool, util::cmpByPair> pathStart;
//  std::map<std::pair<uint64_t, bool>, bool, util::cmpByPair> pathEnd;

  pufg::Graph semiCG;

  void
  processContigSeq(uint64_t& contigId, std::string& contigSeq,
                   spp::sparse_hash_map<std::string, uint64_t>& seq2newid,
                   uint64_t& idCntr);
  bool is_start(uint64_t& nodeId);
  bool is_end(uint64_t& nodeId);
  void update_start(uint64_t& newId, bool newOri);
  void update_end(uint64_t& newId, bool newOri);
  // bool isCornerCase(pufg::Node& n, pufg::edgetuple, bool mergeIn);
  bool isCornerCase(pufg::Node& n, bool mergeIn);
  void mergeIn(pufg::Node& n);
  void mergeOut(pufg::Node& n);
  void eraseFromOldList(uint64_t nodeId);

public:
  GFAConverter(const char* gfaFileName, size_t input_k);
  void parseFile();
  void buildGraph();
  void randomWalk();
  void writeFile(const char* gfaFileName);
};

#endif
