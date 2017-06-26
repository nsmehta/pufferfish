#include "GFAConverter.hpp"
#include "cereal/archives/binary.hpp"
#include "xxhash.h"
#include <algorithm>

GFAConverter::GFAConverter(const char* gfaFileName, size_t input_k) {
		filename_ = gfaFileName;
  std::cerr << "Reading GFA file " << gfaFileName << "\n";
  file.reset(new zstr::ifstream(gfaFileName));
  k = input_k;
}

void GFAConverter::parseFile() {
  std::string ln;
  std::string tag, id, value;
  size_t contig_cnt{0};
  {	
    std::cerr << "Start reading GFA file... \n";
    spp::sparse_hash_map<std::string, uint64_t> seq2newid;
    uint64_t idCntr = 0;
    while (std::getline(*file, ln)) {
      char firstC = ln[0];
      if (firstC != 'S')
        continue;
      stx::string_view lnview(ln);
      std::vector<stx::string_view> splited = util::split(lnview, '\t');
      id = splited[1].to_string();
      value = splited[2].to_string();
        if (util::is_number(id)) {
		  uint64_t contigId = std::stoll(id);
          processContigSeq(contigId, value, seq2newid, idCntr);
        }
        contig_cnt++;
	}

	file.reset(new zstr::ifstream(filename_));
	while (std::getline(*file, ln)) {
 	   char firstC = ln[0];
		   if (firstC != 'P')
				continue;
      stx::string_view lnview(ln);
      std::vector<stx::string_view> splited = util::split(lnview, '\t');
      id = splited[1].to_string();
      //value = splited[2].to_string();
      // A path line
        auto pvalue = splited[2];
        std::vector<std::pair<uint64_t, bool>> contigVec =
            util::explode(pvalue, ',');
        // parse value and add all conitgs to contigVec
        path[id] = contigVec;
	auto& oldId = contigVec[0];	
    auto& newIdList = old2newids[oldId.first];//s.first];
	if (oldId.second) { //Old segment appears in forward orientation
		if (newIdList[0].second) // if first segment matches the old one in forward orientation, it shouldn't be updated from its start point
				newSegs[newIdList[0].first].set_start();
		else // ow, it shouldn't be updated from its end point
				newSegs[newIdList[0].first].set_end();
	}
	else { //Old segment appears in reverse orientation
		auto lastIdx = newIdList.size()-1;
		if (newIdList[lastIdx].second) newSegs[newIdList[lastIdx].first].set_end();
		else newSegs[newIdList[lastIdx].first].set_start();
	}

	oldId = contigVec[contigVec.size()-1];	
    newIdList = old2newids[oldId.first];//s.first];
	if (oldId.second) { //Old segment appears in forward orientation
		auto lastIdx = newIdList.size()-1;
		if (newIdList[lastIdx].second) // if first segment matches the old one in forward orientation, it shouldn't be updated from its start point
				newSegs[newIdList[lastIdx].first].set_end();
		else // ow, it shouldn't be updated from its end point
				newSegs[newIdList[lastIdx].first].set_start();
	}
	else { //Old segment appears in reverse orientation
		if (newIdList[0].second) newSegs[newIdList[0].first].set_start();
		else newSegs[newIdList[0].first].set_end();
	}
    }
    std::cerr << "Done Reading\nStart updating pathStart and pathEnd...\n";
  }

  std::cerr << "Done updating pathStart and pathEnd based on the newIds\n";
  std::cerr << "Total # of Contigs : " << contig_cnt
            << "\tTotal # of numerical Contigs : " << old2newids.size()
            << " \tTotal # of contigs after spliting :" << newSegs.size() //new2seqAoldids.size()
            << "\n";
}

void GFAConverter::buildGraph() {
  std::cerr << "Start building the graph...\n";
  // build the graph
  for (auto& kv : path) {
	//std::string pathId = kv.first;
    auto& contigVec = kv.second;

    std::pair<uint64_t, bool> prev;
    bool pathStart = true;
    for (auto& idOri : contigVec) {
      std::vector<std::pair<uint64_t, bool>> newIdList =
          old2newids[idOri.first];
      if (!idOri.second) {
        std::reverse(newIdList.begin(), newIdList.end()); // As newIdList is a
                                                          // copy it won't
                                                          // affect the actual
                                                          // vector in
                                                          // old2newids (tested)
      }
      if (pathStart) {
        uint64_t id = newIdList[0].first;
        bool ori = (idOri.second ? newIdList[0].second : !newIdList[0].second);
        prev = std::make_pair(id, ori);
        pathStart = false;
      }
      for (size_t i = 1; i < newIdList.size(); i++) {
        uint64_t id = newIdList[i].first;
        bool ori = (idOri.second ? newIdList[i].second : !newIdList[i].second);
        semiCG.addEdge(prev.first, prev.second, id, ori);
        prev = std::make_pair(id, ori);
      }
    }
  }
  std::cerr << "Done Building the graph\n";
}

void GFAConverter::processContigSeq(
    uint64_t& contigId, std::string& contigSeq,
    spp::sparse_hash_map<std::string, uint64_t>& seq2newid,
    uint64_t& idCntr) {
  /**
   * Divide every segment into 3 pieces
   * first k
   * middle one without the first and last nucleotide
   * last k
   * special case is segment with length = k+1 that doesn't have the middle case
   * ATTENTION : The order we insert these 3 pieces into the vector matters. So
   *keep it the way it is
  **/
  // std::string prefix = "00";
  std::vector<std::string> seqParts;
  if (util::isRevcomp(contigSeq)) {
    for (size_t i = 0; i <= contigSeq.size() - k; i++) {
      seqParts.push_back(contigSeq.substr(i, k));
    }
  } else {
    seqParts.push_back(contigSeq.substr(0, k));
    if (contigSeq.size() > k + 1) {
      seqParts.push_back(contigSeq.substr(1, contigSeq.size() - 2));
    }
    seqParts.push_back(contigSeq.substr(contigSeq.size() - k));
  }

  // for each part, search in hash_map whether it already existed or not and
  // assign the proper (new/retrieved) contig id to it
  // int cntr = 1;
  for (std::string seq : seqParts) {
    uint64_t newId;
    bool plus = true;
    if (seq2newid.find(seq) != seq2newid.end()) {
      newId = seq2newid[seq];
    } else if (seq2newid.find(util::revcomp(seq)) != seq2newid.end()) {
      plus = false;
      newId = seq2newid[util::revcomp(seq)];
    } else {
      // newId = prefix + std::to_string(cntr) + contigId;
      newId = idCntr;
      idCntr++;
      seq2newid[seq] = newId;
    }

	if (newId < newSegs.size()) { // newSeg already created. Just add the new oldid
      newSegs[newId].add_oldId(contigId, plus);
	} else {
		SegInfo si(seq);
		si.add_oldId(contigId, plus);
		newSegs.push_back(si);
      	if (seq.size() == k)
        	ksizeContig.push_back(newId);
	}
/*    if (new2seqAoldids.find(newId) != new2seqAoldids.end()) {
      new2seqAoldids[newId].second.emplace_back(contigId, plus);
    } else {
      std::vector<std::pair<uint64_t, bool>> ids;
      ids.emplace_back(contigId, plus);
      new2seqAoldids[newId] = std::make_pair(seq, ids);
      if (seq.size() == k)
        ksizeContig.push_back(newId);
    }*/
    old2newids[contigId].emplace_back(newId, plus);
  }
}

void GFAConverter::randomWalk() {
  spp::sparse_hash_map<uint64_t, pufg::Node>& nodes = semiCG.getVertices();
  std::cerr << "# of contigs with length = 31 : " << ksizeContig.size() << "\n";
  std::cerr << "\nStart merging .. \n";
  for (auto& v : ksizeContig) {
    uint64_t curId = v;
    pufg::Node& curNode = nodes[curId];
    // I strongly rely on TwoPaCo here for not having a case of possible merging
    // for in & out nodes both while none of in/out nodes are of size k!!
    if (curNode.getRealIndeg() == 1 and !isCornerCase(curNode, true)) {
      mergeIn(curNode);
    } else if (curNode.getRealOutdeg() == 1 and !isCornerCase(curNode, false)) {
      mergeOut(curNode);
    }
    // otherwise it is complex and you should keep the node and not merge it
    // with any left or right neighbors
  }
  std::cerr << "Done merging \n";
}

void GFAConverter::eraseFromOldList(uint64_t nodeId) {
//  if (new2seqAoldids.contains(nodeId)) {
  if (newSegs[nodeId].is_valid()) {
  //  auto& seqAoldids = new2seqAoldids[nodeId];
    std::vector<std::pair<uint64_t, bool>>& oldids = newSegs[nodeId].get_oldIds();// seqAoldids.second;
    for (auto& idOri : oldids) {
      uint64_t& id = idOri.first;
      auto& newids = old2newids[id];
      newids.erase(std::remove_if(
                       newids.begin(), newids.end(),
                       [&nodeId](std::pair<uint64_t, bool>& newid) -> bool {
                         return newid.first == nodeId;
                       }),
                   newids.end());
    }
  }
}
void GFAConverter::mergeIn(pufg::Node& n) {
  uint64_t id = n.getId();
  pufg::edgetuple& edge = n.getOnlyRealIn();
	if (!newSegs[id].is_valid()) {
    std::cerr << "[mergeIn] NO; the id " << id
              << " was not in new2seqAoldids!\n";
  }
	if (!newSegs[edge.contigId].is_valid()) {
    std::cerr << "[mergeIn] NO; the edge.contigId " << edge.contigId
              << " was not in new2seqAoldids!\n";
  }
  std::string& tobeMerged = newSegs[id].get_seq();// new2seqAoldids[id].first;
  std::string& seq = newSegs[edge.contigId].get_seq();// new2seqAoldids[edge.contigId].first;
  if (edge.baseSign() != edge.neighborSign()) {
    tobeMerged = util::revcomp(tobeMerged);
    //				if (tobeMerged.substr(tobeMerged.size()-(k-1)) != seq.substr(0,
    //k-1)) std::cerr << "1 " << id << " " << edge.contigId << " " << seq <<
    //"\n" << tobeMerged << "\n";
    seq = tobeMerged.substr(0, tobeMerged.size() - (k - 1)) + seq;
  } else {
    //			if (tobeMerged.substr(0, k-1) != seq.substr(seq.size() - (k-1)))
    //std::cerr << "2 " << seq << "\n" << tobeMerged << "\n";
    seq += tobeMerged.substr(k - 1);
  }
  eraseFromOldList(id);
  newSegs[id].disregard();
  semiCG.removeNode(id);
  if (newSegs[id].is_start())
    update_start(edge.contigId, edge.baseSign() == edge.neighborSign());
  if (newSegs[id].is_end())
    update_end(edge.contigId, edge.baseSign() == edge.neighborSign());
}

void GFAConverter::mergeOut(pufg::Node& n) {
  uint64_t id = n.getId();
  pufg::edgetuple& edge = n.getOnlyRealOut();
 //  if (!new2seqAoldids.contains(id)) {
	if (!newSegs[id].is_valid()) {
    std::cerr << "[mergeOut] NO; the id " << id
              << " was not in new2seqAoldids!\n";
  }
//  if (!new2seqAoldids.contains(edge.contigId)) {
	if (!newSegs[edge.contigId].is_valid()) {
    std::cerr << "[mergeOut] NO; the edge.contigId " << edge.contigId
              << " was not in new2seqAoldids!\n";
  } 
  std::string& tobeMerged =  newSegs[id].get_seq();//new2seqAoldids[id].first;
  std::string& seq = newSegs[edge.contigId].get_seq();// new2seqAoldids[edge.contigId].first;
  if (edge.baseSign() != edge.neighborSign()) {
    tobeMerged = util::revcomp(tobeMerged);
    //				if (tobeMerged.substr(0, k-1) != seq.substr(seq.size() - (k-1)))
    //std::cerr << id << " " << edge.contigId << " " << "3 " << seq << "\n" <<
    //tobeMerged << "\n";
    seq += tobeMerged.substr(k - 1);
  } else {
    //			if (tobeMerged.substr(tobeMerged.size()-(k-1)) != seq.substr(0,
    //k-1)) std::cerr << "4 " << seq << "\n" << tobeMerged << "\n";
    seq = tobeMerged.substr(0, tobeMerged.size() - (k - 1)) + seq;
  }
  //   		if (edge.contigId == "00125208939" or edge.contigId == "00225208939" or
  //   edge.contigId == "00325208939")
  //				std::cerr << id << " is merged out with " << edge.contigId
  //<< " seq:" << seq << "\n";
  //		if (id == "00125208939" or id == "00225208939" or id ==
  //"00325208939")
  //				std::cerr << id << " is merged out with " << edge.contigId
  //<< " seq:" << seq << "\n";
  eraseFromOldList(id);
  newSegs[id].disregard();
  //new2seqAoldids.erase(id);
  semiCG.removeNode(id);
//  if (is_start(id))
  if (newSegs[id].is_start())
    update_start(edge.contigId, edge.baseSign() == edge.neighborSign());
//  if (is_end(id))
  if (newSegs[id].is_end())
    update_end(edge.contigId, edge.baseSign() == edge.neighborSign());
  //		return edge.contigId;
}

/**
 * Check 3 cases
 * 1. node is not corner node
 * 2. the neighbor we want to merge the node with is not corner node
 * 3. the current node is not the same as the node it will be merged with
 */
bool GFAConverter::isCornerCase(pufg::Node& n, bool mergeIn) {
  if (mergeIn) {
		  if (newSegs[n.getId()].is_start()) return true;
/*    if (is_start(n.getId()))
      return true;*/
    pufg::edgetuple& edge = n.getOnlyRealIn();
    if (edge.contigId == n.getId())
      return true;
    pufg::Node& neighbor = semiCG.getVertices()[edge.contigId];
    // if ( (edge.baseSign and edge.neighborSign) or (!edge.baseSign and !
    // edge.neighborSign) )
    if (edge.baseSign() == edge.neighborSign()) {
//      return is_end(edge.contigId) or neighbor.getRealOutdeg() != 1;
      return newSegs[edge.contigId].is_end() or neighbor.getRealOutdeg() != 1;
    } else
//      return is_start(edge.contigId) or neighbor.getRealIndeg() != 1;
      return newSegs[edge.contigId].is_start() or neighbor.getRealIndeg() != 1;
  } else { // Merge out case
		if (newSegs[n.getId()].is_end()) return true;
/*    if (is_end(n.getId()))
      return true;*/
    pufg::edgetuple& edge = n.getOnlyRealOut();
    if (edge.contigId == n.getId())
      return true;
    pufg::Node& neighbor = semiCG.getVertices()[edge.contigId];
    // if ( (edge.baseSign and edge.neighborSign) or (!edge.baseSign and !
    // edge.neighborSign) ) {
    if (edge.baseSign() == edge.neighborSign()) {
//      return is_start(edge.contigId) or neighbor.getRealIndeg() != 1;
      return newSegs[edge.contigId].is_start() or neighbor.getRealIndeg() != 1;

    } else {
//      return is_end(edge.contigId) or neighbor.getRealOutdeg() != 1;
      return newSegs[edge.contigId].is_end() or neighbor.getRealOutdeg() != 1;
    }
  }
  return false;
}

/*bool GFAConverter::is_start(uint64_t& nodeId) {
  if (pathStart.find(std::make_pair(nodeId, true)) != pathStart.end() or
      pathEnd.find(std::make_pair(nodeId, false)) != pathEnd.end())
    return true;
  return false;
}

bool GFAConverter::is_end(uint64_t& nodeId) {
  if (pathStart.find(std::make_pair(nodeId, false)) != pathStart.end() or
      pathEnd.find(std::make_pair(nodeId, true)) != pathEnd.end())
    return true;
  return false;
}*/

void GFAConverter::update_start(uint64_t& newId, bool newOri) {
		if (newOri) newSegs[newId].set_start();
		else newSegs[newId].set_end();
//  pathStart[std::make_pair(newId, newOri)] = true;
//  pathEnd[std::make_pair(newId, !newOri)] = true;
}

void GFAConverter::update_end(uint64_t& newId, bool newOri) {
		if (newOri) newSegs[newId].set_end();
		else newSegs[newId].set_start();
//  pathStart[std::make_pair(newId, !newOri)] = true;
//  pathEnd[std::make_pair(newId, newOri)] = true;
}

void GFAConverter::writeFile(const char* gfaFileName) {

  std::ofstream gfa_file(gfaFileName);
  uint32_t contigCntr = 0;
  /*for (auto& kv : new2seqAoldids) {
    gfa_file << "S"
             << "\t" << kv.first << "\t" << (kv.second).first << "\n";
    contigCntr++;
  }*/
  for (size_t i = 0; i < newSegs.size(); i++) {
  	if (newSegs[i].is_valid()) {
		gfa_file << "S"
				<< "\t" << i << "\t" << newSegs[i].get_seq() << "\n";
		contigCntr++;
	}
  }
  for (auto& p : path) {
    auto tid = p.first;
    gfa_file << "P"
             << "\t" << tid << "\t";
    auto vec = p.second;
    bool first = true;
    std::pair<uint64_t, bool> prev;
    bool firstPiece = true;
    for (size_t i = 0; i < vec.size(); i++) {
      auto newidVec = old2newids[vec[i].first];
      if (!vec[i].second)
        std::reverse(newidVec.begin(), newidVec.end());
      for (auto& idOri : newidVec) {
        uint64_t id = idOri.first;
        bool ori = vec[i].second ? idOri.second : !idOri.second;
        if (!first) {
          if (firstPiece and prev.first == id and prev.second == ori) {
            firstPiece = false;
            // NOTE: Why do we want this conditional?
            // if (new2seqAoldids[id].first.size() == k)
            continue;
          }
          /*						if (firstPiece and prev.first == id and
             prev.second != ori) {
                                      std::cerr << "found the bug " <<
             prev.first << "\n";
                                  }*/
          gfa_file << ",";
        }
        gfa_file << id << (ori ? "+" : "-");
        prev = std::make_pair(id, ori);
        first = false;
      }
      if (newidVec.size() != 0)
        firstPiece = true;
      // here's the crazy case.
      // If an old contig's list of new ids is completely empty, this means that
      // we are actually jumping over this contig,
      // so we shouldn't be worried about the next new contig with previous one,
      // as prev is not actually the real prev.
      // Say we have the case A- B+ A- in old set of contigs
      // A : k-, k'+, k''+
      // B : k-, k''+
      // Now if we merge all 3 pieces of A into one final contig F, all parts of
      // B are removed.
      // Then if we see F+, F+ we shouldn't remove the second F+ b/c F+ is a
      // newly created contig that has overlap of k-1 with itself.
      // The path is actually F+ followed by F+ and not the case that these two
      // contigs are equal kmers with overlap of k.
      // Adding one "else" and a whole story! It's not even counted as a comment
      // anymore!!!!!! :/
      else
        firstPiece = false;
    }
    gfa_file << "\t*\n";
  }
  std::cerr << "# of contigs written to file : " << contigCntr << "\n";
}
