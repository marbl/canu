
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

// Copyright (c) 2011-2015, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following
// disclaimer in the documentation and/or other materials provided
// with the distribution.
//
// * Neither the name of Pacific Biosciences nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

#include "Alignment.H"
#include "AlnGraphBoost.H"

static int MAX_OFFSET = 10000;

AlnGraphBoost::AlnGraphBoost(const std::string& backbone) {
    // initialize the graph structure with the backbone length + enter/exit
    // vertex
    size_t blen = backbone.length();
    _templateLength = blen;
    _g = G(blen+1);
    for (size_t i = 0; i < blen+1; i++)
        boost::add_edge(i, i+1, _g);

    VtxIter curr, last;
    boost::tie(curr, last) = boost::vertices(_g);
    _enterVtx = *curr++;
    _g[_enterVtx].base = '^';
    _g[_enterVtx].tpos = UINT32_MAX;
    _g[_enterVtx].backbone = true;
    for (size_t i = 0; i < blen; i++, ++curr) {
        VtxDesc v = *curr;
        _g[v].backbone = true;
        _g[v].weight = 0;
        _g[v].base = backbone[i];
        _g[v].tpos = UINT32_MAX;
        _bbMap[v] = v;
    }
    _exitVtx = *curr;
    _g[_exitVtx].base = '$';
    _g[_exitVtx].tpos = UINT32_MAX;
    _g[_exitVtx].backbone = true;
}

AlnGraphBoost::AlnGraphBoost(const size_t blen) {
    _templateLength = blen;
    _g = G(blen+1);
    for (size_t i = 0; i < blen+1; i++)
        boost::add_edge(i, i+1, _g);

    VtxIter curr, last;
    boost::tie(curr, last) = boost::vertices(_g);
    _enterVtx = *curr++;
    _g[_enterVtx].base = '^';
    _g[_enterVtx].tpos = UINT32_MAX;
    _g[_enterVtx].backbone = true;
    for (size_t i = 0; i < blen; i++, ++curr) {
        VtxDesc v = *curr;
        _g[v].backbone = true;
        _g[v].weight = 0;
        _g[v].deleted = false;
        _g[v].base = 'N';
        _g[v].tpos = UINT32_MAX;
        _bbMap[v] = v;
    }
    _exitVtx = *curr;
    _g[_exitVtx].base = '$';
    _g[_exitVtx].tpos = UINT32_MAX;
    _g[_exitVtx].backbone = true;
}

void AlnGraphBoost::addAln(dagAlignment& aln) {
    IndexMap index = boost::get(boost::vertex_index, _g);
    // tracks the position on the backbone
    uint32_t bbPos = aln.start;
    VtxDesc prevVtx = _enterVtx;
    for (size_t i = 0; i < aln.length; i++) {
        char queryBase = aln.qstr[i], targetBase = aln.tstr[i];
        VtxDesc currVtx = index[bbPos];
        _g[_bbMap[currVtx]].tpos = bbPos;
        // match
        if (queryBase == targetBase) {
            if (_g[_bbMap[currVtx]].coverage < UINT8_MAX) _g[_bbMap[currVtx]].coverage++;

            // NOTE: for empty backbones
            _g[_bbMap[currVtx]].base = targetBase;

            if (_g[currVtx].weight < UINT8_MAX) _g[currVtx].weight++;
            if (prevVtx != _enterVtx || bbPos <= MAX_OFFSET || MAX_OFFSET == 0)
                addEdge(prevVtx, currVtx);
            else
                addEdge(_bbMap[bbPos-1], currVtx);
            bbPos++;
            prevVtx = currVtx;
        // query deletion
        } else if (queryBase == '-' && targetBase != '-') {
            if (_g[_bbMap[currVtx]].coverage < UINT8_MAX) _g[_bbMap[currVtx]].coverage++;

            // NOTE: for empty backbones
            _g[_bbMap[currVtx]].base = targetBase;

            bbPos++;
        // query insertion
        } else if (queryBase != '-' && targetBase == '-') {
            // create new node and edge
            VtxDesc newVtx = boost::add_vertex(_g);
            _g[newVtx].base = queryBase;
            if (_g[newVtx].weight < UINT8_MAX) _g[newVtx].weight++;
            _g[newVtx].backbone = false;
            _g[newVtx].deleted = false;
            _bbMap[newVtx] = bbPos;

            if (prevVtx != _enterVtx || bbPos <= MAX_OFFSET || MAX_OFFSET == 0)
               addEdge(prevVtx, newVtx);
            else
               addEdge(_bbMap[bbPos-1], newVtx);
            prevVtx = newVtx;
        }
    }
    if (bbPos + MAX_OFFSET >= _templateLength || MAX_OFFSET == 0)
       addEdge(prevVtx, _exitVtx);
    else
       addEdge(prevVtx, _bbMap[bbPos]);
}

void AlnGraphBoost::addEdge(VtxDesc u, VtxDesc v) {
    // Check if edge exists with prev node.  If it does, increment edge counter,
    // otherwise add a new edge.
    InEdgeIter ii, ie;
    bool edgeExists = false;
    for (boost::tie(ii, ie) = boost::in_edges(v, _g); ii != ie; ++ii) {
        EdgeDesc e = *ii;
        if (boost::source(e , _g) == u) {
            // increment edge count
            if (_g[e].count < UINT8_MAX) _g[e].count++;
            edgeExists = true;
        }
    }
    if (! edgeExists) {
        // add new edge
        std::pair<EdgeDesc, bool> p = boost::add_edge(u, v, _g);
        if (_g[p.first].count < UINT8_MAX) _g[p.first].count++;
    }
}

void AlnGraphBoost::mergeNodes() {
    std::queue<VtxDesc> seedNodes;
    seedNodes.push(_enterVtx);

    while(true) {
        if (seedNodes.size() == 0)
            break;

        VtxDesc u = seedNodes.front();
        seedNodes.pop();
        mergeInNodes(u);
        mergeOutNodes(u);

        OutEdgeIter oi, oe;
        for (boost::tie(oi, oe) = boost::out_edges(u, _g); oi != oe; ++oi) {
            EdgeDesc e = *oi;
            _g[e].visited = true;
            VtxDesc v = boost::target(e, _g);
            InEdgeIter ii, ie;
            int notVisited = 0;
            for (boost::tie(ii, ie) = boost::in_edges(v, _g); ii != ie; ++ii) {
                if (_g[*ii].visited == false)
                    notVisited++;
            }

            // move onto the boost::target node after we visit all incoming edges for
            // the boost::target node
            if (notVisited == 0)
                seedNodes.push(v);
        }
    }
}

void AlnGraphBoost::mergeInNodes(VtxDesc n) {
    std::map<char, std::vector<VtxDesc> > nodeGroups;
    InEdgeIter ii, ie;
    // Group neighboring nodes by base
    for(boost::tie(ii, ie) = boost::in_edges(n, _g); ii != ie; ++ii) {
        VtxDesc inNode = boost::source(*ii, _g);
        if (out_degree(inNode, _g) == 1) {
            nodeGroups[_g[inNode].base].push_back(inNode);
        }
    }

    // iterate over node groups, merge an accumulate information
    for(std::map<char, std::vector<VtxDesc> >::iterator kvp = nodeGroups.begin(); kvp != nodeGroups.end(); ++kvp) {
        std::vector<VtxDesc> nodes = (*kvp).second;
        if (nodes.size() <= 1)
            continue;

        std::vector<VtxDesc>::const_iterator ni = nodes.begin();
        VtxDesc an = *ni++;
        OutEdgeIter anoi, anoe;
        boost::tie(anoi, anoe) = boost::out_edges(an, _g);

        // Accumulate out edge information
        for (; ni != nodes.end(); ++ni) {
            OutEdgeIter oi, oe;
            boost::tie(oi, oe) = boost::out_edges(*ni, _g);
            _g[*anoi].count =  ((uint32_t)_g[*anoi].count + (uint32_t)_g[*oi].count < UINT8_MAX)  ? _g[*anoi].count + _g[*oi].count :  UINT8_MAX;
            _g[an].weight   =  ((uint32_t)_g[an].weight + (uint32_t)_g[*ni].weight < UINT8_MAX)   ? _g[an].weight + _g[*ni].weight  : UINT8_MAX;
        }

        // Accumulate in edge information, merges nodes
        ni = nodes.begin();
        ++ni;
        for (; ni != nodes.end(); ++ni) {
            InEdgeIter ii, ie;
            VtxDesc n = *ni;
            for (boost::tie(ii, ie) = boost::in_edges(n, _g); ii != ie; ++ii) {
                VtxDesc n1 = boost::source(*ii, _g);
                EdgeDesc e;
                bool exists;
                boost::tie(e, exists) = edge(n1, an, _g);
                if (exists) {
                    _g[e].count = ((uint32_t)_g[e].count + (uint32_t)_g[*ii].count < UINT8_MAX) ? _g[e].count + _g[*ii].count : UINT8_MAX;
                } else {
                    std::pair<EdgeDesc, bool> p = boost::add_edge(n1, an, _g);
                    _g[p.first].count = _g[*ii].count;
                    _g[p.first].visited = _g[*ii].visited;
                }
            }
            markForReaper(n);
        }
        mergeInNodes(an);
    }
}

void AlnGraphBoost::mergeOutNodes(VtxDesc n) {
    std::map<char, std::vector<VtxDesc> > nodeGroups;
    OutEdgeIter oi, oe;
    for(boost::tie(oi, oe) = boost::out_edges(n, _g); oi != oe; ++oi) {
        VtxDesc outNode = boost::target(*oi, _g);
        if (in_degree(outNode, _g) == 1) {
            nodeGroups[_g[outNode].base].push_back(outNode);
        }
    }

    for(std::map<char, std::vector<VtxDesc> >::iterator kvp = nodeGroups.begin(); kvp != nodeGroups.end(); ++kvp) {
        std::vector<VtxDesc> nodes = (*kvp).second;
        if (nodes.size() <= 1)
            continue;

        std::vector<VtxDesc>::const_iterator ni = nodes.begin();
        VtxDesc an = *ni++;
        InEdgeIter anii, anie;
        boost::tie(anii, anie) = boost::in_edges(an, _g);

        // Accumulate inner edge information
        for (; ni != nodes.end(); ++ni) {
            InEdgeIter ii, ie;
            boost::tie(ii, ie) = boost::in_edges(*ni, _g);
            _g[*anii].count = ((uint32_t)_g[*anii].count + (uint32_t)_g[*ii].count < UINT8_MAX) ? _g[*anii].count + _g[*ii].count : UINT8_MAX;
            _g[an].weight   = ((uint32_t)_g[an].weight + (uint32_t)_g[*ni].weight < UINT8_MAX)  ? _g[an].weight + _g[*ni].weight  : UINT8_MAX;
        }

        // Accumulate and merge outer edge information
        ni = nodes.begin();
        ++ni;
        for (; ni != nodes.end(); ++ni) {
            OutEdgeIter oi, oe;
            VtxDesc n = *ni;
            for (boost::tie(oi, oe) = boost::out_edges(n, _g); oi != oe; ++oi) {
                VtxDesc n2 = boost::target(*oi, _g);
                EdgeDesc e;
                bool exists;
                boost::tie(e, exists) = edge(an, n2, _g);
                if (exists) {
                    _g[e].count = ((uint32_t)_g[e].count + (uint32_t)_g[*oi].count < UINT8_MAX) ? _g[e].count + _g[*oi].count : UINT8_MAX;
                } else {
                    std::pair<EdgeDesc, bool> p = boost::add_edge(an, n2, _g);
                    _g[p.first].count = _g[*oi].count;
                    _g[p.first].visited = _g[*oi].visited;
                }
            }
            markForReaper(n);
        }
    }
}

void AlnGraphBoost::markForReaper(VtxDesc n) {
    _g[n].deleted = true;
    clear_vertex(n, _g);
    _reaperBag.push_back(n);
}

void AlnGraphBoost::reapNodes() {
    int reapCount = 0;
    std::sort(_reaperBag.begin(), _reaperBag.end());
    std::vector<VtxDesc>::iterator curr = _reaperBag.begin();
    for (; curr != _reaperBag.end(); ++curr) {
        assert(_g[*curr].backbone==false);
        remove_vertex(*curr-reapCount++, _g);
    }
}

const std::string AlnGraphBoost::consensus(uint8_t minWeight) {
    // get the best scoring path
    std::vector<AlnNode> path = bestPath();

    // consensus sequence
    std::string cns;

    // track the longest consensus path meeting minimum weight
    int offs = 0, bestOffs = 0, length = 0, idx = 0;
    bool metWeight = false;
    std::vector<AlnNode>::iterator curr = path.begin();
    for (; curr != path.end(); ++curr) {
        AlnNode n = *curr;
        if (n.base == _g[_enterVtx].base || n.base == _g[_exitVtx].base)
            continue;

        cns += n.base;

        // initial beginning of minimum weight section
        if (!metWeight && n.weight >= minWeight) {
            offs = idx;
            metWeight = true;
        } else if (metWeight && n.weight < minWeight) {
        // concluded minimum weight section, update if longest seen so far
            if ((idx - offs) > length) {
                bestOffs = offs;
                length = idx - offs;
            }
            metWeight = false;
        }
        idx++;
    }

    // include end of sequence
    if (metWeight && (idx - offs) > length) {
        bestOffs = offs;
        length = idx - offs;
    }

    return cns.substr(bestOffs, length);
}

const std::string AlnGraphBoost::consensusNoSplit(uint8_t  minWeight,
                                                  uint32_t *templateToFinal,
                                                  uint32_t  templateLength) {
    // get the best scoring path
    std::vector<AlnNode> path = bestPath();

    // consensus sequence
    std::string cns;

    for (uint32_t ii=0; ii<templateLength; ii++)
      templateToFinal[ii] = UINT32_MAX;

    // track the longest consensus path meeting minimum weight
    int offs = 0, offMax = 0, idx = 0;
    bool metWeight = false;
    std::vector<AlnNode>::iterator curr = path.begin();
    for (; curr != path.end(); ++curr) {
        AlnNode n = *curr;
        if (n.base == _g[_enterVtx].base || n.base == _g[_exitVtx].base)
            continue;
        cns += n.base;
        if (metWeight == false && n.weight >= minWeight) {
           metWeight = true;
           offs=idx;
        }
        if (n.weight >= minWeight && idx > offMax) {
           offMax = idx;
        }
        if (n.tpos != UINT32_MAX) {
            templateToFinal[n.tpos] = idx + 1;
        }
        idx++;
    }

    for (uint32_t ii=0; ii<templateLength; ii++)  //  Adjust templateToFinal map:
      if (templateToFinal[ii] == UINT32_MAX)      //    Base not referenced,
        ;                                         //      leave it alone.

      else if (templateToFinal[ii] < offs)        //    Base trimmed out,
        templateToFinal[ii] = 0;                  //      reset to first base.

      else if (templateToFinal[ii] > offMax)      //    Base trimmed out,
        templateToFinal[ii] = offMax - offs;      //      reset to last base.

      else                                        //    Base referenced,
        templateToFinal[ii] -= offs;              //      trim off offset.

    for (uint32_t ii=templateLength; ii--; )      //  Reset all mappings at the end
      if (templateToFinal[ii] == UINT32_MAX)      //  to the end of the cns.
        templateToFinal[ii] = offMax - offs;
      else
        break;

    return cns.substr(offs, (offMax-offs));
}

void AlnGraphBoost::consensus(std::vector<CnsResult>& seqs, uint8_t minWeight, size_t minLen) {
    seqs.clear();

    // get the best scoring path
    std::vector<AlnNode> path = bestPath();

    // consensus sequence
    std::string cns;

    // track the longest consensus path meeting minimum weight
    int offs = 0, idx = 0;
    bool metWeight = false;
    std::vector<AlnNode>::iterator curr = path.begin();
    for (; curr != path.end(); ++curr) {
        AlnNode n = *curr;
        if (n.base == _g[_enterVtx].base || n.base == _g[_exitVtx].base)
            continue;

        cns += n.base;

        // initial beginning of minimum weight section
        if (!metWeight && n.weight >= minWeight) {
            offs = idx;
            metWeight = true;
        } else if (metWeight && n.weight < minWeight) {
        // concluded minimum weight section, add sequence to supplied vector
            metWeight = false;
            CnsResult result;
            result.range[0] = offs;
            result.range[1] = idx;
            size_t length = idx - offs;
            result.seq = cns.substr(offs, length);
            if (length >= minLen) seqs.push_back(result);
        }
        idx++;
    }

    // include end of sequence
    if (metWeight) {
        size_t length = idx - offs;
        CnsResult result;
        result.range[0] = offs;
        result.range[1] = idx;
        result.seq = cns.substr(offs, length);
        if (length >= minLen) seqs.push_back(result);
    }
}

const std::vector<AlnNode> AlnGraphBoost::bestPath() {
    EdgeIter ei, ee;
    for (boost::tie(ei, ee) = edges(_g); ei != ee; ++ei)
        _g[*ei].visited = false;

    std::map<VtxDesc, EdgeDesc> bestNodeScoreEdge;
    std::map<VtxDesc, int64_t> nodeScore;
    std::queue<VtxDesc> seedNodes;

    // start at the end and make our way backwards
    seedNodes.push(_exitVtx);
    nodeScore[_exitVtx] = 0;

    while (true) {
        if (seedNodes.size() == 0)
            break;

        VtxDesc n = seedNodes.front();
        seedNodes.pop();

        bool bestEdgeFound = false;
        int64_t bestScore = INT64_MIN;
        EdgeDesc bestEdgeD = boost::initialized_value;
        OutEdgeIter oi, oe;
        for(boost::tie(oi, oe) = boost::out_edges(n, _g); oi != oe; ++oi) {
            EdgeDesc outEdgeD = *oi;
            VtxDesc outNodeD = boost::target(outEdgeD, _g);
            int64_t newScore, score = nodeScore[outNodeD];
            AlnNode bbNode = _g[_bbMap[outNodeD]];
            newScore = (uint64_t)_g[outEdgeD].count - round((uint64_t)bbNode.coverage*0.5f) + score;

            if (newScore > bestScore) {
                bestScore = newScore;
                bestEdgeD = outEdgeD;
                bestEdgeFound = true;
            }
        }

        if (bestEdgeFound) {
            nodeScore[n]= bestScore;
            bestNodeScoreEdge[n] = bestEdgeD;
        }

        InEdgeIter ii, ie;
        for (boost::tie(ii, ie) = boost::in_edges(n, _g); ii != ie; ++ii) {
            EdgeDesc inEdge = *ii;
            _g[inEdge].visited = true;
            VtxDesc inNode = boost::source(inEdge, _g);
            int notVisited = 0;
            OutEdgeIter oi, oe;
            for (boost::tie(oi, oe) = boost::out_edges(inNode, _g); oi != oe; ++oi) {
                if (_g[*oi].visited == false)
                    notVisited++;
            }

            // move onto the target node after we visit all incoming edges for
            // the target node
            if (notVisited == 0)
                seedNodes.push(inNode);
        }
    }

    // construct the final best path
    VtxDesc prev = _enterVtx, next;
    std::vector<AlnNode> bpath;
    while (true) {
        bpath.push_back(_g[prev]);
        if (bestNodeScoreEdge.count(prev) == 0) {
            break;
        } else {
            EdgeDesc bestOutEdge = bestNodeScoreEdge[prev];
            _g[prev].bestOutEdge = bestOutEdge;
            next = boost::target(bestOutEdge, _g);
            _g[next].bestInEdge = bestOutEdge;
            prev = next;
        }
    }

    return bpath;
}

bool AlnGraphBoost::danglingNodes() {
    VtxIter curr, last;
    boost::tie(curr, last) = boost::vertices(_g);
    bool found = false;
    for (;curr != last; ++curr) {
        if (_g[*curr].deleted)
            continue;
        if (_g[*curr].base == _g[_enterVtx].base || _g[*curr].base == _g[_exitVtx].base)
            continue;

        int indeg = out_degree(*curr, _g);
        int outdeg = in_degree(*curr, _g);
        if (outdeg > 0 && indeg > 0) continue;

        found = true;
    }
    return found;
}

AlnGraphBoost::~AlnGraphBoost(){}
