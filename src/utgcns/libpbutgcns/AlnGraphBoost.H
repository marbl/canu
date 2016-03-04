
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Sergey Koren beginning on 2015-DEC-28
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
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


#ifndef __GCON_ALNGRAPHBOOST_HPP__
#define __GCON_ALNGRAPHBOOST_HPP__

#include <boost/graph/adjacency_list.hpp>

/// Alignment graph representation and consensus caller.  Based on the original
/// Python implementation, pbdagcon.  This class is modelled after its
/// aligngraph.py component, which accumulates alignment information into a
/// partial-order graph and then calls consensus.  Used to error-correct pacbio
/// on pacbio reads.
///
/// Implemented using the boost graph library.

// forward declaration
//struct Alignment;

// this allows me to forward-declare properties with graph descriptors as
// members types
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS> graphTraits;

/// Graph vertex property. An alignment node, which represents one base position
/// in the alignment graph.
struct AlnNode {
    char base; ///< DNA base: [ACTG]
    int coverage; ///< Number of reads align to this position, but not
                  ///< necessarily match
    int weight; ///< Number of reads that align to this node *with the same base*, but not
                ///< necessarily represented in the target.
    bool backbone; ///< Is this node based on the reference
    bool deleted; ///< mark for removed as part of the merging process
    graphTraits::edge_descriptor bestInEdge; ///< Best scoring in edge
    graphTraits::edge_descriptor bestOutEdge; ///< Best scoring out edge
    AlnNode() {
        base = 'N';
        coverage = 0;
        weight = 0;
        backbone = false;
        deleted = false;
    }
};

/// Graph edge property. Represents an edge between alignment nodes.
struct AlnEdge {
    int count; ///< Number of times this edge was confirmed by an alignment
    bool visited; ///< Tracks a visit during algorithm processing
    AlnEdge() {
        count = 0;
        visited = false;
    }
};

// Boost-related typedefs
// XXX: listS, listS?
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, AlnNode, AlnEdge> G;
typedef boost::graph_traits<G>::vertex_descriptor VtxDesc;
typedef boost::graph_traits<G>::vertex_iterator VtxIter;
typedef boost::graph_traits<G>::edge_descriptor EdgeDesc;
typedef boost::graph_traits<G>::edge_iterator EdgeIter;
typedef boost::graph_traits<G>::in_edge_iterator InEdgeIter;
typedef boost::graph_traits<G>::out_edge_iterator OutEdgeIter;
typedef boost::property_map<G, boost::vertex_index_t>::type IndexMap;

///
/// Simple consensus interface datastructure
///
struct CnsResult {
    int range[2]; ///< Range on the target
    std::string seq; ///< Consensus fragment
};

///
/// Core alignments into consensus algorithm, implemented using the boost graph
/// library.  Takes a set of alignments to a reference and builds a higher
/// accuracy (~ 99.9) consensus sequence from it.  Designed for use in the HGAP
/// pipeline as a long read error correction step.
///
class AlnGraphBoost {
public:
    /// Constructor.  Initialize graph based on the given sequence. Graph is
    /// annotated with the bases from the backbone.
    /// \param backbone the reference sequence.
    AlnGraphBoost(const std::string& backbone);

    /// Constructor.  Initialize graph to a given backbone length. Base
    /// information is filled in as alignments are added.
    /// \param blen length of the reference sequence.
    AlnGraphBoost(const size_t blen);

    /// Add alignment to the graph.
    /// \param Alignment an alignment record (see Alignment.hpp)
    void addAln(dagcon::Alignment& aln);

    /// Adds a new or increments an existing edge between two aligned bases.
    /// \param u the 'from' vertex descriptor
    /// \param v the 'to' vertex descriptor
    void addEdge(VtxDesc u, VtxDesc v);

    /// Collapses degenerate nodes (vertices).  Must be called before
    /// consensus(). Calls mergeInNodes() followed by mergeOutNodes().
    void mergeNodes();

    /// Recursive merge of 'in' nodes.
    /// \param n the base node to merge around.
    void mergeInNodes(VtxDesc n);

    /// Non-recursive merge of 'out' nodes.
    /// \param n the base node to merge around.
    void mergeOutNodes(VtxDesc n);

    /// Mark a given node for removal from graph. Doesn't not modify graph.
    /// \param n the node to remove.
    void markForReaper(VtxDesc n);

    /// Removes the set of nodes that have been marked.  Modifies graph.
    /// Prohibitively expensive when using vecS as the vertex container.
    void reapNodes();

    /// Generates the consensus from the graph.  Must be called after
    /// mergeNodes(). Returns the longest contiguous consensus sequence where
    /// each base meets the minimum weight requirement.
    /// \param minWeight sets the minimum weight for each base in the consensus.
    ///        default = 0
    const std::string consensus(int minWeight=0);

    /// Generates all consensus sequences from a target that meet the minimum
    /// weight requirement.
    void consensus(std::vector<CnsResult>& seqs, int minWeight=0, size_t minLength=500);

    /// Locates the optimal path through the graph.  Called by consensus()
    const std::vector<AlnNode> bestPath();

    /// Locate nodes that are missing either in or out edges.
    bool danglingNodes();

    /// Destructor.
    virtual ~AlnGraphBoost();
private:
    G _g;
    VtxDesc _enterVtx;
    VtxDesc _exitVtx;
    std::map<VtxDesc, VtxDesc> _bbMap;
    std::vector<VtxDesc> _reaperBag;
};

#endif // __GCON_ALNGRAPHBOOST_HPP__
