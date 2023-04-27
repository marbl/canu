
/******************************************************************************
 *
 *  This is a tool for linking haplotype specific k-mers in long-reads.
 *
 *  This software is based on:
 *    'Meryl'                  (https://github.com/marbl/meryl)
 *
 *  This is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "kmers.H"
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

class node;

class edge {
public:
  edge()  { }
  ~edge() { }

  int    readCnt = 0;  //  weight.
  node  *fromNode;
  node  *toNode;
  char   orientFrom;   //  orientation of the "from" node 
  char   orientTo;     //  orientation of the "to"   node

  void   addCount(int val);
  int    getRC(void) {  return(readCnt); }
};


class node {
public:
  node()  { }
  ~node() { }

  kmer    nodeId;   //  node name. uint64
  //  char    groupId = 'U';  //  group ID --- replacing with col
  //  Could be anything. M (Maternal), P (Paternal), U (Unknown), H (Homozygous), A, B, ...
  const char*   col = "#AAAAAA";  //  default is gray
  kmvalu        val = 0;          //  default is 0

  static const char MER_FWD='+';
  static const char MER_REV='-';

public:
  vector<edge*>  edgesFwd;  //  outgoing edges from FWD node
  vector<edge*>  edgesRev;  //  outgoing edges from REV node
  vector<node*>  nodesFwd;  //  outgoing nodes in FWD
  vector<node*>  nodesRev;  //  outgoing nodes in REV

public:
  void setNode(kmer mer, kmvalu val) {
    this->nodeId = mer;
    this->val    = val;
  }

  void    addEdge(node*, char, char);

  void  nodeName(char* kmerstring) {
    this->nodeId.toString(kmerstring);
  }

  kmer    reverseComplement(void) {
    kmer _rmer;
    char kmerstring[65];
    this->nodeName(kmerstring);

    for (int _bufferPos = 0; _bufferPos < this->nodeId.merSize(); _bufferPos++) {
      _rmer.addL(kmerstring[_bufferPos]);
    }
    return _rmer;
  }

};

