
#include "graph.H"

using namespace std;

void
edge::addCount(int val) {
  readCnt += val;
}

/***
 *  Edges are always inserted from a smaller mer (this node) to a larger mer (toNode).
 *  nodesFwd: "to" nodes, from this node in +.
 *  nodesRev: "to" nodes, form this node in -.
 */
void
node::addEdge(node* toNode, char orientFromNode, char orientToNode) {

  //  check if edge exists. add new if not, or increment count if it does.
  // fprintf(stderr, "[ DEBUG ] :: We are in node %s %c -> %s %c.\n", this->nodeName, orientFromNode, toNode->nodeName, orientToNode);

  bool added = false;

  if (orientFromNode == MER_FWD) {
  //  #pragma omp critical (fwdlock)
    {
      auto it = find(nodesFwd.begin(), nodesFwd.end(), toNode);
      // fprintf(stderr, "[ DEBUG ] :: Does edge match orientation? %s -> %s found; %c -> %c ?\n",
      //    this->nodeName, toNode->nodeName, orientFromNode, orientToNode);

      // fprintf(stderr, "[ DEBUG ] :: it vs. nodesFwd.end(): %d vs. %d\n", it, nodesFwd.end());

      if (it != nodesFwd.end()) {
        //  edge exists.
        edge* e = edgesFwd.at(it - nodesFwd.begin());
        // fprintf(stderr, "[ DEBUG ] :: Does edge match orientation at fromNode - %s %c -> %s %c found; %c -> %c ?\n",
        //    this->nodeName, e->orientFrom, toNode->nodeName, e->orientTo, orientFromNode, orientToNode);

        if (e->orientTo == orientToNode) {
          //  orientation of the "to" node also matches.
          e->addCount(1);
          added = true;
          //  fprintf(stderr, "[ DEBUG ] :: Increase edge count for %s %c -> %s %c by 1 : %d\n", 
          //    this->nodeName, orientFromNode, toNode->nodeName, orientToNode, e->readCnt);
        }
      }
      if (!added) {
        //  new edge
        edge* newEdge = new edge();
        newEdge->fromNode = this;
        newEdge->toNode   = toNode;
        newEdge->orientFrom = orientFromNode;
        newEdge->orientTo   = orientToNode;
        newEdge->addCount(1);
        nodesFwd.push_back(toNode);
        edgesFwd.push_back(newEdge);
        // fprintf(stderr, "[ DEBUG ] :: New edge made! I am %s %c -> %s %c\n",
        //    this->nodeName, orientFromNode, toNode->nodeName, orientToNode);
      }
    }
  } else {
  //  #pragma omp critical (revlock)
    {
      auto it = find(nodesRev.begin(), nodesRev.end(), toNode);
      // fprintf(stderr, "[ DEBUG ] :: Does edge match orientation? %s -> %s found; %c -> %c ?\n",
      //    this->nodeName, toNode->nodeName, orientFromNode, orientToNode);

      if (it != nodesRev.end()) {
        //  edge exists.
        edge* e = edgesRev.at(it - nodesRev.begin());
        // fprintf(stderr, "[ DEBUG ] :: Does edge match orientation at fromNode - %s %c -> %s %c found; %c -> %c ?\n",
        //    this->nodeName, e->orientFrom, toNode->nodeName, e->orientTo, orientFromNode, orientToNode);

        if (e->orientTo == orientToNode) {
          //  orientation of the "to" node also matches.
          e->addCount(1);
          added = true;
          // fprintf(stderr, "[ DEBUG ] :: Increase edge count for %s %c -> %s %c by 1 : %d\n", 
          //    this->nodeName, orientFromNode, toNode->nodeName, orientToNode, e->readCnt);
        }
      }
      if (!added) {
        //  new edge
        edge* newEdge = new edge();
        newEdge->fromNode = this;
        newEdge->toNode   = toNode;
        newEdge->orientFrom = orientFromNode;
        newEdge->orientTo   = orientToNode;
        newEdge->addCount(1);
        nodesRev.push_back(toNode);
        edgesRev.push_back(newEdge);
        // fprintf(stderr, "[ DEBUG ] :: New edge made! I am %s %c -> %s %c\n",
        //    this->nodeName, orientFromNode, toNode->nodeName, orientToNode);
      }
    }   
  }
  // fprintf(stderr, "[ DEBUG ] :: Size of nodesFwd vs. nodesRev : %d vs. %d\n", nodesFwd.size(), nodesRev.size());
}
