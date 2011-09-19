// $Id: baseline.cpp,v 1.4 2002/08/14 13:46:18 clark Exp $
// Create a simple baseline





#include <iostream>
#include "simplecorpus.h"
#include "clusters.h"

using namespace std;

int main(int argc, char* argv[])
{
  cerr << "Started ok." << endl;
  if (argc != 4){
    cerr << "Usage -- cluster_baseline corpus  additional clusters" << endl;
    exit(-1);
  }
  SimpleCorpusOne corpus(argv[1],argv[2]); 
  cerr << "Loaded corpus ok." << endl;
  int clusters = atoi(argv[3]);
  assert(clusters > 0);
  Clusters cne(clusters,corpus, 0, 0, false);
  cerr << "Created initial clustering" << endl;
  //  cne.clusterNeyEssen();
  cne.dump(corpus);
  return 0;
}

