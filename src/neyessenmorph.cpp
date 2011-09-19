// $Id: neyessenmorph.cpp,v 1.3 2002/08/14 13:46:18 clark Exp $

// Do the ney essen clustering with morphological information
// and random initialisation






#include <iostream>
#include "simplecorpus.h"
#include "clusters.h"


#include <unistd.h> /* for: getopt() */
#include <stdlib.h>
#include <time.h>

using namespace std;

void printUsage(char * arg0)
{
  cerr << "Usage: " << arg0 << "[-r] -t (use true weight for HMM training) [-s numberStates] [-m frequencyCutoff] [-i iterations] [ -p factor] (use prior weights) \n [-x use extra influence of morphology]"
       << " corpus additional clusters" << endl;
  exit(-1);
}

int main(int argc, char* argv[])
{
  int numberStates = 6;
  char ch;
  char * optstr = "xhs:m:i:rtp:"; 
  while( -1 != (ch=getopt(argc,argv,optstr))) {
    switch(ch) {
    case 'h':
      printUsage(argv[0]);
      break;
    case 'r':
      cerr << "Setting random seed\n";
      srand(unsigned(time(NULL)));
      break;
    case 'p':
      PRIOR_BOOST = atof(optarg);
      cerr << " p option uses prior cluster probabilities " << PRIOR_BOOST << endl;
      break;

    case 't':
      cerr << "Setting use true weight to true\n";
      USE_TRUE_WEIGHT = true;
      break;
    case 'm':
      cerr << "m option " << optarg << endl;
      FREQ_CUTOFF = atoi(optarg);
      cerr << "Frequency cutoff is " << FREQ_CUTOFF << endl;
      break;
    case 's':
      cerr << "s option " << optarg << endl;
      numberStates = atoi(optarg);
      break;
    case 'i':
      cerr << "i option " << optarg << endl;
      MAX_ITERATIONS = atoi(optarg);
      break;
    case 'x':
      cerr << "Using eXtra influence of morphology" << endl;
      FULL_MORPHOLOGY_WEIGHT = true;
      break;
    case '?':
      cerr << "unrecognized option: " << optopt << endl;
      printUsage(argv[0]);
      break;
    }
  }
  // we need three arguments
  if (optind + 3 != argc){
    cerr << "Wrong number of non option arguments\n";
    printUsage(argv[0]);
  }
  SimpleCorpusOne corpus(argv[optind],argv[optind+1]); 
  cerr << "Loaded corpus ok." << endl;
  int clusters = atoi(argv[optind+2]);
  assert(clusters > 0);
  // do randomised initialisation
  Clusters cne(clusters, corpus, numberStates, 255, true);
  cerr << "Created initial clustering" << endl;
  cne.clusterNeyEssen();
  cne.dump(corpus);
  return 0;
}

