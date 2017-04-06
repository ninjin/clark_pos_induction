// $Id: traintwohmm.cpp,v 1.2 2002/08/14 13:46:18 clark Exp $
// Train a twohmm on the data
// with many options
// simple corpus format.
// and then tag a corpus with it


#include <iostream>
#include <unistd.h> /* for: getopt() */
#include <stdlib.h> // rnd
#include <time.h> // srand(time)

#include "simplecorpus.h"
#include "clusters.h"
#include "twohmm.h"


using namespace std;


void printUsage(const char * arg0)
{
  cerr << "Usage:" << arg0
       << " traincorpus outputfile" << endl;
  cerr << "Options:\n" 
       <<   " [-n numberStates]  the number of states or clusters at the word level \n"
       <<    " [-s numberSubStates] number of states in the letter level HMM \n" 
       << "  [-i iterations] \n"
       <<  "[-r]  set the random seed" 
       <<   " [-m mixture] 0 means no HMM, 1.0 means no memorisation (default 0.0)\n" 
       <<   " [-f frequency]  is the minimum frequency for a word to be memorised \n";
  exit(-1);
}

int main(int argc, char* argv[])
{
  char ch;                   /* to hold command line option */
  //  char *optstr = "sf:v:3@"; 
  char* optstr = "n:s:i:rhm:f:";
  int numberStates = 32;
  int numberSubStates = 6;

  while( -1 != (ch=getopt(argc,argv,optstr))) {
    switch(ch) {
    case 'h':
      printUsage(argv[0]);
      break;
      
    case 'r':           /* this is just the flag -@ -- no argument expected */
      cerr << "Setting random seed\n";
      srand(unsigned(time(NULL)));
      break;

    case 'n':
      cerr << "n option " << optarg << endl;
      numberStates = atoi(optarg);
      cerr << "numberStates " << numberStates << endl;
      break;
    case 'm':
      cerr << "m (mixture option) option " << optarg << endl;
      INITIAL_MIXTURE = atof(optarg);
      cerr << "INITIAL_MIXTURE "  << INITIAL_MIXTURE << endl;
      break;
    case 'f':
      cerr << "f (minimum frequency option) option " << optarg << endl;
      MIN_FREQUENCY  = atoi(optarg);
      cerr << "MIN_FREQUENCY  "  << MIN_FREQUENCY  << endl;
      break;
    case 's':
      cerr << "s option " << optarg << endl;
      numberSubStates = atoi(optarg);
      cerr << "numberSubStates " << numberSubStates << endl;
      break;
	
    case 'i':
      cerr << "i option " << optarg << endl;
      MAX_ITERATIONS = atoi(optarg);
      break;

    case '?':
      cerr << "unrecognized option: " << optopt << endl;
      printUsage(argv[0]);
      break;
    }
  }
  
  if (optind + 2 != argc){
    cerr << "Wrong number of arguments " << endl;
    printUsage(argv[0]);
  }
  const char * trainCorpus = argv[optind];
  const char * outputFile = argv[optind + 1];
  cerr << "Started ok." << endl;
  SimpleCorpusOne corpus(trainCorpus,NULL);
  cerr << "Loaded corpus ok." << endl;
  TwoHMM* thmmPtr = trainTwoHMM(corpus,
				numberStates,
				numberSubStates,
				255);
  cerr << "Finished training\nStarting to store\n";
  thmmPtr->store(outputFile);
  cerr << "Finished storing .. exiting\n" << endl;
  return 0;
}

