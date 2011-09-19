

// $Id: neyessen.cpp,v 1.4 2002/08/14 13:46:18 clark Exp $
// Do the Ney Essen clustering on some data in the 
// simple corpus format.





#include <iostream>
#include "simplecorpus.h"
#include "clusters.h"
#include <unistd.h> /* for: getopt() */



using namespace std;

void printUsage()
{
  cerr << "Usage -- cluster_neyessen "
       << "[-p factor ] use prior probabilities [-i iterations] [-mMinCount] [-r](random initialisation) corpus additionalCorpus clusters" << endl;
  exit(-1);
}


int main(int argc, char* argv[])
{
  char ch;                   /* to hold command line option */
  //  char *optstr = "sf:v:i:3@"; 
  char* optstr = "m:i:rhp:";

  /* 
   * the option string will be passed to getopt(3), the format
   * of our string "sf:v:" will allow us to accept -s as a flag,
   * and -f or -v with an argument, the colon suffix tells getopt(3)
   * that we're expecting an argument.  Eg:  optest -s -f this -v8
   *
   * getopt(3) takes our argc, and argv, it also takes
   * the option string we set up earlier.  It will assign
   * the switch character to ch, and -1 when there are no more
   * command line options to parse.
   *
   */
  bool randomInitialization = false;
  while( -1 != (ch=getopt(argc,argv,optstr))) {
    switch(ch) {

    case 'h':
      printUsage();
      break;
    case 'i':

      MAX_ITERATIONS = atoi(optarg);
      cerr << " iterations option " << MAX_ITERATIONS << endl;      
      break;
    case 'p':
      PRIOR_BOOST = atof(optarg);
      cerr << " p option prior factor " << PRIOR_BOOST << endl;
      break;

    case 'r':           /* this is just the flag -@ -- no argument expected */
      randomInitialization = true;
      cerr << "r option = random intialiszation" << endl;
      break;
      
      /* 
       * for -m  we expect an argument, which getopt(3) 
       * will lstore in the buffer optarg 
       * 
       */
    case 'm':
      cerr << "m option " << optarg << endl;
      FREQ_CUTOFF = atoi(optarg);
      cerr << "Frequency cutoff is " << FREQ_CUTOFF << endl;
      break;
	
    case '?':
      cerr << "unrecognized option: " << optopt << endl;
      printUsage();
      break;
    }
  }
  cerr << "Frequency cutoff is " << FREQ_CUTOFF << endl;
  if (optind + 3 != argc){
    cerr << "Wrong number of arguments " << endl;
    printUsage();
  }
  cerr << "Started ok." << endl;
  SimpleCorpusOne corpus(argv[optind],argv[optind + 1]); 
  cerr << "Loaded corpus ok." << endl;
  int clusters = atoi(argv[optind + 2]);
  assert(clusters > 0);
  Clusters cne(clusters, corpus, 0, 0, randomInitialization);
  cerr << "Created initial clustering" << endl;
  cne.clusterNeyEssen();
  cne.dump(corpus);
  return 0;
}

