// $Id: tagtwohmm.cpp,v 1.2 2002/08/14 13:46:18 clark Exp $
// Tag a corpus with a traintwohmm


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
  cerr << "Usage: " << arg0
       << " modelfile testcorpus " << endl;
  exit(-1);
}

int main(int argc, char* argv[])
{
  char ch;                   /* to hold command line option */
  //  char *optstr = "sf:v:3@"; 
  char* optstr = "h";

  while( -1 != (ch=getopt(argc,argv,optstr))) {
    switch(ch) {
    case 'h':
      printUsage(argv[0]);
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
  const char * modelFile = argv[optind];
  const char * testCorpus = argv[optind + 1];
  cerr << "Started ok." << endl;

  SimpleCorpusOne corpus(testCorpus,NULL);
  cerr << "Loaded corpus ok." << endl;
  TwoHMM* thmmPtr = loadTwoHMM(modelFile);
  cerr << "Loaded model " << endl;
  // now tag it
  vector<string *> sentence;
  int index = -1;
  while (index < corpus.numberTokens){
    int next = corpus.findNextOccurrence(SENTENCE_BOUNDARY, index);
    int size = next - (index+1);
    if (size > 0){
      sentence.resize(size);
      for (int i = 0; i < size; i++)
	sentence[i] = corpus.wordArray[corpus.data[index + i + 1]];
      vector<int> tags;
      thmmPtr->tagTestSentence(sentence,tags);
      for (int i = 0; i < size; i++){
	cout << *(sentence[i]) << " " << tags[i] << endl;
      }
    }
    index = next;
    cout << endl;
    //cerr << "Sentence " << i << " length " <<  sentence.size() << endl;
  }
  cerr << "Finished storing .. exiting\n" << endl;
  return 0;
}

