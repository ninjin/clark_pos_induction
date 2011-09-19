#ifndef ASC_CORPUS_H
#define ASC_CORPUS_H 1

/**
 * $Id: simplecorpus.h,v 1.3 2002/08/14 13:46:18 clark Exp $
 */

#include <string>
#include <vector>
#include <map>
using namespace std;

//
// load up the corpus
// this is in a simple one-line per word format
// keep a hash table
// of previously seen ones
// linked to an index
// 

class SimpleCorpusOne
{
public:
  SimpleCorpusOne(const char* filename, const char * additional);
  ~SimpleCorpusOne();
  int numberTypes;
  int numberTokens;
  int boundaryType;
  int* data;
  // array stores them as vector of integers
  vector<string *> wordArray;
  vector<int> countArray;
  // dictionary that looks up words in 
  map<string, int> dictionary;
  void dump() const;
  const string & getNthWord(int i) const
    {
      return *(wordArray[i]);
    }
  const string & getNthDataPoint(int i) const
    {
      int w = data[i];
      return  *(wordArray[w]);
    }
  int lookUpWord(const string &) const;
  // find next occurrence of the word
  int findNextOccurrence(const string &, int start) const;
  int countWords(int f) const;
};

#endif
