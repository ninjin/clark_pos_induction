#include <iostream>
#include <fstream>
#include <algorithm>

#include "simplecorpus.h"

// $Id: simplecorpus.cpp,v 1.4 2002/08/14 13:46:18 clark Exp $


SimpleCorpusOne::
~SimpleCorpusOne()
{
  for (size_t i = 0; i < wordArray.size(); i++)
    delete wordArray[i];
  delete[] data;
}

const int BUFLEN = 1000;

// basic constructor
// loads it from filename

const char * BOUNDARY = "</s>";

SimpleCorpusOne::
SimpleCorpusOne(const char* filename,
		const char* additional)
  : 
  numberTypes(0), 
  numberTokens(0)
{
  ifstream in(filename);
  if (!in)
    {
      cerr << "Couldn't open file." << endl;
      exit(-1);
    }
  char buffer[BUFLEN]; 
  while (!in.eof()){
    in.getline(buffer,BUFLEN);
    numberTokens++;
  }
  in.close();
  ifstream inn(filename);
  //  inn.open(filename);
  data = new int[numberTokens];
  //in2.seekg(0);
  int currentToken = 0;
  int currentType = 0;
  //int pass2 = 0;
 while (!inn.eof()){
    inn.getline(buffer,BUFLEN);
    // pass2++;
    string  word(buffer);
    if (dictionary.find(word) == dictionary.end()){
      // new word
      //  cerr << "New word " << word << endl;
      dictionary[word] = currentType;
      data[currentToken] = currentType;
      currentType++;
      wordArray.push_back(new string(word));
      countArray.push_back(1);
    }
    else {
      //cerr << "old word " << endl;
      int wordIndex = dictionary[word];
      countArray[wordIndex]++;
      data[currentToken] = wordIndex;
    }
    currentToken++;
  }
 // cerr << "Pass2 " << pass2 << endl;
  numberTypes = currentType;
 
  cerr << "read " << numberTypes << " types" << endl;
  cerr << "read " << numberTokens << "  tokens" << endl;
  if (additional){
    cerr << "Starting to read additional file for vocabulary\n";
    ifstream in2(additional);
    if (!in2)
      {
	cerr << "Couldn't open file." << endl;
      }
    else {
      cerr << "Opened file " << additional << endl;
      int newTokens = 0;
      while (!in2.eof()){
	newTokens++;
	in2.getline(buffer,BUFLEN);
	string  word(buffer);
	if (dictionary.find(word) == dictionary.end()){
	  // new word
	  //	  cerr << "New word " << word << endl;
	  dictionary[word] = currentType;
	  currentType++;
	  wordArray.push_back(new string(word));
	  countArray.push_back(0);
	}
      }
      cerr << "Read " << newTokens <<" new tokens\n";
    }
  }
  numberTypes = currentType;
  cerr << "read " << numberTypes << " total types" << endl;
}

/**
 * Return the index or -1 if it doesn;t exist
 */

int 
SimpleCorpusOne::
lookUpWord(const string & word) 
  const
{
  map<string,int>::const_iterator pos = dictionary.find(word);
  if (pos != dictionary.end())
    return pos->second;
  else
    return -1;
}

/**
 * Count the number of types in the corpus
 * whose frequency is greater than or equal to f
 */

int 
SimpleCorpusOne::
countWords(int f) 
  const
{
  int count = 0;
  for (int i = 0; i< numberTypes; i++){
    if (countArray[i] >= f)
      count++;
  }
  return count;
}

/**
 * Find the next occurrence of string in the corpus
 * that is AFTER start.
 * If there is none return the size of the array.
 */

int 
SimpleCorpusOne::
findNextOccurrence(const string & target, 
		   int start)
  const
{
  int w = lookUpWord(target);
  for (int index = start + 1; index < numberTokens; index++){
    if (data[index] == w)
      return index;
  }
  return numberTokens;
}
