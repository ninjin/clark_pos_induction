#ifndef ASC_CLUSTERING_H
#define ASC_CLUSTERING_H 1

/**
 * $Id: clusters.h,v 1.3 2002/08/14 13:46:18 clark Exp $
 */

#include <vector>

#include "matrix.h"
#include "simplecorpus.h"
#include "hmm.h"

// All words that have frequency <= FREQ_CUTOFF are in 
// a separate cluster

extern bool USE_TRUE_WEIGHT;
extern int MAX_ITERATIONS;
extern int FREQ_CUTOFF;
extern double PRIOR_BOOST;
extern bool FULL_MORPHOLOGY_WEIGHT;
// this defines a clustering of a set of words

class Clusters
{
 private:
  int numberClasses;
  int numberTypes;
  int numberTokens;
  int numberStates;
  int alphabetSize;
  int* data;
  int* next;
  const SimpleCorpusOne & corpus;
  
  // a vector of hmms one per class
  vector<HMM *> hmms;

  MatrixInt clusterBigrams;
  vector<int> clusterUnigrams;
  vector<int> first;
  vector<int> classVector;
  vector<int> counts;
  vector<int> sortedWords;  // most frequent zero
  
 public:
  // this could be simplified a lot
  Clusters(int numberClasses_, 
	   const SimpleCorpusOne & corpus,
	   int numberStates_,
	   int alphabetSize_,
	   bool randomised);
  ~Clusters() 
    { 
      delete[] next; 
    };
  void dump(const SimpleCorpusOne & ) const;
  void clusterNeyEssen();
  int reclusterNeyEssen();
  
  bool bestCluster(int, 	const vector<double> & prior);
  double calcChange(int , int , const vector<int> & , const vector<int> &, int) const;
  // optimised variants
  double calcChangeFast(int w, 
			double oldStringProb,
			const vector<int> & wordVector,
			int newCluster, 
			const vector<int> & , 
			const vector<int> & , 
			int doubles,
			const vector<double> & prior) const;
  void moveWord(int w,
		int oldCluster,
		int newCluster,
		int doubles,
		vector<int> & left,
		vector<int> & right);
  bool validateBigramCounts();
};


#endif
