#ifndef ASC_TWOHMM_H
#define ASC_TWOHMM_H 1

/**
 * $Id: twohmm.h,v 1.3 2002/08/14 13:46:18 clark Exp $
 */


/**
 * This is a header file for the class TwoHMM
 * which covers a number of different models
 * all of which have a HMM with states relating to words
 * and (optionally) output functions that are HMMs with states
 * relating to words or letters.
 */


#include "simplecorpus.h"
#include "matrix.h"
#include "hmm.h"

#include <vector>
#include <map>

using namespace std;

extern int MAX_ITERATIONS;
extern double INITIAL_MIXTURE;
extern double HMM_DECREMENT;
extern double TOKEN_ML_DECREMENT;
extern double TYPE_ML_DECREMENT;
extern int CHUNK_SIZE;
extern int MIN_FREQUENCY;
extern const char * SENTENCE_BOUNDARY;

// a model where each of the output distributions is a larger HMM
// so we can acquire the basics
// either train supervised or unsupervised
// supervised is we have the boundaries
// unsupervised we don't


// VARIANT
// each output model is an interpolation of a ML model and a HMM



class TwoHMM
{
  friend TwoHMM* loadTwoHMM(const char * filename);
 private:
  static const int initialState = 0;  // 0
  static const int finalState = 1;    // 1
  int numberStates;   // number of states in this one
  int numberSubStates; // number of states in each outputHMM
  int numberLetters;
  int numberWords;
  // state transitions P(s'|s)
  Matrix p;
  // output probabilities P(l|s) q(state, word)
  Matrix q;
  // mixing coefficients
  // this gives you the amount of weight for the HMMs
  vector<double> mixing;
  // HMMs
  vector<HMM*> outputHMMs;

  vector<const string *> wordArray;
  map<const string,int> dictionary;
 public:
  TwoHMM(int numberStates, int numberSubStates, int numberLetters, int numberWords);
  void fillFromCorpus(const SimpleCorpusOne &);
  void fillFromCorpus(const SimpleCorpusOne &, int f);
  ~TwoHMM();
  double condQ(int i, int w) const;
  void generate(const SimpleCorpusOne &) const;
  void dump() const;
  void generate2(const SimpleCorpusOne &, int) const;
  int randomNextState(int state) const;
  int randomQ(int state) const;
  // randomises the parameters
  // of this model and all the submodels?
  void randomise();
  void initialiseFromClusters(const SimpleCorpusOne & corpus, const char * clusteringFile);
  // normalises the probabilities of the whole thing
  void normalise();
  // storing and loading
  void store(ofstream &) const;
  void store(const char * filename) const;
  
  double& getP(int from, int to){
    return p(from,to);
  }
  void normaliseState(int state);
  // normalise the submodels
  void normaliseStateOutputs(int state);
  void normaliseStateTransitions(int state);
  void perturbQ();
  // generate a long string of characters
  // 0 is space
  //  vector<int> *  generate() const;
  char convertToChar(int i) const {
    return char(i);
  }
  int lookUpWord(const string &) const;
  void convertString(const string & input, vector<int> & output) const;
  double transitionSum(int state) const;
  // probability of this 0-separated string
  double probability(const vector<int> &) const;

  TwoHMM* em(const SimpleCorpusOne &) const;
  double emSingle(TwoHMM & newTwoHMM, 
		  double weight, 
		  const vector<int> & data) const;
  // we have precomputed the output probabilities for each word
  // and we accumulate the posteriors into each one
  double emSentence(TwoHMM*, 
		    const vector<int> &, 
		    const Matrix & probabilities , 
		    Matrix & posteriors ) const;
  // returns log probability of the whole thing
  double calculateScaledBackward(const vector<int> & sentence, 
				 const Matrix & probs,  
				 Matrix & beta) const;
  double calculateScaledForward(const vector<int> & sentence, 
				const Matrix & probs, 
				Matrix & alpha) const;
  double calculateScaledProbabilities(const vector<int> & sentence, const Matrix & probs, Matrix & alpha,  Matrix & beta) const;
  // tag the top level with the sequences.
  // tag a sentence from the test corpus
  // we can start from a given start state so we cna break it up into chunks easily
  void tagTestSentence(const vector<string *> & sentence, 
		       vector<int> & tags) const;

};


TwoHMM* trainTwoHMM(const SimpleCorpusOne &, 
		    int numberStates,
		    int numberSubStates,
		    int numberLetters);



TwoHMM* loadTwoHMM(const char * filename);

#endif

