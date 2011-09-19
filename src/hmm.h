#ifndef ASC_HMM_H
#define ASC_HMM_H 1

/**
 * $Id: hmm.h,v 1.3 2002/08/14 13:46:18 clark Exp $
 */

#include "matrix.h"
#include <vector>
#include <string>
#include <fstream>

using namespace std;

class HMM
{
  friend HMM* loadHMM(ifstream & in );
 private:
  static const int initialState = 0;
  static  const int finalState = 1;
  const int numberStates;
  const int numberLetters;
  static const int terminator = 0;
  

  // state transitions P(s'|s)
  Matrix p;
  // output probabilities P(l|s)
  Matrix q;

 public:
  HMM(int numberStates_, int  numberLetters_);

  // store to stream
  void store(ofstream &) const;
  

  void randomise();
  Matrix* calculateForward(const vector<int> &) const;
  Matrix* calculateBackward(const vector<int> &) const; 
  Matrix* calculateLogForward(const vector<int> &) const;
  Matrix* calculateLogBackward(const vector<int> &) const;
  double& getP(int from, int to){
    return p(from,to);
  }
  double& getQ(int state, int letter){
    return q(state,letter);
  }
  bool validate() const;
  bool validateState(int state) const;
  void normalise();
  void normaliseState(int state);
  void normaliseStateOutputs(int state);
  void normaliseStateTransitions(int state);
  
  vector<int> *  generate() const;
  static void convertString(const string input, vector<int> & output);
  double outputSum(int state) const;
  double transitionSum(int state) const;
  double probability(const vector<int> &) const;
  double treeEntropy() const;
  double entropy() const;
  int randomNextState(int state) const;
  int randomOutput(int state) const;
  void dumpString(const vector<int> &) const;
  // training
  HMM* em(const vector< const vector<int> * > & dataVector) const;
  double emSingle(HMM & newHMM, double weight, const vector<int> & data) const;
};

// using the forward-backward algorithm

HMM* trainHMM(const vector< const vector<int> * > & dataVector, 
	      int numberStates,
	      int numberLetters);


HMM* loadHMM(ifstream & in );

#endif
