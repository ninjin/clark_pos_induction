#include <math.h>
#include <iostream>

/**
 * $Id: hmm.cpp,v 1.5 2002/06/09 16:16:12 clark Exp $
 */

#include "hmm.h"
#include "matrix.h"


const char * HMM_HEADER = "ASC_HMM_HEADER";
const char * HMM_FOOTER = "ASC_HMM_FOOTER";

const int BUFLEN = 1000;

using namespace std;

HMM::
HMM(int numberStates_, 
    int numberLetters_)
  : 
  numberStates(numberStates_), 
  numberLetters(numberLetters_ ), 
  p(numberStates_,numberStates_),
  q(numberStates_,numberLetters_)
{
};

// beta(state, i) 
// prob that we start in state s and generate u_i+1, ...

Matrix* 
HMM::
calculateBackward(const vector<int>& word) 
const 
{
  int l = word.size();
  Matrix * betaPtr = new Matrix(numberStates,l+1);
  // calculate correctly;
  Matrix& beta = *betaPtr;
  beta(finalState,l) = 1.0l;
  for (int i = l-1; i >= 0; i--){
    for (int state = 0; state < numberStates; state++){
      for (int state2 = 0; state2 < numberStates; state2++){
	//	cerr << "state " << state << " state2 " << state2 << endl;
	int letter = word[i];
	assert(letter >=0 && letter < numberLetters);
	beta(state,i) += p(state,state2) * q(state2, letter) * beta(state2,i+1);
	//assert(p(state,state2) <= 1.0l);
	//assert(q(state2,letter) <= 1.0l);
	assert(beta(state,i) <= 1.0001l);
      }
    }
  }
  // beta.dump();
  return betaPtr;
}

// alpha(state, i) 
// prob that we generate the first i of them
// and end in state i
// alpha(0, 0) = 1

Matrix* 
HMM::
calculateForward(const vector<int>& word) const 
{
  int l = word.size();
  Matrix * alphaPtr = new Matrix(numberStates,l+1);
  // calculate correctly;
  Matrix& alpha = *alphaPtr;
  alpha(initialState,0) = 1.0l;
  for (int i = 0; i < l; i++){
    for (int state = 0; state < numberStates; state++){
      for (int state2 = 0; state2 < numberStates; state2++){
	int letter = word[i];
	assert(letter >=0 && letter < numberLetters);
	// start in state mover to state2 and output letter 
	alpha(state2,i+1) +=  alpha(state,i) * p(state,state2) * q(state2,letter);
      }
    }
  }
  return alphaPtr;
}

// calculate the probability of a string

double 
HMM::
probability(const vector<int>& word) 
  const
{
  Matrix* betaPtr = calculateBackward(word);
  double answer = (*betaPtr)(initialState,0);
  delete betaPtr;
  return answer;
}

bool 
HMM::
validate() 
  const
{
  for (int s = 0; s < numberStates; s++)
    if (!validateState(s))
      return false;
  return true;
}

void 
HMM::
normalise() 
{
  for (int s = 0; s < numberStates; s++)
    normaliseState(s);
}

void 
HMM::
normaliseState(int state)
{
  normaliseStateOutputs(state);
  normaliseStateTransitions(state);
}

void 
HMM::
normaliseStateOutputs(int state)
{
  if (state == initialState)
    q(initialState, terminator) = 1.0l;
  if (state == finalState)
    q(finalState, terminator) = 1.0l; // may not be necessary
  double sum = outputSum(state);
  assert(sum >0);
  for (int i = 0; i < numberLetters;i++)
    q(state,i) /= sum;
}

void 
HMM::
normaliseStateTransitions(int state)
{
  if (state == finalState)
    p(finalState,finalState) = 1.0l;
  double sum = transitionSum(state);
  assert(sum > 0);
  for (int i = 0; i < numberStates;i++){
    p(state,i) /= sum;
    if (p(state,i) > 1.0l)
      cerr << "Non-fatal error " << p(state,i) << " should be <= 1 " << endl;
  }
}


bool 
HMM::
validateState(int state) 
  const
{
  if (state != finalState && fabs(transitionSum(state) - 1.0l) > 0.1)
    {
      cerr << "Validation error of state " << state << endl;
      cerr << "Transitions = " << transitionSum(state) << ", and outputs = " << outputSum(state) << endl;
      return false;
    }
  else if (state != initialState &&  fabs(outputSum(state) -1.0l) > 0.1)
    {
      cerr << "Validation error of state " << state << endl;
      q.dump();
      cerr << "Transitions = " << transitionSum(state) << ", and outputs = " << outputSum(state) << endl;
      return false;
    }
  else
    return true;
}

double 
HMM::
transitionSum(int state) 
  const
{
  double sum = 0.0l;
  for (int s = 0; s < numberStates; s++)
    sum += p(state,s);
  return sum;
}

double 
HMM::
outputSum(int state) 
  const
{
  double sum = 0.0l;
  for (int l = 0; l < numberLetters; l++)
    sum += q(state,l);
  return sum;
}

//
// set all the p's and q's to random 
// 

void 
HMM::
randomise() 
{
  q(initialState,terminator) = 1;
  q(finalState,terminator)= 1;
  p(finalState,finalState) = 1;
  for (int state = 0; state < numberStates; state++){
    if (state != initialState && state != finalState){
      for (int letter = 0; letter < numberLetters; letter++){
	if (letter != terminator){
	  q(state,letter) = 1 + double(rand());
	}
      }
    }
  }
  // from state to state2
  for (int state = 0; state < numberStates; state++){
    if (state != finalState){
      for (int state2 = 0; state2 < numberStates; state2++){
	if (state2 != initialState){
	  p(state,state2) = 1 + double(rand());
	}
      }
    }
  }
}


//
// Calculate the entropy of the sequence of transitions
// using Miller and Sullivan
// "Entropy and Combinatorics ..."
// IEEE Trans. Inf. Theory 38 (4) 1992
//


const int MAXIMUM_HMM_SEQUENCE = 1000;
const double EPSILON = 1e-50;

double 
HMM::
treeEntropy() const
{
  // calculate the expected number of visits to each state
  vector<double> en(numberStates,0.0l);
  vector<double> currentStateVector(numberStates,0.0l);
  currentStateVector[initialState] = 1.0l;
  for (int i = 0; i < MAXIMUM_HMM_SEQUENCE; i++){
    // add to en
    double sum = 0;
    for (int state = 0; state < numberStates; state++){
      double v = currentStateVector[state];
      en[state] += v;
      sum += v;
    }
    // cerr << "Iteration " << i << ", sum = " << sum << endl;
    if (sum < EPSILON)
      break;
    // maybe terminate if sum < EPSILON
    // update new states
    vector<double> nextStates(numberStates,0.0l);
    for (int nextState = 0; nextState < numberStates; nextState++)
      for (int state = 0; state < numberStates; state++)
	nextStates[nextState] += currentStateVector[state] * p(state,nextState);
    nextStates[finalState] = 0.0l;
    // change over
    for (int state = 0; state < numberStates; state++)
      currentStateVector[state] = nextStates[state];
  }
  // now calculate the per state entropies
  double entropy = 0.0l;
  for (int state = 0; state < numberStates; state++){
    if (state != initialState && state != finalState){
      double stateEntropy = 0.0l;
      for (int nextState = 0; nextState < numberStates; nextState++){
	double pp = p(state,nextState);
	if (pp > 0)
	  stateEntropy += pp * log(pp);
      }
      for (int letter = 0; letter < numberLetters; letter++){
	double qq = q(state,letter);
	if (qq > 0)
	  stateEntropy += qq * log(qq);
      }
      entropy -= stateEntropy * en[state];
    }
  }
  return entropy;
}


const double ENTROPY_SAMPLE_COUNT = 10000;

double 
HMM::
entropy()
  const
{
  // stochatically estimate it 
  // generate and take avrage log p
  double entropy = 0;
  for (int i = 0; i < ENTROPY_SAMPLE_COUNT; i++){
    vector<int> * randomStringPtr = generate();
    //    dumpString(*randomStringPtr);
    double prob = probability(*randomStringPtr);
    //cerr << " prob = " << prob << endl;
    assert(prob > 0);
    entropy  += log(prob);
    delete randomStringPtr;
  }
  return -1.0l * entropy / double(ENTROPY_SAMPLE_COUNT);
}

//
// return a heap allocated vector
// generated randomly from the HMM
//

vector<int> * 
HMM::
generate() 
  const
{
  int currentState = initialState;
  vector<int> * answer = new vector<int>;
  while (currentState != finalState){
    // choose next one
    currentState = randomNextState(currentState);
    answer->push_back(randomOutput(currentState));
  }
  return answer;
}


int 
HMM::
randomNextState(int state) 
  const
{
  double r = double(rand())/double(RAND_MAX + 1.0l);
  for (int s = 0; s < numberStates; s++)
    {
      r -= p(state,s);
      if (r < 0)
	return s;
    }
  assert(false);
  return(finalState);
}

int 
HMM::
randomOutput(int state) 
  const
{
  double r = double(rand())/double(RAND_MAX + 1.0l);
  for (int s = 0; s < numberLetters; s++)
    {
      r -= q(state,s);
      if (r < 0)
	return s;
    }
  assert(false);
  return(terminator);
}

void 
HMM::
dumpString(const vector<int> & str) 
  const
{
  for (size_t i = 0; i < str.size(); i++)
    cout << " " << str[i];
  cout << endl;
}

/**
 * Do one E-step on the data vector
 * newHMM  is the accumulator for the posterior expectations
 */

double 
HMM::
emSingle(HMM & newHMM, 
	 double weight, 
	 const vector<int> & data) 
  const
{
  Matrix * alphaPtr = calculateForward(data);
  Matrix * betaPtr = calculateBackward(data);
  Matrix & alpha = *alphaPtr;
  Matrix & beta = *betaPtr;
  size_t l = data.size();
  double probability = beta(initialState,0);
  assert(probability < 1.0001l);
  if (probability > 0){
    for (size_t i = 0; i < l; i++){
      int u = data[i];
      // transitions from i to i+1
      for (int state = 0; state < numberStates; state++){
	for (int state2 = 0; state2 < numberStates; state2++){
	  double expectation = alpha(state,i) * p(state,state2) * q(state2,u) * beta(state2,i+1);
	  double increment = expectation * weight / probability;
	  newHMM.p(state,state2) += increment;
	  newHMM.q(state2,u) += increment;
	}
      }
    }
  }
  else {
    cerr << "Non-fatal error -- zero probability. factor = " << weight  << endl;
  }
  delete alphaPtr;
  delete betaPtr;
  return probability;
}

//
// return a new HMM based on a single iteration of the HMM
//

HMM* 
HMM::
em(const vector<const vector<int>*> & dataVector) 
  const
{
  HMM* newHMMPtr = new HMM(numberStates,numberLetters);
  size_t n = dataVector.size();
  double logProb = 0;
  for (size_t w = 0; w<n; w++){
    double prob = emSingle(*newHMMPtr, 1.0l, *dataVector[w]);
    // cout << prob << endl;
    if (prob > 0)
      logProb += log(prob);
  }
  newHMMPtr->normalise();
  cout << "LOGPROB = " << logProb << endl;
  return newHMMPtr;
}

//
// Train a new HMM on this weighted data
//

HMM* 
trainHMM(const vector< const vector<int> * > & dataVector, 
	 int numberStates,
	 int numberLetters)
{
  HMM* HMMPtr = new HMM(numberStates, numberLetters);
  cout << "Randomised" << endl;
  HMMPtr->randomise();
  HMMPtr->normalise();
  for (int iteration = 0; iteration < 10; iteration++){
    cout << "Iter " << iteration << endl;
    HMM* newHMMPtr = HMMPtr->em(dataVector);
    delete HMMPtr;
    HMMPtr = newHMMPtr;
  }
  return HMMPtr;
}



void 
HMM::
convertString(const string input, vector<int> & output)
{
  int l = input.size();
  output.resize(l+1);
  for (int i = 0; i < l; i++){
    int x = input[i];
    if (x < 0){
      x += 256;
    }
    output[i] = x;
  }
  output[l] = terminator;
}



void 
HMM::
store(ofstream & out) 
  const
{
  out << HMM_HEADER << "\n";
  out << numberStates << endl;
  out << numberLetters << endl;
  for (int i = 0; i < numberStates; i++){
    // transitions
    for (int j = 0; j < numberStates; j++){
      out << p(i,j) << endl;
    }
    for (int k = 0; k < numberLetters; k++){
      out << q(i,k) << endl;
    }
  }
  out << HMM_FOOTER << "\n";
}




HMM* 
loadHMM(ifstream & in)
{
  char buffer[BUFLEN];
  in.getline(buffer,BUFLEN);
  if (strcmp(buffer,HMM_HEADER))
    {
      cout << "Load error -hmm header" << endl;
      cout << strlen(buffer) << endl;
      cout << buffer <<  endl;
      cout << strlen(HMM_HEADER) << endl;
      cout << HMM_HEADER << endl;
      exit(-1);
    }
  int numberStates, numberLetters;
  in >> numberStates;
  in >> numberLetters;
  HMM* mPtr = new HMM(numberStates, numberLetters);
  for (int s = 0; s < numberStates; s++){
    for (int s2 = 0; s2 < numberStates; s2++)
      in >> mPtr->p(s,s2);
    for (int l = 0; l < numberLetters; l++)
      in >> mPtr->q(s,l);
  }
  in.ignore();
  in.getline(buffer,BUFLEN);
  if (strcmp(buffer,HMM_FOOTER))
    {
      cout << "Load error -hmm footer" << buffer <<  endl;
      exit(-1);
    }
  return mPtr;
}



