#include <math.h>
#include <iostream>
#include <limits.h>
#include <stdio.h>
#include <map>
#include <string>

/**
 * $Id: twohmm.cpp,v 1.6 2002/08/14 13:46:18 clark Exp $
 */

#include "twohmm.h"

using namespace std;


// constants

const char * TWOHMM_HEADER =  "ASC_TWOHMM_HEADER2";
const char * TWOHMM_FOOTER =  "ASC_TWOHMM_FOOTER2";
const char * SENTENCE_BOUNDARY = "";

const int BUFLEN = 1000;
const double EPSILON = 1.0e-50l;


// Some global parameters

// with 0 it should converge rapidly to a ML model 
double HMM_DECREMENT = 0;

// This is the fraction of a count we remove 
// from all of the ML parameters together; i.e. we remove ML_DECREMENT * counts[w] / numberStates


double TOKEN_ML_DECREMENT = 0.0l;
double TYPE_ML_DECREMENT = 0.0l;

/**
 * this gives the initial weight 
 * if it is 0.0l then no HMM
 * if it is 1.0l then only hMM
 */

double INITIAL_MIXTURE = 0.0l;


// this is the minimum frequency that 
// a word must have in order to be memorised

int MIN_FREQUENCY = 1;

//
//
//  A TwoHMM is a HMM where each output function is itself a HMM
// The TwoHMM has a distinguished start and end state
// which do not output anything
// mixing[i] is the proportion of the weight of state i that goes to
// the HMM model
// if mixing[i] is zero then outputHMMs[i] == 0
// 

/**
 * Destructor
 */

TwoHMM::
~TwoHMM()
{
  for (int i = 1; i < numberStates; i++)
    delete outputHMMs[i];
  // dont delete word array since it is owned by the corpus
}


void 
TwoHMM::
fillFromCorpus(const SimpleCorpusOne & corpus,
	       int f)
{
  int w = 0;
  for (int i = 0; i < corpus.numberTypes; i++){
    const string & word = corpus.getNthWord(i);
    if (corpus.countArray[i] >= MIN_FREQUENCY){
      wordArray[w] = &word;
      dictionary[word] = w;
      w++;
    }
  }
  assert(w == numberWords);
}

/**
 * Train a two level hmm on the data 
 * and return it 
 */

TwoHMM* 
trainTwoHMM(const SimpleCorpusOne & corpus,
	    int numberStates,
	    int numberSubStates,
	    int numberLetters)
{
  // first count the number of words with sufficient counts
  int numberWords = corpus.countWords(MIN_FREQUENCY);
  // first create a new one
  TwoHMM* thmmPtr = new TwoHMM(numberStates, numberSubStates, numberLetters, numberWords);
  thmmPtr->fillFromCorpus(corpus,MIN_FREQUENCY);
  thmmPtr->randomise();
  thmmPtr->normalise();
  for (int iteration = 0; iteration < MAX_ITERATIONS; iteration++){
    cerr << "Starting Iteration " << iteration << endl;
    TwoHMM* newPtr = thmmPtr->em(corpus);
    delete thmmPtr;
    thmmPtr = newPtr;
  }
  return thmmPtr;
}


/**
 * Train the model so as to maximise the likelihood of the 
 * corpus.  
 * We perform one iteration of the EM algorithm,
 * and return the new, correctly normalised model.
 * Do NOT assume that corpus.numberTypes == numberWords
 */

TwoHMM* 
TwoHMM::
em(const SimpleCorpusOne & corpus) 
  const
{
  TwoHMM* newPtr = new TwoHMM(numberStates, numberSubStates, numberLetters, numberWords);
  int numberTypes = corpus.numberTypes;
  // So numberTypes is the number of different words in the corpus
  // numberWords is the number memorised by the model
  for (int i = 0; i < numberWords; i++){
    const string * wordPtr = wordArray[i];
    newPtr->wordArray[i] = wordPtr;
    newPtr->dictionary[*wordPtr] = i;
  }
  int boundary = corpus.lookUpWord(SENTENCE_BOUNDARY);
  // first calculate the probs of each word
  // with respect to the submodels
  cerr << "Caching word probabilities." << endl;
  Matrix probabilities(numberStates, numberTypes);
  vector<int> wordIndices(numberTypes,-1);
  for (int w = 0; w < corpus.numberTypes; w++){
    const string & word = *(corpus.wordArray[w]);
    vector<int> wv;
    convertString(word,wv);
    int wordIndex = lookUpWord(word);
    wordIndices[w] = wordIndex;
    //if (wordIndex == -1){
    // cerr << "New word " << word << endl;
    //}
    //cerr << "word " << w << " length = " << corpus.wordArray[w]->size() <<  endl;
    double wp = 0.0l;
    for (int i = 0; i < numberStates; i++){
      if (i != initialState && i != finalState){
	// mixture of the HMM and the ML model
	double qq = (wordIndex == -1) ? 0 : q(i,wordIndex);
	//	cerr << "qq = " << qq << endl;
	double alpha = mixing[i];
	if (alpha > 0)
	  probabilities(i,w) = 
	    (alpha * outputHMMs[i]->probability(wv)) +
	    (1 - alpha) * qq;
	else
	  probabilities(i,w) = qq;
	wp += probabilities(i,w);
	if (probabilities(i,w) > 1.00001l){
	  cerr << "i = " << i << " w = " << w << " p = " << probabilities(i,w) << endl;
	  assert(false);
	}
      }
    }
    if (wp < 1.0e-100l && w != boundary){
      cerr << "Fatal error -- zero word probability for word " << w << " = |" << corpus.getNthWord(w) << "|" << endl;
      exit(-1);
    }

  }
  cerr << "Done." << endl;
  // this matri stores the posterior expectations that
  // each state generated each word
  Matrix posteriors(numberStates,corpus.numberTypes);
  // accumulate the transitions
  double lp = 0.0l;
  cerr << "Calculating the transitions." << endl;

  // find a sentence starting at index
  int index = -1;
  vector<int> sentence;
  while (index < corpus.numberTokens){
    int next = corpus.findNextOccurrence(SENTENCE_BOUNDARY, index);
    int size = next - (index+1);
    if (size > 0){
      sentence.resize(size);
      for (int i = 0; i < size; i++)
	sentence[i] = corpus.data[index + i + 1];
      lp += emSentence(newPtr, sentence, probabilities, posteriors);
    }
    index = next;
    //cerr << "Sentence " << i << " length " <<  sentence.size() << endl;
  }
  double psum = posteriors.sum();
  assert(finite(psum));
  cerr << "Agreggate posterior mass = " << psum << endl;
  //  posteriors.dump();
  cerr << "TOTAL = " << lp << endl;
  cerr << "Training submodels" << endl;
  for (int i = 0; i < numberStates; i++){
    // total probability mass of the state
    double total = 0.0l;
    double alpha = mixing[i];
    // total sums
    double mlWeight = 0.0l;
    double hmmWeight = 0.0l;
    vector<double> mlWeights(numberTypes,0.0l);
    vector<double> hmmWeights(numberTypes,0.0l);
    int numberZeroes = 0; // just for debugging
    for (int w = 0; w < numberTypes; w++){
      double piw = probabilities(i,w);
      double x = posteriors(i,w);
      total += x;
      int wordIndex = wordIndices[w]; 
      double qq = (wordIndex == -1) ? 0.0l : q(i,wordIndex);
      if (qq == 0)
	numberZeroes++;
      assert(x >= 0);
      // proportion that will have been ML
      if (x > EPSILON){
	assert(piw > 0);
	double MLProp =  (1 - alpha) * qq/ piw;
	//cerr << "MLPROP = " << MLProp << endl;
	assert(MLProp >= 0.0l && MLProp <= 1.0l);
	// y is the posterior mass that will made by the ML
	double y = x * MLProp;
	mlWeights[w] = y;
	mlWeight += y;
	hmmWeights[w] = x - y;    // this holds the expected values of the HMM model
	hmmWeight += x - y;
      }
    }
    // now we have the posteriors
    // we adjust them using the negative Dirichlet priors
    // the same for each of the words
    // these two are just for debugging
    double totalDecrement = 0.0l;
    for (int w = 0; w < numberTypes; w++){
      // decrement is the amount we reduce the ML by 
      double decrement = (TOKEN_ML_DECREMENT * corpus.countArray[w] + TYPE_ML_DECREMENT) / numberStates;
      // qq is the adjusted weight for the ML model
      double MLW = mlWeights[w];
      if (MLW < decrement){
	decrement = MLW;
      }
      totalDecrement += decrement;
      double newQ = MLW - decrement;
      mlWeight -= decrement;
      if (newQ > 0){
	int wordIndex = wordIndices[w];
	assert(wordIndex > -1);
	newPtr->q(i,wordIndex) = newQ;
      }
    }

    cerr << "State " << i 
	 << ": Total mass = " << total 
	 << ", total dec. = " << totalDecrement 
	 << ", number zeroes = " << numberZeroes << endl;
    // this might break when we have full blown parameter extinction
    double newmix  = hmmWeight/(hmmWeight + mlWeight);
    if (i == initialState || i == finalState){
      cerr << "Zeroing initial state mixtures.\n";
      newmix = 0.0l;
    }
    newPtr->mixing[i] = newmix;
    // now we train the model
    cerr << "State " << i << " mix  = " << newmix << endl;
    if (newmix > 0){
      HMM* newHMMPtr = newPtr->outputHMMs[i];
      HMM* oldHMMPtr = outputHMMs[i];
      for (int w = 0; w < numberTypes; w++){
	double ww = hmmWeights[w];
	if (ww > EPSILON){
	  vector<int> wv;
	  convertString(*(corpus.wordArray[w]),wv);
	  oldHMMPtr->emSingle(*newHMMPtr, ww, wv);
	}
      }
    }
  }
  newPtr->normalise();
  return newPtr;
}

/**
 * do the EM algorithm (E-step)  on a single sentence
 * The sentence consists of indices into corpus features
 * 
 */

double 
TwoHMM::
emSentence(TwoHMM* newPtr, 
	   const vector<int> & sentence,
	   const Matrix & probabilities, 
	   Matrix & posteriors) 
  const
{
  int l = sentence.size();
  // we want a separate begin and end state that are unrelated
  Matrix alpha(numberStates, l+2);
  Matrix beta(numberStates, l+2);
  double logProb  = calculateScaledProbabilities(sentence, probabilities, alpha, beta);
  assert(logProb <= 0);
  //cerr << " logprob " << logProb << endl;
  //cerr << "Posteriors sum = " << posteriors.sum() <<  " next length = " << l << endl;
  for (int i = 0; i < l; i++){
    int w = sentence[i];
    // cerr << " I = " << i << ", w = " << w << endl;
    // transitions from i to i+1
    double sum = 0.0l;
    for (int state = 0; state < numberStates; state++){
      for (int state2 = 1; state2 < numberStates; state2++){
	//cerr << "alpha " <<  alpha(state,i) << ","  
	//     << p(state,state2) << "," <<  probabilities(state2,w)  << "," <<  beta(state2,i+1) << endl;
	double expectation = alpha(state,i) * p(state,state2) * probabilities(state2,w) * beta(state2,i+1);
	// since they are scaled we don't need to divide by the prob.
	sum += expectation;
	assert(expectation >= 0);
	newPtr->p(state,state2) += expectation;
	//if (expectation > 0 && state == 2)
	//cerr << " p " << state << "," << state2 << " += " << expectation << " i = " << i << " l = " << l << endl;
	posteriors(state2,w) += expectation;
      }
    }
    //cerr << "SUM " << sum << endl;
    if (sum == 0){
      cerr << "w = " << w << " i " << i << ", l = " << l << endl;
    }
    assert(sum > 0);
    // cerr << "Added " << sum << " to word " << w << endl;
  }
  // also do the final one from l to l+1
  for (int state = 0; state < numberStates; state++){
    //if (state == 2)
    //cerr << "2 -> final " << alpha(state,l) * p(state,finalState) * beta(finalState,l+1) << endl;;
    newPtr->p(state,finalState) += alpha(state,l) * p(state,finalState) * beta(finalState,l+1);
  }
  return logProb;
}


//
// alpha(state,i)
//   is the probability that it will have output the first i symbols ( 0, ..., i-1)
//   and ended up in  state 
//   so alpha(initialState,0) = 1
//      alpha(finalState, l+1) = total prob where we ignore the symbol at l
//   
//
// beta(state,i)
//   is the probability that starting from state  at stage i 
//   it will output symbols i, i+1 , ..., l-1 and then end up in final state
//   so beta(finalState,l+1) = 1
//
// These are scaled so that the sum of all the transitions at step i is unity
// The function returns the log probability of the observation sequence
// 
// sum_state alpha(state,i) * beta(state,i) =  1 
//


// For now divide all of the alphas by the total prob
// check to see if that works correctly

double 
TwoHMM::
calculateScaledProbabilities(const vector<int> & sentence, 
			     const Matrix & probabilities,
			     Matrix & alpha,
			     Matrix & beta) 
  const
{
  int l = sentence.size();
  assert(beta.dim1() == numberStates);
  assert(beta.dim2() == l+2);
  assert(alpha.dim1() == numberStates);
  assert(alpha.dim2() == l+2);
  // the first word is produced by the state in slot 1
  alpha(initialState,0) = 1.0l;

  // keep the array of scalars
  vector<double> alphaScales(l+2,1.0l);
  for (int i = 0; i < l; i++){
    // calculate i+1 in terms of i
    double sum = 0.0l;
    int word = sentence[i];
    //cerr << "Position " << i << " word " << word << endl;
    for (int state = 0; state < numberStates; state++){
      for (int state2 = 0; state2 < numberStates; state2++){
	//	cerr << "state " << state << " state2 " << state2 << endl;

	assert(word >= 0 && word < probabilities.dim2());
	double inc = alpha(state,i) * p(state,state2) * probabilities(state2, word);
	//if (inc > 0)
	// cerr << state << "," << state2 << " = " << alpha(state,i) 
	//     << "," << p(state,state2) << "," << probabilities(state2,word) << endl;
	alpha(state2,i+1) += inc;
	sum += inc;
	//assert(p(state,state2) <= 1.0l);
	//assert(q(state2,letter) <= 1.0l);
	//assert(alpha(state,i) <= 1.0001l);
      }
    }
    // calculate scaling factor
    if (sum <= 0 || sum > 1.0l){
      cerr << " SUM " << sum << endl;
    }
    assert(sum > 0);
    assert(sum < 1.00001l);
    alphaScales[i+1] = sum;
    for (int state = 0; state < numberStates; state++)
      alpha(state,i+1) /= sum;
  }
  // now calculate the alphas for i = l;
  for (int state = 0; state < numberStates; state++){
    //cerr << "finally, state = " << state << " alpha " << alpha(state,l) << " p(s,f) = " << p(state,finalState) << endl;
    alpha(finalState,l+1) += alpha(state,l) * p(state,finalState);
  }

  alphaScales[l+1] = alpha(finalState,l+1);
  assert(alphaScales[l+1]>0);
  alpha(finalState,l+1) = 1.0l;
  // the total prob is the product of all of the alphaScales

  //
  // Now calculate the backward probabilities
  //

  //  beta(finalState,l+1) = 1.0l;
  beta(finalState,l+1) = 1.0l/alphaScales[l+1];
  for (int state = 0; state < numberStates; state++)
    beta(state,l) = p(state,finalState) * beta(finalState,l+1) / alphaScales[l];
  for (int i = l-1; i >= 0; i--){
    // calculate beta(*,i) in terms of beta(*,i+1)
    double sum = 0.0l;
    int word = sentence[i];
    //    cerr << "i = " << i << " w " << word << endl;
    // calculate the values for i in terms of i+1
    for (int state = 0; state < numberStates; state++){
      for (int state2 = 0; state2 < numberStates; state2++){
	//  cerr << "state " << state
	//  	     << " state2 " << state2
	//  	     << " p " << p(state,state2)
	//  	     << " pq " << probabilities(state2, word)
	//  	     << " beta " << beta(state2,i+1)
	//  	     << " alphascales " << alphaScales[i+1] << endl;
 	double inc = p(state,state2) * probabilities(state2, word) * beta(state2,i+1) / alphaScales[i];
	beta(state,i) += inc;
	assert(finite(beta(state,i)));
	sum += inc;
      }
    }
    assert(sum > 0);
  }
  // Calculate the log probability and return it
  double lp = 0.0l;
  for (int i = 0; i < l+2; i++)
    lp += log(alphaScales[i]);
  return lp;
}

void 
TwoHMM::
normaliseStateTransitions(int state)
{
  assert(state >= 0 && state < numberStates);
  double sum = 0.0l;
  for (int i = 0; i < numberStates; i++){
    sum += p(state,i);
  }
  if (sum == 0.0l){
    cerr << "Non-fatal error state failed to normalise" << endl;
  }
  else 
    for (int i = 0; i < numberStates; i++)
      p(state,i) /= sum;
}


void 
TwoHMM::
store(const char * filename) 
  const
{
  ofstream out(filename);
  if (!out){
    cerr << "Couldnt open output file " << filename << endl;
    exit(-1);
  }
  store(out);
}

void 
TwoHMM::
store(ofstream & out) 
  const
{
  out << TWOHMM_HEADER << "\n";
  out << numberStates << endl;
  out << numberSubStates << endl;
  out << numberLetters << endl;
  out << numberWords << endl;
  for (int i = 0 ; i < numberWords; i++){
    out << *(wordArray[i]) << endl;
  }
  for (int i = 0; i < numberStates; i++){
    // transitions
    out << mixing[i] << endl;
    for (int j = 0; j < numberStates; j++){
      out << p(i,j) << endl;
    }
    for (int w = 0; w < numberWords; w++){
      out << q(i,w) << endl;
    }
    if (i > 0 && mixing[i] > 0){
      assert(outputHMMs[i]);
      outputHMMs[i]->store(out);
    }
  }
  out << TWOHMM_FOOTER << "\n";
}

// basic constructor

TwoHMM::
TwoHMM(int numberStates_, 
       int numberSubStates_, 
       int numberLetters_, 
       int numberWords_)
  : 
  numberStates(numberStates_), 
  numberSubStates(numberSubStates_), 
  numberLetters(numberLetters_),
  numberWords(numberWords_),
  p(numberStates,numberStates),
  q(numberStates,numberWords)
{
  mixing.resize(numberStates);
  outputHMMs.resize(numberStates);
  for (int s = 0; s < numberStates; s++){
    if (s != finalState && s != initialState &&  numberSubStates > 0)
      outputHMMs[s] = new HMM(numberSubStates, numberLetters);
    else 
      outputHMMs[s] = 0;
  }
  wordArray.resize(numberWords);
  for (int i = 0; i < numberWords; i++)
    wordArray[i] = 0;
}
 
void 
TwoHMM::
fillFromCorpus(const SimpleCorpusOne & corpus)
{
  assert(numberWords == corpus.numberTypes);
  for (int i = 0; i < numberWords; i++){
    const string * wordPtr = corpus.wordArray[i];
    wordArray[i] = wordPtr;
    dictionary[*wordPtr] = i;
  }
}

void 
TwoHMM::
randomise()
{
  cerr << "random initialisation of two hmm" << endl;
  p(finalState,finalState) = 1;
  for (int state = 0; state < numberStates; state++){
    if (state != finalState){
      for (int state2 = 0; state2 < numberStates; state2++)
	if (state2 != initialState) 
	  p(state,state2) = 1 + double(rand());
    }
  }
  for (int state = 0; state < numberStates; state++){
    if (outputHMMs[state]){
      cerr << "random initialisation of hmm from " << state << endl;
      outputHMMs[state]->randomise();
    }
    if (state != initialState && state != finalState){
      cerr << "random initialisation of ml outputs from " << state << endl;
      for (int w = 0; w < numberWords; w++)
	q(state,w) = 1 + double(rand());
      mixing[state] = INITIAL_MIXTURE;
    }
  }
}

void 
TwoHMM::
normalise()
{
  // normalise outputs
  for (int state = 0; state < numberStates; state++)
    if (outputHMMs[state])
      outputHMMs[state]->normalise();
  // normalise transitions
  for (int state = 0; state < numberStates; state++){
    if (state != finalState){
      double sum = 0.0l;
      for (int state2 = 0; state2 < numberStates; state2++)
	sum += p(state,state2);
      if (sum > 0)
	for (int state2 = 0; state2 < numberStates; state2++)
	  p(state,state2) /= sum;
      else
	cerr << "Non-fatal error -- transition normalisation error state " << state << endl;
    }
    if (state != finalState && state != initialState){
      double sumq = 0.0l;
      for (int w = 0; w < numberWords; w++){
	sumq += q(state,w);
      }
      if (sumq > 0){
	for (int w = 0; w < numberWords; w++)
	  q(state,w) /= sumq;
      }
      else if (mixing[state] < 1){
	cerr << "Non-Fatal error -- unnormalised ML model with mixing coefficient < 1, state = " 
	     << state 
	     << " mixing = " 
	     << mixing[state] << endl;
	//q.dump();
	exit(-1);
      }
    }
  }
}

/**
 * Generate at random
 * and output to cerr
 */

void 
TwoHMM::
generate(const SimpleCorpusOne & corpus) 
  const
{
  int currentState = randomNextState(initialState);
  while (currentState != finalState){
    //cerr << "STATE " << currentState << endl;
    // choose the mixing thing
    double alpha = mixing[currentState];
    //cerr << "alpha " << alpha << endl;
    double r = double(rand())/double(RAND_MAX + 1.0l);
    if (r < alpha){
      for (int x = 0; x < 10; x++){
	assert(outputHMMs[currentState]->validate());
	vector<int>* wPtr = outputHMMs[currentState]->generate();
	for (size_t i = 0; i < wPtr->size(); i++){
	  int c = (*wPtr)[i];
	  if (c)
	    cerr << convertToChar(c);
	}
	delete wPtr;
	cerr << " HMM " ; 
      }
    }
    else {
      // ML model.

      int w = randomQ(currentState);
      //      for (int i = 0; i < numberStates; i++)
      //	cerr << i << ": " << q(i,w) << ", " << endl;

      //assert(w > 0);
      assert(w < numberWords);
      cerr <<  *(corpus.wordArray[w]);
      cerr << " ML ";
    }
    cerr << " " <<  currentState << endl;
    currentState = randomNextState(currentState);
  }
}


/**
 * Generate from the second-level HMM
 */

void 
TwoHMM::
generate2(const SimpleCorpusOne & corpus, int n) 
  const
{
  for (int currentState = 0; currentState < numberStates; currentState++){
    if (currentState != initialState && currentState != finalState){
      cerr << "STATE " << currentState << endl;
      for (int w = 0; w < numberWords; w++){
	if (q(currentState,w) > 0){
	  cerr << condQ(currentState,w) << " ";
	  cerr <<  *(corpus.wordArray[w]);
	  cerr << endl;
	}
      }
      for (int i = 0; i < n; i++){
	vector<int>* wPtr = outputHMMs[currentState]->generate();
	for (size_t i = 0; i < wPtr->size(); i++){
	  int c = (*wPtr)[i];
	  if (c)
	    cerr << convertToChar(c);
	}
	cerr << endl;
      }
    }
  }
}

int 
TwoHMM::
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
TwoHMM::
randomQ(int state) 
  const
{
  double r = double(rand())/double(RAND_MAX + 1.0l);
  //  cerr << "r = " << r << endl;
  for (int w = 0; w < numberWords; w++)
    {
      r -= q(state,w);
      if (r < 0)
	return w;
    }
  assert(false);
  return(0);
}

TwoHMM* 
loadTwoHMM(const char * filename)
{
  ifstream in(filename);
  if (!in){
    cerr << "fatal error couldn't open " << filename << endl;
    exit(-1);
  }
  char buffer[BUFLEN];
  in.getline(buffer,BUFLEN);
  cerr << buffer << endl;
  int numberStates,numberSubStates, numberLetters, numberWords;
  in >> numberStates;  
  in >> numberSubStates;
  in >> numberLetters;
  in >> numberWords;
  in.ignore();
  cerr << "Loading model " << numberStates 
       << " , " << numberSubStates 
       << " , " << numberLetters 
       << " , " << numberWords<< endl;
  TwoHMM * mmPtr = new TwoHMM(numberStates,numberSubStates,numberLetters,numberWords);
  cerr << "allocated ok" << endl;

  cerr << "Loading vocab of " << numberWords << endl;
  for (int i = 0; i < numberWords; i++){
    in.getline(buffer,BUFLEN);
    string * sPtr = new string(buffer);
    mmPtr->wordArray[i] = sPtr;;
    mmPtr->dictionary[*sPtr] = i;
  }
  for (int state = 0; state < numberStates; state++){
    in >> mmPtr->mixing[state];
    if (!in.good()){
      cerr << "load error\n state = " << state <<endl ;
      exit(-1);
    }
    for (int s2 =0; s2 < numberStates; s2++){
      in >> mmPtr->p(state,s2);
      //	cerr << "(" << state << "," << s2 << ")= " << mmPtr->p(state,s2) << endl;
      if (!in.good()){
	cerr << "load error p\n state = " << state  <<  "state 2 = " << s2 <<endl ;
	exit(-1);
      }
    }
    for (int w = 0; w < numberWords; w++){
      in >> mmPtr->q(state,w);
      if (!in.good()){
	cerr << "load error q\n state = " << state  <<  " word = " << w <<endl ;
	exit(-1);
      }
    }
    if (state > 0 &&  mmPtr->mixing[state] > 0){
      // load HMM
      in.ignore();
      mmPtr->outputHMMs[state] = loadHMM(in);
    }
  }
  in.getline(buffer,BUFLEN);
  if (buffer[0] == 0){
    in.getline(buffer,BUFLEN);
  }
  if (strcmp(buffer,TWOHMM_FOOTER)){
    cerr << "FILE FORMAT ERROR:" << endl;
    cerr << "FILE END was:" << buffer << endl;
    cerr << "Should be:" << TWOHMM_FOOTER << endl;
    exit(-1);
  }
  return mmPtr;    
}

double 
TwoHMM::
condQ(int j, int w) 
  const
{
  double sum = 0.0l;
  for (int i = 0; i < numberStates; i++)
    sum += q(i,w);
  if (sum > 0)
    return q(j,w)/sum;
  else
    return 0.0l;
}

/**
 * Use the Viterbi algorithm to tag a sentence.
 * Fill the tag vector with the relevant tags
 */

void 
TwoHMM::
tagTestSentence(const vector<string *> & sentence, 
		vector<int> & tags) 
  const
{ 
  int l = sentence.size();
  tags.resize(l);
  // probs(s,i) stores the prob of the most likely state sequence 
  // that generates the first i of them and ends in 
  // state s
  Matrix probs(numberStates,l+1);
  // previousState stores where the state sequence came from
  MatrixInt previousState(numberStates,l+1);
  probs(initialState,0) = 1.0l;
  for (int i = 0; i < l; i++){
    //cerr << "Word " << i << " of " << l << endl;
    const string & w = *(sentence[i]);
    vector<int> wv;
    convertString(w,wv);
    int wordIndex = lookUpWord(w);
    //if (wordIndex == -1){
    // cerr << "New word " << w << endl;
    //}
    // now calculate the output probabilities
    vector<double> op(numberStates,0.0l);
    for (int s = 0; s < numberStates; s++){
      double alpha = mixing[s];
      // cerr << "alpha " << alpha << endl;
      if (alpha > 0)
	op[s] = alpha * outputHMMs[s]->probability(wv);
      if (wordIndex > -1)
	op[s] += (1 - alpha) * q(s,wordIndex);
      assert(finite(op[s]));
    }
    // scale them now to unity
    double sum = 0.0l;
    
    for (int s = 0; s < numberStates; s++)
      sum += op[s];
    assert(finite(sum));
    if (sum == 0.0l){
      cerr << "Fatal error: this word is not generated by the model." << endl;
      exit(-1);
    }
    //    cerr << "Total output probability is " << sum << endl;
    for (int s = 0; s < numberStates; s++)
      op[s] /= sum;
    // now work out the best path
    // this is the sumo f the column of the trellis, which we just scale by
    double sum2 = 0.0l;
    for (int s = 0; s < numberStates; s++){
      // best path ending in state s
      int bestState = -1;
      double bestScore = -1.0l;
      for (int ps = 0; ps < numberStates; ps++){
	double score = probs(ps,i) * p(ps,s) * op[s];
	if (score > bestScore){
	  bestScore = score;
	  bestState = ps;
	}
      }
      assert(bestState>-1);
      //cerr << "(s, i+1) " << s << "," << i+1 << " (score,state) = " << bestScore << "," << bestState << endl;
      probs(s,i+1) = bestScore;
      sum2 += bestScore;
      previousState(s,i+1) = bestState;
    }
    for (int s = 0; s < numberStates; s++){
      probs(s,i+1) /= sum2;
    }
  }
  int bestState = -1;
  double bestScore = -1.0l;
  for (int state = 0; state < numberStates; state++){
    double score = probs(state,l) * p(state,finalState);
    if (score > bestScore){
      bestScore = score;
      bestState = state;
    }
  }
  assert(bestState > -1);
  assert(bestScore > 0);
  // now trace back.
  int currentState = bestState;
  for (int i = l; i > 0; i--){
    tags[i-1] = currentState;
    assert(probs(currentState,i) > 0);
    currentState = previousState(currentState,i);
  }
  assert(currentState == initialState);
}

void 
TwoHMM::
dump() 
  const
{
  cout << "P\n";
  p.dumpNonZero();
  cout << "Q\n";
  q.dumpNonZero();
}

// train it using an initial hard clustering which you read from
// a clustering file

TwoHMM* trainTwoHMMInitialHard(const SimpleCorpusOne & corpus, 
			       const char * clusteringFile,
			       int numberStates,
			       int numberSubStates,
			       int numberLetters,
			       int numberWords)
{
  // first create a new one
  TwoHMM* thmmPtr = new TwoHMM(numberStates, numberSubStates, numberLetters + 1, numberWords);
  thmmPtr->fillFromCorpus(corpus);
  thmmPtr->initialiseFromClusters(corpus,clusteringFile);
  thmmPtr->perturbQ();
  thmmPtr->normalise();
  for (int iteration = 0; iteration < MAX_ITERATIONS; iteration++){
    cout << "Starting Iteration " << iteration << endl;
    TwoHMM* newPtr = thmmPtr->em(corpus);
    delete thmmPtr;
    thmmPtr = newPtr;
  }
  return thmmPtr;
  
}

void TwoHMM::initialiseFromClusters(const SimpleCorpusOne & corpus, const char * clusteringFile)
{
  ifstream in(clusteringFile);
  if (!in){
    cerr << "Fatal error -- couldn't open clustering file " << clusteringFile << endl;
    exit(-1);
  }
  map<string,int> tagDict;
  int currentTag = 0;
  int lineNum = 0;
  while (!in.eof()){
    char buffer[BUFLEN];
    char wordbuffer[BUFLEN];
    char clusterbuffer[BUFLEN];
    in.getline(buffer,BUFLEN);
    lineNum++;
    float param = 0.0l;
    int matches = sscanf(buffer,"%s %s %f",wordbuffer, clusterbuffer, &param);
    if (matches > 1){
      cout << "Word " << wordbuffer << " cluster " << clusterbuffer << endl;
      // look up tag 
      string tag(clusterbuffer);
      if (tagDict.find(tag) == tagDict.end()){
	// a new tag
	while (currentTag == initialState || currentTag == finalState)
	  currentTag++;
	if (currentTag >= numberStates){
	  cerr << "Too many clusters in cluster file = " << currentTag << endl;
	  exit(-1);
	}
	tagDict[tag] = currentTag;
	currentTag++;
      }
      // this is inefficient 
      int thisTag = tagDict[tag];
      int thisWord = corpus.lookUpWord(string(wordbuffer));
      if (thisWord == -1){
	cerr << "Non-fatal error : unknown word in cluster file " << wordbuffer << endl;
      }
      q(thisTag,thisWord) = corpus.countArray[thisWord];
    }
    else {
      cerr << "Non-fatal error : couldn't scan linenumber " << lineNum << endl;
    }
  }
  cout << "Initialising uniform transitions" << endl;
  p(finalState,finalState) = 1;
  for (int state = 0; state < numberStates; state++){
    if (state != finalState){
      for (int state2 = 0; state2 < numberStates; state2++)
	if (state2 != initialState) 
	  p(state,state2) = 1.0l;
    }
  }
  cout << "Randomising HMM submodels" << endl;
  for (int state = 0; state < numberStates; state++){
    if (state != initialState && state != finalState){
      outputHMMs[state]->randomise();
      mixing[state] = INITIAL_MIXTURE;
    }
  }
}

//
// Make all of the parameters non-zero
//

void 
TwoHMM::
perturbQ()
{
  for (int i = 0; i < numberStates; i++)
    if (i != initialState && i != finalState)
      for (int w = 0; w < numberWords; w++)
	q(i,w) += 0.01l / double(numberWords);
}

void 
TwoHMM::
convertString(const string & input, 
	      vector<int> & output) 
  const
{
  return HMM::convertString(input,output);
}



/**
 * Return the index or -1 if it doesn't exist
 */

int 
TwoHMM::
lookUpWord(const string & word) 
  const
{
  map<string,int>::const_iterator pos = dictionary.find(word);
  if (pos != dictionary.end())
    return pos->second;
  else
    return -1;
}
