#include "clusters.h"
#include <algorithm>
#include <utility>
#include <iostream>
#include <math.h>

/**
 * $Id: clusters.cpp,v 1.5 2002/08/14 13:46:18 clark Exp $
 */

int FREQ_CUTOFF = -1;
int MAX_ITERATIONS = 20;
bool USE_TRUE_WEIGHT = false;
double PRIOR_BOOST = 0;
bool FULL_MORPHOLOGY_WEIGHT = false;

Clusters::
Clusters(int numberClasses_, 
	 const SimpleCorpusOne & corpus_,
	 int numberStates_,
	 int alphabetSize_,
	 bool randomised)
  :
  numberClasses(numberClasses_), 
  numberTypes(corpus_.numberTypes), 
  numberTokens(corpus_.numberTokens),
  numberStates(numberStates_),
  alphabetSize(alphabetSize_),
  data(corpus_.data),
  corpus(corpus_),
  clusterBigrams(numberClasses_,numberClasses_)
{
  classVector.resize(numberTypes);
  counts.resize(numberTypes); 
  sortedWords.resize(numberTypes);
  first.resize(numberTypes);
  clusterUnigrams.resize(numberClasses);
  next = new int[numberTokens];
  for (int i = 0; i < numberTokens; i++)
    next[i] = numberTokens;
  for (int w = 0; w < numberTypes; w++){
    counts[w]=0;
    classVector[w] = numberClasses -1;
  }
  // counts are set
  for (int i = 0; i < numberTokens; i++)
    counts[data[i]]++;
  // now find the most frequent numberClasses -1 of them.
  vector< pair<int,int> > countsTable(numberTypes);
  for (int i = 0; i < numberTypes; i++){
    countsTable[i] = pair<int,int>(counts[i],i);
    //cerr << counts[i] << " " << i << endl;
  }
  
  cerr << "Sorting words" << endl;
  sort(countsTable.begin(),countsTable.end());

  for (int i = 0; i < numberTypes; i++){
    first[i] = -1;
    sortedWords[i] = countsTable[numberTypes - 1 - i].second;
    //cerr << "sort " << i << " " << sortedWords[i] << " , n =" << countsTable[numberTypes - 1 - i].first << endl;
  }

  if (randomised)
    {
      for (int i = 0; i < numberTypes; i++){
	if (counts[i] > FREQ_CUTOFF){
	  int rc = (int) (1.0 * numberClasses *rand()/(RAND_MAX+1.0));
	  classVector[i] = rc;
	}
      }
    }
  else {
    for (int i = 0; i < numberClasses-1; i++){
      classVector[sortedWords[i]]= i;
    }
  }
  
  vector<int> last(numberTypes,0);
  cerr << "Indexing data" << endl;
  for (int i = 0; i < numberTokens-1; i++){
    int w = data[i];
    int w2 = data[i+1];
    if (w2 < 0 || w2 > numberTypes -1){
      cerr << i+1 << " " << w2 << endl;
    }
    assert(w >= 0 && w < numberTypes);
    assert(w2 >= 0 && w2 < numberTypes);
    if (first[w] == -1){
      first[w] = i;
      last[w] = i;
    }
    else
      {
	next[last[w]] = i;
	last[w] = i;
      }
    int c1 = classVector[w];
    int c2 = classVector[w2];
    assert(c1 >= 0 && c1 < numberClasses);
    assert(c2 >= 0 && c2 < numberClasses);
    clusterBigrams(c1,c2)++;
    clusterUnigrams[c1]++;
  }
  cerr << "Finished indexing " << endl;

  // be careful
  clusterUnigrams[classVector[data[numberTokens-1]]]++;
  cerr << "Numberstates " << numberStates << endl;
  if (numberStates > 0){
    cerr << "Starting to do the HMMs" << endl;
    hmms.resize(numberClasses);
    for (int i = 0; i < numberClasses; i++){
      HMM* hmmPtr = new HMM(numberStates, alphabetSize);
      hmmPtr->randomise();
      hmmPtr->normalise();
      hmms[i] = hmmPtr;
    }
  }
}




void 
Clusters::
dump(const SimpleCorpusOne & corpus) 
  const
{
  // dump it as word cluster p(word| cluster)
  vector<double> clusterCounts(numberClasses,0.0l);
  for (int i = 0; i < numberTypes; i++)
    clusterCounts[classVector[i]] += counts[i];
  for (int i = 0; i < numberClasses; i++)
    cerr << i << " " << clusterCounts[i] << endl;
  for (int i = 0; i < numberTypes; i++){
    //  if (classVector[i] < numberClasses-1)
    int w = sortedWords[i];
    int cl = classVector[w];
    cout << *(corpus.wordArray[w]) << " " 
	 << cl << " " 
	 << double(counts[w])/clusterCounts[cl] << endl;
  }
}


void 
Clusters::
clusterNeyEssen()
{
  int i = 0;
  int c = 0;
  validateBigramCounts();
  while ((c = reclusterNeyEssen()) > 0 && i < MAX_ITERATIONS){
    i++;
    cerr << "finished iter " << i << ", changed " << c <<  endl;
    validateBigramCounts();
  } 
}

int 
Clusters::
reclusterNeyEssen()
{
  cerr << "Calculating MLE for prior probabilities" << endl;
  vector<double> prior(numberClasses,0.0l);
  double xxx = 1.0l/(double)numberTypes;
  for (int i = 0; i < numberTypes; i++){
    int c = classVector[i];
    prior[c] += xxx;
  }
  for (int i = 0; i < numberClasses;i++)
    cerr << i << " " << prior[i] << endl;
  if (numberStates > 0){
    cerr << "Training all the HMMs" << endl;
    for (int c = 0; c < numberClasses; c++){
      //      cerr << "Training HMM " << c << endl;
      HMM* hmmPtr = hmms[c];
      HMM* newHmmPtr = new HMM(numberStates, alphabetSize);
      for (int i = 0; i < numberTypes; i++){
	if (classVector[i] == c){
	  // then this word is in the right class
	  // so train it on word i
	  const string & word = *(corpus.wordArray[i]);
	  vector<int> v;
	  hmmPtr->convertString(word,v);
	  // FIXME 
	  double weight = 1.0l;
	  if (USE_TRUE_WEIGHT){
	    weight = corpus.countArray[i];
	  }
	  hmmPtr->emSingle(*newHmmPtr, weight, v);
	}
      }
      newHmmPtr->normalise();
      hmms[c] = newHmmPtr;
      delete hmmPtr;
    }
  }
  int something = 0;
  for (int i = 0; i < numberTypes; i++){
    //    cerr << "Word " << i;
    int w = sortedWords[i];
    //cerr << *(corpus.wordArray[w]) << endl;
    if (counts[w] > FREQ_CUTOFF){
      //cerr << "Doing " << w << endl;
      if (bestCluster(w, prior)){
	something++;
      }
    }
  }
 
  return something;
}


// try to swap word w 
// return best cluster

bool 
Clusters::
bestCluster(int w,
	    const vector<double> & prior)
{
  //  cerr << "Starting to work on word " << w << endl;
  double score = 0;
  int best = classVector[w];
  // number of  occurrences of cluster x immediately before word w
  vector<int> left(numberTypes,0);
  // number of  occurrences of cluster x immediately after word w
  vector<int> right(numberTypes,0);
  int currentPos = first[w];
  // number of occurrences of w w in the corpus
  int doubles = 0;
  // assert it occurs in the corpus
  if (currentPos >= 0){
    while (currentPos < numberTokens){
      if (currentPos > 0)
	// do left
	left[classVector[data[currentPos-1]]]++;
      if (currentPos < numberTokens - 1)
	right[classVector[data[currentPos+1]]]++;
      if (currentPos+1 == next[currentPos]){
	doubles++;
	//      cerr << "Double word!" << endl;
      }
      currentPos = next[currentPos];
    }
  }
  // it might be better to preconvert them and store them in the corpus
  vector<int> ws;
  double oldStringProb = 0.0l;
  if (numberStates > 0){
    hmms[0]->convertString(*(corpus.wordArray[w]), ws);
    oldStringProb = hmms[classVector[w]]->probability(ws); 
    if (FULL_MORPHOLOGY_WEIGHT){
      oldStringProb *= counts[w] + 1;
    } 
  }
  for (int i = 0; i < numberClasses; i++){
    //cerr << "Calculating word " << w << " from class " << classVector[w] << " to " << i << endl;  
    // convert it
    double newScore = calcChangeFast(w,
				     oldStringProb,
				     ws,
				     i,left,right,doubles,prior);
    //    double fullScore = calcChange(w,i,left,right,doubles);
    //cerr << newScore << " " << fullScore << endl;
    //if ((newScore * fullScore) < 0)
    // exit(-1);
    if (newScore > score){
      score = newScore;
      best = i;
      //exit(-1);
    }
  }
  int old = classVector[w];
  if (old != best){ 
    moveWord(w,old,best,doubles,left,right);
    return true;
  }
  else {
    //cerr << "leaving " << w << " in class " << old << endl;
    return false;
  }
}

// calculate the score of changing w to newCluster
// should be zero if it is newCluster
// bigger is better

// this is quadratic in numberClasses
// we can optimise this to linear
// by only doing it if we have one ofthem equal to new or old
// We can precompute the before, and cache it, which should double it.


double 
Clusters::
calcChange(int w, 
	   int newCluster, 
	   const vector<int> & left, 
	   const vector<int> & right,
	   int doubles) 
const
{
  int oldCluster = classVector[w];
  // do this very inefficiently
  // calculate the before and the after
  double before = 0.0l;
  int c = counts[w];
  //  cerr << "w " << w << " old " << oldCluster << " new " << newCluster << endl;
  if (c == clusterUnigrams[oldCluster]){
    //    cerr << "ignoring singleton cluster " << endl;
    return 0;
  }
  for (int g1 = 0; g1 < numberClasses; g1++)
    for (int g2 = 0; g2 < numberClasses; g2++){
      if (g1 == newCluster || g1 == oldCluster
	  || g2 == newCluster || g2 == oldCluster)
	{
    
	  double bg = clusterBigrams(g1,g2);
	  //cerr << bg << " " ;
	  if (bg > 0){
	    double c1 = clusterUnigrams[g1];
	    double c2 = clusterUnigrams[g2];
	    assert(c1 > 0);
	    assert(c2 > 0);
	    double factor = bg * log(bg / (c1 * c2));
	    before +=  factor;
	    //	    cerr << g1 << " " <<  g2 << " f = " << factor << endl;
	    // assert(finite(before));
	  }
	}
    }
  //assert(finite(before));
  double after = 0.0l;
  for (int g1 = 0; g1 < numberClasses; g1++){
    double cg1 = clusterUnigrams[g1];
    if (g1 == oldCluster)
      cg1 -= c;
    if (g1 == newCluster)
      cg1 += c;
    assert(cg1 > 0);
    for (int g2 = 0; g2 < numberClasses; g2++){
      if (g1 == newCluster || g1 == oldCluster ||
	  g2 == newCluster || g2 == oldCluster)
	{
	  
	  double cg2 = clusterUnigrams[g2];
	  if (g2 == oldCluster)
	    cg2 -= c;
	  if (g2 == newCluster)
	    cg2 += c;
	  assert(cg2 > 0);
	  double bg = clusterBigrams(g1,g2);
	  // fiddly bit comes here
	  if (g1 == newCluster){
	    bg += right[g2];
	    if (g2 == newCluster)
	      bg += doubles;
	    else if (g2 == oldCluster)
	      bg -= doubles;
	  }
	  if (g1 == oldCluster){
	    bg -= right[g2];
	    if (g2 == newCluster)
	      bg -= doubles;
	    else if (g2 == oldCluster)
	      bg += doubles;
	  }
	  if (g2 == newCluster){
	    bg += left[g1];
	  }
	  if (g2 == oldCluster){
	    bg -= left[g1];
	  }
	  assert(bg >= 0);
	  if (bg > 0){
	    double factor = bg * log(bg / (cg1 * cg2));
	    //cerr << g1 << " " <<  g2 << " f = " << factor << endl;
	    after += factor;
	  }
	  //else
	    //cerr << g1 << " " <<  g2 << " f = 0.0000000"  << endl;
	}
    }
  }
  //  assert(finite(after));
  cerr << "Slow: before " << before << " after " << after << " diff " << after - before << endl;
  //exit(0);
  return after - before;
}

/**
 * just check that all of te bigram counts are ok
 * and print an error if there is a problem
 */

bool 
Clusters::
validateBigramCounts()
{
  //cerr << "Starting validation ... " << endl;
  MatrixInt bg(numberClasses,numberClasses);
  for (int i = 0; i < numberTokens -1; i++){
    int w1 = data[i];
    int w2 = data[i+1];
    bg(classVector[w1],classVector[w2])++;
  }
  for (int x = 0; x < numberClasses;x++)
    for (int y = 0; y < numberClasses;y++){
      int correct = bg(x,y);
      if (correct != clusterBigrams(x,y)){
      	cerr << x << " " << y << " correct " << bg(x,y) << " estimated " << clusterBigrams(x,y) << endl;
	//exit(-1);
      }
      clusterBigrams(x,y) = correct;
      //  assert(bg(x,y) == clusterBigrams(x,y));
    }
  //cerr << ".. finished validation" << endl;
  return true;
}


//
// this calculates the change in linear time
//
// this should return  identical values
// to calcChange()
// 

double mi(double bg, double u1, double u2)
{
  if (bg > 0)
    return bg * log(bg / (u1 * u2));
  else
    return 0.0l;
}
   

double 
Clusters::
calcChangeFast(int w, 
	       double oldStringProb,
	       const vector<int> & wordVector,
	       int newCluster, 
	       const vector<int> & left, 
	       const vector<int> & right,
	       int doubles,
	       const vector<double> & prior) 
const
{
  int oldCluster = classVector[w];
  double before = 0;
  double after = 0;
  before = log(prior[oldCluster]) * PRIOR_BOOST;
  after = log(prior[newCluster]) * PRIOR_BOOST;
  if (numberStates > 0){
    double newStringProb = hmms[newCluster]->probability(wordVector);
    if (newStringProb == 0.0l){
      after += -1e20;
    }
    else {
      after += log(newStringProb);
    }
    before += log(oldStringProb);
  }
  int c = counts[w];
  if (c > 0){
    double oldUnigramCount = clusterUnigrams[oldCluster];
    double newUnigramCount = clusterUnigrams[newCluster];
  //  cerr << "w " << w << " old " << oldCluster << " new " << newCluster << endl;
    if (c == clusterUnigrams[oldCluster]){
      //    cerr << "ignoring singleton cluster " << endl;
      return 0.0l;
    }
    if (oldCluster == newCluster)
      return 0.0l;
    // do the four stripes
    for (int g1 = 0; g1 < numberClasses; g1++){
      if (g1 != oldCluster && g1 != newCluster){
	const double leftVal = left[g1];
	const double rightVal = right[g1];
	const double unigramCount = clusterUnigrams[g1];
	// first do it for the row
	// (newCluster, g1)
	{
	  
	  double bg = clusterBigrams(newCluster,g1);
	  double xu = newUnigramCount;
	  double yu = unigramCount;
	  before += mi(bg,xu,yu);
	  //cerr << newCluster << " " << g1 << " " << mi(bg + rightVal, xu + c, yu) <<  endl;
	  after += mi(bg + rightVal, xu + c, yu);
	}
	// (oldCluster , g1)
	{
	  double bg = clusterBigrams(oldCluster,g1);
	  double xu = oldUnigramCount;
	  double yu = unigramCount;
	  before += mi(bg,xu,yu);
	  after += mi(bg - rightVal, xu - c, yu);
	  //cerr << oldCluster << " " << g1 << " " << mi(bg - rightVal, xu - c, yu) << endl;
	}
	// (g1, newCluster)
	{
	  double bg = clusterBigrams(g1,newCluster);
	  double xu = unigramCount;
	  double yu = newUnigramCount;
	  before += mi(bg,xu,yu);
	  after += mi(bg + leftVal, xu, yu + c);
	  //cerr << g1 << " " << newCluster << " " <<  mi(bg + leftVal, xu, yu + c) << endl;
	}
	// (g1, oldCluster)
	{
	  double bg = clusterBigrams(g1,oldCluster);
	  double xu = unigramCount;
	  double yu = oldUnigramCount;
	  before += mi(bg,xu,yu);
	  after += mi(bg - leftVal, xu, yu - c);
	  //cerr << g1 << " " << oldCluster << " " << mi(bg - leftVal, xu, yu - c) << endl;
	}
      }
    }
  // We just have to do the four intersection points now
    // remembering not to double count
  // (oldCluster, oldCluster)
  {
    double bg = clusterBigrams(oldCluster,oldCluster);
    before += mi(bg,oldUnigramCount,oldUnigramCount);
    // we add the doubles to avoid double counting
    after += mi(bg + doubles - right[oldCluster] - left[oldCluster], oldUnigramCount - c, oldUnigramCount -c);
    //cerr << oldCluster << " " << oldCluster << " " <<  
    //mi(bg + doubles - right[oldCluster] - left[oldCluster], oldUnigramCount - c, oldUnigramCount -c) << endl;
  }
  // (newCluster, newCluster)
  {
    double bg = clusterBigrams(newCluster,newCluster);
    before += mi(bg,newUnigramCount,newUnigramCount);
    after += mi(bg + right[newCluster] + left[newCluster] + doubles, newUnigramCount + c, newUnigramCount +c);
    //    cerr << newCluster << " " << newCluster << " " 
    //<< mi(bg + right[newCluster] + left[newCluster] + doubles, newUnigramCount + c, newUnigramCount +c) << endl;
  }
  // (oldCluster newCluster)
  {
    double bg = clusterBigrams(oldCluster,newCluster);
    before += mi(bg, oldUnigramCount, newUnigramCount);
    after += mi(bg + left[oldCluster] - right[newCluster] - doubles, oldUnigramCount - c, newUnigramCount + c);
    // cerr << oldCluster << " " << newCluster << " "
    //	 << mi(bg + left[oldCluster] - right[newCluster] - doubles, oldUnigramCount - c, newUnigramCount + c) << endl;
  }
  // (newCluster oldCluster)
  {
    double bg = clusterBigrams(newCluster,oldCluster);
    before += mi(bg, newUnigramCount,oldUnigramCount);
    after += mi(bg + right[oldCluster] - left[newCluster] - doubles, newUnigramCount + c, oldUnigramCount - c);
    //cerr << newCluster << " " << oldCluster << " "
    //	 << mi(bg + right[oldCluster] - left[newCluster] - doubles, newUnigramCount + c, oldUnigramCount - c) << endl;
  }
  // Note that the totals of evrythign must be invariant
  //  assert(finite(before));
  //assert(finite(after));
  //  cerr << "Quick: before " << before << " after " << after << " diff " << after - before << endl;
  //exit(0);
  }
  return after - before;
}


/**
 * Move word w from oldCluster to newCluster
 * and adjust all of the figures correctly
 */

void 
Clusters::moveWord(int w,
		   int oldCluster,
		   int newCluster,
		   int doubles,
		   vector<int> & left,
		   vector<int> & right)
{
  //validateBigramCounts();
  classVector[w] = newCluster;
  if (counts[w] > 0){
    // out of old cluster
    clusterUnigrams[oldCluster] -= counts[w];
    for (int i = 0; i < numberClasses; i++){
      clusterBigrams(i,oldCluster) -= left[i];
      clusterBigrams(oldCluster,i) -= right[i];
    }
    clusterBigrams(oldCluster,oldCluster) += doubles;
    left[oldCluster] -= doubles;
    right[oldCluster] -= doubles;
    

    
    left[newCluster] += doubles;
    right[newCluster] += doubles;
    // into new cluster
    clusterUnigrams[newCluster] += counts[w];
    
    
    for (int i = 0; i < numberClasses; i++){
      clusterBigrams(i,newCluster) += left[i];
      clusterBigrams(newCluster,i) += right[i];
      
    }
    // adjust for doubles
    // new new bigrams will be too low
    clusterBigrams(newCluster,newCluster) -= doubles;
    
    
    // and the crossover
    //  clusterBigrams(newCluster,oldCluster) -= doubles;
  //clusterBigrams(oldCluster,newCluster) -= doubles;
    // old old bigrams will be too low
    // since we ill have double counted them
  }
  //cerr << "Moving " << w << "(" <<  *(corpus.wordArray[w])  <<  ") from " << oldCluster << " to " << newCluster << " count = " << counts[w] 
  //       << " doubles " << doubles << endl;
  //  validateBigramCounts();
}


