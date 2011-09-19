// $Id: perplexity.cpp,v 1.2 2002/08/14 13:46:18 clark Exp $
// Evaluate the perplexity of some test data
// with respect to a class-based bigram model



#include <iostream>
#include <stdio.h>

#include "simplecorpus.h"
#include "clusters.h"
#include "matrix.h"
#include <math.h>


int main(int argc, char* argv[])
{
  if (argc != 4){
    cerr << "Usage: " << argv[0] << " trainingData testData clusterFile" << endl;
    exit(-1);
  }
  char * trainingData = argv[1];
  char * testData = argv[2];
  char * clusterData = argv[3];
  int numberClasses = 0;
  vector<int> classVector;
  map<string,int> dictionary;  // maps words to index
  cerr << "Loading clusters now ...";
  ifstream inc(clusterData);
  if (!inc){
    cerr << "Coulnd't open cluster file" << endl;
    exit(-1);
  }
  const int BUFLEN = 1000;
  char buffer[BUFLEN];
  while (inc){
    // read in a line and parse it
    // remember to handle the blank line correctly
    inc.getline(buffer,BUFLEN);
    // cerr << buffer << endl;
    int cluster = -1;
    float prob = -1.0;
    char wordBuffer[BUFLEN];
    int result = sscanf(buffer,"%s %d %f", wordBuffer, &cluster, &prob);
    string word;
    if (result != 3){
      cerr << "Blank line?" << endl;
      result = sscanf(buffer," %d %f", &cluster, &prob);
      word = "";
    }
    else {
      word = wordBuffer;
    }
    if ((cluster > -1) && (prob > -0.5))
      {
	// ignore the prob
	dictionary[word] = classVector.size();
	classVector.push_back(cluster);
	//cout << "Read " << word << " " << cluster << endl;
	if (cluster  >= numberClasses){
	  numberClasses = cluster + 1;
	}
      }
    else {
      cerr << "eof?\n";
    }
  }
  cerr << "Read " << numberClasses << " clusters.\nStarting to read training data\n" << endl;
  Matrix transitions(numberClasses,numberClasses);
  MatrixInt transitionCounts(numberClasses,numberClasses);
  
  int numberTypes = classVector.size();
  vector<double> wordProbs(numberTypes,0.0l);
  ifstream in(trainingData);
  if (!in){
    cerr << "Coulnd't open cluster file" << endl;
    exit(-1);
  }
  int blankWord = dictionary[string("")];
  // cluster1  is the previosu cluster
  int cluster1 = classVector[blankWord];
  while (in){
    // read in a line and parse it
    // remember to handle the blank line correctly
    in.getline(buffer,BUFLEN);
    string word(buffer);
    if (dictionary.find(word) == dictionary.end()){
      cerr << "Fatal error -- OOV word " << word << endl;
      exit(-1);
    }
    int w = dictionary[word];
    int cluster0 = classVector[w];
    transitionCounts(cluster1,cluster0)++;
    wordProbs[w]++;
    cluster1 = cluster0;
  }
  cerr << "Finished loading .. starting to smooth\n";
  vector<double> beta(numberClasses,0.0l);
  vector<int> g0v(numberClasses,0);
  vector<int> nv(numberClasses,0);
  int g1 = 0;  // number of bigrams with count 1
  int g2 = 0; //
  for (int i = 0; i < numberClasses; i++){
    for (int j = 0; j < numberClasses; j++){
      int c = transitionCounts(i,j);
      nv[i] +=c;
      if (c == 0){
	g0v[i]++;
      }
      if (c == 1){
	g1++;
	beta[j]++;
      }
      if (c == 2){
	g2++;
      }
    }
  }
  // normalise beta
  double betaSum = 0.0l;
  for (int i = 0; i < numberClasses; i++){
    if (beta[i] < 1.0){
      beta[i] = 1;
    }
    betaSum += beta[i];
  }
  for (int i = 0; i < numberClasses; i++){
    beta[i] /= betaSum;
  }
  
  double bInd = (1.0 *  g1) / (g1 + 2.0 * g2);
  cerr << "History  Independent Discount = " << bInd << endl;
  for (int i = 0; i < numberClasses; i++){
    for (int j = 0; j < numberClasses; j++){
      int c = transitionCounts(i,j);
      double value = (c - bInd)/(1.0 * nv[i]);
      if (value < 0)
	value = 0;
      value += (numberClasses - g0v[i]) * bInd * beta[j] / nv[i];
      if (value == 0.0l){
	cerr << i << "," << j << " c=" << c << " error -- zero prob" << endl;
	cerr << "beta j " << beta[j] << endl;
	cerr << "g0v " << g0v[i] << endl;
	cerr << "nvi " << nv[i] << endl;
	exit(-1);
      }
      transitions(i,j) = value;
    }
  }
  // smooth the membership distribution
  for (int i = 0; i < numberClasses;i++){
    double sum = 0;
    double n0 = 0;
    double n1 = 0;
    double n2 = 0;
    double num = 0.0l;
    for (int j = 0; j < numberTypes; j++){
      if (i == classVector[j]){
	num++;
	double c = wordProbs[j];
	if (c < 0.5l){
	  n0++;
	}
	else if (c < 1.5l){
	  n1++;
	}
	else if (c < 2.5l){
	  n2++;
	}
	sum += c;
      }
    }
    if (n0 > 0){
      if (n1 == 0) n1++;
      double b = n1 / (1.0 + n1 + 2 * n2);

      //      cerr << "Word class discount = " << b << endl;
      for (int j = 0; j < numberTypes; j++){
	if (i == classVector[j]){
	  double c = wordProbs[j];


	  if (j == 44902){
	    cerr << "DEBUG before cchoppy " << wordProbs[j] << endl;
	    cerr << "num " << num  << " n0 " << n0 << " n1 " << n1 
		 << " n2 " << n2 << " b " << b << endl;
	  }
	  if (c < 0.5l){
	    wordProbs[j] = (num - n0) * b / (n0 * sum);
	  }
	  else {
	    wordProbs[j] = (c - b)/ sum;
	  }
	  if (j == 44902)
	    cerr << "DEBUG after cchoppy " << wordProbs[j] << endl;
	}
      }
    }
    else {
      cerr << "No smoothing in this class" << endl;
      for (int j = 0; j < numberTypes; j++){
	if (i == classVector[j]){
	  double c = wordProbs[j];
	  wordProbs[j] = c /sum;
	}
      }
    }
  }
  cerr << "Checking normalisation\n";
  for (int i = 0; i < numberClasses; i++){
    double sum = 0.0l;
    for (int j = 0; j < numberClasses; j++){
      sum += transitions(i,j);
    }
    double wsum = 0.0l;
    for (int w = 0; w < numberTypes; w++){
      if (i == classVector[w]){
	wsum += wordProbs[w];
      }
    }
    if (fabs(wsum - 1.0) + fabs(sum -1.0) > 0.000001){
      cerr << "LM normalisation error?" << endl;
      cerr << i << "tsum  = " << sum << ", wsum=" << wsum <<  endl;
    }
  }
  cerr << "Starting to load test data " << endl;
  ifstream inTest(testData);
  if (!inTest){
    cerr << "Couldn't open test data  file" << endl;
    exit(-1);
  }
  cluster1 = classVector[blankWord];
  int count = 0;
  double logProb = 0.0l;
  while (inTest){
    inTest.getline(buffer,BUFLEN);
    string word(buffer);
    if (dictionary.find(word) == dictionary.end()){
      cerr << "Fatal error -- OOV word in test data" << word << endl;
      exit(-1);
    }
    int w = dictionary[word];
    int cluster0 = classVector[w];
    double transitionProb = transitions(cluster1, cluster0);
    double outputProb = wordProbs[w];
    double prob = transitionProb * outputProb ;
    if (prob == 0){
      cerr << "Fatal error \n";
      cerr << "class " << cluster0 << endl;
      cerr << "Word index " << w << endl;
      cerr << word << " " << transitionProb << ", " << outputProb << "=" <<  prob << endl;
      exit(-1);
    }
    logProb += log(prob);
    count++;
    //    cout << word << " " << prob << endl;
    cluster1 = cluster0;
  }
  double perplexity = exp(-1.0 * logProb/count);
  cerr << "Perplexity " << perplexity << endl;
  return 0;
}
