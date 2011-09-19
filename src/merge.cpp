// $Id: merge.cpp,v 1.2 2002/08/14 13:46:18 clark Exp $
// Evaluate the perplexity of some test data
// with respect to a class-based bigram model



#include <iostream>
#include <stdio.h>

#include "simplecorpus.h"
#include "clusters.h"
#include "matrix.h"
#include <math.h>



void printUsage(const char * arg0)
{
  cerr << "Usage " << arg0
       << "[options] trainingData  clusterFile" << endl;
  exit(-1);
}

/**
 * This calculates the change of merging the two classes i and j.
 * We want this to reflect the distance and not the change in perplexity.
 */

double calculateChange(const Matrix & distribution, int i, int j)
{
  // (i,j) represents the probability of going from class i to class j
  // so we can just calculate D(D_i || D_i + D_j) + D(D_j || ....)
  double value = 0.0l;
  for (int y = 0; y < distribution.dim2(); y++){
    double vi = distribution(i,y);
    double vj = distribution(j,y);
    //cerr << "vi,vj" << vi << " " << vj << endl;
    double vk = 0.5l * (vi + vj);
    if (vi > 0)
      value += vi * log(vi /vk);
    if (vj > 0)
      value += vj * log(vj/vk);
  }
  for (int x = 0; x < distribution.dim1(); x++){
    double vi = distribution(x,i);
    double vj = distribution(x,j);
    //cerr << "vi,vj" << vi << " " << vj << endl;
    double vk = 0.5l * (vi + vj);
    if (vi > 0)
      value += vi * log(vi /vk);
    if (vj > 0)
      value += vj * log(vj/vk);
  }
  
  return value;
}


/**
 * This program takes a set of clusters and merges them into a smaller set of clusters
 *  using the same sort of criterion, but not weighted by the counts.
 * I hope this means it will merge very frequent very similar clusters together
 * thus getting a more skew class distribution
 */

int main(int argc, char* argv[])
{
  if (argc != 3){
    printUsage(argv[0]);
  }
  char * trainingData = argv[1];
  char * clusterData = argv[2];
  int numberClasses = 0;
  vector<int> classVector;
  vector<string*> wordList;
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
	wordList.push_back(new string(word));
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

  MatrixInt bigrams(numberClasses,numberClasses);
  
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
    bigrams(cluster1,cluster0)++;
    wordProbs[w]++;
    cluster1 = cluster0;
  }
  cerr << "Loaded bigrams ok." << endl;
  vector<bool> ok(numberClasses,true);
  for (int iteration = 0 ; iteration < (numberClasses/2); iteration++){
    cerr << "Starting iteration " << iteration << endl;
    Matrix distribution(numberClasses,numberClasses);
    for (int x = 0; x < numberClasses; x++){
      if (ok[x]){
	double sum = 0.0l;
	for (int y = 0; y < numberClasses; y++){
	  sum += bigrams(x,y);
	}
	if (sum == 0.0l){
	  cerr << "Blank class " << x << endl;
	  for (int y = 0; y < numberClasses; y++){
	    distribution(x,y) = 1.0 / numberClasses;
	  }
	}
	else {
	  for (int y = 0; y < numberClasses; y++){
	    distribution(x,y) = bigrams(x,y) /sum;
	  }
	}
      }
    }
    
    cerr << "Choosing best pair ... " << endl;
    double bestD = -1.0l;
    int bestX = -1;
    int bestY = -1;
    for (int x = 0; x < numberClasses; x++){
      if (ok[x]){
	for (int y = x+1; y < numberClasses; y++){
	  if (ok[y]){
	    double score = calculateChange(distribution,x,y);
	    //cerr << x << "," << y << " " << score << endl;
	    if (bestD < 0.0l || score < bestD){
	      bestD = score;
	      bestX = x;
	      bestY = y;
	      //cerr << "New winner " << score << " at (" << x << "," << y << ")" << endl;
	    }
	  }
	}
      }
    }
    
    cerr << bestX << " " << bestY << endl;
    for (int i = 0; i < numberTypes; i++){
      if (classVector[i] == bestY)
	classVector[i] = bestX;
    }
    cerr << "adjust bigram counts " << endl;
    int cc = bigrams(bestY,bestY);
    bigrams(bestX,bestX) -= cc;
    for (int i = 0; i < numberClasses; i++){
      bigrams(i,bestX) += bigrams(i,bestY);
      bigrams(bestX,i) += bigrams(bestY,i);
      bigrams(i,bestY) = 0;
      bigrams(bestY,i) = 0;
    }
    ok[bestY] = false;
  }
  vector<int> newClassLabels(numberClasses,0);
  int currentClass = 0;
  for (int i = 0; i < numberClasses; i++){
    if (ok[i]){
      newClassLabels[i] = currentClass;
      currentClass++;
    }
    else {
      newClassLabels[i] = -1;
    }
  }
  vector<double> classSums(numberClasses,0.0l);
  for (int i = 0; i < numberTypes; i++){
    int c = classVector[i];
    if (ok[c]){
      classSums[c] += wordProbs[i];
    }
    else {
      cerr << "Found word " << i << " in defunct class " << c << endl;
    }
  }
  cerr << "Dumping...\n";
  for (int i = 0; i < numberTypes; i++){
    cout << *(wordList[i]) << " " 
	 << newClassLabels[classVector[i]] << " " 
	 << double(wordProbs[i])/classSums[classVector[i]] << endl;
  }
  return 0;
}
