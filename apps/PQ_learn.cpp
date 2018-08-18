/*
 * PQ.cpp
 *
 *  Created on: 16 Aug 2018
 *      Author: oruqimaru
 */
#include "learn/PQ/PQ.h"
#include "learn/PQ/point.h"
#include "learn/PQ/kmeans.h"
#include "load/DataUtil.h"
#include <iostream>
#include <string>
//#include "learn/PQ/PQ.h"

using namespace std;
//using namespace arma;

int main(int argc, char** argv){
	int nt = stoi(argv[1]);
	int tDim = stoi(argv[2]);
	string tInputPath(argv[3]);
	int nd = stoi(argv[4]);
	int dDim = stoi(argv[5]);
	string dInputPath(argv[6]);
	int k = stoi(argv[7]);
	int nSep = stoi(argv[8]);
	int subDim = stoi(argv[9]);
	string outputPath(argv[10]);

	vector<float> points(nt * tDim);
	ReadOneDimensionalPoints(tInputPath,FVEC,points,nt,tDim);
	//PQquantizer(int numData,int seg, int datadim, int nCent, const vector<float>& originData){
	PQquantizer myPQ(nt, nSep, dDim,k,points);
	myPQ.learn();
	myPQ.save_codeBooks(outputPath);


}




