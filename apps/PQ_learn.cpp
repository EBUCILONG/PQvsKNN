/*
 * PQ.cpp
 *
 *  Created on: 16 Aug 2018
 *      Author: oruqimaru
 */
#include <DataUtil.h>
#include <iostream>
#include <string>

using namespace std;
using namespace arma;

int main(int argc, char** argv){
	int nt = stoi(argv[0]);
	int tDim = stoi(argv[1]);
	string tInputPath(argv[2]);
	int nd = stoi(argv[3]);
	int dDim = stoi(argv[4]);
	string dInputPath(argv[5]);
	int k = stoi(argv[6]);
	int nSep = stoi(argv[7]);
	int subDim = stoi(argv[8]);

	vector<float> points(nt * tDim);
	ReadOneDimensionalPoints(tInputPath,FVEC,points,nt,tDim);



}




