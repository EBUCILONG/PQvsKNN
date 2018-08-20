/*
 * PQ.cpp
 *
 *  Created on: 16 Aug 2018
 *      Author: oruqimaru
 */
#include "learn/PQ/PQ.h"
#include "load/DataUtil.h"
#include <iostream>
#include <string>
//#include "learn/PQ/PQ.h"

using namespace std;
//using namespace arma;

int main(int argc, char** argv){
	cout << "----------------------start----------------";
    int nt = stoi(argv[1]);
	int tDim = stoi(argv[2]);
	string tInputPath(argv[3]);
	int nd = stoi(argv[4]);
	int dDim = stoi(argv[5]);
	string dInputPath(argv[6]);
	int k = stoi(argv[7]);
	int nSep = stoi(argv[8]);
	int subDim = stoi(argv[9]);
	string codeBookOutputPath(argv[10]);
	string bucketOutputPath(argv[11]);
    
	vector<float> points(nt * tDim);
    cout << "ready to read";
	ReadOneDimensionalPoints(tInputPath,FVEC,points,nt,tDim);
	cout << "read ready";
    //PQquantizer(int numData,int seg, int datadim, int nCent, const vector<float>& originData){
	PQquantizer myPQ(nt, nSep, dDim,k,points);
	cout << "start quantizer ready";
    myPQ.learn();
    cout << "learn ready ";
    myPQ.save_codeBooks(codeBookOutputPath);
    cout << "code book saved"<< endl;
	myPQ.save_buckets_same_vectorfile(bucketOutputPath);
    cout << "buckets saved" << endl;

}




