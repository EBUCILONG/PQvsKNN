/*
 * distV2V.cpp
 *
 *  Created on: 20 Sep 2018
 *      Author: oruqimaru
 */


#include "load/DataUtil.h"
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include "omp.h"
#include <vector>
#include <random>
#include <string>
#include <math.h>
#include <utility>

int nt = 1000000;
int tDim = 960;

float distan(float* a, float* b){
	float sum;
	for(int i =0; i < tDim; i++){
		float sub = a[i] - b[i];
		float squ = sub*sub;
		sum += squ;
	}
	return sqrt(sum);
}

using namespace std;

int main(int argc, char** argv){
	cout << "----------------------start----------------";
	nt = stoi(argv[1]);
	tDim = stoi(argv[2]);
	string tInputPath(argv[3]);
	string idInputPath(argv[4]);
	int nbucket = stoi(argv[5]);
	string outPath(argv[6]);

	ifstream idfile;
	idfile.open(idInputPath,ios::in);

	vector<int> buckets[60000];
	for(int i= 0; i< nbucket; i++){
		int num = 0;
		num << idfile;
		for(int j = 0; j < num; j++){
			int buf;
			buf << idfile;
			buckets[i].push_back(buf);
		}
	}
	idfile.close();
	float* points = malloc(nt * tDim);

	ReadOneDimensionalPoints(tInputPath,FVEC,points,nt,tDim);
	vector<float> thread_out[10];
	for(int i = 0; i < 10; i++) thread_out[i].reserve(1000000);

	omp_set_num_threads(10);

	#pragma omp parallel
	{
	#pragma omp for
		for (int i = 0; i < nbucket; i++){
			vector<int> idList(buckets[i]);
			int n = idList.size();
			vector<int> randomShu;
			randomShu.reserve(n*n);
			for(int i = 0; i < n*n; i++)
				randomShu.push_back(i);
			random_shuffle(randomShu.begin(), randomShu.end());
			pair<int,int> vecPair[1000];
			for(int i = 0; i < 1000; i++){
				int resulting = randomShu[i];
				int x = idList[resulting / n];
				int y = idList[resulting % n];
				thread_out[omp_get_thread_num()].push_back(distan(&points[x * tDim],&points[y * tDim]));
			}
		}
	}

	ifstream outfile;
	outfile.open(outPath,ios::out);
	for(int i = 0; i < 10; i++){
		int len =thread_out[i].size();
		for(int j = 0; j < len; j++){
			outfile << thread_out[i][j];
			outfile << "\n";
		}
	}
	outfile.close();
}




