/*
 * PQ.h
 *
 *  Created on: 16 Aug 2018
 *      Author: oruqimaru
 */

#ifndef LEARN_PQ_PQ_H_
#define LEARN_PQ_PQ_H_

#include "learn/PQ/kmeans.h"
#include "load/DataUtil.h"
#include <unordered_map>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

using namespace std;
using std::vector;
using std::string;


typedef vector<float> Point;
/*
 * this function is used to learn the code book for PQ
 */
class PQquantizer{
public:
	PQquantizer(int numData,int seg, int datadim, int nCent, const vector<float>& originData){
		nData = numData;
		nCodeBook = seg;
		dataDim = datadim;
		subDim = dataDim / nCodeBook;
		nCentro = nCent;
		codeBooks.resize(0);
		bucketBelong.resize(numData);
		modifiedData = new char* [nCodeBook];
		for(int i = 0; i < datadim * numData; i++){
			modifiedData[i] = new float [subDim];
			for(int j = 0 ; j < subDim; j++){
				modifiedData[i][j] = originData[i*subDim + j];
			}
		}
	}

	~PQquantizer(){
		for(int i = 0; i < nCentro; i++)
			delete [] modifiedData[i];
		delete [] modifiedData;
	}

	void learn(){
		for(int i = 0; i < nCodeBook; i++){
			KMeans newK(subDim,nCentro);
			newK.SetInitMode(KMeans::InitUniform);
			int* label = new int[nData];
			newK.Cluster(modifiedData[i], nData, label);
			for(int j = 0; j < nData; j++)
				bucketBelong[i].push_back(label[i]);
			vector<vector<float>> aimingCBB;
			aimingCBB.resize(nCentro);
			for(int j = 0; j < nCentro; j++){
				aimingCBB[j].resize(dataDim);
				for(int k = 0; k < dataDim; k++)
					aimingCBB[j][k] = newK.m_means[j][k];
			}
			codeBooks.push_back(aimingCBB);
			delete label;
		}
	}

	void save_codeBooks(string outputPath){
		/*
		std::ofstream file_stream(filepath, std::ios_base::out);
	  if (!file_stream) {
	    cout << "Could not open file " << filepath << endl;
	    return;
	  }

	  // Copy all data to file_stream, then append a newline.
	  for (auto &mean : means_) {
	    std::ostream_iterator<float> itr(file_stream, " ");
	    std::copy(mean.data_.begin(), mean.data_.end(), itr);
	    file_stream << endl;
	  }
	  return;
		*/
		for(int i = 0; i < nCodeBook; i++){
			stringstream tag;
			tag << "/No." << i << "CB.fvecs";
			string new_outpath = outputPath + tag.str();
			ofstream file_stream(new_outpath);
			if (!file_stream) {
				    cout << "Could not open file " << new_outpath << endl;
				    return;
			}

			for(int j = 0; j < nCentro; j++){
				for(int k = 0; k < subDim; k++){
					file_stream << codeBooks[i][j][k] << " ";
				}
				file_stream << endl;
			}
		}
	}

	void save_buckets_same_vectorfile(string outputPath){
		unordered_map<string, int> buckets;
		for(int i = 0; i < nData; i++){
			stringstream ss;
			stringstream bucketId;
			ss << outputPath << "/" << "buckets";
			for(int j = 0; j < nCodeBook - 1; j++)
				bucketId << bucketBelong[i][j] << "-";
			bucketId << bucketBelong[i][nCodeBook - 1];
			ss << bucketId.str();
			if(buckets.find(bucketId.str()) != buckets.end())
				buckets[bucketId.str()]++;
			else buckets.insert(std::make_pair(bucketId.str(),1));
			ofstream file_stream(bucketId.str(),ios::app);
			file_stream << i << endl;
		}
		stringstream countOut;
		countOut << outputPath << "/count_result";
		ofstream file_stream(countOut.str());
		vector<pair<string, int>> ordered_by_value(buckets.begin(), buckets.end());
		std::sort(ordered_by_value.begin(), ordered_by_value.end(),[](const pair<string, int>& lhs, const pair<string, int>& rhs){
			return lhs.second < rhs.second;
		});
		for(auto & aim : ordered_by_value){
			file_stream << aim.second << endl;
		}
	}
public:
	float** modifiedData;
	vector<vector<Point>> codeBooks;
	vector<vector<int>> bucketBelong;
	int nData;							//number of origin data
	int nCodeBook;                //how much codeBooks are reuquried
	int subDim;							//what's the dim of vectors in a codeBook
	int dataDim;						//what's the orgin dim of the dataDimension
	int nCentro;							//how many codes in one codeBook
};




#endif /* LEARN_PQ_PQ_H_ */
