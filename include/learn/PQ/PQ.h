/*
 * PQ.h
 *
 *  Created on: 16 Aug 2018
 *      Author: oruqimaru
 */

#ifndef LEARN_PQ_PQ_H_
#define LEARN_PQ_PQ_H_

extern "C" {
#include "learn/PQ/vlfeat/vl/kmeans.h"
}
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
		cout << "coming in" << endl;
        nData = numData;
		nCodeBook = seg;
		dataDim = datadim;
		subDim = dataDim / nCodeBook;
		nCentro = nCent;
		codeBooks.resize(0);
		bucketBelong.resize(numData);
		distTable.resize(numData);
		modifiedData = new float* [nCodeBook];
		for(int i = 0; i < nCodeBook; i++){
			modifiedData[i] = new float [subDim*nData];
			for(int j = 0 ; j < nData; j++){
                for(int k = 0; k < subDim; k++){
                    modifiedData[i][j*subDim+k] = originData[j*datadim + i*subDim + k];
                }
            }
            cout << i << " "; 
		}
        cout << "go away" << endl;
	}

	~PQquantizer(){
		for(int i = 0; i < nCodeBook; i++)
			delete [] modifiedData[i];
		delete [] modifiedData;
	}

	void learn(int maxIter = 100){
		for(int i = 0; i < nCodeBook; i++){
			VlKMeans * newK = vl_kmeans_new(VL_TYPE_FLOAT, VlDistanceL2);
			vl_kmeans_set_algorithm(newK, VlKMeansLloyd);
			vl_kmeans_init_centers_with_rand_data(newK, modifiedData[i], subDim, nData, nCentro);
			vl_kmeans_set_max_num_iterations (newK, 100) ;
			vl_kmeans_refine_centers (newK, modifiedData, nData) ;
			float* centers = vl_kmeans_get_centers(newK);
			vector<vector<float>> aimingCBB;
			aimingCBB.resize(nCentro);
			for(int j = 0; j < nCentro; j++){
				aimingCBB[j].resize(subDim);
				for(int k = 0; k < subDim; k++)
					aimingCBB[j][k] = centers[j*subDim + k];
			}
			codeBooks.push_back(aimingCBB);
			vl_uint32 * label = vl_malloc(sizeof(vl_uint) * nData) ;
			float * distances = vl_malloc(sizeof(float) * nData) ;
			vl_kmeans_quantize(newK, label, distances, modifiedData[i], nData);
			for(int j = 0; j < nData; j++)
				bucketBelong[j].push_back((int)label[j]);
			for(int j = 0; j < nData; j++)
				distTable[j].push_back(distances[j]);
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
        cout << "code book save finished" << endl;
	}

	void save_buckets_same_vectorfile(string outputPath){
        unordered_map<string, int> buckets;
		for(int i = 0; i < nData; i++){
            cout << bucketBelong[i].size() << " ";
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
            cout << bucketId.str() << " "<<ss.str();
			ofstream file_stream(ss.str(),ios::app);
            if(!file_stream) cout << "cannot open dir" << endl;
			file_stream << i << endl;
            cout << " done " << endl;
		}
        cout << "finish recursion" << endl;
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
	vector<vector<float>> distTable;
	vector<vector<Point>> codeBooks;
	vector<vector<int>> bucketBelong;
	int nData;							//number of origin data
	int nCodeBook;                //how much codeBooks are reuquried
	int subDim;							//what's the dim of vectors in a codeBook
	int dataDim;						//what's the orgin dim of the dataDimension
	int nCentro;							//how many codes in one codeBook
};




#endif /* LEARN_PQ_PxQ_H_ */
