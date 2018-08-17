/*
 * PQ.h
 *
 *  Created on: 16 Aug 2018
 *      Author: oruqimaru
 */

#ifndef LEARN_PQ_PQ_H_
#define LEARN_PQ_PQ_H_

#include "learn/PQ/point.h"
#include "learn/PQ/kmeans.h"
#include "load/DataUtil.h"
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

using namespace std;
using std::vector;
using std::string;
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
		modifieData.resize(nCodeBook);
		codeBooks.resize(0);
		bucketBelong.resize(numData);
		for(int i = 0; i < nCodeBook; i++){
			 for(int j = 0; j < nData; j++){
				 vector<float> newPoint;
				 newPoint.reserve(subDim);
				 int bias = j*dataDim + i*subDim;
				 for(int k = 0; k < subDim; k++){
					 newPoint.push_back(originData[bias + k]);
				 }
				 modifieData[i].push_back(Point(newPoint));
			 }
		}
	}

	void learn(){
		for(int i = 0; i < nCodeBook; i++){
			KMeans newK(nCentro);
			newK.init(modifieData[i]);
			newK.run();
			codeBooks.push_back(newK.means_);
			for(int i = 0; i < nData; i++){
					bucketBelong[i].push_back(newK.points_[i].cluster_);
			}
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

			for (auto &mean : codeBooks[i]) {
					file_stream << subDim << " ";
				    std::ostream_iterator<float> itr(file_stream, " ");
				    std::copy(mean.data_.begin(), mean.data_.end(), itr);
				    file_stream << endl;
			}
		}
	}

public:
	vector<vector<Point> > modifieData;
	vector<vector<Point> > codeBooks;
	vector<vector<int> > bucketBelong;
	int nData;							//number of origin data
	int nCodeBook;                //how much codeBooks are reuquried
	int subDim;							//what's the dim of vectors in a codeBook
	int dataDim;						//what's the orgin dim of the dataDimension
	int nCentro;							//how many codes in one codeBook
};




#endif /* LEARN_PQ_PQ_H_ */
