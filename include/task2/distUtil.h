/*
 * distUtil.h
 *
 *  Created on: 21 Aug 2018
 *      Author: oruqimaru
 */

#ifndef INCLUDE_TASK2_DISTUTIL_H_
#define INCLUDE_TASK2_DISTUTIL_H_

#include <vector>
#include <unordered_map>
#include <sstream>

#include "learn/PQ/PQ.h"

using namespace std;

vector<float> centro_to_item(PQquantizer myPQ){
	unordered_map<string, vector<int>> buckets;
	for(int i = 0; i < myPQ.nData; i++){
		stringstream bucketId;
		for(int j = 0; j < myPQ.nCodeBook - 1; j++)
			bucketId << myPQ.bucketBelong[j] << "-";
		bucketId << myPQ.bucketBelong[i][myPQ.nCodeBook - 1];
		if(buckets.find(bucketId.str()) != buckets.end())
			buckets[bucketId.str()].push_back(i);
		else {
			vector<int> vin;
			vin.push_back(i);
			buckets.insert(std::make_pair(bucketId.str(),vin));
		}
	}
	for (unordered_map<int, int>::iterator i = buckets.begin(); i != buckets.end(); i++){


	}

}



#endif /* INCLUDE_TASK2_DISTUTIL_H_ */
