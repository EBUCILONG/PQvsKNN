/*
 * ivfPQ.h
 *
 *  Created on: 21 Aug 2018
 *      Author: oruqimaru
 */

#ifndef INCLUDE_LEARN_PQ_IVFPQ_H_
#define INCLUDE_LEARN_PQ_IVFPQ_H_

#include <utility>
#include <unordered_map>
#include <queue>

#include "learn/PQ/PQ.h"

using namespace std;

class ivfPQquantizer : PQquantizer{
public:
	ivfPQquantizer(int numData,int seg, int datadim, int nCent, const vector<float>& originData){
		PQquantizer(numData, seg, datadim,nCent, originData);

	}



public:

};



#endif /* INCLUDE_LEARN_PQ_IVFPQ_H_ */
