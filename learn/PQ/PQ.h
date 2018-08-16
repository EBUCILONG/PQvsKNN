/*
 * PQ.h
 *
 *  Created on: 16 Aug 2018
 *      Author: oruqimaru
 */

#ifndef LEARN_PQ_PQ_H_
#define LEARN_PQ_PQ_H_

#include <point.h>
#include <DataUtil.h>
#include <vector>
#include <string>

using std::vector;
using std::string;
/*
 * this function is used to learn the code book for PQ
 */
class PQquantizer{
public:
	PQquantizer(){}
	vector<mat> codeBooks;

	void operator() (mat& dataMat, int numBook, int ncen){





	}

};




#endif /* LEARN_PQ_PQ_H_ */
