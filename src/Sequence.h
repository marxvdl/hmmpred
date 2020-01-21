/*
 * Copyright 2014-2020, Marx Gomes van der Linden
 *                      marx.linden@ifb.edu.br
 * 
 * This file is part of HmmPred.
 * 
 * HmmPred is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * HmmPred is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with HmmPred.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <map>
#include <vector>
#include <string>

#include "misc.h"
#include "mappings.h"

using namespace std;

class Sequence {
public:
	vector<byte> intSeq;

	Sequence(string& seq, map<char,byte>& mapping){
		intSeq.resize(seq.size());

		for(uint i=0; i<seq.size(); i++)
			intSeq[i] = mapping[seq[i]];
	}

	string charSeq(map<byte,char>& reverseMapping){
		string seq(intSeq.size(), ' ');
		for(uint i=0; i<intSeq.size(); i++)
			seq[i] = reverseMapping[intSeq[i]];
		return seq;
	}

};

#endif /* SEQUENCE_H_ */
