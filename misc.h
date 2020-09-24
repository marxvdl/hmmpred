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

#ifndef MISC_H_
#define MISC_H_

#include <string>
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

using namespace std;

// #define real double
#define real float

typedef unsigned long int ulint;
typedef unsigned char byte;

const real BITS_PER_NAT = 1.44269504088896340735992468100;

// trim from start
static inline std::string &ltrim(std::string &s) {
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
	s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
	return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
	return ltrim(rtrim(s));
}

inline bool is_valid_probability_distribution(vector<real>& v){
	real sum = 0;
	for(uint i=0; i<v.size(); i++)
		sum += v[i];
	
	return abs( sum - 1 ) < 0.00001;
}

inline string vec2string(vector<byte>& v){
	ostringstream ss;
	for(uint i=0; i<v.size(); i++)
		ss << (int)v[i] << " ";
	return ss.str();
}

inline string vec2string(vector<byte>& v, map<byte,char>& revmap){
	ostringstream ss;
	for(uint i=0; i<v.size(); i++)
		ss << revmap[v[i]];
	return ss.str();
}

inline string vec2string(vector<byte>& v, map<byte,char>& revmap, uint maxSize){
	ostringstream ss;
	for(uint i=0; i<maxSize; i++)
		ss << revmap[v[i]];
	return ss.str();
}

inline string vecvec2string(vector< vector<real> >& v){
	ostringstream ss;
	for(uint i=0; i<v.size(); i++){
		for(uint j=0; j<v[i].size(); j++)
			ss << v[i][j] << " ";
		if(i != v.size() -1)
			ss << "; ";
	}
	return ss.str();
}

template<class T> inline void normalize(vector<T>& v){
	real sum = 0;
	for(uint i=0; i<v.size(); i++)
		sum += v[i];
	if(sum != 0)
		for(uint i=0; i<v.size(); i++)
			v[i] /= sum;
}

template<class T> inline void normalize(vector< vector<T> >& v){
    real sum = 0;
    for(uint i=0; i<v.size(); i++)
        for(uint j=0; j<v[i].size(); j++)
            sum += v[i][j];
    if(sum != 0)
        for(uint i=0; i<v.size(); i++)
            for(uint j=0; j<v[i].size(); j++)
                v[i][j] /= sum;
}

inline uint maxElement(vector<real>& v){
	double max = -1;
	uint kmax = 0;
	for(uint k=0; k<v.size(); k++){
		if(v[k] > max){
			max = v[k];
			kmax = k;
		}
	}
	return kmax;
}

inline int string2int(string s){
	int result;
	stringstream(s) >> result;
	return result;
}

inline string int2string(int n){
	stringstream ss;
	ss << n;
	return ss.str();
}

inline void printVec(vector<byte>& v){
	for(uint i=0; i<v.size(); i++)
		cout << v[i];
	cout << endl;
}

inline void zeroVector(vector<real>& v){
	for(uint i=0; i<v.size(); i++)
		v[i] = 0.0;
}

inline void zeroVector(vector< vector<real> >& v){
	for(uint i=0; i<v.size(); i++)
		for(uint j=0; j<v[i].size(); j++)
			v[i][j] = 0.0;
}


inline real stats_mean(vector<real>& v){
	real sum = 0;
	for(uint i=0; i<v.size(); i++)
		sum += v[i];
	return sum / (real)v.size();
}

inline real stats_variance(vector<real>& v) {
	real mean = stats_mean(v);

	real variance = 0;
	for(uint i=0; i<v.size() ; i++) {
		real s = v[i] - mean;
		variance += s*s;
	}

	variance /= (double)(v.size()-1);
	return variance;
}

inline real stats_stddev(vector<real>& v) {
	return sqrt(stats_variance(v));
}

inline ulint myclock() {
    return (ulint)clock() + (ulint)time(NULL);
}

inline vector<string> extract_filename_extendsion(const string& str){
	vector<string>r(2);
	size_t pos = str.find_last_of('.');
	if(pos == string::npos){
		r[0] = str;
		r[1] = "";
	}
	else{
		r[0] = str.substr(0,pos);
		r[1] = str.substr(pos);
	}
	return r;
	
	
}

#endif /* MISC_H_ */

