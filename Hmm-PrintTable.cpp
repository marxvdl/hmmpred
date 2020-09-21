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

#include "Hmm.h"
#include <cassert>

struct LineT1 {
	string str;
	real probState;
	vector<real> probNextFragments;
	LineT1(string s, real p, uint n) : str(s), probState(p), probNextFragments(n,0) {}
};
bool lineCompare1 (LineT1 i,LineT1 j) { return (i.probState > j.probState); }

/**
 * Prints out the transition probabilities that were learned during
 * the training phase.
 */
void Hmm::printTableTransitions(ofstream& out){
	vector<LineT1> lines;
	
	ulint j=-1;
	for(ulint i=0; i<numberOfStates; i++){
		vector<byte> state_i = indexToSequence(i);
		
		LineT1 line(state_vec2string(state_i, windowSize), probState[i], numberOfSecondarySymbols);
		
		for(ulint suc=0; suc<numberOfSecondarySymbols; suc++){
//  			vector<byte> frag_ij(windowSize, 0);                  //  }
// 			for(int x=0; x<windowSize-1; x++)                     //  }
// 				frag_ij[x] = state_i[x+1];                    //  }<- this block of code...
// 			frag_ij[windowSize-1] = suc;                          //  }
//  			ulint j = sequenceToIndex(windowSize, frag_ij);       //  }
			
			vector<byte> frag_ij = indexToSequence(++j);          //  <- ...should be exactly equivalent to this line
			
			for(int x=0; x<windowSize-1; x++)                     //
				assert(frag_ij[x] == state_i[x+1]);           //  assertions just in case
			assert(frag_ij[windowSize-1] == suc);                 //
			
			line.probNextFragments[suc] = probFragment[j]; 
			
		}
		lines.push_back(line);
	}
	
	doPrintTableTransitions(out, lines);
	
	_printTable_HQ = 0;
}

void Hmm::doPrintTableTransitions(ofstream& out, vector<LineT1>& lines){	
	sort(lines.begin(), lines.end(), lineCompare1);
	
	char buffer[100];
	out << "#	p(state)  ";
	for(uint i=0; i<numberOfSecondarySymbols; i++)
		out << "\t" << mySecRevMapping[i] << "         ";
// 	out << "\t\tH";
	out << "\n";
	
	real sum = 0;
	real count = 0;
	
	for(vector<LineT1>::iterator it = lines.begin(); it!=lines.end(); ++it){
		snprintf(buffer, 100, "\t%10.8f", it->probState);
		out << it->str << buffer;
		
// 		real entropy = 0;
		
		for(uint i=0; i<numberOfSecondarySymbols; i++){
			real prob = (it->probState? it->probNextFragments[i] / it->probState : 0);
// 			if(prob != 0)
// 				entropy -= prob * BITS_PER_NAT * log(prob);
			
			snprintf(buffer, 100, "\t%10.8f",  prob);
			out << buffer;
		}
		
// 		if(entropy == entropy) {
// 			sum += entropy * it->probState;
// 			count += it->probState;
// 		}
// 		snprintf(buffer, 100, "\t%10.8f",  entropy);		
// 		out << "\t   " << buffer;

		out << "\n";
	}
	
// 	out << "#  \t          ";
// 	for(uint i=0; i<numberOfSecondarySymbols; i++){
// 		out << "\t          ";
// 	}
// 	snprintf(buffer, 100, "\t%10.8f",  sum/count);
// 	out << "\t   " << buffer;
// 	
// 	snprintf(buffer, 100, "\t%10.8f",  (log(numberOfSecondarySymbols) * BITS_PER_NAT) - (sum/count));
// 	out << "  " << buffer << "\n";
}

struct LineT2 {
	string str;
	real probFragment;
	vector<real> probEmission;
	LineT2(string s, real p, uint n) : str(s), probFragment(p), probEmission(n,0) {}
};
bool lineCompare2 (LineT2 i,LineT2 j) { return (i.probFragment > j.probFragment); }

/**
 * Prints out the emission probabilities that were learned during
 * the training phase.
 */
void Hmm::printTableEmissions(ofstream& out){
	vector<LineT2> lines;
	
	ulint j=-1;
	for(ulint i=0; i<numberOfFragments; i++){
		vector<byte> frag_i = indexToSequence(i);
		LineT2 line(fragment_vec2string(frag_i, windowSize), probFragment[i], numberOfSecondarySymbols);
		
		line.probEmission.resize(numberOfPrimarySymbols - 1, 0.0);
		
		for(byte obs=0; obs<numberOfPrimarySymbols - 1; obs++)
			line.probEmission[obs] = probFragmentEmitsPrimarySymbol[i][obs];
		
		lines.push_back(line);
	}
	
	doPrintTableEmissions(out, lines);
}

void Hmm::doPrintTableEmissions(ofstream& out, vector<LineT2>& lines){
	sort(lines.begin(), lines.end(), lineCompare2);
	char buffer[100];
	
	out << "#\tp(frag)   ";
	for(byte obs=0; obs<numberOfPrimarySymbols - 1; obs++)
		out << "\t" << myPriRevMapping[obs] << "         ";
	out << "\n";
	
	for(vector<LineT2>::iterator it = lines.begin(); it!=lines.end(); ++it){
		snprintf(buffer, 100, "\t%10.8f", it->probFragment);
		out << it->str << buffer;		
		
		for(byte obs=0; obs<numberOfPrimarySymbols - 1; obs++){
			snprintf(buffer, 100, "\t%10.8f",  it->probEmission[obs]);
			out << buffer;
		}
		
		out << "\n";

	}
}


//
// Auxiliary functions for printTable
//
inline string Hmm::fragment_vec2string(const vector<byte>& vec, uint n){
	string str(n, '*');
	for(int i=0; i<n; i++){
		str[i] = mySecRevMapping[vec[i]];
	}
	return str;
}

inline string Hmm::state_vec2string(const vector<byte>& vec, uint n){
	string str(n-1, '*');
	for(int i=1; i<n; i++){
		str[i-1] = mySecRevMapping[vec[i]];
	}
	return str;
}

/// Old implementations from here below

// /**
//  * Prints out the transition probabilities that were learned during
//  * the training phase for given primary symbol.
//  */
// void Hmm::printTableForPrimarySymbol(ofstream& out, byte p){
// 
// 	real probPriSymbol = 0;
// 	
// 	vector<LineT1> lines;
// 	ulint j=-1;
// 	for(ulint i=0; i<numberOfStates; i++){
// 		vector<byte> state_i = indexToSequence(i);
// 		
// // 		out << state_vec2string(state_i, windowSize) << "   p(i) = "<< probState[i] << "\n";
// 		
// 		LineT1 line(state_vec2string(state_i, windowSize), 0, numberOfSecondarySymbols);		
// 		
// 		for(ulint suc=0; suc<numberOfSecondarySymbols; suc++){
// 			vector<byte> frag_ij = indexToSequence(++j);
// 
// 			out << "\t" << fragment_vec2string(frag_ij, windowSize) << " " 
// 				<< myPriRevMapping[p] 
// 				<< "      p(q|fi) = " << probFragmentEmitsPrimarySymbol[j][p]
// 				<< "      p(fi) = " << probFragment[j]
// 				<< "      " << (probFragmentEmitsPrimarySymbol[j][p] * probFragment[j])
// 				<< "\n";
// 				
// 				real probQ = probFragmentEmitsPrimarySymbol[j][p] * probFragment[j];
// 				
// 				probPriSymbol  += probQ;
// 				line.probState += probQ;
// 				line.probNextFragments[suc] = probQ;
// 				
// 				
// // 			line.probNextFragments[suc] = probFragmentEmitsPrimarySymbol[j][p] * probFragment[j];
// 		}
// 		normalize<real>(line.probNextFragments);
// 		lines.push_back(line);		
// 		
// 		out << "\n";
// 	}
// 	
// 	cout << probPriSymbol << "\n";
// 	_printTable_HQ -= probPriSymbol * log(probPriSymbol) * BITS_PER_NAT;
// 	
// 	for(vector<LineT1>::iterator it=lines.begin(); it != lines.end(); it++){
// 		it->probState /= probPriSymbol;
// 		for(j=0; j<it->probNextFragments.size(); j++){
// 			it->probNextFragments[j] *= it->probState;
// 		}
// // 		normalize<real>(it->probNextFragments);
// 	}
// 	
// 	doPrintTableTransitions(out, lines);
// 	
// }
// 
// void Hmm::calculateEntropyForPrimarySymbol(byte p){
// 	real probPriSymbol = 0;
// 	
// 	ulint j=-1;
// 	for(ulint i=0; i<numberOfStates; i++){		
// 		for(ulint suc=0; suc<numberOfSecondarySymbols; suc++){
// 			vector<byte> frag_ij = indexToSequence(++j);
// 				real probQ = probFragmentEmitsPrimarySymbol[j][p] * probFragment[j];				
// 				probPriSymbol  += probQ;
// 		}
// 	}
// 	
// 	cout << probPriSymbol << "\n";
// 	_printTable_HQ -= probPriSymbol * log(probPriSymbol) * BITS_PER_NAT;
// }

// void Hmm::printTableV2(ofstream& out){
// // 	char buffer[512];
// 	
// 	real const HQ  = _printTable_HQ;
// 	real       HQB = 0;
// 	
// 	ulint j=-1;
// 	for(ulint i=0; i<numberOfStates; i++){
// 		vector<byte> state_i = indexToSequence(i);
// 		
// 		for(ulint suc=0; suc<numberOfSecondarySymbols; suc++){
// 			vector<byte> frag_ij = indexToSequence(++j);
// 			
// // 				real p0 = probFragmentEmitsPrimarySymbol[j][0];
// // 				real p1 = probFragmentEmitsPrimarySymbol[j][1];
// // 				real p2 = probFragmentEmitsPrimarySymbol[j][2];
// // 				
// // 				HQB -= ( p0*log(p0) + p1*log(p1) + p2*log(p2) ) * BITS_PER_NAT * probFragment[j];				
// 				
// 				real tmp = 0;
// 				for(int z=0; z<20; z++){
// 					real p = probFragmentEmitsPrimarySymbol[j][z];
// 					tmp += p*log(p);
// 				}				
// 				HQB -= tmp * BITS_PER_NAT * probFragment[j];
// 		}		
// 	}	
// 	
// 	out << "H(Q|B) = " <<   HQB    << "\n"
// 	    << "H(Q)   = " <<    HQ    << "\n"
// 	    << "I(Q;B) = " << HQ - HQB << "\n";
// 	    
// // 	    out << "zxz\n";
// 	
// }