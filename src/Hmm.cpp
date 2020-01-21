/*
 * Copyright 2014, Marx Gomes van der Linden.
 *                 marx@unb.br
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
#include "mappings.h"

#include <cassert>
#include <cmath>
#include <string>

bool performedReduxAlready = false;

void printVec2(vector<byte>& v){
	cout << "{";
	for(uint i=0; i<v.size(); i++)
		cout << " " << REV_SEC_STRUCT[v[i]];
	cout << " }\n  - - - - ' - - -\n ";
	for(uint i=0; i<v.size(); i++)
		cout << " " << i;
	cout << "\n\n";
}

void printVec3(vector<byte>& v){
	for(uint i=0; i<8 && i<v.size(); i++)
		cout << REV_SEC_STRUCT[v[i]];
	
}


/**
 * Initializes everything that the HMM needs:
 * Allocates memory, counts probabilities from the training file and adjusts parameters.
 */
void Hmm::readTrainingData(DoubleFasta& training, uint myWindowSize, uint nPriSymbols, uint nSecSymbols,
		map<char,byte>& priMapping, map<byte,char>& priRevMapping,
		map<byte,char>& secRevMapping, string& redux, bool resample) {

	windowSize = myWindowSize;
	halfWindow = windowSize / 2;
	numberOfPrimarySymbols = nPriSymbols;
	numberOfSecondarySymbols = nSecSymbols - 1;
	numberOfFragments = (ulint) pow((real)numberOfSecondarySymbols, (int)windowSize);
	numberOfStates = numberOfFragments / numberOfSecondarySymbols;
	mySecRevMapping = secRevMapping;
	myPriRevMapping = priRevMapping;
	myPriMapping    = priMapping;
	
	if(headAndTail) {				
		mainHalfSize = (windowSize % 2)? halfWindow : halfWindow - 1;
		xtraHalfSize = windowSize - mainHalfSize - 1;
		
		numberOfMainHalfStates = (ulint) pow((real)numberOfSecondarySymbols, (int)mainHalfSize);
		numberOfXtraHalfStates = (ulint) pow((real)numberOfSecondarySymbols, (int)xtraHalfSize);
	}

	stats.resize(numberOfReplicas+1, Stats(this, priRevMapping, secRevMapping));

	if(!performedReduxAlready){
		numberOfReducedSymbols = 0;
		if(redux != ""){
			if(redux.length() != numberOfSecondarySymbols){
				stringstream ss;
				ss << "Reduced alphabet should have same size as secondary alphabet (" << nSecSymbols << ")";
				throw ss.str().c_str();
			}

			for(uint i=0; i < redux.length(); i++){
				if(!reducedAlphabet.count(redux[i])){
					reducedAlphabet[redux[i]] = numberOfReducedSymbols;
					rev_reducedAlphabet[numberOfReducedSymbols++] = redux[i];
				}
				reduxMapping[i] = reducedAlphabet[redux[i]];
			}

			performedReduxAlready = true;
		}		
				
		reduxMapping[nSecSymbols-1] = nSecSymbols-1; //the lack of this single line cost me a month+ of lost data :(
	}
	
	if(numberOfReducedSymbols)
		numberOfReducedStates = (ulint) pow((real)numberOfReducedSymbols, (int)windowSize-1);

	reducedSecX = numberOfReducedSymbols;
	secX = nSecSymbols-1;
	priX = nPriSymbols-1;

	real primaryPseudocount  = (real) (nPriSymbols - 1);

	//'''''''''''''''''''''''''''''''''''''''''
	// 1. Init counts and probabilities
	//.........................................
	
	vector<real> probHead;
	vector<real> probTail;

	if(resample){
		zeroVector(probFragment);
		zeroVector(probState);
		
		if(headAndTail){
			zeroVector(probHead);
			zeroVector(probTail);
			zeroVector(probState_head);
			zeroVector(probState_tail);
		}
	}
	else{
		probFragment.resize(numberOfFragments, 0.0);
		probState.resize(numberOfStates, 0.0);
		
		if(headAndTail){
			probHead.resize(numberOfMainHalfStates, 0.0);
			probTail.resize(numberOfMainHalfStates, 0.0);
			probState_head.resize(numberOfStates, 0.0);
			probState_tail.resize(numberOfStates, 0.0);
		}
	}
	
	if(resample){		
		zeroVector(probPriGivenSec);
		zeroVector(probFragmentEmitsPrimarySymbol);
	}
	else{		
		probPriGivenSec.resize(numberOfPrimarySymbols, vector<real>(numberOfSecondarySymbols, 0.0));
		probFragmentEmitsPrimarySymbol.resize(numberOfFragments, vector<real>(numberOfPrimarySymbols, 0.0));		
	}


	//'''''''''''''''''''''''''''''''''''''''''
	// 2. Read from file and count occurrences
	//.........................................
	list<Sequence*>::iterator pri = training.sequences.begin();
	list<Sequence*>::iterator sec = training.secondarySequences.begin();

	list<string>::iterator comm = training.comments.begin();

	int weight;

	//Repeat for every pair of sequences
	while(pri != training.sequences.end()){
		weight = resample? gsl_ran_poisson(rng1, 1.0) : 1;
		vector<byte> priSeq = (*pri)->intSeq;
		vector<byte> secSeq = (*sec)->intSeq;
		uint size = priSeq.size();
			
		assert(secSeq.size() == size);

		//Count amino acids and secondary structure segments
		vector<byte> priFrag(windowSize);
		vector<byte> secFrag(windowSize);
		for(uint x = 0; x <= size-windowSize; x++){
			if(priSeq[x] != priX && secSeq[x] != secX)
				probPriGivenSec[priSeq[x]][secSeq[x]] += weight;

			bool validFragment = true;
			bool validState = true;

			for(uint i=0; i<windowSize; i++){
				priFrag[i] = priSeq[i+x];
				secFrag[i] = secSeq[i+x];

				if(priFrag[i] == priX || secFrag[i] == secX){
					validFragment = false;
					if(i != windowSize - 1)
						validState = false;
					break;
				}
			}

			if(validState){
				ulint state = sequenceToIndex(windowSize-1, secFrag);

				probState[state] += weight;
			}

			if(validFragment){
				ulint seg = sequenceToIndex(secFrag);

				probFragment[seg] += weight;
				probFragmentEmitsPrimarySymbol[seg][priFrag[halfWindow]] += weight;
			}
		}
		
		//Read head and tail
		if(headAndTail){
			vector<byte> headFrag(mainHalfSize);
			vector<byte> tailFrag(mainHalfSize);
				
			bool validHead = true;
			bool validTail = true;
			
			for(uint i=0; i<mainHalfSize; i++){
				headFrag[i] = secSeq[i];
				tailFrag[i] = secSeq[size - mainHalfSize + i];
				
				if(headFrag[i] == secX)
					validHead = false;
				
				if(tailFrag[i] == secX)
					validTail = false;
			}			
			
			if(validHead){
				ulint headIndex = sequenceToIndex(mainHalfSize, headFrag);
				probHead[headIndex] += weight;
			}
			
			if(validTail){
				ulint tailIndex = sequenceToIndex(mainHalfSize, tailFrag);
				probTail[tailIndex] += weight;
			}
		}

		//last state
		bool validState = true;
		vector<byte> lastState(windowSize);
		for(uint i=0; i < windowSize-1; i++){
			lastState[i] = secSeq[i + size-windowSize + 1];

			if(lastState[i] == secX){
				validState = false;
				break;
			}
		}
		if(validState){
			uint state = sequenceToIndex(windowSize-1, lastState);

			probState[state] += weight;
		}

		//last probPriGivenSecs
		for(uint x = size-windowSize + 1; x < size; x++)
			if(priSeq[x] != priX && secSeq[x] != secX)
				probPriGivenSec[priSeq[x]][secSeq[x]] += weight;

		pri++; sec++; comm++;		
	}//
	
	

	//'''''''''''''''''''''''''''''''''''''''''
	// 2. Normalize probabilities
	//.........................................

	normalize(probFragment);
	normalize(probState);

	//probPriGivenSec
	for(byte ss=0; ss<numberOfSecondarySymbols; ss++){
		real total = 0.0;
		for(byte aa=0; aa<numberOfPrimarySymbols; aa++)
			total += probPriGivenSec[aa][ss];

		if(total != 0)
			for(byte aa=0; aa<numberOfPrimarySymbols; aa++)
				probPriGivenSec[aa][ss] /= total;
	}

	//probFragmentEmitsPrimarySymbol
	for(ulint seg=0; seg<numberOfFragments; seg++){

		for(byte aa=0; aa<numberOfPrimarySymbols; aa++){
			byte ss = indexToSequence(seg)[halfWindow];

			probFragmentEmitsPrimarySymbol[seg][aa] += primaryPseudocount * probPriGivenSec[aa][ss];
		    //                                              ^-- here comes the cool Pseudocount trick!
		}

		normalize(probFragmentEmitsPrimarySymbol[seg]);
		probFragmentEmitsPrimarySymbol[seg][priX] = 1.0;
	}
	
// 	//'''''''''''''''''''''''''''''''''''''''''
// 	// 3.5. Create and normalize probState_head and probState_tail
// 	//.........................................

	normalize(probHead);
	normalize(probTail);
	
	if(headAndTail){
		//
		// Head
		//
		vector<real> sumProbsAllStatesThatContain(numberOfMainHalfStates);
		for(ulint i=0; i<numberOfStates; i++){
			vector<byte> bigvec = indexToSequence(windowSize-1, i);
			vector<byte> smallvec(&bigvec[halfWindow], &bigvec[windowSize-1]);
			
			sumProbsAllStatesThatContain[ sequenceToIndex(mainHalfSize, smallvec) ] += probState[i];
		}
		
		for(ulint mainHalf=0; mainHalf < numberOfMainHalfStates; mainHalf++){
			vector<byte> mainHalfVec = indexToSequence_mainHalf(mainHalf);
			vector<byte> tmp(windowSize-1, secX);
			
			//populate main half (the one that is counted from the beginning of sequences)
			for(uint p=0; p<mainHalfSize; p++)
				tmp[windowSize-mainHalfSize+p-1] = mainHalfVec[p];
			
			//populate xtra half (the one used to pad the state)
			for(ulint xtrahalf=0; xtrahalf < numberOfXtraHalfStates; xtrahalf++){
				vector<byte> xtraHalfVec = indexToSequence_xtraHalf(xtrahalf);
								
				for(uint p=0; p<xtraHalfSize; p++)
					tmp[p] = xtraHalfVec[p];
				
				ulint mainHalfIndex = sequenceToIndex(mainHalfSize, mainHalfVec);
				ulint fullStateIndex = sequenceToIndex(windowSize-1, tmp);
				
				real probMainHalf = probHead[mainHalfIndex];
				real probFullState = probState[fullStateIndex];
				
				real sum = sumProbsAllStatesThatContain[mainHalfIndex];
				
				if(sum)
					probState_head[fullStateIndex] = probMainHalf * probFullState / sum;
			}
		}
	
// 		assert(is_valid_probability_distribution(probState_head));
		
// 		real sum = 0;
// 		for(ulint i=0; i<numberOfStates; i++){
// 			vector<byte> bigvec = indexToSequence(windowSize-1, i);
// 			vector<byte> smallvec(&bigvec[halfWindow], &bigvec[windowSize-1]);
// 			
// 	// 		if(sequenceToIndex(halfWindow, smallvec) == 46){
// 				if(!probState_head[i]) continue;
// 				printVec3(bigvec); cout << " "; printVec3(smallvec); cout << " " << probState_head[i] << "\n";
// 				sum += probState_head[i];
// 	// 		}		
// 		}
// 		cout << "sum       " << sum << "\n";
		
		//
		// Tail
		//
		zeroVector(sumProbsAllStatesThatContain);
		for(ulint i=0; i<numberOfStates; i++){
			vector<byte> bigvec = indexToSequence(windowSize-1, i);
			vector<byte> smallvec(&bigvec[0], &bigvec[halfWindow]);
		
			sumProbsAllStatesThatContain[ sequenceToIndex(mainHalfSize, smallvec) ] += probState[i];
		}
		
		for(ulint mainHalf=0; mainHalf < numberOfMainHalfStates; mainHalf++){
			vector<byte> mainHalfVec = indexToSequence_mainHalf(mainHalf);
			vector<byte> tmp(windowSize-1, secX);

			//populate main half (the one that is counted from the end of sequences)
			for(uint p=0; p<mainHalfSize; p++)
				tmp[p] = mainHalfVec[p];
				
			//populate xtra half (the one used to pad the state
			for(ulint xtrahalf=0; xtrahalf < numberOfXtraHalfStates; xtrahalf++){
				vector<byte> xtraHalfVec = indexToSequence_xtraHalf(xtrahalf);
				
				for(uint p=0; p<xtraHalfSize; p++)
					tmp[p + mainHalfSize] = xtraHalfVec[p];
					
					ulint mainHalfIndex = sequenceToIndex(mainHalfSize, mainHalfVec);
					ulint fullStateIndex = sequenceToIndex(windowSize-1, tmp);
					
					real probMainHalf = probTail[mainHalfIndex];
					real probFullState = probState[fullStateIndex];
					
					real sum = sumProbsAllStatesThatContain[mainHalfIndex];
				
					if(sum)
						probState_tail[fullStateIndex] = probMainHalf * probFullState / sum;
			}
		}
				
// 		assert(is_valid_probability_distribution(probState_tail));
		
// 		real sum = 0;
// 		for(ulint i=0; i<numberOfStates; i++){
// 			vector<byte> bigvec = indexToSequence(windowSize-1, i);
// 			vector<byte> smallvec(&bigvec[0], &bigvec[halfWindow]);
// 			
// 	// 		if(sequenceToIndex(halfWindow, smallvec) == 46){
// 					if(!probState_tail[i]) continue;
// 				printVec3(bigvec); cout << " "; printVec3(smallvec); cout << " " << probState_tail[i] << "\n";
// 				sum += probState_tail[i];
// 	// 		}		
// 		}
// 		cout << "sum       " << sum << "\n";
		
	}
	

	

	//'''''''''''''''''''''''''''''''''''''''''
	// 3. Deal with 'X'
	//.........................................
	for(uint seg=0; seg<numberOfFragments; seg++)
		probFragmentEmitsPrimarySymbol[seg][priX] = 1;

	for(uint a1=0; a1<numberOfPrimarySymbols - 1; a1++)
		for(uint s1=0; s1<numberOfSecondarySymbols - 1; s1++)
			probPriGivenSec[priX][s1] += probPriGivenSec[a1][s1];

	
	//'''''''''''''''''''''''''''''''''''''''''
	// 4. Build the {reduced states} -> {full states} mapping.
	//.........................................
	if(numberOfReducedSymbols && !fullStatesOfReducedState.size()){		
		fullStatesOfReducedState.resize(numberOfReducedStates, vector<ulint>());
		
		for(ulint i=0; i < numberOfStates; i++){
			vector<byte> state = indexToSequence(windowSize-1, i);
			vector<byte> reducedState(windowSize-1);
			
			for(ulint j=0; j < windowSize-1; j++)
				reducedState[j] = reduxMapping[state[j]];
			
			ulint ri = sequenceToIndex(windowSize-1, reducedState, numberOfReducedSymbols);
			
			fullStatesOfReducedState[ri].push_back(i);
		}
		
	}
}

/**
 * Converts an index into an array of secondary symbols of size smallHalfWindow.
 * This is based on the frag_to_seq function by Gavin E. Crooks on second-hmm.c.
 */
vector<byte> Hmm::indexToSequence_mainHalf(const ulint f) {
	const uint base = numberOfSecondarySymbols;
	
	if(frag_to_seq_init_mainhalf == false) {		
		mainhalf_frag_seq.resize(numberOfMainHalfStates, vector<byte>(mainHalfSize));
		
		for(ulint t=0; t<numberOfMainHalfStates; t++) {
			ulint s= t;
			for(int j=mainHalfSize-1; j>=0; j--) {
				mainhalf_frag_seq[t][j] =  s % base;
				s /= base;
			}
			assert( t == sequenceToIndex(mainHalfSize, mainhalf_frag_seq[t], base) );
		}
		
		frag_to_seq_init_mainhalf = true;
	}
	
	return mainhalf_frag_seq[f];
}

/**
 * Converts an index into an array of secondary symbols of size smallHalfWindow.
 * This is based on the frag_to_seq function by Gavin E. Crooks on second-hmm.c.
 */
vector<byte> Hmm::indexToSequence_xtraHalf(const ulint f) {
	const uint base = numberOfSecondarySymbols;
	
	if(frag_to_seq_init_xtrahalf == false) {		
		xtrahalf_frag_seq.resize(numberOfXtraHalfStates, vector<byte>(xtraHalfSize));
		
		for(ulint t=0; t<numberOfXtraHalfStates; t++) {
			ulint s= t;
			for(int j=xtraHalfSize-1; j>=0; j--) {
				xtrahalf_frag_seq[t][j] =  s % base;
				s /= base;
			}
			assert( t == sequenceToIndex(xtraHalfSize, xtrahalf_frag_seq[t], base) );
		}
		
		frag_to_seq_init_xtrahalf = true;
	}
	
	return xtrahalf_frag_seq[f];
}

/**
 * Converts an index into an array of secondary symbols of size windowSize or windowSize-1.
 * This is based on the frag_to_seq function by Gavin E. Crooks on second-hmm.c.
 */
vector<byte> Hmm::indexToSequence(const size_t K, const ulint f, uint base) {
	int j,k,s;
	ulint t, types;

	//Initialize
	if(frag_to_seq_init==false) {

		//Allocate space to frag_seq
		frag_seq.reserve(2);
		for(uint a = 0; a <= 2; a++){
			vector< vector<byte> > vec;
			vec.reserve(numberOfFragments);
			for(ulint b=0; b<numberOfFragments; b++){
				vector<byte> vec2(windowSize+1, 0);
				vec.push_back(vec2);
			}
			frag_seq.push_back(vec);
		}
		//
		
		frag_to_seq_init = true;

		for(k=windowSize-1; (uint)k<=windowSize;k++) {
			types = (ulint) pow((real)base, k);
			assert(types <=numberOfFragments);

			for(t=0; t<types; t++) {
				s= t;
				for(j=k-1; j>=0; j--) {
					frag_seq[k-windowSize+1][t][j] =  s % base;
					s /= base;
				}
				assert( t == sequenceToIndex(k, frag_seq[k-windowSize+1][t], base) );
			}
		}

	}

	return frag_seq[K - windowSize + 1][f];
}

/**
 * Converts an array of secondary symbols into an index.
 * This is based on the seq_to_frag function by Gavin E. Crooks on second-hmm.c.
 */
ulint Hmm::sequenceToIndex(const size_t K, vector<byte>& seq, uint base) {
	ulint k, f=0;

	for(k=0;k<K;k++) {
		assert(seq[k] < base);
		f = (f*base) + seq[k];
	}
	return f;
}
