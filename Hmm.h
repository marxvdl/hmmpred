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

#include "fasta.h"
#include "mappings.h"
#include "misc.h"
#include "Stats.h"
#include <cmath>
#include <string>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>


#ifndef HMM_H_
#define HMM_H_

class Stats;
struct LineT1;
struct LineT2;

class Hmm {
protected:	
	bool frag_to_seq_init;
	bool frag_to_seq_init_mainhalf;
	bool frag_to_seq_init_xtrahalf;
	void setSeed(ulint seed){
		//Start random number generator
		gsl_rng_env_setup();
		const gsl_rng_type * rng_type = gsl_rng_default;
		rng1 = gsl_rng_alloc(rng_type);
		gsl_rng_set(rng1, seed);
	}
private:	
	// Bootstrapping replicas
	uint numberOfReplicas;
	
	//Adds "XXX" before the sequence.
	vector<byte> addXXXs(vector<byte>& original, int paddingSymbol){
		uint orlen = original.size();
		vector<byte> capped(orlen + 2*halfWindow);
		for(uint i=0; i<orlen; i++)
			capped[i+halfWindow] = original[i];
		for(uint i=0; i<halfWindow; i++)
			capped[i] = capped[orlen+halfWindow+i] = paddingSymbol;
		return capped;
	}



	//Hmm functions
	real probAij(uint i, vector<byte>& vj);

	//
	// Debug functions
	//
	void printSeq(vector<byte>& seq, map<byte,char>& revMap){
		for(uint i=0; i<seq.size(); i++)
			cout << revMap[seq[i]];
		cout << endl;
	}
    
    //
    //insideOutside method, introduced by AFPA
    //
    void insideOutside(vector<byte>& seq);
    void insideOutsideLeftRightTogether(vector<byte>& seq );
    void insideOutsideV1(vector<byte>& seq);
    void insideOutsideV2(vector<byte>& seq);
    void insideOutsideV3(vector<byte>& seq);
    void insideOutsideV4(vector<byte>& seq);
    void insideOutsideV5(vector<byte>& seq);

	void forwardBackward(vector<byte>& seq);
	vector< vector<real> > probSecondarySymbol;

	vector< vector< vector<real> > > probPartialFragment;
	vector< vector< vector<real> > > probPartialFragment_redux;

	vector< vector<real> > gamma;
	vector< vector<real> > reducedGamma;	
	vector< vector<ulint> > fullStatesOfReducedState;

	map<char, byte> reducedAlphabet;
	map<byte, byte> reduxMapping;

	int secX;
	int priX;

	gsl_rng * rng1;
	byte reducedSecX;	
	ulint numberOfMainHalfStates;
	ulint numberOfXtraHalfStates;
	
	//used in vec2string
	map<byte,char> mySecRevMapping;
	map<byte,char> myPriRevMapping;
	map<char,byte> myPriMapping;

public:
	byte numberOfReducedSymbols;
	map<byte, char> rev_reducedAlphabet;
	
	bool calculatePartials;
	bool headAndTail;
	bool loglikes;
    
    //options added by AFPA
    int predGrammar;
    bool seqAlignment;
    bool burialConstraints;
    

	//Statistics
	vector<Stats> stats;
	void printResampledStats(ostream& out);

	//
	// Probability vectors
	//
	vector< vector<real> > probPriGivenSec;	
	vector<real> probFragment;
	vector<real> probState;	
	vector< vector<real> > probFragmentEmitsPrimarySymbol;	
	vector<real> probState_head;	 // Used only if 
	vector<real> probState_tail;    //    headAndTail == true	
	uint mainHalfSize, xtraHalfSize;
	
	vector< vector<real> > mutualEntropyAAFrag;

	// Data
	uint windowSize;
	uint halfWindow;
		
	ulint numberOfFragments;
	ulint numberOfStates;
	ulint numberOfReducedStates;
	uint numberOfPrimarySymbols;
	uint numberOfSecondarySymbols;

	Hmm(uint nreplicas){
		setSeed(0);
		frag_to_seq_init = false;
		frag_to_seq_init_mainhalf = false;
		frag_to_seq_init_xtrahalf = false;
		headAndTail = false;
		loglikes = false;
		numberOfReplicas = nreplicas;
	}
	
	Hmm(ulint seed, uint nreplicas) {
		setSeed(seed);
		frag_to_seq_init = false;
		frag_to_seq_init_mainhalf = false;
		frag_to_seq_init_xtrahalf = false;
		headAndTail = false;
		loglikes = false;
		numberOfReplicas = nreplicas;
	}

	void readTrainingData(DoubleFasta& training, uint myWindowSize, uint nPriSymbols, uint nSecSymbols,
				map<char,byte>& priMapping,	map<byte,char>& priRevMapping,
				map<byte,char>& secRevMapping, string& redux, bool resample );


	//Returns the predicted sequence
	vector<byte> predict(vector<byte>& originalSeq, vector< vector<real> >* probsref = NULL);
	vector<byte> predict(vector<byte>& originalSeq, vector<byte>& secondarySeq, uint r, vector< vector<real> >* probsref = NULL);
	vector< vector<real> > getProbabilities(vector<byte>& originalSeq, uint r);


	//
	// Index <--> Sequence conversions
	//
	vector< vector< vector<byte>  > > frag_seq;
	vector< vector<byte> >   mainhalf_frag_seq;
	vector< vector<byte> >   xtrahalf_frag_seq;
	
	vector<byte> indexToSequence_mainHalf(const ulint f);
	vector<byte> indexToSequence_xtraHalf(const ulint f);

	vector<byte> indexToSequence(const size_t K, const ulint f, uint base);
	ulint sequenceToIndex(const size_t K, vector<byte>& seq, uint base);

	vector<byte> indexToSequence(const size_t K, const ulint f){
			return indexToSequence(K, f, numberOfSecondarySymbols);
	}
	ulint sequenceToIndex(const size_t K, vector<byte>& seq){
			return sequenceToIndex(K, seq, numberOfSecondarySymbols);
	}

	vector<byte> indexToSequence(const ulint f){
			return indexToSequence(windowSize, f, numberOfSecondarySymbols);
	}
	ulint sequenceToIndex(vector<byte>& seq){
			return sequenceToIndex(windowSize, seq, numberOfSecondarySymbols);
	}
	
	vector<byte> indexToSequence_noCache(const size_t K, const ulint f, uint base);
	
	//
	// Print the training probabilities table
	//
	void printTableTransitions(ofstream& out);
	void printTableEmissions(ofstream& out);	
	
// 	void printTableForPrimarySymbol(ofstream& out, byte p);
// 	void printTableV2(ofstream& out);
// 	void calculateEntropyForPrimarySymbol(byte p);
	
	private: string fragment_vec2string(const vector< byte >& vec, uint n);
	         string    state_vec2string(const vector< byte >& vec, uint n);
		 void doPrintTableTransitions(ofstream& out, vector<LineT1>& lines);
		 void doPrintTableEmissions(ofstream& out, vector<LineT2>& lines);
		 real _printTable_HQ;
};


#endif /* HMM_H_ */
