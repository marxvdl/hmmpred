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

#include "Tests.h"
#include "fasta.h"
#include "Hmm.h"

const string TESTDATA_PATH = "/home/marx/Lab/KDevelop/HmmPred/src/test/";

/**
 * Checks if it produces exaclty the same results
 * as second-hmm's on a standard data set.
 */
void Tests::compareWithSecondHMM_1(){
	string trainingFile = TESTDATA_PATH + "miniseclib.training";
	string evalFile = TESTDATA_PATH + "miniseclib.eval";

	TEST_ASSERT( fileExists(trainingFile) );
	TEST_ASSERT( fileExists(evalFile) );

	DoubleFasta training(trainingFile, AMINO_ACIDS, SEC_STRUCT_CK);
	SingleFasta eval(evalFile, AMINO_ACIDS);

	Hmm hmm(myclock(), 50);
	hmm.loglikes = true;
	string redux = "";
	hmm.readTrainingData(training, 9, AMINO_ACIDS_N, SEC_STRUCT_N, AMINO_ACIDS, REV_AMINO_ACIDS, REV_SEC_STRUCT, redux, false);

	list<Sequence*>::iterator evalIt = eval.sequences.begin();
	uint i = 0;
	while(evalIt != eval.sequences.end()){
			vector<byte> predicted = hmm.predict( (*evalIt)->intSeq );
			cout << vec2string(predicted, REV_SEC_STRUCT) << endl;
			TEST_ASSERT( vec2string(predicted, REV_SEC_STRUCT) == results[i++] );
			++evalIt;
	}
}


void Tests::compareWithSecondHMM_2(){
	string trainingFile = TESTDATA_PATH + "miniseclib.training";
	string evalFile = TESTDATA_PATH + "miniseclib.eval+";

	TEST_ASSERT( fileExists(trainingFile) );
	TEST_ASSERT( fileExists(evalFile) );

	DoubleFasta training(trainingFile, AMINO_ACIDS, SEC_STRUCT_CK);
	DoubleFasta eval(evalFile, AMINO_ACIDS, SEC_STRUCT_CK);

	Hmm hmm(myclock(), 50);
	hmm.loglikes = true;
	string redux = "";
	hmm.readTrainingData(training, 9, AMINO_ACIDS_N, SEC_STRUCT_N, AMINO_ACIDS, REV_AMINO_ACIDS, REV_SEC_STRUCT, redux, false);

	list<Sequence*>::iterator evalIt = eval.sequences.begin();
	uint i = 0;
	while(evalIt != eval.sequences.end()){
			vector<byte> predicted = hmm.predict( (*evalIt)->intSeq );
			cout << vec2string(predicted, REV_SEC_STRUCT) << endl;
			TEST_ASSERT( vec2string(predicted, REV_SEC_STRUCT) == results[i++] );
			++evalIt;
	}
}
