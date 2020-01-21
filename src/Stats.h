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

#ifndef STATS_H_
#define STATS_H_

#include "misc.h"
#include "Hmm.h"

#include <iostream>

class Stats {
public:
	bool on;

	Hmm* hmm;
	map<byte,char> priRevMapping;
	map<byte,char> secRevMapping;

	uint numberOfSequences;
	uint numberOfResidues;

	uint totalCorrectResidues;
	real q3perSequence;
	uint totalCorrectFragments;

	vector< vector<uint> > predictionsPerSecondary;

	real ssLogLikelihood;
	real fragLogLikelihood;

	real numberOfFragments;

	vector<real> partialLogLikelihood;
	vector<real> numberOfPartials;
	
	vector<real> ssLogLikelihood_perResidue;
	vector<real> fragLogLikelihood_perResidue;

	Stats(Hmm* _hmm, map<byte,char>& _priRevMapping, map<byte,char>& _secRevMapping):
		on(false),
		hmm(_hmm),
		priRevMapping(_priRevMapping),
		secRevMapping(_secRevMapping),
		numberOfSequences(0),
		numberOfResidues(0),
		totalCorrectResidues(0),
		q3perSequence(0.0),
		totalCorrectFragments(0),
		predictionsPerSecondary(hmm->numberOfSecondarySymbols, vector<uint>(hmm->numberOfSecondarySymbols, 0) ),
		ssLogLikelihood(0.0),
		fragLogLikelihood(0.0),
		numberOfFragments(0.0),
		partialLogLikelihood(hmm->windowSize-1, 0.0),
		numberOfPartials(hmm->windowSize-1, 0.0)
		{}

	void print(ostream& out, bool calculatePartials);
};

#endif /* STATS_H_ */
