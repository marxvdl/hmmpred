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

#include "Stats.h"
#include "misc.h"

void Stats::print(ostream& out, bool calculatePartials){
	out << "\n"
			<< "# --- Result Statistics ---\n"
			<< "\n"
			<< "# Total Sequence: " << numberOfSequences << "\n"
			<< "# Total Residues: " << numberOfResidues << "\n"
			<< "# Average Q3 per sequence: " << (q3perSequence / (real)numberOfSequences) << "\n"
			<< "# Average Q3 per residue: " << ((real)totalCorrectResidues / (real)numberOfResidues) << "\n"
			<< "# Average Fragment Accuracy: " << ((real)totalCorrectFragments / (real)numberOfFragments) << "\n"
			<< "\n"
			<< "# Observed versus predicted secondary symbols\n"
			<< "# Obs.  Pred->\n"
			;

	for(uint i=0; i<hmm->numberOfSecondarySymbols; i++)
		out << "\t" << secRevMapping[i];
	out << "\n";

	if(hmm->rev_reducedAlphabet.size()){
		for(uint i=0; i<hmm->rev_reducedAlphabet.size(); i++){
			out << hmm->rev_reducedAlphabet[i];
			for(uint j=0; j<hmm->numberOfSecondarySymbols; j++)
				out << "\t" << (int) predictionsPerSecondary[i][j];
			out << "\n";
		}
	}
	else{
		for(uint i=0; i<hmm->numberOfSecondarySymbols; i++){
			out << secRevMapping[i];
			for(uint j=0; j<hmm->numberOfSecondarySymbols; j++)
				out << "\t" << (int) predictionsPerSecondary[i][j];
			out << "\n";
		}
	}

	out << "\n"
			<< "# Average Observation Log Likelihood\n"
			<< "# Per Residue, Per Segment (frag length-1)\n"
			<< (ssLogLikelihood / (real)numberOfResidues)
			<< "\t"
			<< (fragLogLikelihood / numberOfFragments)
			<< "\n";
	
	if(calculatePartials){
		out << "\n# Partial Log Likelihoods:\n";
		for(uint part=0; part < hmm->windowSize-1; part++){
			out << (part+1) << " " << partialLogLikelihood[part] / (real)numberOfPartials[part] << "\n";
		}
	}
	
	if(hmm->loglikes){
		uint len = ssLogLikelihood_perResidue.size();
		
		out << "\n# Residue Log Likelihoods\n";
		for(uint i=0; i<len-1; i++)
			out << ssLogLikelihood_perResidue[i] << " ";
		out << ssLogLikelihood_perResidue[len-1] << "\n"
		
		    << "\n# Fragment Log Likelihoods\n";		
		    for(uint i=0; i<len-1; i++)
			out << fragLogLikelihood_perResidue[i] << " ";
		out << fragLogLikelihood_perResidue[len-1] << " ";
	}
}


