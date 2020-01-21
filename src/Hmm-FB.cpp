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

bool jahfoi = false;

/**
 * Performs the Forward-Backward algorithm.
 * The resulting gamma matrix is stored in the private variable "gamma".
 */
void Hmm::forwardBackward(vector<byte>& seq){
	uint fullLength = seq.size();
	uint windowedLength = fullLength - windowSize + 2;

	vector<real> transProb(numberOfFragments, 0.0);
	
	//'''''''''''''''''''''''''''''''''''''''''
	// 1. Backward (beta)
	//.........................................

	//beta[t][i] is the probability of observing the rest of the sequence, from t+1 on to the end,
	//           given that we are in state i at the t'th symbol. (see Rabiner, 1989)
	vector< vector<real> > beta(fullLength, vector<real>(numberOfStates,0.0));

	// 2.1. Initialization
	if(headAndTail){
		for(ulint i=0; i < numberOfStates; i++)
			beta[seq.size() - windowSize + 1][i] = probState_tail[i];
	}
	else {
		for(ulint i=0; i < numberOfStates; i++)
			beta[seq.size() - windowSize + 1][i] = probState[i];
	}

	// 2.2. Induction
	if(seq.size() > 1)
		for(uint t=0; t < windowedLength-1; t++){
			uint z = seq.size() - windowSize - t;

			for(ulint i=0; i<numberOfFragments; i++)
				transProb[i] =  probFragment[i] * probFragmentEmitsPrimarySymbol[i][ seq[z+halfWindow] ];
			normalize(transProb);
			
			for(ulint prev=0; prev<numberOfStates; prev++){
				beta[z][prev] = 0;

				for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
					ulint trans = (prev * numberOfSecondarySymbols) + ss;
					ulint next = trans % numberOfStates;

					if(probState[next] != 0){
						beta[z][prev] += beta[z+1][next]
						              * transProb[trans]
						              / probState[next];
					}
				}
			}

			normalize(beta[z]);
			
		}

	//'''''''''''''''''''''''''''''''''''''''''
	// 2. Forward (alpha) and gamma
	//.........................................

	//	alpha[t][i] would be the probability of finding state i (secondary fragment) when
	//              observing the t'th symbol of the primary sequence. (also see Rabiner, 1989)
	//
	//  To save RAM, we don't store the entire alpha matrix, as we only need the current
	//  and the previous rows.
	//  (we couldn't do the same with beta because it is calculated backwards).
	//
	vector<real> alphaT(numberOfStates,0.0); //alpha[t]
	vector<real> alphaTminus1(numberOfStates,0.0); //alpha[t-1]

	//2.1. Initialization
	if(headAndTail){
		for(ulint i=0; i < numberOfStates; i++)
			alphaTminus1[i] = probState_head[i];
	}
	else {
		for(ulint i=0; i < numberOfStates; i++)
			alphaTminus1[i] = probState[i];
	}
	normalize(alphaTminus1);

	// gamma
	gamma.clear();
	gamma.resize(windowedLength, vector<real>(numberOfStates,0.0));
	
	for(uint t=0; t < windowedLength; t++){
		
		//2.2. Calculate alpha for this value of t
		if(t > 0){
			for(ulint i=0; i<numberOfFragments; i++)
				transProb[i] = (probFragment[i] * probFragmentEmitsPrimarySymbol[i][ seq[t-1+halfWindow] ]);
			normalize(transProb);

			for(ulint next=0; next<numberOfStates; next++){
				alphaT[next] = 0;

				for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
					ulint trans = next + (ss * numberOfStates);
					ulint prev = trans / numberOfSecondarySymbols;

					if(probState[prev] != 0){
						alphaT[next]
						       += alphaTminus1[prev]
						       * transProb[trans]
						       / probState[prev];
					}
				}
			}

			normalize(alphaT);

			for(ulint i=0; i < numberOfStates; i++)
				alphaTminus1[i] = alphaT[i];
		}
		
		//2.3. Use alphaT and the previoulsy calculated beta[t] to get gamma
		for(ulint i=0; i < numberOfStates; i++){
			if(probState[i] == 0)
				continue;

			if(t==0)
				gamma[t][i] = alphaTminus1[i] * beta[t][i] / probState[i];
			else
				gamma[t][i] = alphaT[i] * beta[t][i] / probState[i];
		}

		normalize(gamma[t]);
	}

	//reduced gamma
	if(numberOfReducedSymbols){
		reducedGamma.clear();
		reducedGamma.resize(windowedLength, vector<real>(numberOfReducedStates,0.0));
		
		for(uint t=0; t < windowedLength; t++)
			for(ulint ri=0; ri < numberOfReducedStates; ri++)
				for(uint j=0; j<fullStatesOfReducedState[ri].size(); j++)
					reducedGamma[t][ri] += gamma[t][ fullStatesOfReducedState[ri][j] ];
	}

	//'''''''''''''''''''''''''''''''''''''''''
	//3. Calculate secondary symbol probabilities
	//.........................................

	probSecondarySymbol.clear();

	//Without reduced mapping
	if(numberOfReducedSymbols == 0){
		probSecondarySymbol.resize(fullLength, vector<real>(numberOfSecondarySymbols, 0.0));
		for(uint a=0; a < windowedLength; a++){
			for(ulint i=0; i < numberOfStates; i++){
				byte sec = indexToSequence(windowSize-1, i)[0];
				probSecondarySymbol[a][sec] += gamma[a][i];
			}
			normalize(probSecondarySymbol[a]);
		}

		for(uint a=windowedLength; a < fullLength; a++){
			for(ulint i=0; i < numberOfStates; i++){
				byte sec = indexToSequence(windowSize-1, i)[a - windowedLength + 1];
				probSecondarySymbol[a][sec] += gamma[windowedLength-1][i];
			}
			normalize(probSecondarySymbol[a]);
		}
	}

	//With reduced mapping
	else{
		probSecondarySymbol.resize(fullLength, vector<real>(numberOfReducedSymbols, 0.0));

		for(uint a=0; a < windowedLength; a++){
			for(ulint i=0; i < numberOfStates; i++){
				byte sec = reduxMapping[ indexToSequence(windowSize-1, i)[0] ];

				probSecondarySymbol[a][sec] += gamma[a][i];
			}
			normalize(probSecondarySymbol[a]);
		}

		for(uint a=windowedLength; a < fullLength; a++){
			for(ulint i=0; i < numberOfStates; i++){
				byte sec = reduxMapping[ indexToSequence(windowSize-1, i)[a - windowedLength + 1] ];

				probSecondarySymbol[a][sec] += gamma[windowedLength-1][i];
			}
			normalize(probSecondarySymbol[a]);
		}
	}
	
	//'''''''''''''''''''''''''''''''''''''''''
	//3. Calculate partial fragment probabilities (for partial log likelihood)
	//.........................................

	if(calculatePartials){
		
		//Without reduced mapping
		probPartialFragment.clear();
		for(uint tam=1; tam<=windowSize; tam++){
			ulint numberOfPartials = (ulint) pow((real)numberOfSecondarySymbols,(int)tam);
			vector< vector<real> > partial(fullLength, vector<real>(numberOfPartials, 0.0));

			for(uint a=0; a < windowedLength; a++){

				for(ulint partialIndex=0; partialIndex < numberOfPartials; partialIndex++){

					ulint numberOfRemaining = (ulint) pow((real)numberOfSecondarySymbols,(int)(windowSize-1-tam));
					ulint baseIndex = partialIndex * numberOfRemaining;

					for(ulint r=0; r<numberOfRemaining; r++){
						ulint fullIndex = baseIndex + r;
						partial[a][partialIndex] += gamma[a][fullIndex];
					}
				}//for i

				normalize(partial[a]);
			}//for a
			
			probPartialFragment.push_back(partial);
		}//for tam
		
		
		//With reduced mapping
		probPartialFragment_redux.clear();
		for(uint tam=1; tam<=windowSize; tam++){

			ulint numberOfPartials = (ulint) pow((real)numberOfReducedSymbols,(int)tam);
			vector< vector<real> > partial(fullLength, vector<real>(numberOfPartials, 0.0));
			
			for(uint a=0; a < windowedLength; a++){
				
				for(ulint partialIndex=0; partialIndex < numberOfPartials; partialIndex++){
					
					ulint numberOfRemaining = (ulint) pow((real)numberOfReducedSymbols,(int)(windowSize-1-tam));
					ulint baseIndex = partialIndex * numberOfRemaining;
					
					for(ulint r=0; r<numberOfRemaining; r++){
						ulint fullIndex = baseIndex + r;
						partial[a][partialIndex] += reducedGamma[a][fullIndex];
					}
				}//for i
				
				normalize(partial[a]);
			}//for a
			
			probPartialFragment_redux.push_back(partial);
		}//for tam
		
	}

}

/**
 * Returns, for the given sequence of primary symbols, a sequence containing
 * the probabilities of all secondary symbols for each position.
 */
vector< vector<real> > Hmm::getProbabilities(vector<byte>& originalSeq, uint r){
	vector<byte> seq = addXXXs( originalSeq, numberOfPrimarySymbols-1 );
	forwardBackward(seq);
	uint fullLength = seq.size();

	if(r != 0)                      //<- this is just to sync getProbabilities(..) with predict(..)
		gsl_ran_poisson(rng1, 1.0); //   when running with a fixed seed

	byte tmax = numberOfReducedSymbols? numberOfReducedSymbols : numberOfSecondarySymbols;

	vector< vector<real> > probs(fullLength - 2*halfWindow, vector<real>(tmax, 0.0));
	for(uint a=halfWindow; a < fullLength - halfWindow; a++)
		probs[a-halfWindow] = probSecondarySymbol[a];
	return probs;
}

/**
 * Returns, for the given sequence of primary symbols, a sequence containing
 * the most probable secondary symbol at each position.
 */
vector<byte> Hmm::predict(vector<byte>& originalSeq, vector< vector<real> >* probsref){
	//TODO: improve this
	vector<byte> empty(0); // <- this is ugly!
	return predict(originalSeq, empty, 0, probsref);
}

vector<byte> Hmm::predict(vector<byte>& originalSeq, vector<byte> &secondarySeq, uint r, vector< vector<real> >* probsref){
	vector<byte> seq = addXXXs( originalSeq, numberOfPrimarySymbols-1 );
	uint fullLength = seq.size();
	forwardBackward(seq);

	vector<uint> totalPartials(windowSize-1, 0);
	uint totalResidues = 0;
	uint correctlyPredicted = 0;

	byte tmax = numberOfReducedSymbols? numberOfReducedSymbols : numberOfSecondarySymbols;

	vector<byte> pred(fullLength - 2*halfWindow);
	
	vector<ulint> observedFragmentIndices(pred.size());
	
	vector<string> datalog(pred.size());
	vector<bool> datalog_validfrag(pred.size(), false);

	int weight = (r==0)? 1 : gsl_ran_poisson(rng1, 1.0);
	
	if(probsref != NULL){
		vector< vector<real> >& probs = *probsref;
		probs.clear();
		probs.resize(fullLength - 2*halfWindow, vector<real>(tmax, 0.0));
		for(uint a=halfWindow; a < fullLength - halfWindow; a++)
			probs[a-halfWindow] = probSecondarySymbol[a];
	}
	
	if(loglikes){
		stats[r].  ssLogLikelihood_perResidue.resize(fullLength - 2*halfWindow);
		stats[r].fragLogLikelihood_perResidue.resize(fullLength - 2*halfWindow);
	}

	if(stats[0].on)
		stats[r].numberOfSequences += weight;
	
	uint col = 0;

	for(uint a=halfWindow; a < fullLength - halfWindow; a++){
		col++;
		
		stringstream logss;
		
		
		real maxProb = -1.0;
		byte maxSymbol = 0;
		for(byte t=0; t <  tmax; t++){
			if(probSecondarySymbol[a][t] > maxProb){
				maxProb = probSecondarySymbol[a][t];
				maxSymbol = t;
			}
		}
		assert(maxProb != -1);
		pred[a-halfWindow] = maxSymbol;

		if(stats[0].on){
			byte pri = originalSeq[a-halfWindow];
			byte obsSec = secondarySeq[a-halfWindow];
			byte predSec = maxSymbol;

			if(pri == priX || obsSec == secX || predSec == secX)
				continue;
			
			totalResidues += weight;

			if(
				(  numberOfReducedSymbols && (reduxMapping[obsSec] == predSec) )
				||
				( !numberOfReducedSymbols && (             obsSec  == predSec) )
			)
					correctlyPredicted += weight;

			stats[r].predictionsPerSecondary[predSec][obsSec] += weight;

			real probSS = numberOfReducedSymbols ?
				probSecondarySymbol[a][reduxMapping[obsSec]] : probSecondarySymbol[a][obsSec];

			//Secondary element log likelihood
			if(probSS != 0.0){
				real ssLL = (real)weight * BITS_PER_NAT * log(probSS);
				stats[r].ssLogLikelihood += ssLL;
				
				if(loglikes){
					stats[r].ssLogLikelihood_perResidue[a-halfWindow] = ssLL;
				}
			}
			
			//Fragment log likelihood
			bool hasX = false;
			vector<byte> observedFragment(windowSize-1);
			for(uint i=0; i<windowSize-1; i++){
				if(a-halfWindow+i >= originalSeq.size())
					observedFragment[i] = secX;
				else
					observedFragment[i] = numberOfReducedSymbols? reduxMapping[secondarySeq[a-halfWindow+i]] : secondarySeq[a-halfWindow+i];

				if(observedFragment[i] == secX){
					hasX = true;
					break;
				}
			}
			
			//xtra info
			if(r == 0){				
				for(uint i=0; i<probSecondarySymbol[a].size(); i++)
					logss << probSecondarySymbol[a][i] << " ";
				
				logss << " " << (numberOfReducedSymbols? REV_DEBUG_000111[obsSec] : REV_DEBUG_01[obsSec]) << " " << REV_DEBUG_01[predSec];
				
				
				if(hasX)
					logss << "\n";
			}
			// /xtra
			
			if(!hasX){				
				stats[r].numberOfFragments += weight;
				
				real gammaValue;
				ulint observedFragmentIndex;
				if(numberOfReducedSymbols){
					observedFragmentIndex = sequenceToIndex(windowSize-1, observedFragment, numberOfReducedSymbols);
					gammaValue = reducedGamma[a][observedFragmentIndex];
				}
				else{
					observedFragmentIndex = sequenceToIndex(windowSize-1, observedFragment, numberOfSecondarySymbols);;
					gammaValue = gamma[a][observedFragmentIndex];
				}
				
				observedFragmentIndices[a-halfWindow] = observedFragmentIndex;
				
				if(gammaValue){
					real fragLL = (real)weight * BITS_PER_NAT * log(gammaValue);
					stats[r].fragLogLikelihood += fragLL;
					
					if(loglikes){
						stats[r].fragLogLikelihood_perResidue[a-halfWindow] = fragLL;
					}
				}
				
				
				datalog_validfrag[a-halfWindow] = true;
			}

			if(calculatePartials){
				//Partial log likelihood
				hasX = false;
				for(uint part=0; part<windowSize-1; part++){

					//Check for X
					for(uint i=0; i<part+1; i++){
						if(observedFragment[i] == secX){
							hasX = true;
							break;
						}
					}

					if(!hasX){
						real partialValue;
						if(numberOfReducedSymbols){
							ulint partialIndex = sequenceToIndex(part+1, observedFragment, numberOfReducedSymbols);
							partialValue = probPartialFragment_redux[part][a][partialIndex];
						}
						else{
							ulint partialIndex = sequenceToIndex(part+1, observedFragment, numberOfSecondarySymbols);
							partialValue = probPartialFragment[part][a][partialIndex];
						}					
						
						if(partialValue)
							stats[r].partialLogLikelihood[part] += (real)weight * BITS_PER_NAT * log(partialValue);

						totalPartials[part] += weight;
					}
				}
			}
		}//if stats
		
		datalog[a-halfWindow] = logss.str();

	}//for a
	
	// log
		for(uint i=0; i<datalog.size()-halfWindow; i++)			
			if(datalog_validfrag[i]){
				vector<byte> predicedFrag(windowSize-1);				

				for(uint j=0; j<windowSize-1; j++)
					predicedFrag[j] = pred[i+j];
				
				ulint predictedIndex = sequenceToIndex(windowSize-1, predicedFrag, numberOfReducedSymbols? numberOfReducedSymbols : numberOfSecondarySymbols);
				
				if(predictedIndex == observedFragmentIndices[i])
					stats[r].totalCorrectFragments += weight;
			}			
			
	if(stats[0].on){
  		for(uint part=0; part<windowSize-1; part++)
			stats[r].numberOfPartials[part] += totalPartials[part];

		stats[r].numberOfResidues += totalResidues;
		stats[r].totalCorrectResidues += correctlyPredicted;
		if(totalResidues != 0)
			stats[r].q3perSequence += weight * (real)correctlyPredicted / (real)totalResidues;
	}

	return pred;
}

void Hmm::printResampledStats(ostream& out){
	vector<real> q3perSeq(numberOfReplicas);
	vector<real> q3perRes(numberOfReplicas, 0.0);
	vector<real> accFragment(numberOfReplicas);
	vector<real> llperSeq(numberOfReplicas);
	vector<real> llperRes(numberOfReplicas);

	for(uint i=1; i<=numberOfReplicas; i++){
		q3perSeq[i-1] = stats[i].q3perSequence / (real)stats[i].numberOfSequences;
		q3perRes[i-1] = (real)stats[i].totalCorrectResidues / (real)stats[i].numberOfResidues;
		accFragment[i-1] = (real)stats[i].totalCorrectFragments / (real)stats[i].numberOfFragments;

		llperRes[i-1] = stats[i].ssLogLikelihood / (real)stats[i].numberOfResidues;
		llperSeq[i-1] = stats[i].fragLogLikelihood / (real)stats[i].numberOfFragments;
	}

	//
	real q3seq_original = stats[0].q3perSequence/(real)stats[0].numberOfSequences;
	real q3seq_bias = stats_mean(q3perSeq) - q3seq_original;

	real q3res_original = (real)stats[0].totalCorrectResidues / (real)stats[0].numberOfResidues;
	real q3res_bias = stats_mean(q3perRes) - q3res_original;
	
	real accfrag_original = (real)stats[0].totalCorrectFragments / (real)stats[0].numberOfFragments;
	real accfrag_bias = stats_mean(accFragment) - accfrag_original;

	//
	real llres_original = stats[0].ssLogLikelihood / (real)stats[0].numberOfResidues;
	real llres_bias = stats_mean(llperRes) - llres_original;

	//
	real llseq_original = stats[0].fragLogLikelihood / (real)stats[0].numberOfFragments;
	real llseq_bias = stats_mean(llperSeq) - llseq_original;

	out 
         << "\n\n# --- Resampled Statistics ---\n"
		 << "\n"
		 << "# Average Q3 per sequence: [ SD, bias corr.]\n"
		 << q3seq_original - q3seq_bias << "\t" << stats_stddev(q3perSeq) << "\t" << q3seq_bias << "\n"
		 << "\n"
	     << "# Average Q3 per residue: [ SD, bias corr.]\n"
		 << q3res_original - q3res_bias << "\t" << stats_stddev(q3perRes) << "\t" << q3res_bias << "\n"
		 << "\n"
	     << "# Average Observation Log Likelihood per Residue [SD, bias correction.]\n"
	     << llres_original - llres_bias << "\t" << stats_stddev(llperRes) << "\t" << llres_bias << "\n"
	     << "\n"
	     << "# Average Observation Log Likelihood per Segment [SD, bias correction.]\n"
	     << llseq_original - llseq_bias << "\t" << stats_stddev(llperSeq) << "\t" << llseq_bias << "\n";

	if(calculatePartials){
		out << "\n# Partial Log Likelihoods [SD, bias correction.]:\n";
		for(uint part=0; part<windowSize-1; part++){

			vector<real> llPart(numberOfReplicas);
			for(uint i=1; i<=numberOfReplicas; i++)
				llPart[i-1] = stats[i].partialLogLikelihood[part] / (real)stats[i].numberOfPartials[part];

			real llpart_original = stats[0].partialLogLikelihood[part] / (real)stats[0].numberOfPartials[part];
			real llpart_bias = stats_mean(llPart) - llpart_original;

			out << (part+1) << " "
			<< llpart_original - llpart_bias << "\t" << stats_stddev(llPart) << "\t" << llpart_bias << "\n";
		}
	}
	
	out  << "\n# Average Fragment Accuracy: [ SD, bias corr.]\n"
	     << accfrag_original - accfrag_bias << "\t" << stats_stddev(accFragment) << "\t" << accfrag_bias << "\n"
	     << "\n"

	     << "\n# --- END HmmPred 1.0 ---\n";

}

