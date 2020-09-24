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

#include "mappings.h"

#include <cassert>
#include <cmath>


/**
 * Performs the Inside-Outside algorithm.
 * The resulting gamma matrix is stored in the private variable "gamma".
 */

//Marx, acho que fica mais fácil se você ler lado a lado com o arquivo Hmm-FB.cpp
//do qual este foi derivado
//Ainda vou fazer as funções isEmittedFromLeft e isEmittedFromRight,
//Que definem a gramática livre de context dependente da estrutura secundária.
//Mas já dá para entender a ideia geral. (FEITO, matrizes isLeftProduced and isRightProduced)

void Hmm::insideOutside(vector<byte>& seq){
   	uint fullLength = seq.size();
	uint windowedLength = fullLength - windowSize + 2;

	vector<real> transProb(numberOfFragments, 0.0);
    
    
    // SScpair is intended to be the pair of secondary structure symbols
    // at the two central positions of the State fragment.
    // Secondary structure symbols (0 for E, 1 for H, 2 for C)
    // are obtained from SecondarySymbol%3,
    // and therefore SecondarySymbols must be constructed accordingly (AFPA)
    
    uint SScpair[numberOfStates];
    ulint number1 = numberOfXtraHalfStates * numberOfSecondarySymbols;
    ulint number2 = numberOfXtraHalfStates;
    ulint number3 = numberOfXtraHalfStates / numberOfSecondarySymbols;
    
    for(uint i = 0; i < numberOfStates; i++){
        uint SScpair1 = i % number1 / number2;
        SScpair1 %= 3;
        uint SScpair2 = i % number2 / number3;
        SScpair2  %=3;
        SScpair[i] = 3 * SScpair1 + SScpair2;
    }
    
    //The Stochastic Context-free grammar is defined by these two 9x9 bool matrices
    //that determine which pairs of States can be produced from the left and/or from the right.
    //Each pair of states acts as a non-terminal that produces other non-terminals
    //while emitting terminals (PrimarySymbols) at its left and/or right side.
    
    bool isLeftProduced[9][9] = {{1, 0, 0, 1, 0, 0, 1, 0, 0},
                                 {0, 0, 0, 1, 0, 0, 0, 0, 0},
                                 {0, 0, 0, 1, 0, 0, 1, 0, 0},
                                 {0, 1, 1, 0, 0, 0, 0, 0, 0},
                                 {0, 1, 1, 0, 1, 0, 0, 1, 0},
                                 {0, 1, 1, 0, 0, 0, 0, 1, 0},
                                 {0, 1, 1, 0, 0, 0, 0, 0, 0},
                                 {0, 1, 1, 1, 0, 1, 0, 0, 0},
                                 {0, 1, 1, 1, 0, 1, 0, 0, 1}};
                                    
    bool isRightProduced[9][9] = {{1, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {1, 0, 0, 1, 0, 0, 1, 0, 0},
                                   {1, 0, 0, 1, 0, 0, 1, 0, 0},
                                   {0, 1, 1, 0, 1, 1, 0, 1, 1},
                                   {0, 1, 1, 0, 1, 0, 0, 1, 0},
                                   {0, 1, 1, 0, 1, 0, 0, 1, 0},
                                   {0, 1, 1, 0, 1, 1, 0, 1, 1},
                                   {0, 0, 0, 1, 0, 1, 0, 0, 1},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 1}};
    
    
    //''''''''''''''''''''''''''''''''''''''''''''''
    // 1. Outside
    //''''''''''''''''''''''''''''''''''''''''''''''
    
    // gammaLeft[t][v][i] X gammaRight[t][v][j] should be the
    // probability of state i at position t and state j at position v
    // given the outside sequence, i.e., excluding the sequence fragment from t to v
    
 //   vector< vector< vector<real> > > gammaLeft(windowedLength,
 //                                               vector< vector<real> > (windowedLength, vector<real> (numberOfStates, 1.0)));
 //   vector< vector< vector<real> > > gammaRight(windowedLength,
 //                                               vector< vector<real> > (windowedLength, vector<real> (numberOfStates, 1.0)));
    
    vector< vector < vector < vector<real> > > > gammaLeftRight(windowedLength,
                                                                vector< vector < vector<real> > > (windowedLength, vector < vector<real> > (numberOfStates, vector <real> (numberOfStates, 0.0))));
    //
    //Talvez tenha que colocar em uma mesma matriz bem maior gamaLeftRight
    //para normalizar corretamente !!!
    //
    
    vector<real> transLeftProb(numberOfFragments, 0.0);
    vector<real> transRightProb(numberOfFragments, 0.0);
    
    // 1.1 Initialization
    if(headAndTail)
        for(ulint i=0; i < numberOfStates; i++)
            for(ulint j=0; j < numberOfStates; j++)
                gammaLeftRight[0][windowedLength-1][i][j] = probState_head[i] * probState_tail[j];
    else
        for(ulint i=0; i < numberOfStates; i++)
            for(ulint j=0; j < numberOfStates; j++)
                gammaLeftRight[0][windowedLength-1][i][j] = 1.0;
        
    // 1.2 Induction
    for(int t=0; t < windowedLength; t++){
        if(t > 0) {
            for(ulint i=0; i<numberOfFragments; i++)
                    transLeftProb[i] = probFragment[i] *
                        probFragmentEmitsPrimarySymbol[i][ seq[t-1+halfWindow] ];
            normalize(transLeftProb);
        }
        for (int v=windowedLength-1; v >= t; v--) {
 //           assert(v >= 0);
//            if(t==0 && v==windowedLength-1) continue;
            if(v < windowedLength - 1) {
                for(ulint i=0; i<numberOfFragments; i++)
                    transRightProb[i] = probFragment[i] *
                            probFragmentEmitsPrimarySymbol[i][ seq[v+1+halfWindow] ];
                    normalize(transRightProb);
            }
            for(ulint nextleft=0; nextleft<numberOfStates; nextleft++){
                for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                    ulint transleft = nextleft + (ss * numberOfStates);
                    ulint prevleft = transleft / numberOfSecondarySymbols;
                    
                    for(ulint prevright=0; prevright<numberOfStates; prevright++){
                        for(uint tt=0; tt<numberOfSecondarySymbols; tt++){
                            ulint transright = (prevright * numberOfSecondarySymbols) + tt;
                            ulint nextright = transright % numberOfStates;
                                                                                   
                            if( t > 0 && v < windowedLength - 1 &&
                                isLeftProduced[SScpair[nextleft]][SScpair[prevright]] &&
                                probState[prevleft] != 0 &&
                                isRightProduced[SScpair[nextleft]][SScpair[prevright]] &&
                                probState[nextright] != 0){
                                real probleftright = gammaLeftRight[t-1][v+1][prevleft][nextright] * transLeftProb[transleft] / probState[prevleft] * transRightProb[transright] / probState[nextright];
                                gammaLeftRight[t][v][nextleft][prevright] += probleftright;
                            }
                                                                                                
                            else if(t > 0 &&
                                isLeftProduced[SScpair[nextleft]][SScpair[prevright]] &&
                                        probState[prevleft] != 0){
                                    real probleft = gammaLeftRight[t-1][v][prevleft][prevright] * transLeftProb[transleft] / probState[prevleft];
                                    gammaLeftRight[t][v][nextleft][prevright] += probleft;
                                }
            
                            else if(v < windowedLength - 1 &&
                                isRightProduced[SScpair[nextleft]][SScpair[prevright]] &&
                                        probState[nextright] != 0){
                                    real probright = gammaLeftRight[t][v+1][nextleft][nextright] * transRightProb[transright] / probState[nextright];
                                gammaLeftRight[t][v][nextleft][prevright] += probright;
                            }
                            
                   
                                else ;
                        }
                    }
                    
                }
            }
        
            normalize(gammaLeftRight[t][v]);
            
    //        normalize(gammaLeft[t][v]);
    //        normalize(gammaRight[t][v]);
        }
    }
  
    //''''''''''''''''''''''''''''''''''''''''''''''
    // 2. Inside and gamma
    //''''''''''''''''''''''''''''''''''''''''''''''
    

    // gamma
    gamma.clear();
    gamma.resize(windowedLength, vector<real>(numberOfStates,0.0));
    
    for(uint t=0; t < windowedLength; t++){
        
        //Inside
        for(ulint i=0; i<numberOfFragments; i++)
                transLeftProb[i] = probFragment[i] *
                    probFragmentEmitsPrimarySymbol[i][ seq[t+halfWindow] ];
        normalize(transLeftProb);
        
        //gamma
        for(ulint next=0; next<numberOfStates; next++)
            for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                ulint trans = next + (ss * numberOfStates);
                ulint prev = trans / numberOfSecondarySymbols;
                
                if(probState[prev] == 0||probState[next] == 0)
                    continue;
         
    //            gamma[t][next] = gammaLeft[t][t][prev] * transLeftProb[trans] / probState[prev]
    //                            * gammaRight[t][t][next] / probState[next];
    
                gamma[t][next] = gammaLeftRight[t][t][prev][next] /probState[prev] / probState[next];
                
            }
        
        normalize(gamma[t]);
    }
    
    
  //'''''''''''''''''''''''''''''''''''''''''''''''''
  //Daqui em diante está igual ao Hmm-FB.cpp
  //'''''''''''''''''''''''''''''''''''''''''''''''''
    
    
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

