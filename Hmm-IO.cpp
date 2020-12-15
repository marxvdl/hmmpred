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
    
/*
 
 // this grammar allows adjacent E and H seconsary structure symbols
 
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
                                   {1, 0, 0, 0, 0, 0, 1, 0, 0},
                                   {0, 1, 1, 0, 1, 1, 0, 1, 1},
                                   {0, 0, 0, 0, 1, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 1, 0, 0, 1, 0},
                                   {0, 1, 1, 0, 1, 1, 0, 1, 1},
                                   {0, 0, 0, 1, 0, 1, 0, 0, 1},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 1}};
  */
 
    // this grammar does not allow E and H secondary structure symbols
    
    bool isLeftProduced[9][9] = {{1, 0, 0, 0, 0, 0, 1, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 1, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 1, 0, 1, 0, 0, 1, 0},
                                 {0, 0, 1, 0, 0, 0, 0, 1, 0},
                                 {0, 0, 1, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 1, 0, 0, 1, 0, 0, 0},
                                 {0, 0, 1, 0, 0, 1, 0, 0, 1}};
                                    
    bool isRightProduced[9][9] = {{1, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {1, 0, 0, 0, 0, 0, 1, 0, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 1, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 1, 0, 0, 1, 0},
                                   {0, 0, 1, 0, 1, 1, 0, 1, 1},
                                   {0, 0, 0, 0, 0, 1, 0, 0, 1},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 1}};
    
    
    //''''''''''''''''''''''''''''''''''''''''''''''
    // 1. Outside
    //''''''''''''''''''''''''''''''''''''''''''''''
    
    // gammaLeft[t][v][i] X gammaRight[t][v][j] should be the
    // probability of state i at position t and state j at position v
    // given the outside sequence, i.e., excluding the sequence fragment from t to v
    
    vector< vector< vector<real> > > gammaLeft(windowedLength,
                                                vector< vector<real> > (windowedLength, vector<real> (numberOfStates, 0.0)));
    vector< vector< vector<real> > > gammaRight(windowedLength,
                                                vector< vector<real> > (windowedLength, vector<real> (numberOfStates, 0.0)));
    
 //   vector< vector < vector < vector<real> > > > gammaLeftRight(windowedLength,
 //                                                               vector< vector < vector<real> > > (windowedLength, vector < vector<real> > (numberOfStates, vector <real> (numberOfStates, 0.0))));
  
    
 //   ulint numberOfDoubleStates = numberOfStates * numberOfStates;
 //   vector < vector < vector<real> > > gammaLeftRight(windowedLength, vector < vector<real> > (windowedLength, vector <real> (numberOfDoubleStates, 0.0)));
    
    
    //
    //Talvez tenha que colocar em uma mesma matriz bem maior gamaLeftRight
    //para normalizar corretamente !!!
    //
    
    vector<real> transLeftProb(numberOfFragments, 0.0);
    vector<real> transRightProb(numberOfFragments, 0.0);
    
    //
    // 1.1 Initialization
    //
    // Note that only pair of states that could actually be generated by the grammar
    // will have non zero probability at the beginning and end of the chain !!!!
    //
    
    if(headAndTail){
        for(ulint i=0; i < numberOfStates; i++)
            for(ulint j=0; j < numberOfStates; j++)
                if(isLeftProduced[SScpair[i]][SScpair[j]] ||
                   isRightProduced[SScpair[i]][SScpair[j]])
                {
                 //       real prod = sqrt(probState_head[i] * probState_tail[j]);
                    real prod = probState_head[i] * probState_tail[j];
                    gammaLeft[0][windowedLength-1][i] += prod;
                    gammaRight[0][windowedLength-1][j] += prod;
                }
    }
    else{
        for(ulint i=0; i < numberOfStates; i++)
            for(ulint j=0; j < numberOfStates; j++)
                if(isLeftProduced[SScpair[i]][SScpair[j]] ||
                   isRightProduced[SScpair[i]][SScpair[j]])
                {
                        gammaLeft[0][windowedLength-1][i] = 1.0;
                        gammaRight[0][windowedLength-1][i] = 1.0;
                }
    }
    
    normalize(gammaLeft[0][windowedLength-1]);
    normalize(gammaRight[0][windowedLength-1]);
    
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
            if(t==0 && v==windowedLength-1) continue;
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
                            
                     //       ulint doubleindex = nextleft * numberOfStates + prevright;
                                                                                   
                            if( t > 0 && v < windowedLength - 1 &&
                                isLeftProduced[SScpair[nextleft]][SScpair[prevright]] &&
                                probState[prevleft] != 0 &&
                                isRightProduced[SScpair[nextleft]][SScpair[prevright]] &&
                                probState[nextright] != 0)
                            {
                            //    ulint parentdoubleindex = prevleft*numberOfStates + nextright;
                                real probleftright = gammaLeft[t-1][v+1][prevleft] *
                                    gammaRight[t-1][v+1][nextright] *
                                    transLeftProb[transleft] / probState[prevleft] *
                                    transRightProb[transright] / probState[nextright];
                                gammaLeft[t][v][nextleft] += probleftright;
                                gammaRight[t][v][prevright] += probleftright;
                            
                            //    gammaLeft[t][v][nextleft] += sqrt(probleftright);
                            //    gammaRight[t][v][prevright] += sqrt(probleftright);
                            }
                                                                                                
                            else if(t > 0 && probState[prevleft] != 0 &&
                                    isLeftProduced[SScpair[nextleft]][SScpair[prevright]])
                            {
                            //    ulint parentdoubleindex = prevleft*numberOfStates + prevright;
                                real probleft = gammaLeft[t-1][v][prevleft] *
                                    gammaRight[t-1][v][prevright] *
                                    transLeftProb[transleft] / probState[prevleft];
                                gammaLeft[t][v][nextleft] += probleft;
                                gammaRight[t][v][prevright] += probleft;
                                
                             //   gammaLeft[t][v][nextleft] += sqrt(probleft);
                             //   gammaRight[t][v][prevright] += sqrt(probleft);
                                
                              
                                
                            }
            
                            else if(v < windowedLength - 1 && probState[nextright] != 0 &&
                                    isRightProduced[SScpair[nextleft]][SScpair[prevright]])
                            {
                             //   ulint parentdoubleindex = nextleft*numberOfStates + nextright;
                                real probright = gammaLeft[t][v+1][nextleft] *
                                    gammaRight[t][v+1][nextright] *
                                    transRightProb[transright] / probState[nextright];
                                gammaLeft[t][v][nextleft] += probright;
                                gammaRight[t][v][prevright] += probright;
                                
                             //   gammaLeft[t][v][nextleft] += sqrt(probright);
                             //   gammaRight[t][v][prevright] += sqrt(probright);
                            }
                            
                            else ;
                        }
                    }
                    
                }
            }
        
    //       normalize(gammaLeftRight[t][v]);
            
            normalize(gammaLeft[t][v]);
            normalize(gammaRight[t][v]);
            
    //        normalize(gammaLeft[t][v],gammaRight[t][v]);
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
        
        for(ulint prev=0; prev<numberOfStates; prev++){
            
            if(probState[prev] == 0) continue;
            real probStateEmitsSeq = 0.0;
            
            for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                ulint trans = (prev * numberOfSecondarySymbols) + ss;
                probStateEmitsSeq += transLeftProb[trans];
            }
                
            probStateEmitsSeq /= probState[prev];
                        
            // ulint doubleindex = prev * numberOfStates + prev;
            gamma[t][prev] = gammaLeft[t][t][prev] * gammaRight[t][t][prev] *
                probStateEmitsSeq;
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



/**
 * Performs the Inside-Outside algorithm.
 * The resulting gamma matrix is stored in the private variable "gamma".
 */

//Marx, acho que fica mais fácil se você ler lado a lado com o arquivo Hmm-FB.cpp
//do qual este foi derivado
//Ainda vou fazer as funções isEmittedFromLeft e isEmittedFromRight,
//Que definem a gramática livre de context dependente da estrutura secundária.
//Mas já dá para entender a ideia geral. (FEITO, matrizes isLeftProduced and isRightProduced)

void Hmm::insideOutsideV1(vector<byte>& seq){
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
    
   // C-doubleState-C emits at left, H-doubleState-H and E-doubleState-E emit at left and right
    
   /*
    

    bool isLeftProduced[9][9] = {{0, 0, 0, 0, 0, 0, 1, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 1, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 1, 0, 0, 0, 0, 1, 0},
                                 {0, 0, 1, 0, 0, 0, 0, 1, 0},
                                 {0, 0, 1, 0, 0, 1, 1, 1, 0},
                                 {0, 0, 1, 0, 0, 1, 1, 1, 0},
                                 {0, 0, 1, 0, 0, 1, 1, 1, 0}};
    
    
    bool isLeftRightProduced[9][9] ={{1, 0, 0, 0, 0, 0, 1, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {1, 0, 0, 0, 0, 0, 1, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 1, 0, 0, 1, 0},
                                     {0, 0, 0, 0, 1, 0, 0, 1, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0}};
    
                                    
    bool isRightProduced[9][9] = {{0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {1, 0, 0, 0, 0, 0, 1, 0, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 1, 0, 0, 1, 0},
                                   {0, 0, 1, 0, 1, 1, 0, 1, 1},
                                   {0, 0, 1, 0, 0, 1, 0, 0, 1},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0}};
    
    */
  
    
  
  //  E-doubleState-E, H-doubleState-H and C-doubleSttate-C emit at left and right
   
    bool isLeftProduced[9][9] = {{0, 0, 0, 0, 0, 0, 1, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 1, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 1, 0, 0, 0, 0, 1, 0},
                                 {0, 0, 1, 0, 0, 0, 0, 1, 0},
                                 {0, 0, 1, 0, 0, 1, 0, 0, 0},
                                 {0, 0, 1, 0, 0, 1, 0, 0, 0},
                                 {0, 0, 1, 0, 0, 1, 0, 0, 0}};
    
    
    bool isLeftRightProduced[9][9] ={{1, 0, 0, 0, 0, 0, 1, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {1, 0, 0, 0, 0, 0, 1, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 1, 0, 0, 1, 0},
                                     {0, 0, 0, 0, 1, 0, 0, 1, 0},
                                     {0, 0, 1, 0, 0, 1, 0, 0, 1},
                                     {0, 0, 1, 0, 0, 1, 0, 0, 1},
                                     {0, 0, 1, 0, 0, 1, 0, 0, 1}};
    
                                    
    bool isRightProduced[9][9] = {{0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {1, 0, 0, 0, 0, 0, 1, 0, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 1, 0, 0, 1, 0},
                                   {0, 0, 1, 0, 1, 1, 0, 1, 1},
                                   {0, 0, 1, 0, 0, 1, 0, 0, 1},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0}};
    
    
    /*
    
    bool isLeftProduced[9][9] = {{1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1}};
    
    
    bool isLeftRightProduced[9][9] ={{0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0}};
    
                                    
    bool isRightProduced[9][9] = {{1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1}};
    
    */
   
    
    /*
    bool isLeftProduced[9][9] = {{1, 0, 1, 0, 0, 0, 1, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 1, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 1, 0, 1, 1, 0, 1, 0},
                                 {0, 0, 1, 0, 0, 0, 1, 1, 0},
                                 {0, 0, 1, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 1, 0, 0, 1, 0, 0, 0},
                                 {0, 0, 1, 0, 0, 1, 1, 0, 1}};
                                    
    bool isRightProduced[9][9] = {{1, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {1, 0, 0, 0, 0, 0, 1, 1, 1},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 1, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 1, 0, 0, 1, 1},
                                   {1, 0, 1, 0, 1, 1, 0, 1, 1},
                                   {0, 0, 0, 0, 1, 1, 0, 0, 1},
                                   {0, 0, 0, 0, 0, 0, 0, 1, 1}};
    */
     
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
    
 //   vector< vector < vector < vector<real> > > > gammaLeftRight(windowedLength,
 //                                                               vector< vector < vector<real> > > (windowedLength, vector < vector<real> > (numberOfStates, vector <real> (numberOfStates, 0.0))));
  
    
    ulint numberOfDoubleStates = numberOfStates * numberOfStates;
    vector < vector < vector<real> > > gammaLeftRight(windowedLength, vector < vector<real> > (windowedLength, vector <real> (numberOfDoubleStates, 0.0)));
    
    
    //
    //Talvez tenha que colocar em uma mesma matriz bem maior gamaLeftRight
    //para normalizar corretamente !!!
    //
    
    vector<real> transLeftProb(numberOfFragments, 0.0);
    vector<real> transRightProb(numberOfFragments, 0.0);
    vector<real> probDoubleState(numberOfDoubleStates, 0.0);
    
    for(ulint i=0; i < numberOfStates; i++)
        for(ulint j=0; j < numberOfStates; j++)
            if(isLeftProduced[SScpair[i]][SScpair[j]] ||
               isRightProduced[SScpair[i]][SScpair[j]] ||
               isLeftRightProduced[SScpair[i]][SScpair[j]])
            {
                ulint doubleindex = i*numberOfStates +j;
                probDoubleState[doubleindex] = probState[i]*probState[j];
            }
    normalize(probDoubleState);
    
    
    //
    // 1.1 Initialization
    //
    // Note that only pair of states that could actually be generated by the grammar
    // will have non zero probability at the beginning and end of the chain !!!!
    //
    
    if(headAndTail){
        for(ulint i=0; i < numberOfStates; i++)
            for(ulint j=0; j < numberOfStates; j++)
   //             if(isLeftProduced[SScpair[i]][SScpair[j]] ||
   //                isRightProduced[SScpair[i]][SScpair[j]] ||
   //                isLeftRightProduced[SScpair[i]][SScpair[j]])
                        gammaLeftRight[0][windowedLength-1][i*numberOfStates + j] =
                            probState_head[i] * probState_tail[j];
    }
    else{
        for(ulint i=0; i < numberOfStates; i++)
            for(ulint j=0; j < numberOfStates; j++)
                if(isLeftProduced[SScpair[i]][SScpair[j]] &&
                   isRightProduced[SScpair[i]][SScpair[j]])
                        gammaLeftRight[0][windowedLength-1][i*numberOfStates + j] = 1.0;
    }
    normalize(gammaLeftRight[0][windowedLength-1]);
    
    // 1.2 Induction
//    for(int t=0; t < windowedLength; t++){
    
    for(int u = windowedLength - 2; u >= 0; u-- ){
        for(int v = u; v < windowedLength; v++){
            
            int t = v - u;
        
          
            if(t > 0) {
                for(ulint i=0; i<numberOfFragments; i++)
                    transLeftProb[i] = probFragment[i] *
                        probFragmentEmitsPrimarySymbol[i][ seq[t-1+halfWindow] ];
                normalize(transLeftProb);
            }
    
 //       for (int v=windowedLength-2; v >= t; v--) {
 //           assert(v >= 0);
        
 //       for(uint u=0; u < windowedLength-1-t; u++) {
 //           uint v = windowedLength - 2 - u;
 //           uint v = seq.size() - windowSize - u;
 //           if(t==0 && v==windowedLength-1)
 //               continue;
            
            if(v < windowedLength - 1) {
                for(ulint i=0; i<numberOfFragments; i++)
                    transRightProb[i] = probFragment[i] *
                            probFragmentEmitsPrimarySymbol[i][ seq[v+halfWindow] ];
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
                            
                            ulint doubleindex = nextleft * numberOfStates + prevright;
                                                                                   
                            if(isLeftRightProduced[SScpair[nextleft]][SScpair[prevright]])
                            {
                                ulint parentdoubleindex = prevleft*numberOfStates + nextright;
                                if(t > 0 && v < (windowedLength - 1) &&
                                   probState[prevleft] != 0 && probState[nextright] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                                 //   ulint parentdoubleindex = prevleft*numberOfStates + nextright;
                                    real probleftright = gammaLeftRight[t-1][v+1][parentdoubleindex]*
                                        transLeftProb[transleft] * transRightProb[transright] /
                                        probState[prevleft] / probState[nextright];
                                        //probDoubleState[parentdoubleindex];
                                    gammaLeftRight[t][v][doubleindex] += probleftright;
                                }
                            }
                                                                                                
                            if(isLeftProduced[SScpair[nextleft]][SScpair[prevright]]) //||
                             //  isRightProduced[SScpair[nextleft]][SScpair[prevleft]] )
                             //  v == windowedLength - 1)
                            {
                                ulint parentdoubleindex = prevleft*numberOfStates + prevright;
                                if( t > 0 && probState[prevleft] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                              //      ulint parentdoubleindex = prevleft*numberOfStates + prevright;
                                    real probleft = gammaLeftRight[t-1][v][parentdoubleindex] *
                                        transLeftProb[transleft] /
                                        probState[prevleft];
                                        //probDoubleState[parentdoubleindex];
                                    gammaLeftRight[t][v][doubleindex] += probleft;
                                }
                            }
                            
                            if(isRightProduced[SScpair[nextleft]][SScpair[prevright]]) // ||
                            //   isLeftProduced[SScpair[prevright]][SScpair[nextright]])
                            //   t == 0)
                            {
                                ulint parentdoubleindex = nextleft*numberOfStates + nextright;
                                if(v < windowedLength - 1 && probState[nextright] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                                 //   ulint parentdoubleindex = nextleft*numberOfStates + nextright;
                                    real probright = gammaLeftRight[t][v+1][parentdoubleindex] *
                                        transRightProb[transright] /
                                        probState[nextright];
                                        //probDoubleState[parentdoubleindex];
                                    gammaLeftRight[t][v][doubleindex] += probright;
                                }
                            }
                            
                            if(//isLeftProduced[SScpair[nextleft]][SScpair[prevright]] ||
                               isRightProduced[SScpair[prevright]][SScpair[nextright]] )
                             //  v == windowedLength - 1)
                            {
                                ulint parentdoubleindex = prevleft*numberOfStates + prevright;
                                if( t > 0 && probState[prevleft] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                              //      ulint parentdoubleindex = prevleft*numberOfStates + prevright;
                                    real probleft = gammaLeftRight[t-1][v][parentdoubleindex] *
                                        transLeftProb[transleft] /
                                        probState[prevleft];
                                        //probDoubleState[parentdoubleindex];
                                    gammaLeftRight[t][v][doubleindex] += probleft;
                                }
                            }
                            
                            if(//isRightProduced[SScpair[nextleft]][SScpair[prevright]] ||
                               isLeftProduced[SScpair[prevleft]][SScpair[nextleft]])
                            //   t == 0)
                            {
                                ulint parentdoubleindex = nextleft*numberOfStates + nextright;
                                if(v < windowedLength - 1 && probState[nextright] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                                 //   ulint parentdoubleindex = nextleft*numberOfStates + nextright;
                                    real probright = gammaLeftRight[t][v+1][parentdoubleindex] *
                                        transRightProb[transright] /
                                        probState[nextright];
                                        //probDoubleState[parentdoubleindex];
                                    gammaLeftRight[t][v][doubleindex] += probright;
                                }
                            }
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
 
        
        for(ulint i=0; i<numberOfStates; i++){
            if(probState[i] == 0)
                continue;
    
            gamma[t][i] = gammaLeftRight[t][t][i*numberOfStates  + i] / probState[i];
        }

 /*       for(ulint i=0; i < numberOfFragments; i++)
            transProb[i] = probFragment[i] *
                probFragmentEmitsPrimarySymbol[i][ seq[t+halfWindow] ];
            normalize(transProb);
        
         for(ulint prev=0; prev<numberOfStates; prev++){
            
            if(probState[prev] == 0) continue;
            real probStateEmitsSeq = 0.0;
            
            for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                ulint trans =  prev * numberOfSecondarySymbols + ss;
                ulint next = trans % numberOfStates;
                
                if(probState[next] == 0) continue;
                
                ulint doubleindex = prev*numberOfStates + next;
                
                probStateEmitsSeq += transProb[trans] *
                                    gammaLeftRight[t][t+1][doubleindex] /
                                    probState[next];
            }
             
            probStateEmitsSeq /= probState[prev];
                  
            gamma[t][prev] = probStateEmitsSeq;
            
        }
  */
    
 
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

void Hmm::insideOutsideV2(vector<byte>& seq){
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
    
   // C-doubleState-C emits at left, H-doubleState-H and E-doubleState-E emit at left and right
    
   /*
    

    bool isLeftProduced[9][9] = {{0, 0, 0, 0, 0, 0, 1, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 1, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 1, 0, 0, 0, 0, 1, 0},
                                 {0, 0, 1, 0, 0, 0, 0, 1, 0},
                                 {0, 0, 1, 0, 0, 1, 1, 1, 0},
                                 {0, 0, 1, 0, 0, 1, 1, 1, 0},
                                 {0, 0, 1, 0, 0, 1, 1, 1, 0}};
    
    
    bool isLeftRightProduced[9][9] ={{1, 0, 0, 0, 0, 0, 1, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {1, 0, 0, 0, 0, 0, 1, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 1, 0, 0, 1, 0},
                                     {0, 0, 0, 0, 1, 0, 0, 1, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0}};
    
                                    
    bool isRightProduced[9][9] = {{0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {1, 0, 0, 0, 0, 0, 1, 0, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 1, 0, 0, 1, 0},
                                   {0, 0, 1, 0, 1, 1, 0, 1, 1},
                                   {0, 0, 1, 0, 0, 1, 0, 0, 1},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0}};
    
    */
  
    
  
  //  E-doubleState-E, H-doubleState-H and C-doubleSttate-C emit at left and right
   
    bool isLeftProduced[9][9] = {{0, 0, 0, 0, 0, 0, 1, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 1, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 1, 0, 0, 0, 0, 1, 0},
                                 {0, 0, 1, 0, 0, 0, 0, 1, 0},
                                 {0, 0, 1, 0, 0, 1, 0, 0, 0},
                                 {0, 0, 1, 0, 0, 1, 0, 0, 0},
                                 {0, 0, 1, 0, 0, 1, 0, 0, 0}};
    
    
    bool isLeftRightProduced[9][9] ={{1, 0, 0, 0, 0, 0, 1, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {1, 0, 0, 0, 0, 0, 1, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 1, 0, 0, 1, 0},
                                     {0, 0, 0, 0, 1, 0, 0, 1, 0},
                                     {0, 0, 1, 0, 0, 1, 0, 0, 1},
                                     {0, 0, 1, 0, 0, 1, 0, 0, 1},
                                     {0, 0, 1, 0, 0, 1, 0, 0, 1}};
    
                                    
    bool isRightProduced[9][9] = {{0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {1, 0, 0, 0, 0, 0, 1, 0, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 1, 0, 0, 1, 0},
                                   {0, 0, 1, 0, 1, 1, 0, 1, 1},
                                   {0, 0, 1, 0, 0, 1, 0, 0, 1},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0}};
    
    
    /*
    
    bool isLeftProduced[9][9] = {{1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1}};
    
    
    bool isLeftRightProduced[9][9] ={{0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0}};
    
                                    
    bool isRightProduced[9][9] = {{1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1}};
    
    */
   
    
    /*
    bool isLeftProduced[9][9] = {{1, 0, 1, 0, 0, 0, 1, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 1, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 1, 0, 1, 1, 0, 1, 0},
                                 {0, 0, 1, 0, 0, 0, 1, 1, 0},
                                 {0, 0, 1, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 1, 0, 0, 1, 0, 0, 0},
                                 {0, 0, 1, 0, 0, 1, 1, 0, 1}};
                                    
    bool isRightProduced[9][9] = {{1, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {1, 0, 0, 0, 0, 0, 1, 1, 1},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 1, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 1, 0, 0, 1, 1},
                                   {1, 0, 1, 0, 1, 1, 0, 1, 1},
                                   {0, 0, 0, 0, 1, 1, 0, 0, 1},
                                   {0, 0, 0, 0, 0, 0, 0, 1, 1}};
    */
     
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
    
 //   vector< vector < vector < vector<real> > > > gammaLeftRight(windowedLength,
 //                                                               vector< vector < vector<real> > > (windowedLength, vector < vector<real> > (numberOfStates, vector <real> (numberOfStates, 0.0))));
  
    
    ulint numberOfDoubleStates = numberOfStates * numberOfStates;
    vector < vector < vector<real> > > gammaLeftRight(windowedLength, vector < vector<real> > (windowedLength, vector <real> (numberOfDoubleStates, 0.0)));
    
    
    //
    //Talvez tenha que colocar em uma mesma matriz bem maior gamaLeftRight
    //para normalizar corretamente !!!
    //
    
    vector<real> transLeftProb(numberOfFragments, 0.0);
    vector<real> transRightProb(numberOfFragments, 0.0);
    vector<real> probDoubleState(numberOfDoubleStates, 0.0);
    
    for(ulint i=0; i < numberOfStates; i++)
        for(ulint j=0; j < numberOfStates; j++)
            if(isLeftProduced[SScpair[i]][SScpair[j]] ||
               isRightProduced[SScpair[i]][SScpair[j]] ||
               isLeftRightProduced[SScpair[i]][SScpair[j]])
            {
                ulint doubleindex = i*numberOfStates +j;
                probDoubleState[doubleindex] = probState[i]*probState[j];
            }
    normalize(probDoubleState);
    
    
    //
    // 1.1 Initialization
    //
    // Note that only pair of states that could actually be generated by the grammar
    // will have non zero probability at the beginning and end of the chain !!!!
    //
    
    if(headAndTail){
        for(ulint i=0; i < numberOfStates; i++)
            for(ulint j=0; j < numberOfStates; j++)
   //             if(isLeftProduced[SScpair[i]][SScpair[j]] ||
   //                isRightProduced[SScpair[i]][SScpair[j]] ||
   //                isLeftRightProduced[SScpair[i]][SScpair[j]])
                        gammaLeftRight[0][windowedLength-1][i*numberOfStates + j] =
                            probState_head[i] * probState_tail[j];
    }
    else{
        for(ulint i=0; i < numberOfStates; i++)
            for(ulint j=0; j < numberOfStates; j++)
                if(isLeftProduced[SScpair[i]][SScpair[j]] &&
                   isRightProduced[SScpair[i]][SScpair[j]])
                        gammaLeftRight[0][windowedLength-1][i*numberOfStates + j] = 1.0;
    }
    normalize(gammaLeftRight[0][windowedLength-1]);
    
    // 1.2 Induction
    for(int t=0; t < windowedLength; t++){
    
  //  for(int u = windowedLength - 2; u >= 0; u-- ){
  //      for(int v = u; v < windowedLength; v++){
            
  //          int t = v - u;
        
          
            if(t > 0) {
                for(ulint i=0; i<numberOfFragments; i++)
                    transLeftProb[i] = probFragment[i] *
                        probFragmentEmitsPrimarySymbol[i][ seq[t-1+halfWindow] ];
                normalize(transLeftProb);
            }
    
        for (int v=windowedLength-1; v >= t; v--) {
 //           assert(v >= 0);
        
 //       for(uint u=0; u < windowedLength-1-t; u++) {
 //           uint v = windowedLength - 2 - u;
 //           uint v = seq.size() - windowSize - u;
            if(t==0 && v==windowedLength-1)
                continue;
            
            if(v < windowedLength - 1) {
                for(ulint i=0; i<numberOfFragments; i++)
                    transRightProb[i] = probFragment[i] *
                            probFragmentEmitsPrimarySymbol[i][ seq[v+halfWindow] ];
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
                            
                            ulint doubleindex = nextleft * numberOfStates + prevright;
                                                                                   
                            if(isLeftRightProduced[SScpair[nextleft]][SScpair[prevright]])
                            {
                                ulint parentdoubleindex = prevleft*numberOfStates + nextright;
                                if(t > 0 && v < (windowedLength - 1) &&
                                   probState[prevleft] != 0 && probState[nextright] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                                 //   ulint parentdoubleindex = prevleft*numberOfStates + nextright;
                                    real probleftright = gammaLeftRight[t-1][v+1][parentdoubleindex]*
                                        transLeftProb[transleft] * transRightProb[transright] /
                                        probState[prevleft] / probState[nextright];
                                        //probDoubleState[parentdoubleindex];
                                    gammaLeftRight[t][v][doubleindex] += probleftright;
                                }
                            }
                                                                                                
                            if(isLeftProduced[SScpair[nextleft]][SScpair[prevright]]) //||
                             //  isRightProduced[SScpair[nextleft]][SScpair[prevleft]] )
                             //  v == windowedLength - 1)
                            {
                                ulint parentdoubleindex = prevleft*numberOfStates + prevright;
                                if( t > 0 && probState[prevleft] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                              //      ulint parentdoubleindex = prevleft*numberOfStates + prevright;
                                    real probleft = gammaLeftRight[t-1][v][parentdoubleindex] *
                                        transLeftProb[transleft] /
                                        probState[prevleft];
                                        //probDoubleState[parentdoubleindex];
                                    gammaLeftRight[t][v][doubleindex] += probleft;
                                }
                            }
                            
                            if(isRightProduced[SScpair[nextleft]][SScpair[prevright]]) // ||
                            //   isLeftProduced[SScpair[prevright]][SScpair[nextright]])
                            //   t == 0)
                            {
                                ulint parentdoubleindex = nextleft*numberOfStates + nextright;
                                if(v < windowedLength - 1 && probState[nextright] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                                 //   ulint parentdoubleindex = nextleft*numberOfStates + nextright;
                                    real probright = gammaLeftRight[t][v+1][parentdoubleindex] *
                                        transRightProb[transright] /
                                        probState[nextright];
                                        //probDoubleState[parentdoubleindex];
                                    gammaLeftRight[t][v][doubleindex] += probright;
                                }
                            }
                            
                            if(//isLeftProduced[SScpair[nextleft]][SScpair[prevright]] ||
                               isRightProduced[SScpair[prevright]][SScpair[nextright]] )
                             //  v == windowedLength - 1)
                            {
                                ulint parentdoubleindex = prevleft*numberOfStates + prevright;
                                if( t > 0 && probState[prevleft] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                              //      ulint parentdoubleindex = prevleft*numberOfStates + prevright;
                                    real probleft = gammaLeftRight[t-1][v][parentdoubleindex] *
                                        transLeftProb[transleft] /
                                        probState[prevleft];
                                        //probDoubleState[parentdoubleindex];
                                    gammaLeftRight[t][v][doubleindex] += probleft;
                                }
                            }
                            
                            if(//isRightProduced[SScpair[nextleft]][SScpair[prevright]] ||
                               isLeftProduced[SScpair[prevleft]][SScpair[nextleft]])
                            //   t == 0)
                            {
                                ulint parentdoubleindex = nextleft*numberOfStates + nextright;
                                if(v < windowedLength - 1 && probState[nextright] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                                 //   ulint parentdoubleindex = nextleft*numberOfStates + nextright;
                                    real probright = gammaLeftRight[t][v+1][parentdoubleindex] *
                                        transRightProb[transright] /
                                        probState[nextright];
                                        //probDoubleState[parentdoubleindex];
                                    gammaLeftRight[t][v][doubleindex] += probright;
                                }
                            }
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
 
        
        for(ulint i=0; i<numberOfStates; i++){
            if(probState[i] == 0)
                continue;
    
            gamma[t][i] = gammaLeftRight[t][t][i*numberOfStates  + i] / probState[i];
        }

 /*       for(ulint i=0; i < numberOfFragments; i++)
            transProb[i] = probFragment[i] *
                probFragmentEmitsPrimarySymbol[i][ seq[t+halfWindow] ];
            normalize(transProb);
        
         for(ulint prev=0; prev<numberOfStates; prev++){
            
            if(probState[prev] == 0) continue;
            real probStateEmitsSeq = 0.0;
            
            for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                ulint trans =  prev * numberOfSecondarySymbols + ss;
                ulint next = trans % numberOfStates;
                
                if(probState[next] == 0) continue;
                
                ulint doubleindex = prev*numberOfStates + next;
                
                probStateEmitsSeq += transProb[trans] *
                                    gammaLeftRight[t][t+1][doubleindex] /
                                    probState[next];
            }
             
            probStateEmitsSeq /= probState[prev];
                  
            gamma[t][prev] = probStateEmitsSeq;
            
        }
  */
    
 
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


void Hmm::insideOutsideV3(vector<byte>& seq){
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
    
   // C-doubleState-C emits at left, H-doubleState-H and E-doubleState-E emit at left and right
    
   /*
    

    bool isLeftProduced[9][9] = {{0, 0, 0, 0, 0, 0, 1, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 1, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 1, 0, 0, 0, 0, 1, 0},
                                 {0, 0, 1, 0, 0, 0, 0, 1, 0},
                                 {0, 0, 1, 0, 0, 1, 1, 1, 0},
                                 {0, 0, 1, 0, 0, 1, 1, 1, 0},
                                 {0, 0, 1, 0, 0, 1, 1, 1, 0}};
    
    
    bool isLeftRightProduced[9][9] ={{1, 0, 0, 0, 0, 0, 1, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {1, 0, 0, 0, 0, 0, 1, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 1, 0, 0, 1, 0},
                                     {0, 0, 0, 0, 1, 0, 0, 1, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0}};
    
                                    
    bool isRightProduced[9][9] = {{0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {1, 0, 0, 0, 0, 0, 1, 0, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 1, 0, 0, 1, 0},
                                   {0, 0, 1, 0, 1, 1, 0, 1, 1},
                                   {0, 0, 1, 0, 0, 1, 0, 0, 1},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0}};
    
    */
  
    
  
  //  E-doubleState-E, H-doubleState-H and C-doubleSttate-C emit at left and right
   
    bool isLeftProduced[9][9] = {{0, 0, 0, 0, 0, 0, 1, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 1, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 1, 0, 0, 0, 0, 1, 0},
                                 {0, 0, 1, 0, 0, 0, 0, 1, 0},
                                 {0, 0, 1, 0, 0, 1, 0, 0, 0},
                                 {0, 0, 1, 0, 0, 1, 0, 0, 0},
                                 {0, 0, 1, 0, 0, 1, 0, 0, 0}};
    
    
    bool isLeftRightProduced[9][9] ={{1, 0, 0, 0, 0, 0, 1, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {1, 0, 0, 0, 0, 0, 1, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 1, 0, 0, 1, 0},
                                     {0, 0, 0, 0, 1, 0, 0, 1, 0},
                                     {0, 0, 1, 0, 0, 1, 0, 0, 1},
                                     {0, 0, 1, 0, 0, 1, 0, 0, 1},
                                     {0, 0, 1, 0, 0, 1, 0, 0, 1}};
    
                                    
    bool isRightProduced[9][9] = {{0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {1, 0, 0, 0, 0, 0, 1, 0, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 1, 0, 0, 1, 0},
                                   {0, 0, 1, 0, 1, 1, 0, 1, 1},
                                   {0, 0, 1, 0, 0, 1, 0, 0, 1},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0}};
    
    
    /*
    
    bool isLeftProduced[9][9] = {{1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1}};
    
    
    bool isLeftRightProduced[9][9] ={{0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0}};
    
                                    
    bool isRightProduced[9][9] = {{1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1}};
    
    */
   
    
    /*
    bool isLeftProduced[9][9] = {{1, 0, 1, 0, 0, 0, 1, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 1, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 1, 0, 1, 1, 0, 1, 0},
                                 {0, 0, 1, 0, 0, 0, 1, 1, 0},
                                 {0, 0, 1, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 1, 0, 0, 1, 0, 0, 0},
                                 {0, 0, 1, 0, 0, 1, 1, 0, 1}};
                                    
    bool isRightProduced[9][9] = {{1, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {1, 0, 0, 0, 0, 0, 1, 1, 1},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 1, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 1, 0, 0, 1, 1},
                                   {1, 0, 1, 0, 1, 1, 0, 1, 1},
                                   {0, 0, 0, 0, 1, 1, 0, 0, 1},
                                   {0, 0, 0, 0, 0, 0, 0, 1, 1}};
    */
     
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
    
 //   vector< vector < vector < vector<real> > > > gammaLeftRight(windowedLength,
 //                                                               vector< vector < vector<real> > > (windowedLength, vector < vector<real> > (numberOfStates, vector <real> (numberOfStates, 0.0))));
  
    
    ulint numberOfDoubleStates = numberOfStates * numberOfStates;
    vector < vector < vector<real> > > alpha(windowedLength, vector < vector<real> > (windowedLength, vector <real> (numberOfDoubleStates, 0.0)));
    vector < vector < vector<real> > > beta(windowedLength, vector < vector<real> > (windowedLength, vector <real> (numberOfDoubleStates, 0.0)));
    
    
    //
    //Talvez tenha que colocar em uma mesma matriz bem maior gamaLeftRight
    //para normalizar corretamente !!!
    //
    
    vector<real> transLeftProb(numberOfFragments, 0.0);
    vector<real> transRightProb(numberOfFragments, 0.0);
    vector<real> probDoubleState(numberOfDoubleStates, 0.0);
    
    for(ulint i=0; i < numberOfStates; i++)
        for(ulint j=0; j < numberOfStates; j++)
            if(isLeftProduced[SScpair[i]][SScpair[j]] ||
               isRightProduced[SScpair[i]][SScpair[j]] ||
               isLeftRightProduced[SScpair[i]][SScpair[j]])
            {
                ulint doubleindex = i*numberOfStates +j;
                probDoubleState[doubleindex] = probState[i]*probState[j];
            }
    normalize(probDoubleState);
    
    
    
    //
    // 1. Inside
    //
    
    //
    // 1.1 Initialization
    //
    for(int t=0; t < windowedLength-1; t++){
        for(ulint i=0; i<numberOfFragments; i++){
            ulint doubleindex = i / numberOfSecondarySymbols * numberOfStates + i % numberOfStates;
            alpha[t][t+1][doubleindex] = probFragment[i] *
                probFragmentEmitsPrimarySymbol[i][ seq[t+halfWindow] ];
        }
        normalize(alpha[t][t+1]);
    }
    
    //
    // 1.2 Induction
    //
    for(int v=2; v < windowedLength; v++){
        
        if(v < windowedLength - 1) {
            for(ulint i=0; i<numberOfFragments; i++)
                transRightProb[i] = probFragment[i] *
                        probFragmentEmitsPrimarySymbol[i][ seq[v+halfWindow] ];
                normalize(transRightProb);
        }
    
        for (int t=v-2; t >= 0; t--) {
            if(t > 0) {
                for(ulint i=0; i<numberOfFragments; i++)
                    transLeftProb[i] = probFragment[i] *
                        probFragmentEmitsPrimarySymbol[i][ seq[t-1+halfWindow] ];
                normalize(transLeftProb);
            }
            
            for(ulint nextleft=0; nextleft<numberOfStates; nextleft++){
                for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                    ulint transleft = nextleft + (ss * numberOfStates);
                    ulint prevleft = transleft / numberOfSecondarySymbols;
                    
                    for(ulint prevright=0; prevright<numberOfStates; prevright++){
                        for(uint tt=0; tt<numberOfSecondarySymbols; tt++){
                            ulint transright = (prevright * numberOfSecondarySymbols) + tt;
                            ulint nextright = transright % numberOfStates;
                            
                            ulint doubleindex = nextleft * numberOfStates + prevright;
                                                                                   
                            if(isLeftRightProduced[SScpair[nextleft]][SScpair[prevright]])
                            {
                                ulint parentdoubleindex = prevleft*numberOfStates + nextright;
                                if( (v - t) > 2 &&
                                   probState[nextleft] != 0 && probState[prevright] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                                 //   ulint parentdoubleindex = prevleft*numberOfStates + nextright;
                                    real probleftright = alpha[t+1][v-1][doubleindex]*
                                        transLeftProb[transleft] * transRightProb[transright] /
                                        probState[nextleft] / probState[prevright];
                                        //probDoubleState[parentdoubleindex];
                                    alpha[t][v][parentdoubleindex] += probleftright;
                                }
                            }
                                                                                                
                            if(isLeftProduced[SScpair[nextleft]][SScpair[prevright]]) //||
                             //  isRightProduced[SScpair[nextleft]][SScpair[prevleft]] )
                             //  v == windowedLength - 1)
                            {
                                ulint parentdoubleindex = prevleft*numberOfStates + prevright;
                                if( probState[nextleft] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                              //      ulint parentdoubleindex = prevleft*numberOfStates + prevright;
                                    real probleft = alpha[t+1][v][doubleindex] *
                                        transLeftProb[transleft] /
                                        probState[nextleft];
                                        //probDoubleState[parentdoubleindex];
                                    alpha[t][v][parentdoubleindex] += probleft;
                                }
                            }
                            
                            if(isRightProduced[SScpair[nextleft]][SScpair[prevright]]) // ||
                            //   isLeftProduced[SScpair[prevright]][SScpair[nextright]])
                            //   t == 0)
                            {
                                ulint parentdoubleindex = nextleft*numberOfStates + nextright;
                                if(probState[prevright] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                                 //   ulint parentdoubleindex = nextleft*numberOfStates + nextright;
                                    real probright = alpha[t][v-1][doubleindex] *
                                        transRightProb[transright] /
                                        probState[prevright];
                                        //probDoubleState[parentdoubleindex];
                                    alpha[t][v][parentdoubleindex] += probright;
                                }
                            }
                            
                       
                        }
                    }
                    
                }
            }
        
            normalize(alpha[t][v]);
            
            
            

        }
    }
    
    //
    //2. Outside
    //
    
    
    //
    // 2.1 Initialization
    //
    // Note that only pair of states that could actually be generated by the grammar
    // will have non zero probability at the beginning and end of the chain !!!!
    //
    
    if(headAndTail){
        for(ulint i=0; i < numberOfStates; i++)
            for(ulint j=0; j < numberOfStates; j++)
   //             if(isLeftProduced[SScpair[i]][SScpair[j]] ||
   //                isRightProduced[SScpair[i]][SScpair[j]] ||
   //                isLeftRightProduced[SScpair[i]][SScpair[j]])
                        beta[0][windowedLength-1][i*numberOfStates + j] =
                            probState_head[i] * probState_tail[j];
    }
    else{
        for(ulint i=0; i < numberOfStates; i++)
            for(ulint j=0; j < numberOfStates; j++)
                if(isLeftProduced[SScpair[i]][SScpair[j]] &&
                   isRightProduced[SScpair[i]][SScpair[j]])
                        beta[0][windowedLength-1][i*numberOfStates + j] = 1.0;
    }
    normalize(beta[0][windowedLength-1]);
    
    // 2.2 Induction
    for(int t=0; t < windowedLength-1; t++){
    
  //  for(int u = windowedLength - 2; u >= 0; u-- ){
  //      for(int v = u; v < windowedLength; v++){
            
  //          int t = v - u;
        
          
            if(t > 0) {
                for(ulint i=0; i<numberOfFragments; i++)
                    transLeftProb[i] = probFragment[i] *
                        probFragmentEmitsPrimarySymbol[i][ seq[t-1+halfWindow] ];
                normalize(transLeftProb);
            }
    
        for (int v=windowedLength-1; v > t; v--) {
 //           assert(v >= 0);
        
 //       for(uint u=0; u < windowedLength-1-t; u++) {
 //           uint v = windowedLength - 2 - u;
 //           uint v = seq.size() - windowSize - u;
            if(t==0 && v==windowedLength-1)
                continue;
            
            if(v < windowedLength - 1) {
                for(ulint i=0; i<numberOfFragments; i++)
                    transRightProb[i] = probFragment[i] *
                            probFragmentEmitsPrimarySymbol[i][ seq[v+halfWindow] ];
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
                            
                            ulint doubleindex = nextleft * numberOfStates + prevright;
                                                                                   
                            if(isLeftRightProduced[SScpair[nextleft]][SScpair[prevright]])
                            {
                                ulint parentdoubleindex = prevleft*numberOfStates + nextright;
                                if(t > 0 && v < (windowedLength - 1) &&
                                   probState[prevleft] != 0 && probState[nextright] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                                 //   ulint parentdoubleindex = prevleft*numberOfStates + nextright;
                                    real probleftright = beta[t-1][v+1][parentdoubleindex]*
                                        transLeftProb[transleft] * transRightProb[transright] /
                                        probState[prevleft] / probState[nextright];
                                        //probDoubleState[parentdoubleindex];
                                    beta[t][v][doubleindex] += probleftright;
                                }
                            }
                                                                                                
                            if(isLeftProduced[SScpair[nextleft]][SScpair[prevright]]) //||
                             //  isRightProduced[SScpair[nextleft]][SScpair[prevleft]] )
                             //  v == windowedLength - 1)
                            {
                                ulint parentdoubleindex = prevleft*numberOfStates + prevright;
                                if( t > 0 && probState[prevleft] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                              //      ulint parentdoubleindex = prevleft*numberOfStates + prevright;
                                    real probleft = beta[t-1][v][parentdoubleindex] *
                                        transLeftProb[transleft] /
                                        probState[prevleft];
                                        //probDoubleState[parentdoubleindex];
                                    beta[t][v][doubleindex] += probleft;
                                }
                            }
                            
                            if(isRightProduced[SScpair[nextleft]][SScpair[prevright]]) // ||
                            //   isLeftProduced[SScpair[prevright]][SScpair[nextright]])
                            //   t == 0)
                            {
                                ulint parentdoubleindex = nextleft*numberOfStates + nextright;
                                if(v < windowedLength - 1 && probState[nextright] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                                 //   ulint parentdoubleindex = nextleft*numberOfStates + nextright;
                                    real probright = beta[t][v+1][parentdoubleindex] *
                                        transRightProb[transright] /
                                        probState[nextright];
                                        //probDoubleState[parentdoubleindex];
                                    beta[t][v][doubleindex] += probright;
                                }
                            }
                            
                        }
                    }
                    
                }
            }
        
            if ( v > t+1 )
                normalize(beta[t][v]);
            
            else if (v == t+1){
                for(ulint i=0; i<numberOfFragments; i++){
                    if(probFragment[i] == 0)
                        continue;
                    ulint prev = i / numberOfSecondarySymbols;
                    ulint next = i % numberOfStates;
                    ulint doubleindex0 = prev * numberOfStates + next;
                    real transRightProb0 = probFragment[i] / probState[next];
                    real transLeftProb0 = probFragment[i] / probState[prev];
 
                    // left side
                    
                    for(ulint j = 0; j < numberOfStates; j++){
                        ulint doubleindex1 = j*numberOfStates + prev;
                                               
                        if(isRightProduced[SScpair[j]][SScpair[prev]]){
                            ulint doubleindex2 = j*numberOfStates + next;
                            for(int k = 0; k < t; k++){
                                real deltabeta = beta[k][v][doubleindex2] *
                                                    alpha[k][t][doubleindex1] * transRightProb0 ;
                                beta[t][v][doubleindex0] += deltabeta;
                            }
                        }
                        
                        if(isLeftRightProduced[SScpair[j]][SScpair[prev]]){
                            for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                                ulint transleft = j + (ss * numberOfStates);
                                ulint prevleft = transleft / numberOfSecondarySymbols;
                                
                                if(probState[prevleft] == 0)
                                    continue;
                                    
                                ulint doubleindex2 = prevleft*numberOfStates + next;
                                ulint doubleindex3 = prevleft*numberOfStates + j;
                                for(int k = 1; k < t; k++){
                                    real deltabeta = beta[k-1][v][doubleindex2] *
                                                        alpha[k-1][k][doubleindex3] *
                                                        alpha[k][t][doubleindex1] *
                                                        transRightProb0 / probState[prevleft] ;
                                    beta[t][v][doubleindex0] += deltabeta;
                                
                                }
                            }
                        }
                    }
                    
                    // right side
                    
                    for(ulint j = 0; j < numberOfStates; j++){
                        ulint doubleindex1 = next*numberOfStates + j;
                                               
                        if(isLeftProduced[SScpair[next]][SScpair[j]]){
                            ulint doubleindex2 = prev*numberOfStates + j;
                            for(int k = windowedLength - 1; k > v; k--){
                                
                                real deltabeta = beta[t][k][doubleindex2] *
                                                    alpha[v][k][doubleindex1] * transLeftProb0 ;
                                beta[t][v][doubleindex0] += deltabeta;
                            }
                        }
                        
                        if(isLeftRightProduced[SScpair[next]][SScpair[j]]){
                            for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                                ulint transright = j*numberOfSecondarySymbols + ss;
                                ulint nextright = transright % numberOfStates;
                                if(probState[nextright] == 0)
                                    continue;
                                ulint doubleindex2 = prev*numberOfStates + nextright;
                                ulint doubleindex3 = j*numberOfStates + nextright;
                                for(int k = windowedLength - 2; k > v; k--){
                                    real deltabeta = beta[t][k+1][doubleindex2] *
                                                        alpha[k][k+1][doubleindex3] *
                                                        alpha[v][k][doubleindex1] *
                                                        transLeftProb0 / probState[nextright] ;
                                    beta[t][v][doubleindex0]+= deltabeta;
                                }
                            }
                        }
                    }
                }
                
                normalize(beta[t][v]);

            }
            
            else cout << "this shouldn't be happening";
        
        }
    }
  
    //''''''''''''''''''''''''''''''''''''''''''''''
    // 3. gamma
    //''''''''''''''''''''''''''''''''''''''''''''''
    

    // gamma
    gamma.clear();
    gamma.resize(windowedLength, vector<real>(numberOfStates,0.0));
    
    for(uint t=0; t < windowedLength-1; t++){
        
        for(ulint i=0; i<numberOfStates; i++){
            for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                ulint trans = i*numberOfSecondarySymbols + ss;
                if(probFragment[trans] == 0)
                    continue;
                ulint next = trans % numberOfStates;
                ulint doubleindex = i*numberOfStates + next;
                real deltagamma = beta[t][t+1][doubleindex] *
                                    alpha[t][t+1][doubleindex] / probFragment[trans];
                gamma[t][i] += deltagamma;
                if(t == windowedLength - 2)
                    gamma[windowedLength - 1][next] += deltagamma;
            }
        }
        normalize(gamma[t]);
    }
    normalize(gamma[windowedLength - 1]);
            

 
 
 //       normalize(gamma[t]);
 //   }
    
    
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


void Hmm::insideOutsideV4(vector<byte>& seq){
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
    
   // C-doubleState-C emits at left, H-doubleState-H and E-doubleState-E emit at left and right
    
   /*

    bool isLeftProduced[9][9] = {{0, 0, 0, 0, 0, 0, 1, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 1, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 1, 0, 0, 0, 0, 1, 0},
                                 {0, 0, 1, 0, 0, 0, 0, 1, 0},
                                 {0, 0, 1, 0, 0, 1, 1, 1, 0},
                                 {0, 0, 1, 0, 0, 1, 1, 1, 0},
                                 {0, 0, 1, 0, 0, 1, 1, 1, 0}};
    
    
    bool isLeftRightProduced[9][9] ={{1, 0, 0, 0, 0, 0, 1, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {1, 0, 0, 0, 0, 0, 1, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 1, 0, 0, 1, 0},
                                     {0, 0, 0, 0, 1, 0, 0, 1, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0}};
    
                                    
    bool isRightProduced[9][9] = {{0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {1, 0, 0, 0, 0, 0, 1, 0, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 1, 0, 0, 1, 0},
                                   {0, 0, 1, 0, 1, 1, 0, 1, 1},
                                   {0, 0, 1, 0, 0, 1, 0, 0, 1},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0}};
    
    
 */
 /*
  
  //  E-doubleState-E, H-doubleState-H and C-doubleSttate-C emit at left and right
   
    bool isLeftProduced[9][9] = {{0, 0, 0, 0, 0, 0, 1, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 1, 0, 0},
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 {0, 0, 1, 0, 0, 0, 0, 1, 0},
                                 {0, 0, 1, 0, 0, 0, 0, 1, 0},
                                 {0, 0, 1, 0, 0, 1, 0, 0, 0},
                                 {0, 0, 1, 0, 0, 1, 0, 0, 0},
                                 {0, 0, 1, 0, 0, 1, 0, 0, 0}};
    
    
    bool isLeftRightProduced[9][9] ={{1, 0, 0, 0, 0, 0, 1, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {1, 0, 0, 0, 0, 0, 1, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 1, 0, 0, 1, 0},
                                     {0, 0, 0, 0, 1, 0, 0, 1, 0},
                                     {0, 0, 1, 0, 0, 1, 0, 0, 1},
                                     {0, 0, 1, 0, 0, 1, 0, 0, 1},
                                     {0, 0, 1, 0, 0, 1, 0, 0, 1}};
    
                                    
    bool isRightProduced[9][9] = {{0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {1, 0, 0, 0, 0, 0, 1, 0, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                   {0, 0, 0, 0, 1, 0, 0, 1, 0},
                                   {0, 0, 1, 0, 1, 1, 0, 1, 1},
                                   {0, 0, 1, 0, 0, 1, 0, 0, 1},
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0}};
  
  */
    
 //   /*
    
    // Everybody is produced from left only.
    // This should provide the same result as the ForwardBacward of the HMM
    
 /*
    bool isLeftProduced[9][9] = {{1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1}};
 */
  // /*
    
    bool isLeftProduced[9][9] ={{0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0}};
  //  */
    
    bool isLeftRightProduced[9][9] = {{1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                 {1, 1, 1, 1, 1, 1, 1, 1, 1}};
    
    /*
    bool isLeftRightProduced[9][9] ={{0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0}};
     */
    
 // /*
    bool isRightProduced[9][9] ={{0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0}};
    
  // */
    /*
    bool isRightProduced[9][9] = {{1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1},
                                   {1, 1, 1, 1, 1, 1, 1, 1, 1}};
    
    */
    
    ulint numberOfDoubleStates = numberOfStates * numberOfStates;
    vector < vector < vector<real> > > alpha(windowedLength, vector < vector<real> > (windowedLength+1, vector <real> (numberOfDoubleStates, 0.0)));
    vector < vector < vector<real> > > beta(windowedLength, vector < vector<real> > (windowedLength+1, vector <real> (numberOfDoubleStates, 0.0)));
    
    
    
    vector<real> transLeftProb(numberOfFragments, 0.0);
    vector<real> transRightProb(numberOfFragments, 0.0);
    vector<real> probDoubleState(numberOfDoubleStates, 0.0);
    
    for(ulint i=0; i < numberOfStates; i++)
        for(ulint j=0; j < numberOfStates; j++)
            if(isLeftProduced[SScpair[i]][SScpair[j]] ||
               isRightProduced[SScpair[i]][SScpair[j]] ||
               isLeftRightProduced[SScpair[i]][SScpair[j]])
            {
                ulint doubleindex = i*numberOfStates +j;
                probDoubleState[doubleindex] = probState[i]*probState[j];
            }
    normalize(probDoubleState);
    
    
    
    //
    // 1. Inside
    //
    
    //
    // 1.1 Initialization
    //
    
    for(int t=0; t < windowedLength; t++){
        for(ulint i=0; i<numberOfStates; i++){
            ulint doubleindex = i * numberOfStates + i;
            alpha[t][t][doubleindex] = 1.0;
        }
        normalize(alpha[t][t]);
    }
    
    for(int t=0; t < windowedLength-1; t++){
        for(ulint i=0; i<numberOfFragments; i++){
            ulint doubleindex = i / numberOfSecondarySymbols * numberOfStates + i % numberOfStates;
            alpha[t][t+1][doubleindex] = probFragment[i] *
                probFragmentEmitsPrimarySymbol[i][ seq[t+halfWindow] ];
        }
        normalize(alpha[t][t+1]);
    }
    
    //
    // 1.2 Induction
    //
    for(int v=2; v < windowedLength; v++){
        
   
        /*
        if(v <= windowedLength - 1) {
            for(ulint i=0; i<numberOfFragments; i++)
                transRightProb[i] = probFragment[i] *
                        probFragmentEmitsPrimarySymbol[i][ seq[v+halfWindow] ];
                normalize(transRightProb);
        }
    */
      
        for (int t=v-2; t >= 0; t--) {
        
            /*
            if(t > 0) {
                for(ulint i=0; i<numberOfFragments; i++)
                    transLeftProb[i] = probFragment[i] *
                        probFragmentEmitsPrimarySymbol[i][ seq[t-1+halfWindow] ];
                normalize(transLeftProb);
            }
             */
            
            for(ulint nextleft=0; nextleft<numberOfStates; nextleft++){
                for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                    ulint transleft = nextleft + (ss * numberOfStates);
                    ulint prevleft = transleft / numberOfSecondarySymbols;
                    ulint doubleindexleft = prevleft*numberOfStates + nextleft;
                    
                    for(ulint prevright=0; prevright<numberOfStates; prevright++){
                        for(uint tt=0; tt<numberOfSecondarySymbols; tt++){
                            ulint transright = (prevright * numberOfSecondarySymbols) + tt;
                            ulint nextright = transright % numberOfStates;
                            ulint doubleindexright = prevright * numberOfStates + nextright;
                            
                            ulint doubleindex = nextleft * numberOfStates + prevright;
                                                                                   
                            if(isLeftRightProduced[SScpair[nextleft]][SScpair[prevright]])
                            {
                                ulint parentdoubleindex = prevleft*numberOfStates + nextright;
                                if( (v - t) > 2 &&
                                   probState[prevleft] != 0 && probState[nextright] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                                 //   ulint parentdoubleindex = prevleft*numberOfStates + nextright;
                                    real probleftright = alpha[t+1][v-1][doubleindex]*
                                        alpha[t][t+1][doubleindexleft] *
                                        alpha[v-1][v][doubleindexright] /
                                        // transLeftProb[transleft] * transRightProb[transright] /
                                        probState[prevleft] / probState[nextright];
                                        //probDoubleState[parentdoubleindex];
                                    alpha[t][v][parentdoubleindex] += probleftright;
                                }
                            }
                                                                                                
                            if(isLeftProduced[SScpair[nextleft]][SScpair[prevright]]) //||
                             //  isRightProduced[SScpair[nextleft]][SScpair[prevleft]] )
                             //  v == windowedLength - 1)
                            {
                                ulint parentdoubleindex = prevleft*numberOfStates + prevright;
                                if( probState[prevleft] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                              //      ulint parentdoubleindex = prevleft*numberOfStates + prevright;
                                    real probleft = alpha[t+1][v][doubleindex] *
                                        alpha[t][t+1][doubleindexleft] /
                                      //  transLeftProb[transleft] /
                                        probState[prevleft];
                                        //probDoubleState[parentdoubleindex];
                                    alpha[t][v][parentdoubleindex] += probleft;
                                }
                            }
                            
                            if(isRightProduced[SScpair[nextleft]][SScpair[prevright]]) // ||
                            //   isLeftProduced[SScpair[prevright]][SScpair[nextright]])
                            //   t == 0)
                            {
                                ulint parentdoubleindex = nextleft*numberOfStates + nextright;
                                if(probState[nextright] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                                 //   ulint parentdoubleindex = nextleft*numberOfStates + nextright;
                                    real probright = alpha[t][v-1][doubleindex] *
                                        alpha[v-1][v][doubleindexright] /
                                  //      transRightProb[transright] /
                                        probState[nextright];
                                        //probDoubleState[parentdoubleindex];
                                    alpha[t][v][parentdoubleindex] += probright;
                                }
                            }
                            
                       
                        }
                    }
                    
                }
            }
        
            normalize(alpha[t][v]);
            
            
            

        }
    }
    
    //
    //2. Outside
    //
    
    
    //
    // 2.1 Initialization
    //
    // Note that only pair of states that could actually be generated by the grammar
    // will have non zero probability at the beginning and end of the chain !!!!
    //
    
    if(headAndTail){
        for(ulint i=0; i < numberOfStates; i++)
            for(ulint j=0; j < numberOfStates; j++)
                if( isLeftProduced[SScpair[i]][SScpair[j]] ||
                    isRightProduced[SScpair[i]][SScpair[j]] ||
                   isLeftRightProduced[SScpair[i]][SScpair[j]])
                        beta[0][windowedLength-1][i*numberOfStates + j] =
                            probState_head[i] * probState_tail[j];
    }
    else{
        for(ulint i=0; i < numberOfStates; i++)
            for(ulint j=0; j < numberOfStates; j++)
                if(isLeftProduced[SScpair[i]][SScpair[j]] &&
                   isRightProduced[SScpair[i]][SScpair[j]])
                        beta[0][windowedLength-1][i*numberOfStates + j] = 1.0;
    }
    normalize(beta[0][windowedLength-1]);
    
    // 2.2 Induction
    for(int t=0; t < windowedLength-1; t++){
        int v;
        if(t > 0) {
            for(ulint i=0; i<numberOfFragments; i++)
                transLeftProb[i] = probFragment[i] *
                        probFragmentEmitsPrimarySymbol[i][ seq[t-1+halfWindow] ];
            normalize(transLeftProb);
        }
    
        for (v=windowedLength-1; v > (t+1); v--) {
 
            if(t==0 && v==windowedLength-1)
                continue;
            
            if(v < windowedLength - 1) {
                for(ulint i=0; i<numberOfFragments; i++)
                    transRightProb[i] = probFragment[i] *
                            probFragmentEmitsPrimarySymbol[i][ seq[v+halfWindow] ];
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
                            
                            ulint doubleindex = nextleft * numberOfStates + prevright;
                                                                                   
                            if(isLeftRightProduced[SScpair[nextleft]][SScpair[prevright]])
                            {
                                ulint parentdoubleindex = prevleft*numberOfStates + nextright;
                                if(t > 0 && v < (windowedLength - 1) &&
                                   probState[prevleft] != 0 && probState[nextright] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                                 //   ulint parentdoubleindex = prevleft*numberOfStates + nextright;
                                    real probleftright = beta[t-1][v+1][parentdoubleindex]*
                                        transLeftProb[transleft] * transRightProb[transright] /
                                        probState[prevleft] / probState[nextright];
                                        //probDoubleState[parentdoubleindex];
                                    beta[t][v][doubleindex] += probleftright;
                                }
                            }
                                                                                                
                            if(isLeftProduced[SScpair[nextleft]][SScpair[prevright]]) //||
                             //  isRightProduced[SScpair[nextleft]][SScpair[prevleft]] )
                             //  v == windowedLength - 1)
                            {
                                  ulint parentdoubleindex = prevleft*numberOfStates + prevright;
                                if( t > 0 && probState[prevleft] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                              //      ulint parentdoubleindex = prevleft*numberOfStates + prevright;
                                    real probleft = beta[t-1][v][parentdoubleindex] *
                                        transLeftProb[transleft] /
                                        probState[prevleft];
                                        //probDoubleState[parentdoubleindex];
                                    beta[t][v][doubleindex] += probleft;
                                }
                            }
                            
                            if(isRightProduced[SScpair[nextleft]][SScpair[prevright]]) // ||
                            //   isLeftProduced[SScpair[prevright]][SScpair[nextright]])
                            //   t == 0)
                            {
                                ulint parentdoubleindex = nextleft*numberOfStates + nextright;
                                if(v <= windowedLength - 1 && probState[nextright] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                                 //   ulint parentdoubleindex = nextleft*numberOfStates + nextright;
                                    real probright = beta[t][v+1][parentdoubleindex] *
                                        transRightProb[transright] /
                                        probState[nextright];
                                        //probDoubleState[parentdoubleindex];
                                    beta[t][v][doubleindex] += probright;
                                }
                            }
                            
                        }
                    }
                    
                }
            }
        
 //           if ( v > t+1 )
 
            normalize(beta[t][v]);
        }
            
 //           else
 
 //           if (v == t+1){
        
                
        for(ulint i=0; i<numberOfFragments; i++){
            if(probFragment[i] == 0)
                continue;
            ulint prev = i / numberOfSecondarySymbols;
            ulint next = i % numberOfStates;
            ulint doubleindex0 = prev * numberOfStates + next;
            real transRightProb0 = probFragment[i] / probState[next];
            real transLeftProb0 = probFragment[i] / probState[prev];
            
            
            // left and right sides not producing central fragment
            
            if(isLeftRightProduced[SScpair[prev]][SScpair[next]])
            {
                for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                    ulint transleft = prev + (ss * numberOfStates);
                    ulint prevleft = transleft / numberOfSecondarySymbols;
                    for(uint tt=0; tt<numberOfSecondarySymbols; tt++){
                        ulint transright = (next * numberOfSecondarySymbols) + tt;
                        ulint nextright = transright % numberOfStates;
                        if(t > 0 && v < (windowedLength - 1) &&
                           probState[prevleft] != 0 && probState[nextright] != 0)
                        {
                            ulint parentdoubleindex = prevleft*numberOfStates + nextright;
                            ulint leftdoubleindex = prevleft*numberOfStates + prev;
                            ulint rightdoubleindex = next*numberOfStates + nextright;
                            real probleftright = beta[t-1][v+1][parentdoubleindex]*
                                alpha[t-1][t][leftdoubleindex] * alpha[v][v+1][rightdoubleindex] /
                                probState[prevleft] / probState[nextright];
                                //probDoubleState[parentdoubleindex];
                            beta[t][v][doubleindex0] += probleftright;
                        }
                    }
                }
            }
                                                                                        
            //left side not producing central fragment
            
            if(isLeftProduced[SScpair[prev]][SScpair[next]])
            {
                for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                    ulint transleft = prev + (ss * numberOfStates);
                    ulint prevleft = transleft / numberOfSecondarySymbols;
                    if( t > 0 && probState[prevleft] != 0)
                    {
                        ulint parentdoubleindex = prevleft*numberOfStates + next;
                        ulint leftdoubleindex =prevleft*numberOfStates + prev;
                        real probleft = beta[t-1][v][parentdoubleindex] *
                                    alpha[t-1][t][leftdoubleindex] / probState[prevleft];
                        beta[t][v][doubleindex0] += probleft;
                    }
                }
            }
            
            //right side not producing central fragment
                    
            if(isRightProduced[SScpair[prev]][SScpair[next]])
            {
                for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                    ulint transright = next * numberOfSecondarySymbols + ss;
                    ulint nextright = transright % numberOfStates;
                    if( v < (windowedLength - 1) && probState[nextright] != 0)
                    {
                        ulint parentdoubleindex = prev*numberOfStates + nextright;
                        ulint rightdoubleindex = next*numberOfStates + nextright;
                        real probright = beta[t][v+1][parentdoubleindex] *
                                    alpha[v][v+1][rightdoubleindex] / probState[nextright];
                        beta[t][v][doubleindex0] += probright;
                    }
                }
            }
            
 
            // left side producing central fragment
                    
                    
            for(ulint j = 0; j < numberOfStates; j++){
                ulint doubleindex1 = j*numberOfStates + prev;

                if(isRightProduced[SScpair[j]][SScpair[prev]]){
                    ulint doubleindex2 = j*numberOfStates + next;
                        for(int w = 0; w < t; w++){
                            real deltabeta = beta[w][v][doubleindex2] *
                                                    alpha[w][t][doubleindex1] * transRightProb0 ;
                            beta[t][v][doubleindex0] += deltabeta;
                        }
                }
                        
                if(isLeftRightProduced[SScpair[j]][SScpair[prev]]){
                    for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                        ulint transleft = j + (ss * numberOfStates);
                        ulint prevleft = transleft / numberOfSecondarySymbols;
                        if(probState[prevleft] == 0)
                            continue;
                                    
                                
                        ulint doubleindex2 = prevleft*numberOfStates + next;
                        ulint doubleindex3 = prevleft*numberOfStates + j;
                            
                        for(int w = 1; w < t; w++){
                            real deltabeta = beta[w-1][v][doubleindex2] *
                                                        alpha[w-1][w][doubleindex3] *
                                                        alpha[w][t][doubleindex1] *
                                                        transRightProb0 / probState[prevleft] ;
                            beta[t][v][doubleindex0] += deltabeta;
                        }
                            
                    }
                        
                }
                    
            }
                    
            // right side producing central fragment
                    
                    
            for(ulint j = 0; j < numberOfStates; j++){
            
                ulint doubleindex1 = next*numberOfStates + j;
                                               
                if(isLeftProduced[SScpair[next]][SScpair[j]]){
                    ulint doubleindex2 = prev*numberOfStates + j;
                    for(int w = windowedLength - 1 ; w > v; w--){
                        real deltabeta = beta[t][w][doubleindex2] *
                                            alpha[v][w][doubleindex1] * transLeftProb0 ;
                        
                        beta[t][v][doubleindex0] += deltabeta;
                            
                    }
                }
                        
                if(isLeftRightProduced[SScpair[next]][SScpair[j]]){
                    for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                        ulint transright = j*numberOfSecondarySymbols + ss;
                        ulint nextright = transright % numberOfStates;
                        if(probState[nextright] == 0)
                            continue;
                        ulint doubleindex2 = prev*numberOfStates + nextright;
                        ulint doubleindex3 = j*numberOfStates + nextright;
                            for(int w = windowedLength - 2; w > v; w--){
                                real deltabeta = beta[t][w+1][doubleindex2] *
                                                    alpha[w][w+1][doubleindex3] *
                                                    alpha[v][w][doubleindex1] *
                                                    transLeftProb0 / probState[nextright] ;
                                beta[t][v][doubleindex0]+= deltabeta;
                            }
                        }
                }
                    
            }
                
        }
                
        normalize(beta[t][v]);
            
            //else cout << "this shouldn't be happening";
    }
  
    //''''''''''''''''''''''''''''''''''''''''''''''
    // 3. gamma
    //''''''''''''''''''''''''''''''''''''''''''''''
    

    // gamma
    gamma.clear();
    gamma.resize(windowedLength, vector<real>(numberOfStates,0.0));
    
    for(uint t=0; t < windowedLength-1; t++){
        
        for(ulint i=0; i<numberOfStates; i++){
            for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                ulint trans = i*numberOfSecondarySymbols + ss;
                if(probFragment[trans] == 0)
                    continue;
                ulint next = trans % numberOfStates;
                ulint doubleindex = i*numberOfStates + next;
                real deltagamma = beta[t][t+1][doubleindex] *
                                    probFragmentEmitsPrimarySymbol[trans][seq[t+halfWindow]];
                              //      alpha[t][t+1][doubleindex];
                                      // /  probState[next];
                // / probFragment[trans];
                gamma[t][i] += deltagamma; // / probFragment[trans];
                if(t == windowedLength - 2)
                    gamma[windowedLength - 1][next] += deltagamma; // /probFragment[trans];
            }
        }
        normalize(gamma[t]);
    }
    normalize(gamma[windowedLength - 1]);
            

 
 
 //       normalize(gamma[t]);
 //   }
    
    
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
