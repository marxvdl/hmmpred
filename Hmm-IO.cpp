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

#include <unordered_map>


/**
 * Performs the Inside-Outside algorithm.
 * The resulting gamma matrix is stored in the private variable "gamma".
 */


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
          
            for(ulint nextLeft=0; nextLeft<numberOfStates; nextLeft++){
                for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                    ulint transLeft = nextLeft + (ss * numberOfStates);
                    ulint prevLeft = transLeft / numberOfSecondarySymbols;
                    
                    for(ulint prevRight=0; prevRight<numberOfStates; prevRight++){
                        for(uint tt=0; tt<numberOfSecondarySymbols; tt++){
                            ulint transRight = (prevRight * numberOfSecondarySymbols) + tt;
                            ulint nextRight = transRight % numberOfStates;
                            
                            ulint doubleindex = nextLeft * numberOfStates + prevRight;
                                                                                   
                            if(isLeftRightProduced[SScpair[nextLeft]][SScpair[prevRight]])
                            {
                                ulint parentdoubleindex = prevLeft*numberOfStates + nextRight;
                                if(t > 0 && v < (windowedLength - 1) &&
                                   probState[prevLeft] != 0 && probState[nextRight] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                                 //   ulint parentdoubleindex = prevLeft*numberOfStates + nextRight;
                                    real probleftright = gammaLeftRight[t-1][v+1][parentdoubleindex]*
                                        transLeftProb[transLeft] * transRightProb[transRight] /
                                        probState[prevLeft] / probState[nextRight];
                                        //probDoubleState[parentdoubleindex];
                                    gammaLeftRight[t][v][doubleindex] += probleftright;
                                }
                            }
                                                                                                
                            if(isLeftProduced[SScpair[nextLeft]][SScpair[prevRight]]) //||
                             //  isRightProduced[SScpair[nextLeft]][SScpair[prevLeft]] )
                             //  v == windowedLength - 1)
                            {
                                ulint parentdoubleindex = prevLeft*numberOfStates + prevRight;
                                if( t > 0 && probState[prevLeft] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                              //      ulint parentdoubleindex = prevLeft*numberOfStates + prevRight;
                                    real probleft = gammaLeftRight[t-1][v][parentdoubleindex] *
                                        transLeftProb[transLeft] /
                                        probState[prevLeft];
                                        //probDoubleState[parentdoubleindex];
                                    gammaLeftRight[t][v][doubleindex] += probleft;
                                }
                            }
                            
                            if(isRightProduced[SScpair[nextLeft]][SScpair[prevRight]]) // ||
                            //   isLeftProduced[SScpair[prevRight]][SScpair[nextRight]])
                            //   t == 0)
                            {
                                ulint parentdoubleindex = nextLeft*numberOfStates + nextRight;
                                if(v < windowedLength - 1 && probState[nextRight] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                                 //   ulint parentdoubleindex = nextLeft*numberOfStates + nextRight;
                                    real probright = gammaLeftRight[t][v+1][parentdoubleindex] *
                                        transRightProb[transRight] /
                                        probState[nextRight];
                                        //probDoubleState[parentdoubleindex];
                                    gammaLeftRight[t][v][doubleindex] += probright;
                                }
                            }
                            
                            if(//isLeftProduced[SScpair[nextLeft]][SScpair[prevRight]] ||
                               isRightProduced[SScpair[prevRight]][SScpair[nextRight]] )
                             //  v == windowedLength - 1)
                            {
                                ulint parentdoubleindex = prevLeft*numberOfStates + prevRight;
                                if( t > 0 && probState[prevLeft] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                              //      ulint parentdoubleindex = prevLeft*numberOfStates + prevRight;
                                    real probleft = gammaLeftRight[t-1][v][parentdoubleindex] *
                                        transLeftProb[transLeft] /
                                        probState[prevLeft];
                                        //probDoubleState[parentdoubleindex];
                                    gammaLeftRight[t][v][doubleindex] += probleft;
                                }
                            }
                            
                            if(//isRightProduced[SScpair[nextLeft]][SScpair[prevRight]] ||
                               isLeftProduced[SScpair[prevLeft]][SScpair[nextLeft]])
                            //   t == 0)
                            {
                                ulint parentdoubleindex = nextLeft*numberOfStates + nextRight;
                                if(v < windowedLength - 1 && probState[nextRight] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                                 //   ulint parentdoubleindex = nextLeft*numberOfStates + nextRight;
                                    real probright = gammaLeftRight[t][v+1][parentdoubleindex] *
                                        transRightProb[transRight] /
                                        probState[nextRight];
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
            
            for(ulint nextLeft=0; nextLeft<numberOfStates; nextLeft++){
                for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                    ulint transLeft = nextLeft + (ss * numberOfStates);
                    ulint prevLeft = transLeft / numberOfSecondarySymbols;
                    
                    for(ulint prevRight=0; prevRight<numberOfStates; prevRight++){
                        for(uint tt=0; tt<numberOfSecondarySymbols; tt++){
                            ulint transRight = (prevRight * numberOfSecondarySymbols) + tt;
                            ulint nextRight = transRight % numberOfStates;
                            
                            ulint doubleindex = nextLeft * numberOfStates + prevRight;
                                                                                   
                            if(isLeftRightProduced[SScpair[nextLeft]][SScpair[prevRight]])
                            {
                                ulint parentdoubleindex = prevLeft*numberOfStates + nextRight;
                                if( (v - t) > 2 &&
                                   probState[nextLeft] != 0 && probState[prevRight] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                                 //   ulint parentdoubleindex = prevLeft*numberOfStates + nextRight;
                                    real probleftright = alpha[t+1][v-1][doubleindex]*
                                        transLeftProb[transLeft] * transRightProb[transRight] /
                                        probState[nextLeft] / probState[prevRight];
                                        //probDoubleState[parentdoubleindex];
                                    alpha[t][v][parentdoubleindex] += probleftright;
                                }
                            }
                                                                                                
                            if(isLeftProduced[SScpair[nextLeft]][SScpair[prevRight]]) //||
                             //  isRightProduced[SScpair[nextLeft]][SScpair[prevLeft]] )
                             //  v == windowedLength - 1)
                            {
                                ulint parentdoubleindex = prevLeft*numberOfStates + prevRight;
                                if( probState[nextLeft] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                              //      ulint parentdoubleindex = prevLeft*numberOfStates + prevRight;
                                    real probleft = alpha[t+1][v][doubleindex] *
                                        transLeftProb[transLeft] /
                                        probState[nextLeft];
                                        //probDoubleState[parentdoubleindex];
                                    alpha[t][v][parentdoubleindex] += probleft;
                                }
                            }
                            
                            if(isRightProduced[SScpair[nextLeft]][SScpair[prevRight]]) // ||
                            //   isLeftProduced[SScpair[prevRight]][SScpair[nextRight]])
                            //   t == 0)
                            {
                                ulint parentdoubleindex = nextLeft*numberOfStates + nextRight;
                                if(probState[prevRight] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                                 //   ulint parentdoubleindex = nextLeft*numberOfStates + nextRight;
                                    real probright = alpha[t][v-1][doubleindex] *
                                        transRightProb[transRight] /
                                        probState[prevRight];
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
          
            
            for(ulint nextLeft=0; nextLeft<numberOfStates; nextLeft++){
                for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                    ulint transLeft = nextLeft + (ss * numberOfStates);
                    ulint prevLeft = transLeft / numberOfSecondarySymbols;
                    
                    for(ulint prevRight=0; prevRight<numberOfStates; prevRight++){
                        for(uint tt=0; tt<numberOfSecondarySymbols; tt++){
                            ulint transRight = (prevRight * numberOfSecondarySymbols) + tt;
                            ulint nextRight = transRight % numberOfStates;
                            
                            ulint doubleindex = nextLeft * numberOfStates + prevRight;
                                                                                   
                            if(isLeftRightProduced[SScpair[nextLeft]][SScpair[prevRight]])
                            {
                                ulint parentdoubleindex = prevLeft*numberOfStates + nextRight;
                                if(t > 0 && v < (windowedLength - 1) &&
                                   probState[prevLeft] != 0 && probState[nextRight] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                                 //   ulint parentdoubleindex = prevLeft*numberOfStates + nextRight;
                                    real probleftright = beta[t-1][v+1][parentdoubleindex]*
                                        transLeftProb[transLeft] * transRightProb[transRight] /
                                        probState[prevLeft] / probState[nextRight];
                                        //probDoubleState[parentdoubleindex];
                                    beta[t][v][doubleindex] += probleftright;
                                }
                            }
                                                                                                
                            if(isLeftProduced[SScpair[nextLeft]][SScpair[prevRight]]) //||
                             //  isRightProduced[SScpair[nextLeft]][SScpair[prevLeft]] )
                             //  v == windowedLength - 1)
                            {
                                ulint parentdoubleindex = prevLeft*numberOfStates + prevRight;
                                if( t > 0 && probState[prevLeft] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                              //      ulint parentdoubleindex = prevLeft*numberOfStates + prevRight;
                                    real probleft = beta[t-1][v][parentdoubleindex] *
                                        transLeftProb[transLeft] /
                                        probState[prevLeft];
                                        //probDoubleState[parentdoubleindex];
                                    beta[t][v][doubleindex] += probleft;
                                }
                            }
                            
                            if(isRightProduced[SScpair[nextLeft]][SScpair[prevRight]]) // ||
                            //   isLeftProduced[SScpair[prevRight]][SScpair[nextRight]])
                            //   t == 0)
                            {
                                ulint parentdoubleindex = nextLeft*numberOfStates + nextRight;
                                if(v < windowedLength - 1 && probState[nextRight] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                                 //   ulint parentdoubleindex = nextLeft*numberOfStates + nextRight;
                                    real probright = beta[t][v+1][parentdoubleindex] *
                                        transRightProb[transRight] /
                                        probState[nextRight];
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
                                ulint transLeft = j + (ss * numberOfStates);
                                ulint prevLeft = transLeft / numberOfSecondarySymbols;
                                
                                if(probState[prevLeft] == 0)
                                    continue;
                                    
                                ulint doubleindex2 = prevLeft*numberOfStates + next;
                                ulint doubleindex3 = prevLeft*numberOfStates + j;
                                for(int k = 1; k < t; k++){
                                    real deltabeta = beta[k-1][v][doubleindex2] *
                                                        alpha[k-1][k][doubleindex3] *
                                                        alpha[k][t][doubleindex1] *
                                                        transRightProb0 / probState[prevLeft] ;
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
                                ulint transRight = j*numberOfSecondarySymbols + ss;
                                ulint nextRight = transRight % numberOfStates;
                                if(probState[nextRight] == 0)
                                    continue;
                                ulint doubleindex2 = prev*numberOfStates + nextRight;
                                ulint doubleindex3 = j*numberOfStates + nextRight;
                                for(int k = windowedLength - 2; k > v; k--){
                                    real deltabeta = beta[t][k+1][doubleindex2] *
                                                        alpha[k][k+1][doubleindex3] *
                                                        alpha[v][k][doubleindex1] *
                                                        transLeftProb0 / probState[nextRight] ;
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


void Hmm::insideOutsideV2(vector<byte>& seq){
       uint fullLength = (uint) seq.size();
    uint windowedLength = fullLength - windowSize + 2;

//    vector<real> transProb(numberOfFragments, 0.0);
    
    
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
    
   // /*

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
    
    
  // */
 
  
   //  E-doubleState-E, H-doubleState-H and C-doubleSttate-C emit at left and right
   
    /*
    
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
    
    // Everybody is produced from left only, or right only, or leftright only.
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
    /*
    
    bool isLeftProduced[9][9] ={{0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0}};
      */
    
     /*
    
    bool isLeftRightProduced[9][9] = {{1, 1, 1, 1, 1, 1, 1, 1, 1},
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
    
     /*
    bool isRightProduced[9][9] ={{0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0}};
    
       */
    
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
    vector<real> transProb(numberOfFragments, 0.0);
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
            
            for(ulint nextLeft=0; nextLeft<numberOfStates; nextLeft++){
                for(ulint prevRight=0; prevRight<numberOfStates; prevRight++){
                    
                    ulint doubleindex = nextLeft * numberOfStates + prevRight;

                    if(isLeftRightProduced[SScpair[nextLeft]][SScpair[prevRight]])
                    {
                        for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                            ulint transLeft = nextLeft + (ss * numberOfStates);
                            ulint prevLeft = transLeft / numberOfSecondarySymbols;
                            ulint doubleindexleft = prevLeft*numberOfStates + nextLeft;
                                        
                            for(uint tt=0; tt<numberOfSecondarySymbols; tt++){
                                ulint transRight = (prevRight * numberOfSecondarySymbols) + tt;
                                ulint nextRight = transRight % numberOfStates;
                                ulint doubleindexright = prevRight * numberOfStates + nextRight;
                                
                                ulint parentdoubleindex = prevLeft*numberOfStates + nextRight;
                                if( (v - t) > 2 &&
                                   probState[nextLeft] != 0 && probState[prevRight] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                                 //   ulint parentdoubleindex = prevLeft*numberOfStates + nextRight;
                                    real probleftright = alpha[t+1][v-1][doubleindex]*
                                        alpha[t][t+1][doubleindexleft] *
                                        alpha[v-1][v][doubleindexright] /
                                        // transLeftProb[transLeft] * transRightProb[transRight] /
                                        probState[nextLeft] / probState[prevRight];
                                        //probDoubleState[parentdoubleindex];
                                    alpha[t][v][parentdoubleindex] += probleftright;
                                }
                            }
                        }
                    }
                                                                                                
                    if(isLeftProduced[SScpair[nextLeft]][SScpair[prevRight]]) //||
                             //  isRightProduced[SScpair[nextLeft]][SScpair[prevLeft]] )
                             //  v == windowedLength - 1)
                    {
                        for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                            ulint transLeft = nextLeft + (ss * numberOfStates);
                            ulint prevLeft = transLeft / numberOfSecondarySymbols;
                            ulint doubleindexleft = prevLeft*numberOfStates + nextLeft;
                      
                            ulint parentdoubleindex = prevLeft*numberOfStates + prevRight;
                            if( probState[nextLeft] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                            {
                              //      ulint parentdoubleindex = prevLeft*numberOfStates + prevRight;
                                real probleft = alpha[t+1][v][doubleindex] *
                                        alpha[t][t+1][doubleindexleft] /
                                      //  transLeftProb[transLeft] /
                                        probState[nextLeft];
                                        //probDoubleState[parentdoubleindex];
                                alpha[t][v][parentdoubleindex] += probleft;
                            }
                        }
                    }
                            
                    if(isRightProduced[SScpair[nextLeft]][SScpair[prevRight]]) // ||
                            //   isLeftProduced[SScpair[prevRight]][SScpair[nextRight]])
                            //   t == 0)
                    {
                        for(uint tt=0; tt<numberOfSecondarySymbols; tt++){
                            ulint transRight = (prevRight * numberOfSecondarySymbols) + tt;
                            ulint nextRight = transRight % numberOfStates;
                            ulint doubleindexright = prevRight * numberOfStates + nextRight;

                            ulint parentdoubleindex = nextLeft*numberOfStates + nextRight;
                            if(probState[prevRight] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                            {
                                 //   ulint parentdoubleindex = nextLeft*numberOfStates + nextRight;
                                real probright = alpha[t][v-1][doubleindex] *
                                        alpha[v-1][v][doubleindexright] /
                                  //      transRightProb[transRight] /
                                probState[prevRight];
                                        //probDoubleState[parentdoubleindex];
                                alpha[t][v][parentdoubleindex] += probright;
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
        
        /*
        if(t > 0) {
            for(ulint i=0; i<numberOfFragments; i++)
                transLeftProb[i] = probFragment[i] *
                        probFragmentEmitsPrimarySymbol[i][ seq[t-1+halfWindow] ];
            normalize(transLeftProb);
        }
         */
    
        for (v=windowedLength-1; v > t; v--) {
 
            if(t==0 && v==windowedLength-1)
                continue;
            
  
            /*
            if(v < windowedLength - 1) {
                for(ulint i=0; i<numberOfFragments; i++)
                    transRightProb[i] = probFragment[i] *
                            probFragmentEmitsPrimarySymbol[i][ seq[v+halfWindow] ];
                normalize(transRightProb);
            }
            */
            
            for(ulint nextLeft=0; nextLeft<numberOfStates; nextLeft++){
                for(ulint prevRight=0; prevRight<numberOfStates; prevRight++){
                    
                    ulint doubleindex = nextLeft * numberOfStates + prevRight;
                    beta[t][v][doubleindex]=0;
                       
                    if(isLeftRightProduced[SScpair[nextLeft]][SScpair[prevRight]]){
                        for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                            ulint transLeft = nextLeft + (ss * numberOfStates);
                            ulint prevLeft = transLeft / numberOfSecondarySymbols;
                            ulint doubleindexleft = prevLeft * numberOfStates + nextLeft;
                                    
                            for(uint tt=0; tt<numberOfSecondarySymbols; tt++){
                                ulint transRight = (prevRight * numberOfSecondarySymbols) + tt;
                                ulint nextRight = transRight % numberOfStates;
                                ulint doubleindexright = prevRight * numberOfStates + nextRight;
                                 
                                ulint parentdoubleindex = prevLeft*numberOfStates + nextRight;
                                if(t > 0 && v < (windowedLength - 1) &&
                                   probState[prevLeft] != 0 && probState[nextRight] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                                {
                                 //   ulint parentdoubleindex = prevLeft*numberOfStates + nextRight;
                                    real probleftright = beta[t-1][v+1][parentdoubleindex]*
                                        alpha[t-1][t][doubleindexleft]*
                                        alpha[v][v+1][doubleindexright]/
                                      //  transLeftProb[transLeft] * transRightProb[transRight] /
                                        probState[prevLeft] / probState[nextRight];
                                        //probDoubleState[parentdoubleindex];
                                    beta[t][v][doubleindex] += probleftright;
                                }
                            }
                        }
                    }
                                                                                                
                    if(isLeftProduced[SScpair[nextLeft]][SScpair[prevRight]]) //||
                             //  isRightProduced[SScpair[nextLeft]][SScpair[prevLeft]] )
                             //  v == windowedLength - 1)
                    {
                        for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                            ulint transLeft = nextLeft + (ss * numberOfStates);
                            ulint prevLeft = transLeft / numberOfSecondarySymbols;
                            ulint doubleindexleft = prevLeft * numberOfStates + nextLeft;
                                  
                            ulint parentdoubleindex = prevLeft*numberOfStates + prevRight;
                            if( t > 0 && probState[prevLeft] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                            {
                              //      ulint parentdoubleindex = prevLeft*numberOfStates + prevRight;
                                real probleft = beta[t-1][v][parentdoubleindex] *
                                        alpha[t-1][t][doubleindexleft] /
                                       // transLeftProb[transLeft] /
                                        probState[prevLeft];
                                        //probDoubleState[parentdoubleindex];
                                beta[t][v][doubleindex] += probleft;
                            }
                        }
                    }
                            
                    if(isRightProduced[SScpair[nextLeft]][SScpair[prevRight]]) // ||
                            //   isLeftProduced[SScpair[prevRight]][SScpair[nextRight]])
                            //   t == 0)
                    {
                        for(uint tt=0; tt<numberOfSecondarySymbols; tt++){
                            ulint transRight = (prevRight * numberOfSecondarySymbols) + tt;
                            ulint nextRight = transRight % numberOfStates;
                            ulint doubleindexright = prevRight * numberOfStates + nextRight;
                             
                            ulint parentdoubleindex = nextLeft*numberOfStates + nextRight;
                            if(v <= windowedLength - 1 && probState[nextRight] != 0)
                                   //probDoubleState[parentdoubleindex] != 0)
                            {
                                 //   ulint parentdoubleindex = nextLeft*numberOfStates + nextRight;
                                real probright = beta[t][v+1][parentdoubleindex] *
                                        alpha[v][v+1][doubleindexright] /
                                     //   transRightProb[transRight] /
                                        probState[nextRight];
                                        //probDoubleState[parentdoubleindex];
                                beta[t][v][doubleindex] += probright;
                            }
                        }
                    }
                    
                }
            }

 //           if ( v > t+1 )
 
  //          normalize(beta[t][v]);
  //      }
            
 //           else
 
            
            if (v == t+1){
                
        //        for(ulint i=0; i<numberOfFragments; i++)
        //            transProb[i] = probFragment[i] *
        //                    probFragmentEmitsPrimarySymbol[i][ seq[t+halfWindow] ];
        //        normalize(transProb);
        
                for(ulint i=0; i<numberOfFragments; i++){
                    if(probFragment[i] == 0)
                        continue;
                    ulint prev = i / numberOfSecondarySymbols;
                    ulint next = i % numberOfStates;
                    ulint doubleindex0 = prev * numberOfStates + next;
                    real transRightProb0 = probFragment[i] / probState[next];
                    real transLeftProb0 = probFragment[i] / probState[prev];
                
            
  
            
  
            
            // left side producing central fragment
                    
                    
                    for(ulint j = 0; j < numberOfStates; j++){
                        ulint doubleindex1 = j*numberOfStates + prev;

                        if(isRightProduced[SScpair[j]][SScpair[prev]] && probState[j] != 0){
                            ulint doubleindex2 = j*numberOfStates + next;
                            for(int w = 0; w < t; w++){
                                real deltabeta = beta[w][v][doubleindex2] *
                                                    alpha[w][t][doubleindex1] / probState[j]; // *
                                                        //transProb[i]
                                                        //alpha[t][v][doubleindex0]/
                                                      //  probFragment[i]/
                                                       // probState[next];
                                                         //transRightProb0 ;
                                beta[t][v][doubleindex0] += deltabeta;
                            }
                        }
                        
                        if(isLeftRightProduced[SScpair[j]][SScpair[prev]] && probState[j] != 0){
                            for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                                    ulint transLeft = j + (ss * numberOfStates);
                                ulint prevLeft = transLeft / numberOfSecondarySymbols;
                                if(probState[prevLeft] == 0)
                                    continue;
                                    
                                
                                ulint doubleindex2 = prevLeft*numberOfStates + next;
                                ulint doubleindex3 = prevLeft*numberOfStates + j;
                            
                                for(int w = 1; w < t; w++){
                                    real deltabeta = beta[w-1][v][doubleindex2] *
                                                        alpha[w-1][w][doubleindex3] / probState[prevLeft] *
                                    alpha[w][t][doubleindex1] / probState[j];
                                                        //transProb[i]
                                                        //alpha[t][v][doubleindex0] /
                                                        //probState[next] /
                                                       // probState[prevLeft] / probState[j];
                                                     //   transRightProb0 / probState[prevLeft] ;
                                    beta[t][v][doubleindex0] += deltabeta;
                                }
                            
                            }
                        
                        }
                    
                    }
                    
            // right side producing central fragment
                    
                    
                    for(ulint j = 0; j < numberOfStates; j++){
            
                        ulint doubleindex1 = next*numberOfStates + j;
                                               
                        if(isLeftProduced[SScpair[next]][SScpair[j]] && probState[j] != 0){
                            ulint doubleindex2 = prev*numberOfStates + j;
                            for(int w = windowedLength - 1 ; w > v; w--){
                                real deltabeta = beta[t][w][doubleindex2] *
                                alpha[v][w][doubleindex1] / probState[j]; // *
                                                //alpha[t][v][doubleindex0] /
   //                                             probFragment[i] /
   //                                             probState[prev];
                                                //transLeftProb0 ;
                        
                                beta[t][v][doubleindex0] += deltabeta;
                            
                            }
                        }
                        
                        if(isLeftRightProduced[SScpair[next]][SScpair[j]] && probState[j] != 0 ){
                            for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                                ulint transRight = j*numberOfSecondarySymbols + ss;
                                ulint nextRight = transRight % numberOfStates;
                                if(probState[nextRight] == 0)
                                    continue;
                                ulint doubleindex2 = prev*numberOfStates + nextRight;
                                ulint doubleindex3 = j*numberOfStates + nextRight;
                                for(int w = windowedLength - 2; w > v; w--){
                                    real deltabeta = beta[t][w+1][doubleindex2] *
                                                        alpha[w][w+1][doubleindex3] *
                                                        alpha[v][w][doubleindex1] /
                                                //* alpha[t][v][doubleindex0] / probState[prev] /
                                                        probState[nextRight] / probState[j];
                                                       // transLeftProb0 / probState[nextRight] ;
                                    beta[t][v][doubleindex0]+= deltabeta;
                                }
                            }
                        }
                    
                    }
                
                }
            }
             
            normalize(beta[t][v]);

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
         //       real deltagamma = beta[t][t+1][doubleindex];
                real deltagamma = beta[t][t+1][doubleindex] * alpha[t][t+1][doubleindex] ;
                 // / probState[next];
                //* probFragmentEmitsPrimarySymbol[trans][seq[t+halfWindow]]; // / probState[i];
                                    // / probState[next] / probState[i];
                              //      alpha[t][t+1][doubleindex];
                                      // /  probState[next];
                // / probFragment[trans];
                gamma[t][i] += deltagamma / probState[i] /probState[next];
                 // probFragment[trans];
                if(t == windowedLength - 2)
                    gamma[windowedLength - 1][next] += deltagamma / probState[i] / probState[next];
                 // /probFragment[trans] / probState[next];
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


void Hmm::insideOutsideV3(vector<byte>& seq){
       
    uint fullLength = (uint) seq.size();
    uint windowedLength = fullLength - windowSize + 2;
    
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
    
    //The Stochastic Context-free grammar is now defined by a 9x9 byte matrix
    //that determines the unique type of emisson for each non-terminal symbol aWu,
    //depending on the central pairs of secondary structure of States â and û.
    //The indexes, from 0 to 8, correspond to pairs EE, EH, EC, HE, ..., CC,
    //in this order on both dimensions. The value in each cell indicates the type of
    //emission, 3 (binary 11) for left/right emission, 2 (binary 10) for
    //left emission, 1 (binay 01) for right emission, and 0 for combinations
    //not expected to appear during any derivation of the grammar.
    
    
    // Double productions for the same secondary structures: XE-W-EX, XH-W-HX, XC-W-CX.
    // Otherwise first left and then right, E > H > C, no EH or HE.
    
    //  /*
    byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 2, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {1, 0, 0, 0, 0, 0, 3, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 2, 0, 3, 0, 0, 2, 0},
                              {0, 0, 2, 0, 1, 0, 0, 3, 0},
                              {0, 0, 3, 0, 1, 1, 0, 1, 1},
                              {0, 0, 2, 0, 0, 3, 0, 0, 1},
                              {0, 0, 2, 0, 0, 2, 0, 0, 3}};
     //  */
    
    // No double production for XC-W-CX. Always FIRST LEFT and then at right.
    
     /*
    byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 2, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {1, 0, 0, 0, 0, 0, 2, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 2, 0, 3, 0, 0, 2, 0},
                              {0, 0, 2, 0, 1, 0, 0, 2, 0},
                              {0, 0, 3, 0, 1, 1, 0, 1, 1},
                              {0, 0, 2, 0, 0, 3, 0, 0, 1},
                              {0, 0, 2, 0, 0, 2, 2, 2, 0}};
      */
    
    // No double production for XC-W-CX. Always FIRST RIGHT and then at left.
    
     /*
    byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 2, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {1, 0, 0, 0, 0, 0, 1, 0, 1},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 2, 0, 3, 0, 0, 2, 0},
                              {0, 0, 2, 0, 1, 0, 0, 1, 1},
                              {0, 0, 3, 0, 1, 1, 0, 1, 1},
                              {0, 0, 2, 0, 0, 3, 0, 0, 1},
                              {0, 0, 2, 0, 0, 2, 0, 0, 0}};
      */
    
    // No double production for XH-W-HX. Always FIRST LEFT and then right.
    
     /*
    byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 2, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {1, 0, 0, 0, 0, 0, 3, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 2, 0, 0, 2, 0, 0, 0},
                              {0, 0, 2, 0, 1, 1, 0, 3, 0},
                              {0, 0, 3, 0, 1, 1, 0, 1, 1},
                              {0, 0, 2, 0, 0, 2, 0, 0, 1},
                              {0, 0, 2, 0, 0, 2, 0, 0, 3}};
     */
    
    // No double production for XH-W-HX. Always FIRST RIGHT and then left.
    
     /*
    byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 2, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {1, 0, 0, 0, 0, 0, 3, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 2, 0, 0, 0, 0, 2, 0},
                              {0, 0, 2, 0, 0, 0, 0, 3, 0},
                              {0, 0, 3, 0, 1, 1, 0, 1, 1},
                              {0, 0, 2, 0, 1, 1, 0, 2, 1},
                              {0, 0, 2, 0, 0, 2, 0, 0, 3}};
     */
    
    
    // No double production for XC-W-CX nor XH-W-HX . Always FIRST LEFT and then right.
    
    
     /*
    byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 2, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {1, 0, 0, 0, 0, 0, 2, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 2, 0, 0, 2, 0, 0, 0},
                              {0, 0, 2, 0, 1, 1, 0, 2, 0},
                              {0, 0, 3, 0, 1, 1, 0, 1, 1},
                              {0, 0, 2, 0, 0, 2, 0, 0, 1},
                              {0, 0, 2, 0, 0, 2, 2, 2, 0}};
      */
    
    // No double production for XC-W-CX nor XH-W-HX . Always FIRST RIGHT and then left.
    
    
     /*
    byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 2, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {1, 0, 0, 0, 0, 0, 1, 0, 1},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 2, 0, 0, 0, 0, 2, 0},
                              {0, 0, 2, 0, 0, 0, 0, 1, 1},
                              {0, 0, 3, 0, 1, 1, 0, 1, 1},
                              {0, 0, 2, 0, 1, 1, 0, 2, 1},
                              {0, 0, 2, 0, 0, 2, 0, 0, 0}};
     */
    
    
    // No H anywhere !!! Double E and C
    
    /*
   byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 2, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {1, 0, 0, 0, 0, 0, 3, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 3, 0, 0, 0, 0, 0, 1},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 2, 0, 0, 0, 0, 0, 3}};
     */
    
    
    
     /*
   // Everybody emits at left, equivalent to the HMM

   byte emissionType[9][9] ={{2, 2, 2, 2, 2, 2, 2, 2, 2},
                             {2, 2, 2, 2, 2, 2, 2, 2, 2},
                             {2, 2, 2, 2, 2, 2, 2, 2, 2},
                             {2, 2, 2, 2, 2, 2, 2, 2, 2},
                             {2, 2, 2, 2, 2, 2, 2, 2, 2},
                             {2, 2, 2, 2, 2, 2, 2, 2, 2},
                             {2, 2, 2, 2, 2, 2, 2, 2, 2},
                             {2, 2, 2, 2, 2, 2, 2, 2, 2},
                             {2, 2, 2, 2, 2, 2, 2, 2, 2}};
       */
    
   
    
      /*
    // Everybody emits at right, also equivalent to the HMM
   
   byte emissionType[9][9] ={{1, 1, 1, 1, 1, 1, 1, 1, 1},
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
    //Everybody emits at left/right, again, also equivalent to the HMM
    
    byte emissionType[9][9] ={{3, 3, 3, 3, 3, 3, 3, 3, 3},
                              {3, 3, 3, 3, 3, 3, 3, 3, 3},
                              {3, 3, 3, 3, 3, 3, 3, 3, 3},
                              {3, 3, 3, 3, 3, 3, 3, 3, 3},
                              {3, 3, 3, 3, 3, 3, 3, 3, 3},
                              {3, 3, 3, 3, 3, 3, 3, 3, 3},
                              {3, 3, 3, 3, 3, 3, 3, 3, 3},
                              {3, 3, 3, 3, 3, 3, 3, 3, 3},
                              {3, 3, 3, 3, 3, 3, 3, 3, 3}};
        */
    
    
    
    ulint numberOfDoubleStates = numberOfStates * numberOfStates;
//    int numberOfSecondarySymbols2 = numberOfSecondarySymbols*numberOfSecondarySymbols;
    
    // map<ulint, ulint> INDEX;
    unordered_map<ulint, ulint> betaINDEX;
    
    ulint maxIndex = 0;
    for(ulint i = 0; i < numberOfStates; i++){
        for(ulint j = 0; j < numberOfStates; j++){
            if(emissionType[SScpair[i]][SScpair[j]]) {
                ulint doubleindex = i*numberOfStates +j;
                betaINDEX[doubleindex] = maxIndex;
                maxIndex ++;
            }
        }
    }
    
    vector < vector < vector<real> > > alpha(windowedLength, vector < vector<real> > (windowedLength+1, vector <real> (numberOfDoubleStates, 0.0)));
  //  vector < vector < vector<real> > > beta(windowedLength, vector < vector<real> > (windowedLength+1, vector <real> (numberOfDoubleStates, 0.0)));
    
    
 //   vector < vector < vector<real> > > alpha(windowedLength, vector < vector<real> > (windowedLength+1, vector <real> (maxIndex, 0.0)));
    vector < vector < vector<real> > > beta(windowedLength, vector < vector<real> > (windowedLength+1, vector <real> (maxIndex, 0.0)));
    
    

    
    //
    // 1. Inside
    //
    
    //
    // 1.1 Initialization
    //
    
  /*
    for(int t=0; t < windowedLength; t++){
        for(ulint i=0; i<numberOfStates; i++){
            ulint doubleindex = i * numberOfStates + i;
            alpha[t][t][doubleindex] = 1.0;
        }
        normalize(alpha[t][t]);
    }
   */
    
    for(int t=0; t < windowedLength-1; t++){
        for(ulint i=0; i<numberOfFragments; i++){
            ulint left = i / numberOfSecondarySymbols;
            ulint right = i % numberOfStates;
            ulint doubleindex = left * numberOfStates + right;
            alpha[t][t+1][doubleindex] = probFragment[i] *
                        probFragmentEmitsPrimarySymbol[i][ seq[t+halfWindow] ];
        }
        normalize(alpha[t][t+1]);
    }
    
    //
    // 1.2 Induction
    //
    for(int v=2; v < windowedLength; v++){
      
        for (int t=v-2; t >= 0; t--) {
            
            for(ulint prevLeft=0; prevLeft<numberOfStates; prevLeft++){
          //      if(probState[prevLeft] == 0) continue;
                for(ulint nextRight=0; nextRight<numberOfStates; nextRight++){
           //         if(probState[nextRight] == 0) continue;
                    
                    if(! emissionType[SScpair[prevLeft]][SScpair[nextRight]])
                        continue;
                    
                    ulint parentdoubleindex = prevLeft * numberOfStates + nextRight;
  //                  ulint indexOfParentdoubleindex = INDEX.at(parentdoubleindex);
                    
                    switch(emissionType[SScpair[prevLeft]][SScpair[nextRight]])
                    {
                        
                        case 1:
                            for(uint tt=0; tt<numberOfSecondarySymbols; tt++){
                                
                                ulint transRight = tt * numberOfStates + nextRight;
                                ulint prevRight = transRight / numberOfSecondarySymbols;
                                ulint doubleindexright = prevRight * numberOfStates + nextRight;
                                ulint doubleindex = prevLeft*numberOfStates + prevRight;
                                
                                if(probState[prevRight] != 0)
                                {
                                    real probright = alpha[t][v-1][doubleindex] *
                                            alpha[v-1][v][doubleindexright] /
                                            probState[prevRight];
                                alpha[t][v][parentdoubleindex] += probright;
                                }
                            }
                            break;
                        
                        case 2:
                            for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                                
                                ulint transLeft = prevLeft * numberOfSecondarySymbols + ss;
                                ulint nextLeft = transLeft % numberOfStates;
                                ulint doubleindexleft = prevLeft*numberOfStates + nextLeft;
                                ulint doubleindex = nextLeft*numberOfStates + nextRight;
                                
                                if( probState[nextLeft] != 0)
                                {
                                    real probleft = alpha[t+1][v][doubleindex] *
                                        alpha[t][t+1][doubleindexleft] /
                                        probState[nextLeft];
                                    alpha[t][v][parentdoubleindex] += probleft;
                                }
                            }
                            break;
                            
                        case 3:
                            for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                                
                                ulint transLeft = prevLeft * numberOfSecondarySymbols + ss;
                                ulint nextLeft = transLeft % numberOfStates;
                                if(probState[nextLeft] == 0) continue;
                                ulint doubleindexleft = prevLeft*numberOfStates + nextLeft;
                                        
                                for(uint tt=0; tt<numberOfSecondarySymbols; tt++){
                                    ulint transRight = tt * numberOfStates + nextRight;
                                    ulint prevRight = transRight / numberOfSecondarySymbols;
                                    if(probState[prevRight] == 0) continue;
                                    
                                    ulint doubleindexright = prevRight * numberOfStates + nextRight;
                                    ulint doubleindex = nextLeft*numberOfStates + prevRight;
                                    
                                    if( (v - t) > 2 &&
                                       probState[nextLeft] != 0 && probState[prevRight] != 0)
                                    {
                                        real probleftright = alpha[t+1][v-1][doubleindex]
                                                            * alpha[t][t+1][doubleindexleft]
                                                            / probState[nextLeft]
                                                            * alpha[v-1][v][doubleindexright]
                                                            / probState[prevRight];
                                    alpha[t][v][parentdoubleindex] += probleftright;
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
    // Maybe not necessary anymore in version 5.
    //
    
    if(headAndTail){
        for(ulint i=0; i < numberOfStates; i++){
            for(ulint j=0; j < numberOfStates; j++){
                if(emissionType[SScpair[i]][SScpair[j]] > 0){
                    ulint doubleindex = i*numberOfStates + j;
                    beta[0][windowedLength-1][betaINDEX.at(doubleindex)] =
                            probState_head[i] * probState_tail[j];
                }
            }
        }
    }
    else{
        for(ulint i=0; i < numberOfStates; i++){
            for(ulint j=0; j < numberOfStates; j++){
                if(emissionType[SScpair[i]][SScpair[j]] > 0){
                    ulint doubleindex = i*numberOfStates + j;
                    beta[0][windowedLength-1][betaINDEX.at(doubleindex)] = 1.0;
                }
            }
        }
    }
    normalize(beta[0][windowedLength-1]);
    
    // 2.2 Induction
    
    for(int t=0; t < windowedLength-1; t++){
    
        for (int v=windowedLength-1; v > t; v--){
 
     //       if(t==0 && v==windowedLength-1)
     //           continue;
            
            for(ulint prevLeft=0; prevLeft<numberOfStates; prevLeft++){
           //     if(probState[prevLeft] == 0) continue;
                
                for(ulint nextRight=0; nextRight<numberOfStates; nextRight++){
               //     if(probState[nextRight] == 0) continue;
                    
                    if(emissionType[SScpair[prevLeft]][SScpair[nextRight]] == 0)
                        continue;
                    
                    ulint parentdoubleindex = prevLeft * numberOfStates + nextRight;
                    ulint indexOfParentdoubleindex = betaINDEX.at(parentdoubleindex);
                    
                    if(emissionType[SScpair[prevLeft]][SScpair[nextRight]] == 1
                            && v < windowedLength - 1)
                    {
                        for(uint tt=0; tt<numberOfSecondarySymbols; tt++)
                        {
                            ulint transRight = tt * numberOfStates + nextRight;
                            ulint prevRight = transRight /  numberOfSecondarySymbols;
                            ulint doubleindexright = prevRight * numberOfStates + nextRight;
                            
                            if(emissionType[SScpair[prevLeft]][SScpair[prevRight]] == 0)
                                continue;
                            ulint doubleindex = prevLeft*numberOfStates + prevRight;
                            if(probState[nextRight] != 0)
                            {
                                real probright = beta[t][v+1][indexOfParentdoubleindex] *
                                                alpha[v][v+1][doubleindexright] /
                                                probState[nextRight];
                                beta[t][v][betaINDEX.at(doubleindex)] += probright;
                            }
                        }
                    }
                    
                    else if(emissionType[SScpair[prevLeft]][SScpair[nextRight]] == 2
                            && t > 0)
                    {
                        for(uint ss=0; ss<numberOfSecondarySymbols; ss++)
                        {
                            ulint transLeft = prevLeft * numberOfSecondarySymbols + ss;
                            ulint nextLeft = transLeft % numberOfStates;
                            ulint doubleindexleft = prevLeft * numberOfStates + nextLeft;
                            
                            if(emissionType[SScpair[nextLeft]][SScpair[nextRight]] == 0)
                                continue;
                            ulint doubleindex = nextLeft*numberOfStates + nextRight;
                            if(probState[prevLeft] != 0)
                            {
                                real probleft = beta[t-1][v][indexOfParentdoubleindex] *
                                                alpha[t-1][t][doubleindexleft] /
                                                probState[prevLeft];
                                beta[t][v][betaINDEX.at(doubleindex)] += probleft;
                            }
                        }
                    }
                    
                    else if(emissionType[SScpair[prevLeft]][SScpair[nextRight]] == 3
                       && t > 0 && v < windowedLength -1)
                    {
                        for(uint ss=0; ss<numberOfSecondarySymbols; ss++)
                        {
                            ulint transLeft = prevLeft * numberOfSecondarySymbols + ss;
                            ulint nextLeft = transLeft % numberOfStates;
                            ulint doubleindexleft = prevLeft * numberOfStates + nextLeft;
                                    
                            for(uint tt=0; tt<numberOfSecondarySymbols; tt++)
                            {
                                ulint transRight = tt * numberOfStates + nextRight;
                                ulint prevRight = transRight / numberOfSecondarySymbols;
                                ulint doubleindexright = prevRight * numberOfStates + nextRight;
                                
                                if(emissionType[SScpair[nextLeft]][SScpair[prevRight]] == 0)
                                    continue;
                                ulint doubleindex = nextLeft*numberOfStates + prevRight;
                                
                                if(probState[prevLeft] != 0 && probState[nextRight] != 0)
                                {
                                    real probleftright =  alpha[t-1][t][doubleindexleft] *
                                                        alpha[v][v+1][doubleindexright] *
                                                        beta[t-1][v+1][indexOfParentdoubleindex] /
                                                        probState[prevLeft] /
                                                        probState[nextRight]; // / probsum;
                                    beta[t][v][betaINDEX.at(doubleindex)] += probleftright;
                                }
                            }
                        }
                    }
                }
            }
            
            if (v == t+1){
                
                for(ulint trans=0; trans < numberOfFragments; trans++){
                    if(probFragment[trans] == 0)
                        continue;
                    ulint prev = trans / numberOfSecondarySymbols;
                    ulint next = trans % numberOfStates;
                    if(emissionType[SScpair[prev]][SScpair[next]] == 0)
                        continue;
                    ulint doubleindex0 = prev * numberOfStates + next;
                    ulint indexOfDoubleindex0 = betaINDEX.at(doubleindex0);
                
                    //
                    // possible production with other possibly long fragments at left of t
                    // current short fragment emits at right
                    //
                    
                    for(ulint prevLeft = 0; prevLeft < numberOfStates; prevLeft++){
                        if(probState[prevLeft] == 0
                           || emissionType[SScpair[prevLeft]][SScpair[next]] == 0)
                                continue;
                        ulint doubleindex2 = prevLeft*numberOfStates + next;
                        ulint indexOfDoubleindex2 = betaINDEX.at(doubleindex2);
                    
                        if(emissionType[SScpair[prevLeft]][SScpair[next]] == 1) {
                            ulint doubleindex1 = prevLeft*numberOfStates + prev;
                                
                            real deltabeta = 0.0;
                            for(int w = 0; w < t; w++){
                                deltabeta += beta[w][v][indexOfDoubleindex2] *
                                                alpha[w][t][doubleindex1];
                            }
                                    
                            beta[t][v][indexOfDoubleindex0] += deltabeta/ probState[prevLeft];
                        }
                            
                        if(emissionType[SScpair[prevLeft]][SScpair[next]] == 3){
                        
                            for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                                ulint transLeft = prevLeft * numberOfSecondarySymbols + ss;
                                ulint nextLeft = transLeft % numberOfStates;
                                if(probState[nextLeft] == 0) continue;
                                    
                                ulint doubleindex1 = nextLeft * numberOfStates + prev;
                                ulint doubleindex3 = prevLeft * numberOfStates + nextLeft;
                                
                                real deltabeta = 0.0;
                                for(int w = 0; w < t; w++){
                                    deltabeta += beta[w][v][indexOfDoubleindex2]
                                                    * alpha[w][w+1][doubleindex3]
                                                    * alpha[w+1][t][doubleindex1];
                                }
                                beta[t][v][indexOfDoubleindex0] += deltabeta
                                                                    / probState[prevLeft]
                                                                    / probState[nextLeft];
                            }
                        }
                    }
                        
                    //
                    // possible production with other possibly long fragments at right of v
                    // current short fragment emits at left
                    //
                        
                    for(ulint nextRight = 0; nextRight < numberOfStates; nextRight++){
                        if(probState[nextRight] == 0
                           || emissionType[SScpair[prev]][SScpair[nextRight]] == 0)
                                continue;
                        ulint doubleindex2 = prev*numberOfStates + nextRight;
                        ulint indexOfDoubleindex2 = betaINDEX.at(doubleindex2);
                            
                        if(emissionType[SScpair[prev]][SScpair[nextRight]] == 2) {
                                
                            ulint doubleindex1 = next*numberOfStates + nextRight;
                        
                            real deltabeta = 0.0;
                            for(int w = windowedLength - 1; w > v; w--){
                                deltabeta += beta[t][w][indexOfDoubleindex2]
                                            * alpha[v][w][doubleindex1];
                            }
                              
                            beta[t][v][indexOfDoubleindex0] += deltabeta / probState[nextRight];
                        }
                                
                        if(emissionType[SScpair[prev]][SScpair[nextRight]] == 3){
  
                            for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                                ulint transRight = nextRight+(numberOfStates * ss);
                                ulint prevRight = transRight / numberOfSecondarySymbols;
                                if (probState[prevRight] == 0) continue;
                                    
                                ulint doubleindex1 = next * numberOfStates + prevRight;
                                ulint doubleindex3 = prevRight * numberOfStates + nextRight;
                                
                                real deltabeta = 0.0;
                                for(int w = windowedLength - 1; w > v; w--){
                                    deltabeta += beta[t][w][indexOfDoubleindex2]
                                                    * alpha[w-1][w][doubleindex3]
                                                    * alpha[v][w-1][doubleindex1];
                                }
                                beta[t][v][indexOfDoubleindex0] += deltabeta
                                                                    / probState[nextRight]
                                                                    / probState[prevRight];
                            }
                        }
                    }
                
                }
            }
             
            normalize(beta[t][v]);

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
            
            if(probState[i] == 0)
                continue;
            real deltagamma = 0.0;
            for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                
                ulint trans = i*numberOfSecondarySymbols + ss;
                ulint next = trans % numberOfStates;
                if(probState[next] == 0)
                    continue;
                if(emissionType[SScpair[i]][SScpair[next]] == 0)
                    continue;
                ulint doubleindex = i*numberOfStates + next;
                //ulint indexOfDoubleindex = betaINDEX.at(doubleindex);
        
                deltagamma += beta[t][t+1][betaINDEX.at(doubleindex)]
                                    * alpha[t][t+1][doubleindex]
                                    / probState[next];
                if(t == windowedLength - 2)
                    gamma[windowedLength - 1][next] += deltagamma / probState[i];
            }
            gamma[t][i] += deltagamma / probState[i];
        }
        normalize(gamma[t]);
    }
    normalize(gamma[windowedLength - 1]);
            
    
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
       
    uint fullLength = (uint) seq.size();
    uint windowedLength = fullLength - windowSize + 2;
        
    // SScpair is intended to be the pair of secondary structure symbols
    // at the two central positions of the State fragment.
    // Secondary structure symbols (0 for E, 1 for H, 2 for C)
    // are obtained from SecondarySymbol%3,
    // and therefore SecondarySymbols must be constructed accordingly (AFPA)
    
  
    
    uint SScpair[numberOfStates];
    
    /*
    ulint number1 = numberOfXtraHalfStates * numberOfSecondarySymbols;
    ulint number2 = numberOfXtraHalfStates;
    ulint number3 = numberOfXtraHalfStates / numberOfSecondarySymbols;
    */
    
    ulint number3 = 1.0;
    for(uint i = 1; i < windowSize - halfWindow - 1; i++)
        number3 *= numberOfSecondarySymbols;
    ulint number2 = number3 * numberOfSecondarySymbols;
    ulint number1 = number2 * numberOfSecondarySymbols;
    
    
    for(uint i = 0; i < numberOfStates; i++){
        uint SScpair1 = i % number1 / number2;
        SScpair1 %= 3;
        uint SScpair2 = i % number2 / number3;
        SScpair2  %=3;
        SScpair[i] = 3 * SScpair1 + SScpair2;
    }
    
    //The Stochastic Context-free grammar is now defined by a 9x9 byte matrix
    //that determines the unique type of emisson for each non-terminal symbol aWu,
    //depending on the central pairs of secondary structure of States â and û.
    //The indexes, from 0 to 8, correspond to pairs EE, EH, EC, HE, ..., CC,
    //in this order on both dimensions. The value in each cell indicates the type of
    //emission, 3 (binary 11) for left/right emission, 2 (binary 10) for
    //left emission, 1 (binay 01) for right emission, and 0 for combinations
    //not expected to appear during any derivation of the grammar.
    
    
    // Double productions for the same secondary structures: XE-W-EX, XH-W-HX, XC-W-CX.
    // Otherwise first left and then right, E > H > C, no EH or HE.
    
      /*
    byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 2, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {1, 0, 0, 0, 0, 0, 3, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 2, 0, 3, 0, 0, 2, 0},
                              {0, 0, 2, 0, 1, 0, 0, 3, 0},
                              {0, 0, 3, 0, 1, 1, 0, 1, 1},
                              {0, 0, 2, 0, 0, 3, 0, 0, 1},
                              {0, 0, 2, 0, 0, 2, 0, 0, 3}};
       */
    
    // No double production for XC-W-CX. Always FIRST LEFT and then at right.
    
     /*
    byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 2, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {1, 0, 0, 0, 0, 0, 2, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 2, 0, 3, 0, 0, 2, 0},
                              {0, 0, 2, 0, 1, 0, 0, 2, 0},
                              {0, 0, 3, 0, 1, 1, 0, 1, 1},
                              {0, 0, 2, 0, 0, 3, 0, 0, 1},
                              {0, 0, 2, 0, 0, 2, 2, 2, 0}};
      */
    
    // No double production for XC-W-CX. Always FIRST RIGHT and then at left.
    
     /*
    byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 2, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {1, 0, 0, 0, 0, 0, 1, 0, 1},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 2, 0, 3, 0, 0, 2, 0},
                              {0, 0, 2, 0, 1, 0, 0, 1, 1},
                              {0, 0, 3, 0, 1, 1, 0, 1, 1},
                              {0, 0, 2, 0, 0, 3, 0, 0, 1},
                              {0, 0, 2, 0, 0, 2, 0, 0, 0}};
      */
    
    // No double production for XH-W-HX. Always FIRST LEFT and then right.
    
     /*
    byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 2, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {1, 0, 0, 0, 0, 0, 3, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 2, 0, 0, 2, 0, 0, 0},
                              {0, 0, 2, 0, 1, 1, 0, 3, 0},
                              {0, 0, 3, 0, 1, 1, 0, 1, 1},
                              {0, 0, 2, 0, 0, 2, 0, 0, 1},
                              {0, 0, 2, 0, 0, 2, 0, 0, 3}};
     */
    
    // No double production for XH-W-HX. Always FIRST RIGHT and then left.
    
     /*
    byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 2, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {1, 0, 0, 0, 0, 0, 3, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 2, 0, 0, 0, 0, 2, 0},
                              {0, 0, 2, 0, 0, 0, 0, 3, 0},
                              {0, 0, 3, 0, 1, 1, 0, 1, 1},
                              {0, 0, 2, 0, 1, 1, 0, 2, 1},
                              {0, 0, 2, 0, 0, 2, 0, 0, 3}};
     */
    
    
    // No double production for XC-W-CX nor XH-W-HX . Always FIRST LEFT and then right.
    
    
     /*
    byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 2, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {1, 0, 0, 0, 0, 0, 2, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 2, 0, 0, 2, 0, 0, 0},
                              {0, 0, 2, 0, 1, 1, 0, 2, 0},
                              {0, 0, 3, 0, 1, 1, 0, 1, 1},
                              {0, 0, 2, 0, 0, 2, 0, 0, 1},
                              {0, 0, 2, 0, 0, 2, 2, 2, 0}};
      */
    
    // No double production for XC-W-CX nor XH-W-HX . Always FIRST RIGHT and then left.
    
    
     /*
    byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 2, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {1, 0, 0, 0, 0, 0, 1, 0, 1},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 2, 0, 0, 0, 0, 2, 0},
                              {0, 0, 2, 0, 0, 0, 0, 1, 1},
                              {0, 0, 3, 0, 1, 1, 0, 1, 1},
                              {0, 0, 2, 0, 1, 1, 0, 2, 1},
                              {0, 0, 2, 0, 0, 2, 0, 0, 0}};
     */
    
    
    // No H anywhere !!! Double E and C
    
    /*
   byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 2, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {1, 0, 0, 0, 0, 0, 3, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 3, 0, 0, 0, 0, 0, 1},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 2, 0, 0, 0, 0, 0, 3}};
     */
    
    
    
     /*
   // Everybody emits at left, equivalent to the HMM

   byte emissionType[9][9] ={{2, 2, 2, 2, 2, 2, 2, 2, 2},
                             {2, 2, 2, 2, 2, 2, 2, 2, 2},
                             {2, 2, 2, 2, 2, 2, 2, 2, 2},
                             {2, 2, 2, 2, 2, 2, 2, 2, 2},
                             {2, 2, 2, 2, 2, 2, 2, 2, 2},
                             {2, 2, 2, 2, 2, 2, 2, 2, 2},
                             {2, 2, 2, 2, 2, 2, 2, 2, 2},
                             {2, 2, 2, 2, 2, 2, 2, 2, 2},
                             {2, 2, 2, 2, 2, 2, 2, 2, 2}};
       */
    
   
    
      /*
    // Everybody emits at right, also equivalent to the HMM
   
   byte emissionType[9][9] ={{1, 1, 1, 1, 1, 1, 1, 1, 1},
                             {1, 1, 1, 1, 1, 1, 1, 1, 1},
                             {1, 1, 1, 1, 1, 1, 1, 1, 1},
                             {1, 1, 1, 1, 1, 1, 1, 1, 1},
                             {1, 1, 1, 1, 1, 1, 1, 1, 1},
                             {1, 1, 1, 1, 1, 1, 1, 1, 1},
                             {1, 1, 1, 1, 1, 1, 1, 1, 1},
                             {1, 1, 1, 1, 1, 1, 1, 1, 1},
                             {1, 1, 1, 1, 1, 1, 1, 1, 1}};
      */
    
    
    
    //  /*
    //Everybody emits at left/right, again, also equivalent to the HMM
    
    byte emissionType[9][9] ={{3, 3, 3, 3, 3, 3, 3, 3, 3},
                              {3, 3, 3, 3, 3, 3, 3, 3, 3},
                              {3, 3, 3, 3, 3, 3, 3, 3, 3},
                              {3, 3, 3, 3, 3, 3, 3, 3, 3},
                              {3, 3, 3, 3, 3, 3, 3, 3, 3},
                              {3, 3, 3, 3, 3, 3, 3, 3, 3},
                              {3, 3, 3, 3, 3, 3, 3, 3, 3},
                              {3, 3, 3, 3, 3, 3, 3, 3, 3},
                              {3, 3, 3, 3, 3, 3, 3, 3, 3}};
    //   */
    
    
    
    ulint numberOfDoubleStates = numberOfStates * numberOfStates;
 //   int numberOfSecondarySymbols2 = numberOfSecondarySymbols*numberOfSecondarySymbols;
    
    vector < vector < vector<real> > > alpha(windowedLength, vector < vector<real> > (windowedLength, vector <real> (numberOfDoubleStates, 0.0)));
    vector < vector < vector<real> > > beta(windowedLength, vector < vector<real> > (windowedLength, vector <real> (numberOfDoubleStates, 0.0)));
    
    
    //
    // 1. Inside
    //
    
    //
    // 1.1 Initialization
    //
    
  /*
    for(int t=0; t < windowedLength; t++){
        for(ulint i=0; i<numberOfStates; i++){
            ulint doubleindex = i * numberOfStates + i;
            alpha[t][t][doubleindex] = 1.0;
        }
        normalize(alpha[t][t]);
    }
   */
    
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
        
        for (int t=v-2; t >= 0; t--) {
        
            for(ulint prevLeft=0; prevLeft<numberOfStates; prevLeft++){
                if(probState[prevLeft] == 0) continue;
                for(ulint nextRight=0; nextRight<numberOfStates; nextRight++){
                    if(probState[nextRight] == 0) continue;
                    
                    byte emission = emissionType[SScpair[prevLeft]][SScpair[nextRight]];
                    if(emission == 0) continue;
                    
                    ulint doubleIndex = prevLeft * numberOfStates + nextRight;
                    
                        
                    if ( emission == 1 ){
                        for(uint tt=0; tt<numberOfSecondarySymbols; tt++){
                                
                            ulint transRight = tt * numberOfStates + nextRight;
                            ulint prevRight = transRight / numberOfSecondarySymbols;
                            if(probState[prevRight] == 0) continue;
                            ulint doubleIndexRight = prevRight * numberOfStates + nextRight;
                            ulint childDoubleIndex = prevLeft*numberOfStates + prevRight;
                            
                              //  if(probState[prevRight] != 0){
                            real deltaAlpha = alpha[t][v-1][childDoubleIndex] *
                                                alpha[v-1][v][doubleIndexRight] /
                                                probState[prevRight];
                            alpha[t][v][doubleIndex] += deltaAlpha;
                             //   }
                        }
                    }
                    
                    else if ( emission == 2 ){
                        for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                                
                            ulint transLeft = prevLeft * numberOfSecondarySymbols + ss;
                            ulint nextLeft = transLeft % numberOfStates;
                            if(probState[nextLeft] == 0) continue;
                            ulint doubleIndexLeft = prevLeft*numberOfStates + nextLeft;
                            ulint childDoubleIndex = nextLeft*numberOfStates + nextRight;
                                
                           // if( probState[nextLeft] != 0) {
                            real deltaAlpha = alpha[t+1][v][childDoubleIndex] *
                                    alpha[t][t+1][doubleIndexLeft] /
                                    probState[nextLeft];
                            alpha[t][v][doubleIndex] += deltaAlpha;
                           // }
                        }
                    }
                            
                    else if (emission == 3) { //} && t < (v - 2)) {
                        
                        for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                                
                            ulint transLeft = prevLeft * numberOfSecondarySymbols + ss;
                            ulint nextLeft = transLeft % numberOfStates;
                            if(probState[nextLeft] == 0) continue;
                            ulint doubleIndexLeft = prevLeft*numberOfStates + nextLeft;
                                        
                            for(uint tt=0; tt<numberOfSecondarySymbols; tt++){
                                ulint transRight = tt * numberOfStates + nextRight;
                                ulint prevRight = transRight / numberOfSecondarySymbols;
                                if(probState[prevRight] == 0) continue;
                                ulint doubleIndexRight = prevRight * numberOfStates + nextRight;
                                
                                real deltaAlpha = 0.0;
                                if( (v - t) > 2){ //} &&
                              //     probState[nextLeft] != 0 && probState[prevRight] != 0) {
                                    ulint childDoubleIndex = nextLeft*numberOfStates + prevRight;
                                    deltaAlpha = alpha[t+1][v-1][childDoubleIndex]
                                                    * alpha[t][t+1][doubleIndexLeft]
                                                    / probState[nextLeft]
                                                    * alpha[v-1][v][doubleIndexRight]
                                                    / probState[prevRight];
                                }
                                else if(nextLeft == prevRight) {
                                    ulint childDoubleIndex = nextLeft*numberOfStates + nextRight;
                                    deltaAlpha =  alpha[t+1][v][childDoubleIndex]
                                                    * alpha[t][t+1][doubleIndexLeft]
                                                    / probState[nextLeft];
                                }
                                alpha[t][v][doubleIndex] += deltaAlpha;
                             //  }
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
                if(emissionType[SScpair[i]][SScpair[j]] > 0)
                        beta[0][windowedLength-1][i*numberOfStates + j] =
                            probState_head[i] * probState_tail[j];
    }
    else{
        for(ulint i=0; i < numberOfStates; i++)
            for(ulint j=0; j < numberOfStates; j++)
                if(emissionType[SScpair[i]][SScpair[j]] > 0)
                        beta[0][windowedLength-1][i*numberOfStates + j] = 1.0;
    }
    normalize(beta[0][windowedLength-1]);
    
    // 2.2 Induction
    
    for(int t=0; t < windowedLength-1; t++){
        
        for (int v=windowedLength-1; v > t; v--){
 
            if(t==0 && v==windowedLength-1)
                continue;
            
            for(ulint prevLeft=0; prevLeft<numberOfStates; prevLeft++){
                if(probState[prevLeft] == 0) continue;
                
                for(ulint nextRight=0; nextRight<numberOfStates; nextRight++){
                    if(probState[nextRight] == 0) continue;
                    
                    real emission = emissionType[SScpair[prevLeft]][SScpair[nextRight]];
                    if(emission == 0) continue;
                    
                    ulint parentDoubleIndex = prevLeft * numberOfStates + nextRight;
                    
                    if(emission == 1
                            && v < windowedLength - 1)
                    {
                        for(uint tt=0; tt<numberOfSecondarySymbols; tt++)
                        {
                            ulint transRight = tt * numberOfStates + nextRight;
                            ulint prevRight = transRight /  numberOfSecondarySymbols;
                            if(probState[prevRight] == 0) continue;
                            ulint doubleIndexRight = prevRight * numberOfStates + nextRight;
                             
                            ulint doubleIndex = prevLeft*numberOfStates + prevRight;
                         //   if(probState[nextRight] != 0)
                         //   {
                                real deltaBeta = beta[t][v+1][parentDoubleIndex] *
                                                alpha[v][v+1][doubleIndexRight] /
                                                probState[nextRight];
                                beta[t][v][doubleIndex] += deltaBeta;
                          //  }
                        }
                    }
                    
                    else if(emission == 2
                            && t > 0)
                    {
                        for(uint ss=0; ss<numberOfSecondarySymbols; ss++)
                        {
                            ulint transLeft = prevLeft * numberOfSecondarySymbols + ss;
                            ulint nextLeft = transLeft % numberOfStates;
                            if(probState[nextLeft] == 0) continue;
                            ulint doubleIndexLeft = prevLeft * numberOfStates + nextLeft;
                                  
                            ulint doubleIndex = nextLeft*numberOfStates + nextRight;
                          //  if(probState[prevLeft] != 0)
                          //  {
                                real deltaBeta = beta[t-1][v][parentDoubleIndex] *
                                                alpha[t-1][t][doubleIndexLeft] /
                                                probState[prevLeft];
                                beta[t][v][doubleIndex] += deltaBeta;
                           // }
                        }
                    }
                    
                    else if(emission == 3
                       && t > 0 && v < windowedLength -1)
                    {
                        for(uint ss=0; ss<numberOfSecondarySymbols; ss++)
                        {
                            ulint transLeft = prevLeft * numberOfSecondarySymbols + ss;
                            ulint nextLeft = transLeft % numberOfStates;
                            if(probState[nextLeft] == 0) continue;
                            ulint doubleIndexLeft = prevLeft * numberOfStates + nextLeft;
                                    
                            for(uint tt=0; tt<numberOfSecondarySymbols; tt++)
                            {
                                ulint transRight = tt * numberOfStates + nextRight;
                                ulint prevRight = transRight / numberOfSecondarySymbols;
                                if(probState[prevRight] == 0) continue;
                                ulint doubleIndexRight = prevRight * numberOfStates + nextRight;
                                 
                                ulint doubleIndex = nextLeft*numberOfStates + prevRight;
                                
                             //   if(probState[prevLeft] != 0 && probState[nextRight] != 0)
                             //   {
                                    
                                real deltaBeta =  alpha[t-1][t][doubleIndexLeft] *
                                                    alpha[v][v+1][doubleIndexRight] *
                                                    beta[t-1][v+1][parentDoubleIndex] /
                                                    probState[prevLeft] /
                                                    probState[nextRight];
                                if(v == t+1 && nextLeft == prevRight){
                                    deltaBeta += alpha[t-1][t][doubleIndexLeft] *
                                                    beta[t-1][v][parentDoubleIndex] /
                                                    probState[prevLeft];
                                }
                                
                                beta[t][v][doubleIndex] += deltaBeta;
                               // }
                                
                            }
                        }
                    }
                }
            }
            
            if (v == t+1){
                
                for(ulint trans=0; trans < numberOfFragments; trans++){
                    if(probFragment[trans] == 0)
                        continue;
                    ulint prev = trans / numberOfSecondarySymbols;
                    ulint next = trans % numberOfStates;
                    ulint doubleindex0 = prev * numberOfStates + next;
            
                    //
                    // possible production with other possibly long fragments at left of t
                    // current short fragment emits at right
                    //
                    
                    for(ulint prevLeft = 0; prevLeft < numberOfStates; prevLeft++){
                        if(probState[prevLeft] == 0)
                            continue;
                        ulint doubleindex2 = prevLeft*numberOfStates + next;
                            
                        if(emissionType[SScpair[prevLeft]][SScpair[next]] == 1) {
                            
                            ulint doubleindex1 = prevLeft*numberOfStates + prev;
                            
                            real deltabeta = 0.0;
                            for(int w = 0; w < t; w++){
                                deltabeta += beta[w][v][doubleindex2] *
                                                alpha[w][t][doubleindex1];
                            }
                            
                            beta[t][v][doubleindex0] += deltabeta / probState[prevLeft];
                        }
                            
                        if(emissionType[SScpair[prevLeft]][SScpair[next]] == 3){
                            
                            for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                                ulint transLeft = prevLeft * numberOfSecondarySymbols + ss;
                                ulint nextLeft = transLeft % numberOfStates;
                                if(probState[nextLeft] == 0) continue;
                                
                                ulint doubleindex1 = nextLeft * numberOfStates + prev;
                                ulint doubleindex3 = prevLeft * numberOfStates + nextLeft;
                            
                                real deltabeta = 0.0;
                                for(int w = 0; w < t; w++){
                                    deltabeta += beta[w][v][doubleindex2]
                                                    * alpha[w][w+1][doubleindex3]
                                                    * alpha[w+1][t][doubleindex1];
                                }
                                beta[t][v][doubleindex0] += deltabeta
                                                            / probState[prevLeft]
                                                            / probState[nextLeft];
                            }
                        }
                    }
                
                    //
                    // possible production with other possibly long fragments at right of v
                    // current short fragment emits at left
                    //
                        
                    for(ulint nextRight = 0; nextRight < numberOfStates; nextRight++){
                        if(probState[nextRight] == 0)
                            continue;
                        ulint doubleindex2 = prev*numberOfStates + nextRight;
                                
                        if(emissionType[SScpair[prev]][SScpair[nextRight]] == 2) {
                                
                            ulint doubleindex1 = next*numberOfStates + nextRight;
                            
                            real deltabeta = 0.0;
                            for(int w = windowedLength - 1; w > v; w--){
                                deltabeta += beta[t][w][doubleindex2]
                                                * alpha[v][w][doubleindex1];
                            }
                            
                            beta[t][v][doubleindex0] += deltabeta / probState[nextRight];;
                        }
                                
                        if(emissionType[SScpair[prev]][SScpair[nextRight]] == 3){
                                
                            for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                                ulint transRight = nextRight+(numberOfStates * ss);
                                ulint prevRight = transRight / numberOfSecondarySymbols;
                                if (probState[prevRight] == 0) continue;
                                    
                                ulint doubleindex1 = next * numberOfStates + prevRight;
                                ulint doubleindex3 = prevRight * numberOfStates + nextRight;
                                
                                real deltabeta = 0.0;
                                for(int w = windowedLength - 1; w > v; w--){
                                    deltabeta += beta[t][w][doubleindex2]
                                                        * alpha[w-1][w][doubleindex3]
                                                       // / probState[nextRight]
                                                        * alpha[v][w-1][doubleindex1];
                                                        // / probState[prevRight];
                                }
                                beta[t][v][doubleindex0] += deltabeta
                                                            / probState[nextRight]
                                                            / probState[prevRight];
                                if((next == prevRight) && (v < windowedLength - 1)){
                                    beta[t][v][doubleindex0] += beta[t][v+1][doubleindex2]
                                                                    * alpha[v][v+1][doubleindex3]
                                                                    / probState[nextRight];
                                }
                            }
                        }
                    }
                
                }
            }
             
            normalize(beta[t][v]);

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
            if(probState[i] == 0)
                continue;
            real deltagamma = 0.0;
            for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                
                ulint trans = i*numberOfSecondarySymbols + ss;
                ulint next = trans % numberOfStates;
                if(probState[next] == 0)
                    continue;
                ulint doubleindex = i*numberOfStates + next;

                deltagamma += beta[t][t+1][doubleindex] * alpha[t][t+1][doubleindex]
                                / probState[next];
                if(t == windowedLength - 2)
                    gamma[windowedLength - 1][next] += deltagamma / probState[i];
                
            }
            gamma[t][i] += deltagamma / probState[i];
        }
        normalize(gamma[t]);
    }
    normalize(gamma[windowedLength - 1]);
            
    
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




void Hmm::insideOutsideV5(vector<byte>& seq){
       
    uint fullLength = (uint) seq.size();
    uint windowedLength = fullLength - windowSize + 2;
    
    // SScpair is intended to be the pair of secondary structure symbols
    // at the two central positions of the State fragment.
    // Secondary structure symbols (0 for E, 1 for H, 2 for C)
    // are obtained from SecondarySymbol%3,
    // and therefore SecondarySymbols must be constructed accordingly (AFPA)
    
 //   uint SScpair[numberOfStates];
    uint LeftPairSS[numberOfStates];
    uint RightPairSS[numberOfStates];
    
 //   ulint number1 = numberOfXtraHalfStates * numberOfSecondarySymbols;
 //   ulint number2 = numberOfXtraHalfStates;
 //   ulint number3 = numberOfXtraHalfStates / numberOfSecondarySymbols;
    
  //  if(windowSize % 2 == 0) {
  //      number1 /= numberOfSecondarySymbols;
  //      number2 /= numberOfSecondarySymbols;
  //      number3 /= numberOfSecondarySymbols;
  //  }

    ulint number1 = 1.0;
    for (uint i=1; i < (windowSize - halfWindow - 1); i++)
        number1 *= numberOfSecondarySymbols;
    ulint number2 = number1 * numberOfSecondarySymbols;
    ulint number3 = number2 * numberOfSecondarySymbols;
    ulint number4 = number3 * numberOfSecondarySymbols;
    
//    for(uint i = 0; i < numberOfStates; i++){
//        uint SScpair1 = i % number1 / number2;
//        SScpair1 %= 3;
//        uint SScpair2 = i % number2 / number3;
//        SScpair2  %=3;
//        SScpair[i] = 3 * SScpair1 + SScpair2;
//    }
    
    for(uint i = 0; i < numberOfStates; i++){
        uint LeftPairSS1 = i % number3 / number2;
        LeftPairSS1 %= 3;
        uint LeftPairSS2 = i % number2 / number1;
        LeftPairSS2  %= 3;
        LeftPairSS[i] = 3 * LeftPairSS1 + LeftPairSS2;
    //    if(windowSize % 2){
            RightPairSS[i] = LeftPairSS[i];
      //  }
      //  else{
      //      uint RightPairSS1 = i % number4 / number3;
      //      RightPairSS1 %= 3;
      //      uint RightPairSS2 = i % number3 / number2;
      //      RightPairSS2 %= 3;
      //      RightPairSS[i] = 3 * RightPairSS1 + RightPairSS2;
      //  }
    }
    
    
    //The Stochastic Context-free grammar is now defined by a 9x9 byte matrix
    //that determines the unique type of emisson for each non-terminal symbol aWu,
    //depending on the central pairs of secondary structure of States â and û.
    //The indexes, from 0 to 8, correspond to pairs EE, EH, EC, HE, ..., CC,
    //in this order on both dimensions. The value in each cell indicates the type of
    //emission, 3 (binary 11) for left/right emission, 2 (binary 10) for
    //left emission, 1 (binay 01) for right emission, and 0 for combinations
    //not expected to appear during any derivation of the grammar.
    
    
    // Double productions for the same secondary structures: XE-W-EX, XH-W-HX, XC-W-CX.
    // Otherwise first left and then right, E > H > C, no EH or HE.
    
      /*
    byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 2, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {1, 0, 0, 0, 0, 0, 3, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 2, 0, 3, 0, 0, 2, 0},
                              {0, 0, 2, 0, 1, 0, 0, 3, 0},
                              {0, 0, 3, 0, 1, 1, 0, 1, 1},
                              {0, 0, 2, 0, 0, 3, 0, 0, 1},
                              {0, 0, 2, 0, 0, 2, 0, 0, 3}};
       */
    
       /*
     byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 3, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 2, 0, 3, 0, 0, 2, 0},
                               {0, 0, 2, 0, 1, 0, 0, 3, 0},
                               {0, 0, 3, 0, 1, 1, 0, 1, 1},
                               {0, 0, 2, 0, 0, 3, 0, 0, 1},
                               {0, 0, 2, 0, 0, 2, 0, 0, 3}};
        */
    
    // No double production for XC-W-CX. Always FIRST LEFT and then at right.
    
     /*
    byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 2, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {1, 0, 0, 0, 0, 0, 2, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 2, 0, 3, 0, 0, 2, 0},
                              {0, 0, 2, 0, 1, 0, 0, 2, 0},
                              {0, 0, 3, 0, 1, 1, 0, 1, 1},
                              {0, 0, 2, 0, 0, 3, 0, 0, 1},
                              {0, 0, 2, 0, 0, 2, 2, 2, 0}};
      */
    
    // No double production for XC-W-CX. Always FIRST RIGHT and then at left.
    
     /*
    byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 2, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {1, 0, 0, 0, 0, 0, 1, 0, 1},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 2, 0, 3, 0, 0, 2, 0},
                              {0, 0, 2, 0, 1, 0, 0, 1, 1},
                              {0, 0, 3, 0, 1, 1, 0, 1, 1},
                              {0, 0, 2, 0, 0, 3, 0, 0, 1},
                              {0, 0, 2, 0, 0, 2, 0, 0, 0}};
      */
    
    // No double production for XH-W-HX. Always FIRST LEFT and then right.
    
     /*
    byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 2, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {1, 0, 0, 0, 0, 0, 3, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 2, 0, 0, 2, 0, 0, 0},
                              {0, 0, 2, 0, 1, 1, 0, 3, 0},
                              {0, 0, 3, 0, 1, 1, 0, 1, 1},
                              {0, 0, 2, 0, 0, 2, 0, 0, 1},
                              {0, 0, 2, 0, 0, 2, 0, 0, 3}};
     */
    
    // No double production for XH-W-HX. Always FIRST RIGHT and then left.
    
     /*
    byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 2, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {1, 0, 0, 0, 0, 0, 3, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 2, 0, 0, 0, 0, 2, 0},
                              {0, 0, 2, 0, 0, 0, 0, 3, 0},
                              {0, 0, 3, 0, 1, 1, 0, 1, 1},
                              {0, 0, 2, 0, 1, 1, 0, 2, 1},
                              {0, 0, 2, 0, 0, 2, 0, 0, 3}};
     */
    
    
    // No double production for XC-W-CX nor XH-W-HX . Always FIRST LEFT and then right.
    
    
     /*
    byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 2, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {1, 0, 0, 0, 0, 0, 2, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 2, 0, 0, 2, 0, 0, 0},
                              {0, 0, 2, 0, 1, 1, 0, 2, 0},
                              {0, 0, 3, 0, 1, 1, 0, 1, 1},
                              {0, 0, 2, 0, 0, 2, 0, 0, 1},
                              {0, 0, 2, 0, 0, 2, 2, 2, 0}};
      */
    
    // No double production for XC-W-CX nor XH-W-HX . Always FIRST RIGHT and then left.
    
    
     /*
    byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 2, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {1, 0, 0, 0, 0, 0, 1, 0, 1},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 2, 0, 0, 0, 0, 2, 0},
                              {0, 0, 2, 0, 0, 0, 0, 1, 1},
                              {0, 0, 3, 0, 1, 1, 0, 1, 1},
                              {0, 0, 2, 0, 1, 1, 0, 2, 1},
                              {0, 0, 2, 0, 0, 2, 0, 0, 0}};
     */
    
    
    // No H anywhere !!! Double E and C
    
    /*
   byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 2, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {1, 0, 0, 0, 0, 0, 3, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 3, 0, 0, 0, 0, 0, 1},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 2, 0, 0, 0, 0, 0, 3}};
     */
    
    
    
     /*
   // Everybody emits at left, equivalent to the HMM

   byte emissionType[9][9] ={{2, 2, 2, 2, 2, 2, 2, 2, 2},
                             {2, 2, 2, 2, 2, 2, 2, 2, 2},
                             {2, 2, 2, 2, 2, 2, 2, 2, 2},
                             {2, 2, 2, 2, 2, 2, 2, 2, 2},
                             {2, 2, 2, 2, 2, 2, 2, 2, 2},
                             {2, 2, 2, 2, 2, 2, 2, 2, 2},
                             {2, 2, 2, 2, 2, 2, 2, 2, 2},
                             {2, 2, 2, 2, 2, 2, 2, 2, 2},
                             {2, 2, 2, 2, 2, 2, 2, 2, 2}};
       */
    
   
    
      /*
    // Everybody emits at right, also equivalent to the HMM
   
   byte emissionType[9][9] ={{1, 1, 1, 1, 1, 1, 1, 1, 1},
                             {1, 1, 1, 1, 1, 1, 1, 1, 1},
                             {1, 1, 1, 1, 1, 1, 1, 1, 1},
                             {1, 1, 1, 1, 1, 1, 1, 1, 1},
                             {1, 1, 1, 1, 1, 1, 1, 1, 1},
                             {1, 1, 1, 1, 1, 1, 1, 1, 1},
                             {1, 1, 1, 1, 1, 1, 1, 1, 1},
                             {1, 1, 1, 1, 1, 1, 1, 1, 1},
                             {1, 1, 1, 1, 1, 1, 1, 1, 1}};
      */
    
    
    
   //   /*
    //Everybody emits at left/right, again, also equivalent to the HMM
    
    byte emissionType[9][9] ={{3, 3, 3, 3, 3, 3, 3, 3, 3},
                              {3, 3, 3, 3, 3, 3, 3, 3, 3},
                              {3, 3, 3, 3, 3, 3, 3, 3, 3},
                              {3, 3, 3, 3, 3, 3, 3, 3, 3},
                              {3, 3, 3, 3, 3, 3, 3, 3, 3},
                              {3, 3, 3, 3, 3, 3, 3, 3, 3},
                              {3, 3, 3, 3, 3, 3, 3, 3, 3},
                              {3, 3, 3, 3, 3, 3, 3, 3, 3},
                              {3, 3, 3, 3, 3, 3, 3, 3, 3}};
     //   */
    
    
    
    ulint numberOfDoubleStates = numberOfStates * numberOfStates;
//    int numberOfSecondarySymbols2 = numberOfSecondarySymbols*numberOfSecondarySymbols;
    
    // map<ulint, ulint> INDEX;
    unordered_map<ulint, ulint> betaIndex;
    
    ulint maxIndex = 0;
    for(ulint i = 0; i < numberOfStates; i++){
        for(ulint j = 0; j < numberOfStates; j++){
            if(emissionType[LeftPairSS[i]][RightPairSS[j]]) {
                ulint doubleIndex = i*numberOfStates +j;
                betaIndex[doubleIndex] = maxIndex;
                maxIndex ++;
            }
        }
    }
    
    vector < vector < vector<real> > > alpha(windowedLength, vector < vector<real> > (windowedLength, vector <real> (numberOfDoubleStates, 0.0)));
  //  vector < vector < vector<real> > > beta(windowedLength, vector < vector<real> > (windowedLength+1, vector <real> (numberOfDoubleStates, 0.0)));
    
    
 //   vector < vector < vector<real> > > alpha(windowedLength, vector < vector<real> > (windowedLength+1, vector <real> (maxIndex, 0.0)));
    vector < vector < vector<real> > > beta(windowedLength, vector < vector<real> > (windowedLength, vector <real> (maxIndex, 0.0)));
    
 //   vector < vector <real> > alphaTTplusOne (windowedLength,
 //                                            vector <real> (numberOfFragments, 0.0));
    vector < vector <real> > betaTTPlusOne (windowedLength,
                                            vector <real> (numberOfFragments, 0.0));
 
 //   for(ulint i = 0; i< numberOfStates; i++){
 //       if(SScpair[i] == 1 || SScpair[i] == 3)
 //           probState[i] = 0;
 //   }
    
    //
    //Removes probabilities for states containing any 'EH' or 'HE' pair of secondary structures,
    //which are inconsistent with the grammar but might appear rarely in the training data bank.
    //
    
  // /*
    for(ulint i = 0; i < numberOfStates; i++){
        if(probState[i] == 0) continue;
        ulint n1=numberOfSecondarySymbols;
        for(uint t = 2; t < windowSize; t++){
            n1 *= numberOfSecondarySymbols;
            ulint n2 = n1 / numberOfSecondarySymbols;
            ulint n3 = n2 / numberOfSecondarySymbols;
            ulint a = i % n1 / n2 % 3;
            ulint b = i % n2 / n3 % 3;
            if((a == 0 && b == 1) || (a ==1 &&  b==0)){
                probState[i] = 0;
                continue;
            }
        }
        normalize(probState);
    }
  //  */
    
    //
    // 1. Inside
    //
    
    //
    // 1.1 Initialization
    //
    
  /*
    for(int t=0; t < windowedLength; t++){
        for(ulint i=0; i<numberOfStates; i++){
            ulint doubleindex = i * numberOfStates + i;
            alpha[t][t][doubleindex] = 1.0;
        }
        normalize(alpha[t][t]);
    }
   */
    
    for(int t=0; t < windowedLength-1; t++){
        for(ulint i=0; i<numberOfFragments; i++){
            ulint left = i / numberOfSecondarySymbols;
            ulint right = i % numberOfStates;
            ulint doubleIndex = left * numberOfStates + right;
            alpha[t][t+1][doubleIndex] = probFragment[i] *
                        probFragmentEmitsPrimarySymbol[i][ seq[t+halfWindow] ];
     //       alphaTTplusOne[t][i] = alpha[t][t+1][doubleIndex];
        }
        normalize(alpha[t][t+1]);
     //   normalize(alphaTTplusOne[t]);
    }
    
    
    
    //
    // 1.2 Induction
    //
    for(int v=2; v < windowedLength; v++){
      
        for (int t=v-2; t >= 0; t--) {
            
            for(ulint nextLeft=0; nextLeft<numberOfStates; nextLeft++){
                if(probState[nextLeft] == 0) continue;
                for(ulint prevRight=0; prevRight<numberOfStates; prevRight++){
                    if(probState[prevRight] == 0) continue;
                    
                    byte emission = emissionType[LeftPairSS[nextLeft]][RightPairSS[prevRight]];
                    if(emission == 0) continue;
                    
                    ulint parentDoubleIndex = nextLeft * numberOfStates + prevRight;
                
                    if( emission == 1){
                        
                        for(uint tt=0; tt<numberOfSecondarySymbols; tt++){
                                
                            ulint transRight = prevRight * numberOfSecondarySymbols + tt;
                            ulint nextRight = transRight % numberOfStates;
                            if(probState[nextRight] == 0) continue;
                            ulint doubleIndexRight = prevRight * numberOfStates + nextRight;
                            ulint doubleIndex = nextLeft*numberOfStates + nextRight;
                            
                         //   if(probState[prevRight] != 0) {
                            real deltaAlpha = alpha[t][v-1][parentDoubleIndex] *
                                            alpha[v-1][v][doubleIndexRight] /
                                            probState[prevRight];
                            alpha[t][v][doubleIndex] += deltaAlpha;
                           // }
                        }
                    }
            
                    else if(emission ==2 ) {

                        for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                                
                            ulint transLeft = ss * numberOfStates + nextLeft;
                            ulint prevLeft = transLeft / numberOfSecondarySymbols;
                            if(probState[prevLeft] == 0) continue;
                            ulint doubleIndexLeft = prevLeft*numberOfStates + nextLeft;
                            ulint doubleIndex = prevLeft*numberOfStates + prevRight;
                                
                           // if( probState[nextLeft] != 0) {
                            real deltaAlpha = alpha[t+1][v][parentDoubleIndex] *
                                    alpha[t][t+1][doubleIndexLeft] /
                                    probState[nextLeft];
                            alpha[t][v][doubleIndex] += deltaAlpha;
                           // }
                        }
                    }
                    
                    else if (emission == 3 && t < (v - 2)) {
                            
                        for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                                
                            ulint transLeft = ss * numberOfStates + nextLeft;
                            ulint prevLeft = transLeft / numberOfSecondarySymbols;
                            if(probState[prevLeft] == 0) continue;
                            ulint doubleIndexLeft = prevLeft*numberOfStates + nextLeft;
                                        
                            for(uint tt=0; tt<numberOfSecondarySymbols; tt++){
                                ulint transRight = prevRight * numberOfSecondarySymbols + tt;
                                ulint nextRight = transRight % numberOfStates;
                                if(probState[nextRight] == 0) continue;
                                ulint doubleIndexRight = prevRight * numberOfStates + nextRight;
                                
                                ulint doubleIndex = prevLeft*numberOfStates + nextRight;
                                    
                              //  if( (v - t) > 2){ //} &&
                              //     probState[nextLeft] != 0 && probState[prevRight] != 0) {
                                        real deltaAlpha = alpha[t+1][v-1][parentDoubleIndex]
                                                            * alpha[t][t+1][doubleIndexLeft]
                                                            / probState[nextLeft]
                                                            * alpha[v-1][v][doubleIndexRight]
                                                            / probState[prevRight];
                                        alpha[t][v][doubleIndex] += deltaAlpha;
                             //  }
                            }
                        }
                    }
                    
              //      if(emission != 3) continue;
                
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
    // Maybe not necessary anymore in version 5.
    //
    
    if(headAndTail){
        for(ulint i=0; i < numberOfStates; i++){
            for(ulint j=0; j < numberOfStates; j++){
                if(emissionType[LeftPairSS[i]][RightPairSS[j]]){
                    ulint doubleIndex = i*numberOfStates + j;
                    beta[0][windowedLength-1][betaIndex.at(doubleIndex)] =
                            probState_head[i] * probState_tail[j];
                }
            }
        }
    }
    else{
        for(ulint i=0; i < numberOfStates; i++){
            for(ulint j=0; j < numberOfStates; j++){
                if(emissionType[LeftPairSS[i]][RightPairSS[j]]){
                    ulint doubleIndex = i*numberOfStates + j;
                    beta[0][windowedLength-1][betaIndex.at(doubleIndex)] = 1.0;
                }
            }
        }
    }
    normalize(beta[0][windowedLength-1]);
    
    // 2.2 Induction
    
    for(int t=0; t < windowedLength-1; t++){
    
        for (int v=windowedLength-1; v > t; v--){
 
            if(t==0 && v==windowedLength-1)
                continue;
            
            for(ulint prevLeft=0; prevLeft<numberOfStates; prevLeft++){
                if(probState[prevLeft] == 0) continue;
                
                for(ulint nextRight=0; nextRight<numberOfStates; nextRight++){
                    if(probState[nextRight] == 0) continue;
                    
                    real emission = emissionType[LeftPairSS[prevLeft]][RightPairSS[nextRight]];
                    if(emission == 0) continue;
                    
                    ulint parentDoubleIndex = prevLeft * numberOfStates + nextRight;
                    ulint indexOfParentDoubleIndex = betaIndex.at(parentDoubleIndex);
                    
                    if(emission == 1
                            && v < windowedLength - 1)
                    {
                        for(uint tt=0; tt<numberOfSecondarySymbols; tt++)
                        {
                            ulint transRight = tt * numberOfStates + nextRight;
                            ulint prevRight = transRight /  numberOfSecondarySymbols;
                            if(probState[prevRight] == 0) continue;
                            ulint doubleIndexRight = prevRight * numberOfStates + nextRight;
                            
                            if(emissionType[LeftPairSS[prevLeft]][RightPairSS[prevRight]] == 0)
                                continue;
                            ulint doubleIndex = prevLeft*numberOfStates + prevRight;
                          //  if(probState[nextRight] != 0)
                          //  {
                                real deltaBeta = beta[t][v+1][indexOfParentDoubleIndex] *
                                                alpha[v][v+1][doubleIndexRight] /
                                                probState[nextRight];
                                beta[t][v][betaIndex.at(doubleIndex)] += deltaBeta;
                           // }
                        }
                    }
                    
                    else if(emission == 2
                            && t > 0)
                    {
                        for(uint ss=0; ss<numberOfSecondarySymbols; ss++)
                        {
                            ulint transLeft = prevLeft * numberOfSecondarySymbols + ss;
                            ulint nextLeft = transLeft % numberOfStates;
                            if(probState[nextLeft] == 0) continue;
                            ulint doubleIndexLeft = prevLeft * numberOfStates + nextLeft;
                            
                            if(emissionType[LeftPairSS[nextLeft]][RightPairSS[nextRight]] == 0)
                                continue;
                            ulint doubleIndex = nextLeft*numberOfStates + nextRight;
                           // if(probState[prevLeft] != 0)
                           // {
                                real deltaBeta = beta[t-1][v][indexOfParentDoubleIndex] *
                                                alpha[t-1][t][doubleIndexLeft] /
                                                probState[prevLeft];
                                beta[t][v][betaIndex.at(doubleIndex)] += deltaBeta;
                           // }
                        }
                    }
                    
                    else if(emission == 3
                       && t > 0 && v < windowedLength -1)
                    {
                        for(uint ss=0; ss<numberOfSecondarySymbols; ss++)
                        {
                            ulint transLeft = prevLeft * numberOfSecondarySymbols + ss;
                            ulint nextLeft = transLeft % numberOfStates;
                            if(probState[nextLeft] == 0) continue;
                            ulint doubleIndexLeft = prevLeft * numberOfStates + nextLeft;
                                    
                            for(uint tt=0; tt<numberOfSecondarySymbols; tt++)
                            {
                                ulint transRight = tt * numberOfStates + nextRight;
                                ulint prevRight = transRight / numberOfSecondarySymbols;
                                if(probState[prevRight] == 0) continue;
                                ulint doubleIndexRight = prevRight * numberOfStates + nextRight;
                                
                                if(emissionType[LeftPairSS[nextLeft]][RightPairSS[prevRight]] == 0)
                                    continue;
                                ulint doubleIndex = nextLeft*numberOfStates + prevRight;
                                
                               // if(probState[prevLeft] != 0 && probState[nextRight] != 0)
                               // {
                                    real deltaBeta =  alpha[t-1][t][doubleIndexLeft] *
                                                        alpha[v][v+1][doubleIndexRight] *
                                                        beta[t-1][v+1][indexOfParentDoubleIndex] /
                                                        probState[prevLeft] /
                                                        probState[nextRight]; // / probsum;
                                    beta[t][v][betaIndex.at(doubleIndex)] += deltaBeta;
                               // }
                            }
                        }
                    }
                    // this line prevents single probabilities to beadded more than once
                    // if(emission != 3) continue;
                }
            }
          
            normalize(beta[t][v]);

        }
    }
    
    //2.3 Termination
    
    for(int t=0; t < windowedLength-1; t++){
    
        for(ulint trans=0; trans < numberOfFragments; trans++){
            if(probFragment[trans] == 0)
                continue;
            ulint prev = trans / numberOfSecondarySymbols;
            ulint next = trans % numberOfStates;
      
            ulint doubleIndex0 = prev * numberOfStates + next;
     //       ulint indexOfDoubleIndex0 = betaINDEX.at(doubleindex0);
            
            if(emissionType[LeftPairSS[prev]][RightPairSS[next]])
                betaTTPlusOne[t][trans] = beta[t][t+1][betaIndex.at(doubleIndex0)];
   
            //
            // possible production with other possibly long fragments at left of t
            // current short fragment emits at right
            //
        
            for(ulint prevLeft = 0; prevLeft < numberOfStates; prevLeft++){
                if(probState[prevLeft] == 0) continue;
                
                real emission = emissionType[LeftPairSS[prevLeft]][RightPairSS[next]];
                if (emission == 0) continue;
                    
                ulint doubleIndex2 = prevLeft*numberOfStates + next;
                ulint indexOfDoubleIndex2 = betaIndex.at(doubleIndex2);
        
                if(emission == 1) {
                    ulint doubleIndex1 = prevLeft*numberOfStates + prev;
                    
                    real deltaBeta = 0.0;
                    for(int w = 0; w < t; w++){
                        deltaBeta += beta[w][t+1][indexOfDoubleIndex2] *
                                    alpha[w][t][doubleIndex1];
                    }
                        
                    //  beta[t][v][indexOfDoubleindex0] += deltabeta/ probState[prevLeft];
                    betaTTPlusOne[t][trans] += deltaBeta / probState[prevLeft];
                }
                
                if(emission == 3){
            
                    for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                        ulint transLeft = prevLeft * numberOfSecondarySymbols + ss;
                        ulint nextLeft = transLeft % numberOfStates;
                        if(probState[nextLeft] == 0) continue;
                        
                        ulint doubleIndex1 = nextLeft * numberOfStates + prev;
                        ulint doubleIndex3 = prevLeft * numberOfStates + nextLeft;
                    
                        real deltaBeta = 0.0;
                        for(int w = 0; w < t-1; w++){
                            deltaBeta += beta[w][t+1][indexOfDoubleIndex2]
                                        * alpha[w][w+1][doubleIndex3]
                                        * alpha[w+1][t][doubleIndex1];
                        }
                        betaTTPlusOne[t][trans] += deltaBeta
                                                    / probState[prevLeft]
                                                    / probState[nextLeft];
                    }
                }
            }
            
            //
            // possible production with other possibly long fragments at right of v
            // current short fragment emits at left
            //
            
            for(ulint nextRight = 0; nextRight < numberOfStates; nextRight++){
                if(probState[nextRight] == 0) continue;
                
                real emission = emissionType[LeftPairSS[prev]][RightPairSS[nextRight]];
                if (emission == 0) continue;
                
                ulint doubleIndex2 = prev*numberOfStates + nextRight;
                ulint indexOfDoubleIndex2 = betaIndex.at(doubleIndex2);
                
                if(emission == 2) {
                    
                    ulint doubleIndex1 = next*numberOfStates + nextRight;
            
                    real deltaBeta = 0.0;
                    for(int w = windowedLength - 1; w > t+1; w--){
                        deltaBeta += beta[t][w][indexOfDoubleIndex2]
                                    * alpha[t+1][w][doubleIndex1];
                    }
                    betaTTPlusOne[t][trans] += deltaBeta / probState[nextRight];
                }
                    
                if(emission == 3){

                    for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                        ulint transRight = nextRight+(numberOfStates * ss);
                        ulint prevRight = transRight / numberOfSecondarySymbols;
                        if (probState[prevRight] == 0) continue;
                        
                        ulint doubleIndex1 = next * numberOfStates + prevRight;
                        ulint doubleIndex3 = prevRight * numberOfStates + nextRight;
                    
                        real deltaBeta = 0.0;
                        for(int w = windowedLength - 1; w > t+2; w--){
                            deltaBeta += beta[t][w][indexOfDoubleIndex2]
                                        * alpha[w-1][w][doubleIndex3]
                                        * alpha[t+1][w-1][doubleIndex1];
                        }
                        betaTTPlusOne[t][trans] += deltaBeta
                                                    / probState[nextRight]
                                                    / probState[prevRight];
                    }
                }
            }
        }
        normalize(betaTTPlusOne[t]);
    }
 
  
    //''''''''''''''''''''''''''''''''''''''''''''''
    // 3. gamma
    //''''''''''''''''''''''''''''''''''''''''''''''
    

    // gamma
    gamma.clear();
    gamma.resize(windowedLength, vector<real>(numberOfStates,0.0));
    
    for(uint t=0; t < windowedLength-1; t++){
        
        for(ulint i=0; i<numberOfStates; i++){
            
            if(probState[i] == 0)
                continue;
            real deltagamma = 0.0;
            for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                
                ulint trans = i*numberOfSecondarySymbols + ss;
                ulint next = trans % numberOfStates;
                if(probState[next] == 0)
                    continue;
       //         if(emissionType[SScpair[i]][SScpair[next]] == 0)
       //             continue;
                ulint doubleindex = i*numberOfStates + next;
                //ulint indexOfDoubleindex = betaINDEX.at(doubleindex);
        
                deltagamma += betaTTPlusOne[t][trans]
                                    * alpha[t][t+1][doubleindex]
                                    / probState[next];
                if(t == windowedLength - 2)
                    gamma[windowedLength - 1][next] += deltagamma / probState[i];
            }
            gamma[t][i] += deltagamma / probState[i];
        }
        normalize(gamma[t]);
    }
    normalize(gamma[windowedLength - 1]);
            
    
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
