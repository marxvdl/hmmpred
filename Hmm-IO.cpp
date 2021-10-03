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


void Hmm::insideOutsideV1(vector<byte>& seq){
       
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



void Hmm::insideOutsideV2(vector<byte>& seq){
       
    uint fullLength = (uint) seq.size();
    uint windowedLength = fullLength - windowSize + 2;
        
    // SScpair is intended to be the central pair of secondary structure symbols
    // at the two central positions of a State fragment.
    // In terms of the indexing used in this code they should
    // correspond to positions (halfWindow - 1) and halfWindow.
    // For left states, therefore, the residue identity of the first position has already been emitted
    // while for the second position (halfWindow) it is waiting to be emitted.
    // For right states the second position (halfWindow) has already been emitted while the first position
    // is waiting to be emitted.
    // Secondary structure symbols (0 for E, 1 for H, 2 for C)
    // are obtained from SecondarySymbol%3,
    // and therefore SecondarySymbols must be constructed accordingly (AFPA)
    
    uint SScPair[numberOfStates];
    
    
     ulint number1 = numberOfXtraHalfStates * numberOfSecondarySymbols;
     ulint number2 = numberOfXtraHalfStates;
     ulint number3 = numberOfXtraHalfStates / numberOfSecondarySymbols;
    

/*
    // no need to use the variable numberOfXtraHalfStates, which is only defined when -h is used in the command line

    ulint number3 = 1.0;
    for(uint i = 1; i < windowSize - halfWindow - 1; i++)
        number3 *= numberOfSecondarySymbols;
    ulint number2 = number3 * numberOfSecondarySymbols;
    ulint number1 = number2 * numberOfSecondarySymbols;
   */
    for(uint i = 0; i < numberOfStates; i++){
        uint SScPair1 = i % number1 / number2 % 3;
        uint SScPair2 = i % number2 / number3 % 3;
        SScPair[i] = 3 * SScPair1 + SScPair2;
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
    
    // exclusively double emission for XE-W-EX
    
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
    
    
    // No H anywhere !!! Double E and C and otherwise first LEFT and then RIGHT
    
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
    
    // No H anywhere !!! Double E and C EXCLUSEVELY !!!!
    
    /*
   byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 3, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 3, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 3}};
      */
    
    // No H anywhere !!! C from RIGHT and E from LEFT EXCLUSEVELY !!!!
    
    /*
   byte emissionType[9][9] ={{0, 0, 0, 0, 0, 0, 0, 0, 2},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 2},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 2},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 2}};
      */
    
     /*
    byte emissionType[9][9] ={{0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {1, 1, 1, 0, 0, 0, 1, 0, 1}};
       */
    
   
    // No H anywhere !!! Double E exclusively C first LEFT and then RIGHT !!!!
    
     /*
   byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 2, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 3, 0, 0, 0, 1, 0, 1},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 2, 0, 0}};
      */
    
  
    // No H anywhere !!! Double E exclusively C first RIGHT and then LEFT !!!!
    
    /*
   byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 2, 0, 0, 0, 1, 0, 1},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 3, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 2, 0, 0, 0, 0, 0, 0}};
      */
    
    // No H anywhere !!! Double E exclusively C LEFT exclusively !!!!
    
     /*
   byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 2, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 3, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 2, 0, 0}};
      */
    
  
    // No H anywhere !!! Double E exclusively C exclusively RIGHT !!!!
    
    /*
   byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 1, 0, 1},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 3, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0}};
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
     byte emissionType[9][9] ={{0, 0, 0, 0, 0, 0, 0, 0, 2},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 2},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 2},
                               {0, 0, 0, 0, 0, 0, 0, 0, 2},
                               {0, 0, 0, 0, 0, 0, 0, 0, 2},
                               {0, 0, 0, 0, 0, 0, 0, 0, 2},
                               {0, 0, 0, 0, 0, 0, 0, 0, 2}};
      */
    

  
    // Alternating sides for C and nonC (E or H) fragments beginning at RIGHT, should also be equivalent to the HMM. NOT WORKING !!!

    /*
   byte emissionType[9][9] ={{0, 0, 0, 0, 0, 0, 2, 2, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 2, 0, 0, 2, 1, 1, 1},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 2, 2, 0},
                             {0, 0, 2, 0, 0, 2, 1, 1, 1},
                             {1, 0, 1, 0, 1, 1, 2, 2, 0},
                             {1, 0, 1, 0, 1, 1, 2, 2, 0},
                             {0, 0, 2, 0, 0, 2, 0, 0, 0}};
     */
    
    // Alternating sides for C and nonC (E or H) fragments beginning at LEFT, should also be equivalent to the HMM. NOT WORKING !!!

    /*
   byte emissionType[9][9] ={{0, 0, 2, 0, 0, 2, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {1, 0, 1, 0, 1, 1, 2, 2, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 2, 0, 0, 2, 0, 0, 0},
                             {1, 0, 1, 0, 1, 1, 2, 2, 0},
                             {0, 0, 2, 0, 0, 2, 1, 1, 1},
                             {0, 0, 2, 0, 0, 2, 1, 1, 1},
                             {0, 0, 0, 0, 0, 0, 2, 2, 0}};
      */
    
    // Alternating sides for C and nonC (only E and no H) fragments beginning at RIGHT, also equivalent to the HMM ?

    /*
   byte emissionType[9][9] ={{0, 0, 0, 0, 0, 0, 2, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 2, 0, 0, 0, 1, 0, 1},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {1, 0, 1, 0, 0, 0, 2, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 2, 0, 0, 0, 0, 0, 0}};
     */
    
    // Alternating sides for C and nonC (only E and no H) fragments beginning at LEFT, also equivalent to the HMM ?

    /*
   byte emissionType[9][9] ={{0, 0, 2, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {1, 0, 1, 0, 0, 0, 2, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 2, 0, 0, 0, 1, 0, 1},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 2, 0, 0}};
     */
    
    
      // Everybody emits at right, also equivalent to the HMM
    
    //  /*
     byte emissionType[9][9] ={{0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {1, 0, 1, 0, 1, 1, 1, 1, 1}};
    //  */
    
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
    
    
    
    //Everybody emits at left/right, again, also equivalent to the HMM
    
    /*
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
    
      /*
    //Everybody emits at left/right, again, but no HE or HE
    
    byte emissionType[9][9] ={{3, 0, 3, 0, 3, 3, 3, 3, 3},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {3, 0, 3, 0, 3, 3, 3, 3, 3},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {3, 0, 3, 0, 3, 3, 3, 3, 3},
                              {3, 0, 3, 0, 3, 3, 3, 3, 3},
                              {3, 0, 3, 0, 3, 3, 3, 3, 3},
                              {3, 0, 3, 0, 3, 3, 3, 3, 3},
                              {3, 0, 3, 0, 3, 3, 3, 3, 3}};
       */
    
   
    
    ulint numberOfDoubleStates = numberOfStates * numberOfStates;
    
    vector < vector < vector<real> > > alpha(windowedLength, vector < vector<real> > (windowedLength, vector <real> (numberOfDoubleStates, 0.0)));
    vector < vector < vector<real> > > beta(windowedLength, vector < vector<real> > (windowedLength, vector <real> (numberOfDoubleStates, 0.0)));
    vector < vector < vector<real> > > probTrans(windowedLength, vector < vector<real> > (windowedLength, vector <real> (numberOfFragments, 0.0)));
    
    vector < real > probDoubleState(numberOfDoubleStates, 0.0);
    for(ulint left = 0; left < numberOfStates; left ++){
        if(probState[left] == 0)
            continue;
        for(ulint right = 0; right < numberOfStates; right ++){
            if(probState[right] == 0)
                continue;
            byte emission = emissionType[SScPair[left]][SScPair[right]];
            if(emission == 0)
                continue;
            ulint doubleIndex = left * numberOfStates + right;
            probDoubleState[doubleIndex] = probState[left] * probState[right];
        }
    }
//    normalize(probDoubleState);
    
    
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
   /*
    for(int t=0; t < windowedLength-1; t++){
        for(ulint i=0; i<numberOfFragments; i++){
            ulint doubleIndex = i / numberOfSecondarySymbols * numberOfStates + i % numberOfStates;
            alpha[t][t+1][doubleIndex] = probFragment[i] *
                probFragmentEmitsPrimarySymbol[i][ seq[t+halfWindow] ];
        }
        normalize(alpha[t][t+1]);
    }
    */
    
    for(int t=0; t < windowedLength-1; t++){
        for(ulint trans=0; trans<numberOfFragments; trans++){
            
            ulint left = trans / numberOfSecondarySymbols;
            ulint right = trans % numberOfStates;
            if(probState[left] == 0 || probState[right] == 0)
                continue;
            ulint doubleIndex = left * numberOfStates + right;
     //       ulint doubleIndex = trans / numberOfSecondarySymbols * numberOfStates + trans % numberOfStates;
            probTrans[t][t+1][trans] = probFragment[trans] *
                probFragmentEmitsPrimarySymbol[trans][ seq[t+halfWindow] ];
           
            alpha[t][t+1][doubleIndex] = probFragmentEmitsPrimarySymbol[trans][ seq[t+halfWindow] ];
                                            // * probFragment[trans];
                                            // / probState[left] / probState[right];
        }
    //    normalize(probTrans[t][t+1]);
    //    normalize(alpha[t][t+1]);
    }
    
    //
    // 1.2 Induction
    //
    for(int v=2; v < windowedLength; v++){
        
        for (int t=v-2; t >= 0; t--) {
        
            for(ulint prevLeft=0; prevLeft<numberOfStates; prevLeft++){
                if(probState[prevLeft] == 0) continue;
                real probPrevLeft = (t > 0 ? probState[prevLeft] : probState_head[prevLeft]);
                if(probPrevLeft == 0) continue;
                
                for(ulint nextRight=0; nextRight<numberOfStates; nextRight++){
                    if(probState[nextRight] == 0) continue;
                    real probNextRight = (v < windowedLength - 1 ? probState[nextRight] : probState_tail[nextRight]);
                    if(probNextRight == 0) continue;
                    
                    byte emission = emissionType[SScPair[prevLeft]][SScPair[nextRight]];
                    if(emission == 0) continue;
                    
                    ulint doubleIndex = prevLeft * numberOfStates + nextRight;
                    
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
                            
                            /*
                            if (t == v-2 &&
                                ( nextLeft != prevRight
                                  || alpha[t+1][t+2][nextLeft*numberOfStates+nextRight]
                                  || alpha[t][t+1][prevLeft*numberOfStates+prevRight] ) )
                                    continue;
                            */
                            
                            if ( emission == 1 ) { // && ss == 0){
                                
                                ulint childDoubleIndex = prevLeft*numberOfStates + prevRight;
                            
                                if(probFragment[transRight]){
                                    real deltaAlpha = alpha[t][v-1][childDoubleIndex]
                                                    * alpha[v-1][v][doubleIndexRight]
                                                 // * probTrans[v-1][v][transRight]
                                                    // * probState[nextRight]
                                                    // * probState[prevLeft]
                                                    // / probFragment[transRight]
                                                    // ` * probState[nextRight]
                                                    // / probState[prevRight]
                                
                                                  //  / probPrevLeft;
                                                  * probFragment[transRight] / probNextRight;
                                                
                                    alpha[t][v][doubleIndex] += deltaAlpha; // * probState[nextRight]; // * (t == 2 ? probState[nextRight] : 1);
                                }
                            }
        
                            else if ( emission == 2 ) { // } && tt == 0){
                                ulint childDoubleIndex = nextLeft*numberOfStates + nextRight;
                                
                                if(probFragment[transLeft]){
                                    real deltaAlpha = alpha[t+1][v][childDoubleIndex]
                                                        * alpha[t][t+1][doubleIndexLeft]
                                                        // / probTrans[t][t+1][transLeft];
                                                        // * probState[prevLeft]
                                                        // * probState[nextRight]
                                                        // * probFragment[prevRight]
                                                        // * probState[prevLeft]
                                                        // / probState[nextLeft]
                                                        // / probPrevLeft;
                                                       // / probNextRight;
                                                       * probFragment[transLeft] / probPrevLeft;
                                                     
                                alpha[t][v][doubleIndex] += deltaAlpha; //  * (t == 2 && prevLeft > 0 ?
                                                                        // probState[prevLeft] : 1);
                               }
                                /*
                                else {
                                    real deltaAlpha = alpha[t+1][v][childDoubleIndex]
                                                        * alpha[t][t+1][doubleIndexLeft];
                                                    // * probTrans[t-1][t][transLeft]
                                                    // / probState[nextLeft];
                                                    // / probState[prevLeft];
                                    alpha[t][v][doubleIndex] += deltaAlpha;
                                }
                                 */
                            }
                            
                            else if (emission == 3) {
                                ulint childDoubleIndex = nextLeft*numberOfStates + prevRight;
                                
                                real deltaAlpha = 0.0;
                                
                                if(t == v - 2 && nextLeft == prevRight) {
                                    deltaAlpha =  probTrans[t+1][v][transRight]
                                                    * probTrans[t][t+1][transLeft]
                                                    / probState[nextLeft]; // / probState[nextRight];
                                }
                                if( (v - t) > 2){
                                    //ulint childDoubleIndex = nextLeft*numberOfStates + prevRight;
                                    deltaAlpha = alpha[t+1][v-1][childDoubleIndex]
                                                    * probTrans[t][t+1][transLeft]
                                                    / probState[nextLeft]
                                                    * probTrans[v-1][v][transRight]
                                                    / probState[prevRight];
                                }
                                
                                 /*
                                else if(t == 0 && v == windowedLength - 1) {
                                
                                    deltaAlpha =  (alpha[t][v-1][childDoubleIndex]
                                                    * probTrans[v-1][v][transRight]
                                                    / probState[prevRight])
                                                    +
                                                    (alpha[t+1][v][childDoubleIndex]
                                                    * probTrans[t][t+1][transLeft]
                                                     / probState[nextLeft]); //; // / probState[nextRight];
                                   
                                }
                                 */
                                alpha[t][v][doubleIndex] += deltaAlpha;

                            }
                        }
                    }
                
                }
            }
        
//             normalize(alpha[t][v]);
            
        }
    }
    
    
    //
    //2. Outside
    //
    
    
    //
    // 2.1 Initialization
    //
    // Note that only pair of states that could actually be generated by the grammar
    // will have non zero probability at the beginning and end of the chain.
    //
    
    if(headAndTail){
        for(ulint i=0; i < numberOfStates; i++){
            if(probState[i] == 0) // || SScPair[i] == 8)
                continue;
            for(ulint j=0; j < numberOfStates; j++){
                if(probState[j] == 0) // || SScPair[j] == 8)
                    continue;
                if(emissionType[SScPair[i]][SScPair[j]]){
                        beta[0][windowedLength-1][i*numberOfStates + j] =
                            //probState[i] * probState[j];
                            probState_head[i] * probState_tail[j];
                }
            }
        }
    }
    else{
        for(ulint i=0; i < numberOfStates; i++)
         
            for(ulint j=0; j < numberOfStates; j++)
                if(emissionType[SScPair[i]][SScPair[j]] > 0)
                    beta[0][windowedLength-1][i*numberOfStates + j] = 1.0;
    }
    
    
    normalize(beta[0][windowedLength-1]);
    
    // 2.2 Induction
    
    for(int t=0; t < windowedLength-1; t++){
        
        for (int v=windowedLength-1; v > t; v--){
 
            if(t==0 && v==windowedLength-1)
                continue;
            
            for(ulint prevLeft=0; prevLeft<numberOfStates; prevLeft++){
                // if(probState[prevLeft] == 0) continue;
                real probPrevLeft = (t == 0 ? probState_head[prevLeft] : probState[prevLeft]);
                real probPrevLeftTminus1 = (t == 1 ? probState_head[prevLeft] : probState[prevLeft]);
                // if(probPrevLeft == 0 || probPrevLeftTminus1 == 0) continue;
                
                for(ulint nextRight=0; nextRight<numberOfStates; nextRight++){
                    // if(probState[nextRight] == 0) continue;
                    real probNextRight = (v == windowedLength - 1 ? probState_tail[nextRight] : probState[nextRight]);
                    real probNextRightVplus1 = (v == windowedLength - 2 ? probState_tail[nextRight] : probState[nextRight]);
                    // if(probNextRight == 0 || probNextRightTplus1 == 0) continue;
                    
                    real probPrevLeftNextRight = probPrevLeft * probNextRight;
                    
                    real emission = emissionType[SScPair[prevLeft]][SScPair[nextRight]];
                    if(emission == 0) continue;
                    
                    ulint parentDoubleIndex = prevLeft * numberOfStates + nextRight;
                    
                    for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                        ulint transLeft = prevLeft * numberOfSecondarySymbols + ss;
                        // if(probFragment[transLeft] == 0) continue;
                        // if(probTrans[t-1][t][transLeft] == 0) continue;
                        ulint nextLeft = transLeft % numberOfStates;
                        if(probState[nextLeft] == 0) continue;
                        ulint doubleIndexLeft = prevLeft * numberOfStates + nextLeft;
                        
                         /*
                        if(0 < t && t==2 && v == t+1){
                           
                            beta[t-1][t][doubleIndexLeft] *= probState[nextLeft];
                        }
                         */
                    
                        for(uint tt=0; tt<numberOfSecondarySymbols; tt++){
                            ulint transRight = tt * numberOfStates + nextRight;
                            // if(probFragment[transRight] == 0) continue;
                            // if(probTrans[v][v+1][transRight] == 0) continue;
                            ulint prevRight = transRight /  numberOfSecondarySymbols;
                            if(probState[prevRight] == 0 && emission != 2) continue;
                            ulint doubleIndexRight = prevRight * numberOfStates + nextRight;
                            
                            if(v == t+1 && nextLeft != prevRight)
                                continue;
                            
                            /*
                            if(0 < t && t==2 && v == t+1 && nextLeft != 2 && ss <= 5 && (tt <= 5)){
                               
                                beta[t-1][t][doubleIndexLeft] *= probState[nextLeft];
                            }
                            */
                            
                             /*
                            if(emission == 1 && v == windowedLength - 1 && probNextRight){
                                
                                // if(probTrans[v-1][v][transRight] == 0) continue;
                                ulint doubleIndex = prevLeft*numberOfStates + prevRight;
                            
                                beta[t][v-1][doubleIndex] += beta[t][v][parentDoubleIndex]; // / probNextRight;
                                
                                beta[v-1][v][doubleIndexRight] += beta[t][v][parentDoubleIndex]
                                                                * alpha[t][v-1][doubleIndex]
                                                                / probPrevLeft;
                                                                // * probPrevLeftNextRight;
                                                                // / probState[prevLeft];
                            }
                             */
                            
                    
                            if(emission == 1 && v < windowedLength-1) { //} && v > t && probPrevLeft && probNextRightVplus1){
                                              //  && (v > t + 1 || nextLeft == prevRight)){ // && beta[t][v+1][parentDoubleIndex]){
                                
                                if(probTrans[v][v+1][transRight] == 0) continue;
                             
                                ulint doubleIndex = prevLeft*numberOfStates + prevRight;
                                //if(probDoubleState[doubleIndex] == 0)
                                //    continue;
                                
                               /*
                               if(0 < t && t==2 && v == t+1){
                                  // && emissionType[SScPair[prevLeft]][SScPair[nextRight]] == 2
                                  //  && emissionType[SScPair[nextLeft]][SScPair[nextRight]] == 2){
                                  
                                   beta[t][t+1][doubleIndexLeft] *= probState[nextLeft];
                                }
                                */
                                
                                real deltaBeta = beta[t][v+1][parentDoubleIndex]
                                                    * alpha[v][v+1][doubleIndexRight]
                                                    // * probTrans[v][v+1][transRight]
                                                     // / probNextRightVplus1;
                                                    // * probPrevLeftNextRight;
                                                     * probFragment[transRight]
                                                    / probState[nextRight];
                                                      //  / probPrevLeft;
                                                    // / probState[prevLeft];
                                // deltaBeta /= (v == windowedLength - 2) ? probState_tail[nextRight] : probState[nextRight];
                                // if(v > t+1)
                                
                                /*
                                if(t == 1 && v == t+1){
                                    // deltaBetaTminus1 *= probState[nextLeft];
                                 deltaBeta *= probState[nextLeft];
                                }
                                */
                                // if(t != 1 || (t == 1 && v != 2))
                                
                                // if(v > t+1 || alpha[t][v][doubleIndex])
                                beta[t][v][doubleIndex] += deltaBeta;  // * (t <= 2 && v <= 71
                                                                          //&& emissionType[SScPair[prevLe ft]][SScPair[prevRight]] == 2
                                                                          //   ? probFragment[transLeft] : 1);
                         
                                //this line corresponds to the termination step that was previously in separate loop
                                if(v > 0)
                                    beta[v][v+1][doubleIndexRight] += beta[t][v+1][parentDoubleIndex]
                                                                        * alpha[t][v][doubleIndex]
                                                                     // / probPrevLeft;
                                                                    // / probState[prevLeft];
                                                                    // / probState[prevRight]
                                                                     * probFragment[transRight] / probState[nextRight];
                                                                 //  * probFragment[transRight]; // / probNextRightVplus1;
                                                                // * probPrevLeftNextRight;
                                                                 // * (v == 1 ? probState[nextRight] : 1);
                                // if(t == 1 && v == t+1){
                                    // deltaBetaTminus1 *= probState[nextLeft];
                                   // beta[t][t+1][doubleIndexLeft] *= probState[nextLeft];
                                // }
                            
                            
                            }
                            
                             /*
                            if(0 < t && t==1 && v == t+1){
                               // && emissionType[SScPair[prevLeft]][SScPair[nextRight]] == 2
                               //  && emissionType[SScPair[nextLeft]][SScPair[nextRight]] == 2){
                               
                                beta[t][t+1][doubleIndexLeft] *= probState[nextLeft];
                            }
                             */
                            
                             /*
                            else if(emission == 2 && t == 0 && probPrevLeft && v > t+1){
                                
                                
                                // if(probTrans[t][t+1][transLeft] == 0) continue;
                                ulint doubleIndex = nextLeft*numberOfStates + nextRight;
                                
                                beta[t+1][v][doubleIndex] += beta[t][v][parentDoubleIndex];
                                                                // / probPrevLeft;
                                
                                beta[t][t+1][doubleIndexLeft] += beta[t][v][parentDoubleIndex]
                                                                    * alpha[t+1][v][doubleIndex]
                                                                    / probNextRight;
                                                                // * probPrevLeftNextRight;
                                                                // / probState[nextRight];
                            }
                             */
                            
                            else if(emission == 2 && t > 0 && tt==0) { //&& v > t && probPrevLeftTminus1 && probNextRight){
                                                  //  && (v > t + 1 || nextLeft == prevRight)){ // && beta[t-1][v][parentDoubleIndex]){
                    
                                // if(probTrans[t-1][t][transLeft] == 0) continue;
                                  
                                ulint doubleIndex = nextLeft*numberOfStates + nextRight;
                                // if(probDoubleState[doubleIndex] == 0)
                                //    continue;
                                
                                /*
                                if(t <=3 && v == t+1){
                                    // deltaBetaTminus1 *= probState[nextLeft];
                                    beta[t-1][t][doubleIndexLeft] *= probState[nextLeft];
                                }
                                */

                                real deltaBeta = beta[t-1][v][parentDoubleIndex]
                                                    * alpha[t-1][t][doubleIndexLeft]
                                                     // * probTrans[t-1][t][transLeft]
                                                     // /  probPrevLeftTminus1
                                                    // * probPrevLeftNextRight;
                                                    * probFragment[transLeft]
                                                    / probState[prevLeft];
                                                    // / probPrevLeft
                                                    // / probNextRight;
                                                    // / probState[nextRight];
                               
                                // deltaBeta /= (t == 1) ? probState_head[prevLeft] : probState[prevLeft];
                                // if(v > t+1)
                                    
                                // if(t == 1 && v == t+1 && probFragment[transRight])
                                //    deltaBeta *= probFragment[transRight];
                                
                                // if(v > t) // +1 || alpha[t][v][doubleIndex])
                                    beta[t][v][doubleIndex] += deltaBeta; // * (v >= 70
                                                                      //  && emissionType[SScPair[nextLeft]][SScPair[nextRight]] == 1 ? probState[prevRight] : 1);
                            
                                //this line corresponds to the termination step that was previously in separate loop
                                real deltaBetaTminus1 = beta[t-1][v][parentDoubleIndex]
                                                                * alpha[t][v][doubleIndex]
                                                                * probFragment[transLeft] / probState[prevLeft];
                                                                // / probNextRight;
                                                                 // / probState[nextRight];
                                // if(t != 1)
                                beta[t-1][t][doubleIndexLeft] += deltaBetaTminus1;
                                                                 // / probState[nextLeft]
                                                                 //    * probFragment[transLeft] / probState[prevLeft];
                                                                // * ( t >0 && t <= 3 && (prevLeft == 12 || prevLeft == 15 ) ? probState[nextLeft] : 1);
                                
                                 /*
                                if(t == 2 && v == t+1){
                                    // deltaBetaTminus1 *= probState[nextLeft];
                                    beta[t-1][t][doubleIndexLeft] *= probState[nextLeft];
                                }
                                 */
                                
                                // beta[t-1][t][doubleIndexLeft] += deltaBetaTminus1;
                              
                                
                            }
                            
                            /*
                           if(0 < t && t==2 && v == t+1){
                              
                               beta[t-1][t][doubleIndexLeft] *= probState[nextLeft];
                           }
                            */
                            
                           
                    
                            else if(emission == 3) {
                                //   && t > 0 && v < windowedLength -1)
                            
                                
                                // if(emissionType[SScPair[nextLeft]][SScPair[prevRight]] == 0)
                                //      continue;
                                
                                ulint doubleIndex = nextLeft*numberOfStates + prevRight;
                                ulint doubleIndex1 = nextLeft*numberOfStates + nextRight;
                                ulint doubleIndex2 = prevLeft*numberOfStates + prevRight;
                                
                                // real deltaBeta = 0.0;
                                if(t > 0 && v < windowedLength - 1 // && v > t + 1
                                   && beta[t-1][v+1][parentDoubleIndex]){
                                    real deltaBeta  =   beta[t-1][v+1][parentDoubleIndex]
                                                        * probTrans[t-1][t][transLeft]
                                                        * probTrans[v][v+1][transRight]
                                                        / probState[prevLeft]
                                                        / probState[nextRight];
                                    beta[t][v][doubleIndex] += deltaBeta;
                                    
                                    //these lines corresponds to the termination step
                                    
                                    /*
                                    beta[t-1][t][doubleIndexLeft] += beta[t-1][v+1][parentDoubleIndex]
                                                                        * alpha[t-1][v+1][parentDoubleIndex]
                                                                        / probTrans[t-1][t][transLeft]
                                                                        * probState[nextLeft]
                                                                        / probState[prevRight]
                                                                        / probState[nextRight];
                                                                        // / probFragment[transRight];
                                    beta[v][v+1][doubleIndexRight] += beta[t-1][v+1][parentDoubleIndex]
                                                                        * alpha[t-1][v+1][parentDoubleIndex]
                                                                        / probTrans[v][v+1][transRight]
                                                                        * probState[prevRight]
                                                                        / probState[prevLeft]
                                                                        / probState[nextLeft];
                                                                        // / probFragment[transLeft];
                                    */
                                    
                                    // /*
                                    beta[t-1][t][doubleIndexLeft] += beta[t-1][v+1][parentDoubleIndex]
                                                                        // * alpha[t][v+1][doubleIndex1]
                                                                        // * alpha[v][v+1][doubleIndexRight]
                                                                        * probTrans[v][v+1][transRight]
                                                                        * alpha[t][v][doubleIndex]
                                                                        / probState[prevRight]
                                                                        / probState[nextRight];
                                                                        // * probFragment[transLeft]
                                                                        // / probState[prevLeft];
                                    beta[v][v+1][doubleIndexRight] += beta[t-1][v+1][parentDoubleIndex]
                                                                        // * alpha[t-1][v][doubleIndex2]
                                                                        * probTrans[t-1][t][transLeft]
                                                                        * alpha[t][v][doubleIndex]
                                                                        / probState[prevLeft]
                                                                        / probState[nextLeft];
                                                                        // * probFragment[transRight]
                                                                        // / probState[nextRight];
                                    // */
                                }
                                
                                /*
                                else if(t == 0 && v == windowedLength - 2){
                                    
                                    ulint parentDoubleIndex1 = nextLeft * numberOfStates + nextRight;
                                    real deltaBeta = beta[t][v+1][parentDoubleIndex1]
                                                        * probTrans[v][v+1][transRight]
                                                        / probState[nextRight];
                                    beta[t][v][doubleIndex] += deltaBeta;
                                    
                                    //this line corresponds to the termination step that was previously in separate loop
                                    beta[v][v+1][doubleIndexRight] += beta[t][v+1][parentDoubleIndex1]
                                                                        * alpha[t][v+1][parentDoubleIndex1]
                                                                        / probTrans[v][v+1][transRight]
                                                                        * probState[prevRight]
                                                                        / probState[nextLeft];
                                    
                                }
                                
                                else if(t == 1 && v == windowedLength - 1){
                                    
                                    ulint parentDoubleIndex2 = prevLeft * numberOfStates + prevRight;
                                    real deltaBeta = beta[t-1][v][parentDoubleIndex2]
                                                        * probTrans[t-1][t][transLeft]
                                                        / probState[prevLeft];
                                    beta[t][v][doubleIndex] += deltaBeta;
                                    
                                    //this line corresponds to the termination step that was previously in separate loop
                                    beta[t-1][t][doubleIndexLeft] += beta[t-1][v][parentDoubleIndex2]
                                                                        * alpha[t-1][v][parentDoubleIndex2]
                                                                        / probTrans[t-1][t][transLeft]
                                                                        * probState[nextLeft]
                                                                        / probState[prevRight];
                                    
                                }
                                 */
                                
                                /*
                                if(v == t+1 && nextLeft == prevRight){
                                    
                                    if(t > 0 && beta[t-1][v][parentDoubleIndex]){
                                        beta[t][t+1][doubleIndexRight] += beta[t-1][v][parentDoubleIndex]
                                                                        * alpha[t-1][t][doubleIndexLeft]
                                                                        / probState[prevLeft];
                                    
                                        //this line corresponds to the termination step
                                        beta[t-1][t][doubleIndexLeft] += beta[t-1][v][parentDoubleIndex]
                                                                            * alpha[t][v][doubleIndexRight]
                                                                            / probState[nextRight];
                                    }
                                    if(v < windowedLength-1 && beta[t][v+1][parentDoubleIndex]){
                                        beta[t][t+1][doubleIndexLeft] += beta[t][v+1][parentDoubleIndex]
                                                                        * alpha[v][v+1][doubleIndexRight]
                                                                        / probState[nextRight];
                                    
                                        //this line corresponds to the termination step
                                        beta[v][v+1][doubleIndexRight] += beta[t][v+1][parentDoubleIndex]
                                                                            * alpha[t][t+1][doubleIndexLeft]
                                                                            / probState[prevLeft];
                                    }
                                    
                                }
                                 */
                                
                                /*
                                if(v == t+1 && t > 0 && nextLeft == prevRight
                                                    && beta[t-1][v][parentDoubleIndex]){
                                    beta[t][v][doubleIndexRight] +=  beta[t-1][v][parentDoubleIndex]
                                                                * alpha[t-1][t][doubleIndexLeft]
                                                                / probState[prevLeft];
                                    
                                    //this line corresponds to the termination step
                                    beta[t-1][t][doubleIndexLeft] += beta[t-1][v][parentDoubleIndex]
                                                                        * alpha[t][v][doubleIndex]
                                                                        / probState[nextRight];
                                }
                                */
                                 
                         //       beta[t][v][doubleIndex] += deltaBeta;
                                
                            }
                        }
                        /*
                        if(0 < t && t==1 && v == t+1 && ss==0){
                           // && emissionType[SScPair[prevLeft]][SScPair[nextRight]] == 2
                           //  && emissionType[SScPair[nextLeft]][SScPair[nextRight]] == 2){
                           
                            beta[t][t+1][doubleIndexLeft] *= probState[nextLeft];
                        }
                        */
                    }
                }
            }
            
//            if(v > t+1)
//                normalize(beta[t][v]);

        }
    }
    
//     for(int t = 0; t < windowedLength - 1; t++)
//        normalize(beta[t][t+1]);
  
    //''''''''''''''''''''''''''''''''''''''''''''''
    // 3. gamma
    //''''''''''''''''''''''''''''''''''''''''''''''
    

    // gamma
    gamma.clear();
    gamma.resize(windowedLength, vector<real>(numberOfStates,0.0));
    
    for(uint t=0; t < windowedLength-1; t++){
        
        for(ulint i=0; i<numberOfStates; i++){
            if(probState[i] == 0) continue;
            // real probI = (t == 0 ? probState_head[i] : probState[i]);
            // if(probI == 0)
            //    continue;
            real deltagamma = 0.0;
            for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                
                ulint trans = i*numberOfSecondarySymbols + ss;
                if(probFragment[trans] == 0)
                    continue;
                ulint next = trans % numberOfStates;
                if(probState[next] == 0) continue;
                // real probNext = (t == windowedLength - 2 ? probState_tail[next] : probState[next]);
                // if(probNext == 0)
                //    continue;
                ulint doubleindex = i*numberOfStates + next;

                deltagamma += beta[t][t+1][doubleindex] * alpha[t][t+1][doubleindex];
                                  // / probFragment[trans];
                                // /  probState[next];
                                 // / probNext;
                if(t == windowedLength - 2)
                    gamma[windowedLength - 1][next] += deltagamma; // / probState[i];
                
            }
            gamma[t][i] += deltagamma; // / probState[i];
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
                //byte sec = indexToSequence(windowSize-1, i)[0];
                //probSecondarySymbol[a][sec] += gamma[a][i];
                 byte sec = indexToSequence(windowSize-1, i)[halfWindow];
                 probSecondarySymbol[a+halfWindow][sec] += gamma[a][i];
            }
            // normalize(probSecondarySymbol[a]);
            normalize(probSecondarySymbol[a+halfWindow]);
        }

        /*
        for(uint a=windowedLength; a < fullLength; a++){
            for(ulint i=0; i < numberOfStates; i++){
                byte sec = indexToSequence(windowSize-1, i)[a - windowedLength + 1];
                probSecondarySymbol[a][sec] += gamma[windowedLength-1][i];
            }
            normalize(probSecondarySymbol[a]);
        }
         */
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
        
    // SScpair is intended to be the central pair of secondary structure symbols
    // at the two central positions of a State fragment.
    // In terms of the indexing used in this code they should
    // correspond to positions (halfWindow - 1) and halfWindow.
    // For left states, therefore, the residue identity of the first position has already been emitted
    // while for the second position (halfWindow) it is waiting to be emitted.
    // For right states the second position (halfWindow) has already been emitted while the first position
    // is waiting to be emitted.
    // Secondary structure symbols (0 for E, 1 for H, 2 for C)
    // are obtained from SecondarySymbol%3,
    // and therefore SecondarySymbols must be constructed accordingly (AFPA)
    
    uint SScPair[numberOfStates];
    
    
     ulint number1 = numberOfXtraHalfStates * numberOfSecondarySymbols;
     ulint number2 = numberOfXtraHalfStates;
     ulint number3 = numberOfXtraHalfStates / numberOfSecondarySymbols;
    

/*
    // no need to use the variable numberOfXtraHalfStates, which is only defined when -h is used in the command line

    ulint number3 = 1.0;
    for(uint i = 1; i < windowSize - halfWindow - 1; i++)
        number3 *= numberOfSecondarySymbols;
    ulint number2 = number3 * numberOfSecondarySymbols;
    ulint number1 = number2 * numberOfSecondarySymbols;
   */
    for(uint i = 0; i < numberOfStates; i++){
        uint SScPair1 = i % number1 / number2 % 3;
        uint SScPair2 = i % number2 / number3 % 3;
        SScPair[i] = 3 * SScPair1 + SScPair2;
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
    
    // exclusively double emission for XE-W-EX
    
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
    
    
    // No H anywhere !!! Double E and C and otherwise first LEFT and then RIGHT
    
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
    
    // No H anywhere !!! Double E and C EXCLUSEVELY !!!!
    
    /*
   byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 3, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 3, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 3}};
      */
    
    // No H anywhere !!! C from RIGHT and E from LEFT EXCLUSEVELY !!!!
    
    /*
   byte emissionType[9][9] ={{0, 0, 0, 0, 0, 0, 0, 0, 2},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 2},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 2},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 2}};
      */
    
     /*
    byte emissionType[9][9] ={{0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {1, 1, 1, 0, 0, 0, 1, 0, 1}};
       */
    
   
    // No H anywhere !!! Double E exclusively C first LEFT and then RIGHT !!!!
    
     /*
   byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 2, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 3, 0, 0, 0, 1, 0, 1},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 2, 0, 0}};
      */
    
  
    // No H anywhere !!! Double E exclusively C first RIGHT and then LEFT !!!!
    
    /*
   byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 2, 0, 0, 0, 1, 0, 1},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 3, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 2, 0, 0, 0, 0, 0, 0}};
      */
    
    // No H anywhere !!! Double E exclusively C LEFT exclusively !!!!
    
     /*
   byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 2, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 3, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 2, 0, 0}};
      */
    
  
    // No H anywhere !!! Double E exclusively C exclusively RIGHT !!!!
    
    /*
   byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 1, 0, 1},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 3, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0}};
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
    
     // /*
     byte emissionType[9][9] ={{0, 0, 0, 0, 0, 0, 0, 0, 2},
                               {0, 0, 0, 0, 0, 0, 0, 0, 2},
                               {0, 0, 0, 0, 0, 0, 0, 0, 2},
                               {0, 0, 0, 0, 0, 0, 0, 0, 2},
                               {0, 0, 0, 0, 0, 0, 0, 0, 2},
                               {0, 0, 0, 0, 0, 0, 0, 0, 2},
                               {0, 0, 0, 0, 0, 0, 0, 0, 2},
                               {0, 0, 0, 0, 0, 0, 0, 0, 2},
                               {0, 0, 0, 0, 0, 0, 0, 0, 2}};
     //  */
    

  
    // Alternating sides for C and nonC (E or H) fragments beginning at RIGHT, should also be equivalent to the HMM !!!

    /*
   byte emissionType[9][9] ={{0, 0, 0, 0, 0, 0, 2, 2, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 2, 0, 0, 2, 1, 1, 1},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 2, 2, 0},
                             {0, 0, 2, 0, 0, 2, 1, 1, 1},
                             {1, 0, 1, 0, 1, 1, 2, 2, 0},
                             {1, 0, 1, 0, 1, 1, 2, 2, 0},
                             {0, 0, 2, 0, 0, 2, 0, 0, 0}};
     */
    
    // Alternating sides for C and nonC (E or H) fragments beginning at LEFT, should also be equivalent to the HMM. NOT WORKING !!!

    /*
   byte emissionType[9][9] ={{0, 0, 2, 0, 0, 2, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {1, 0, 1, 0, 1, 1, 2, 2, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 2, 0, 0, 2, 0, 0, 0},
                             {1, 0, 1, 0, 1, 1, 2, 2, 0},
                             {0, 0, 2, 0, 0, 2, 1, 1, 1},
                             {0, 0, 2, 0, 0, 2, 1, 1, 1},
                             {0, 0, 0, 0, 0, 0, 2, 2, 0}};
      */
    
    // Alternating sides for C and nonC (only E and no H) fragments beginning at RIGHT, also equivalent to the HMM ?

    /*
   byte emissionType[9][9] ={{0, 0, 0, 0, 0, 0, 2, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 2, 0, 0, 0, 1, 0, 1},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {1, 0, 1, 0, 0, 0, 2, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 2, 0, 0, 0, 0, 0, 0}};
     */
    
    // Alternating sides for C and nonC (only E and no H) fragments beginning at LEFT, also equivalent to the HMM ?

    /*
   byte emissionType[9][9] ={{0, 0, 2, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {1, 0, 1, 0, 0, 0, 2, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 2, 0, 0, 0, 1, 0, 1},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 2, 0, 0}};
     */
    
    
      // Everybody emits at right, also equivalent to the HMM
    
      /*
     byte emissionType[9][9] ={{0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {1, 1, 1, 1, 1, 1, 1, 1, 1}};
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
    
    
    
    //Everybody emits at left/right, again, also equivalent to the HMM
    
     /*
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
    
      /*
    //Everybody emits at left/right, again, but no HE or HE
    
    byte emissionType[9][9] ={{3, 0, 3, 0, 3, 3, 3, 3, 3},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {3, 0, 3, 0, 3, 3, 3, 3, 3},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {3, 0, 3, 0, 3, 3, 3, 3, 3},
                              {3, 0, 3, 0, 3, 3, 3, 3, 3},
                              {3, 0, 3, 0, 3, 3, 3, 3, 3},
                              {3, 0, 3, 0, 3, 3, 3, 3, 3},
                              {3, 0, 3, 0, 3, 3, 3, 3, 3}};
       */
    
   
    
    ulint numberOfDoubleStates = numberOfStates * numberOfStates;
    
    vector < vector < vector<real> > > alpha(windowedLength, vector < vector<real> > (windowedLength, vector <real> (numberOfDoubleStates, 0.0)));
    vector < vector < vector<real> > > beta(windowedLength, vector < vector<real> > (windowedLength, vector <real> (numberOfDoubleStates, 0.0)));
    vector < vector < vector<real> > > probTrans(windowedLength, vector < vector<real> > (windowedLength, vector <real> (numberOfFragments, 0.0)));
    
    vector < real > probDoubleState(numberOfDoubleStates, 0.0);
    for(ulint left = 0; left < numberOfStates; left ++){
        if(probState[left] == 0)
            continue;
        for(ulint right = 0; right < numberOfStates; right ++){
            if(probState[right] == 0)
                continue;
            byte emission = emissionType[SScPair[left]][SScPair[right]];
            if(emission == 0)
                continue;
            ulint doubleIndex = left * numberOfStates + right;
            probDoubleState[doubleIndex] = probState[left] * probState[right];
        }
    }
//    normalize(probDoubleState);
    
    
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
   /*
    for(int t=0; t < windowedLength-1; t++){
        for(ulint i=0; i<numberOfFragments; i++){
            ulint doubleIndex = i / numberOfSecondarySymbols * numberOfStates + i % numberOfStates;
            alpha[t][t+1][doubleIndex] = probFragment[i] *
                probFragmentEmitsPrimarySymbol[i][ seq[t+halfWindow] ];
        }
        normalize(alpha[t][t+1]);
    }
    */
    
    for(int t=0; t < windowedLength-1; t++){
        for(ulint trans=0; trans<numberOfFragments; trans++){
            
            ulint left = trans / numberOfSecondarySymbols;
            ulint right = trans % numberOfStates;
            if(probState[left] == 0 || probState[right] == 0)
                continue;
            ulint doubleIndex = left * numberOfStates + right;
     //       ulint doubleIndex = trans / numberOfSecondarySymbols * numberOfStates + trans % numberOfStates;
            probTrans[t][t+1][trans] = probFragment[trans] *
                probFragmentEmitsPrimarySymbol[trans][ seq[t+halfWindow] ];
           
            alpha[t][t+1][doubleIndex] = probFragmentEmitsPrimarySymbol[trans][ seq[t+halfWindow] ];
                                            // * probFragment[trans]
                                            // / probState[left] / probState[right];
        }
    //    normalize(probTrans[t][t+1]);
    //    normalize(alpha[t][t+1]);
    }
    
    //
    // 1.2 Induction
    //
    for(int v=2; v < windowedLength; v++){
        
        for (int t=v-2; t >= 0; t--) {
        
            for(ulint prevLeft=0; prevLeft<numberOfStates; prevLeft++){
                if(probState[prevLeft] == 0) continue;
                real probPrevLeft = (t > 0 ? probState[prevLeft] : probState_head[prevLeft]);
                if(probPrevLeft == 0) continue;
                
                for(ulint nextRight=0; nextRight<numberOfStates; nextRight++){
                    if(probState[nextRight] == 0) continue;
                    real probNextRight = (v < windowedLength - 1 ? probState[nextRight] : probState_tail[nextRight]);
                    if(probNextRight == 0) continue;
                    
                    byte emission = emissionType[SScPair[prevLeft]][SScPair[nextRight]];
                    if(emission == 0) continue;
                    
                    ulint doubleIndex = prevLeft * numberOfStates + nextRight;
                    
                    if ( emission == 1 ) {
                        
                        for(uint tt=0; tt<numberOfSecondarySymbols; tt++){
                                
                            ulint transRight = tt * numberOfStates + nextRight;
                            ulint prevRight = transRight / numberOfSecondarySymbols;
                            if(probState[prevRight] == 0) continue;
                            ulint doubleIndexRight = prevRight * numberOfStates + nextRight;
                            
                            ulint childDoubleIndex = prevLeft*numberOfStates + prevRight;
                    
                            if(probFragment[transRight]){
                                real deltaAlpha = alpha[t][v-1][childDoubleIndex]
                                                        * alpha[v-1][v][doubleIndexRight]
                                                        * probFragment[transRight] / probState[nextRight]; // / probPrevLeft;
                                
                                if(v > t + 2) {
                                    alpha[t][v][doubleIndex] += deltaAlpha;
                                }
                                else {                             // i.e. if (v == t + 2)
                                    for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                                                    
                                        ulint transLeft =  prevLeft * numberOfSecondarySymbols + ss;
                                        ulint nextLeft = transLeft % numberOfStates;
                                        if(prevRight == nextLeft) {
                                            alpha[t][v][doubleIndex] += deltaAlpha
                                                                        * probFragment[transLeft] / probState[nextLeft];
                                        }
                                    }
                                }
                                /*
                            if(probFragment[transRight]){
                                real deltaAlpha = alpha[t][v-1][childDoubleIndex]
                                                    * alpha[v-1][v][doubleIndexRight]
                                * probFragment[transRight] / probState[nextRight]; // / probNextRight;
                                
                                if(t == v-2) {
                                    ulint transMiddle = prevLeft*numberOfSecondarySymbols +
                                                        prevRight%numberOfSecondarySymbols;
                                    deltaAlpha *= probFragment[transMiddle] / probState[prevRight];
                                }
                                alpha[t][v][doubleIndex] += deltaAlpha;
                                 */
                            }
                        }
                    }
                            
                    else if ( emission == 2 ) {
                        
                        for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                                
                            ulint transLeft = prevLeft * numberOfSecondarySymbols + ss;
                            ulint nextLeft = transLeft % numberOfStates;
                            if(probState[nextLeft] == 0) continue;
                            ulint doubleIndexLeft = prevLeft*numberOfStates + nextLeft;
                        
                            ulint childDoubleIndex = nextLeft*numberOfStates + nextRight;
                                
                            if(probFragment[transLeft]){
                                real deltaAlpha = alpha[t+1][v][childDoubleIndex]
                                                        * alpha[t][t+1][doubleIndexLeft]
                                                        * probFragment[transLeft] / probState[prevLeft]; // / probPrevLeft;
                                
                                if(v > t + 2) {
                                    alpha[t][v][doubleIndex] += deltaAlpha;
                                }
                                else {                             // i.e. if (v == t + 2)
                                    for(uint tt=0; tt<numberOfSecondarySymbols; tt++){
                                                    
                                        ulint transRight = tt * numberOfStates + nextRight;
                                        ulint prevRight = transRight / numberOfSecondarySymbols;
                                        if(prevRight == nextLeft) {
                                            alpha[t][v][doubleIndex] += deltaAlpha
                                                                        * probFragment[transRight] / probState[prevRight];
                                        }
                                    }
                                }
                                // if(t == v-2) deltaAlpha /= probState[nextLeft];
                                //  * (t == 2 && prevLeft > 0 ?
                                                                        // probState[prevLeft] : 1);
                            }
                        }
                    }
                            
                    else if (emission == 3) {
                                
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
                                
                                ulint childDoubleIndex = nextLeft*numberOfStates + prevRight;
                                
                                real deltaAlpha = 0.0;
                                
                                // /*
                                if(v - t == 2 && nextLeft == prevRight) {
                                    deltaAlpha = alpha[t][t+1][doubleIndexLeft]
                                                    * probFragment[transLeft]
                                                    / probState[prevLeft]
                                                    
                                                    * alpha[v-1][v][doubleIndexRight]
                                                    * probFragment[transRight]
                                                    / probState[nextRight]
                                                  
                                                    / probState[prevRight]; // / probState[nextRight];
                                }
                                // */
                                
                                // /*
                                if(v - t == 3 && alpha[t+1][v-1][childDoubleIndex]) {
                            
                                    ulint transMiddle = nextLeft*numberOfSecondarySymbols +
                                                        (prevRight%numberOfSecondarySymbols);
                                    ulint transMiddle2 = nextLeft * numberOfSecondarySymbols
                                                            / numberOfStates * numberOfStates + prevRight;
                                  //  if(transMiddle == transMiddle2) { //} && probFragment[transMiddle]) {
                                        deltaAlpha = alpha[t+1][v-1][childDoubleIndex]
                                      
                                        * alpha[t][t+1][doubleIndexLeft]
                                        * probFragment[transLeft]
                                        / probState[prevLeft]
                                      
                                        * alpha[v-1][v][doubleIndexRight]
                                        * probFragment[transRight]
                                        / probState[nextRight]
                                        
                                        * probFragment[transMiddle2]
                                                    / probState[nextLeft]
                                                    / probState[prevRight];
                                    // }
                                 
                                
                                }
                                // */
                                
                                if( (v - t) > 3){
                                    //ulint childDoubleIndex = nextLeft*numberOfStates + prevRight;
                                    deltaAlpha = alpha[t+1][v-1][childDoubleIndex]
                                                  
                                                    * alpha[t][t+1][doubleIndexLeft]
                                                    * probFragment[transLeft]
                                                    / probState[prevLeft]
                                                 
                                                    * alpha[v-1][v][doubleIndexRight]
                                                    * probFragment[transRight]
                                                    / probState[nextRight];
                                }
                                
                                 /*
                                else if(t == 0 && v == windowedLength - 1) {
                                
                                    deltaAlpha =  (alpha[t][v-1][childDoubleIndex]
                                                    * probTrans[v-1][v][transRight]
                                                    / probState[prevRight])
                                                    +
                                                    (alpha[t+1][v][childDoubleIndex]
                                                    * probTrans[t][t+1][transLeft]
                                                     / probState[nextLeft]); //; // / probState[nextRight];
                                   
                                }
                                 */
                                alpha[t][v][doubleIndex] += deltaAlpha;
                            
                            }
                        }
                    }
                
                }
            }
        
//             normalize(alpha[t][v]);
            
        }
    }
    
    
    //
    //2. Outside
    //
    
    
    //
    // 2.1 Initialization
    //
    // Note that only pair of states that could actually be generated by the grammar
    // will have non zero probability at the beginning and end of the chain.
    //
    
    if(headAndTail){
        for(ulint i=0; i < numberOfStates; i++){
            if(probState[i] == 0) // || SScPair[i] == 8)
                continue;
            for(ulint j=0; j < numberOfStates; j++){
                if(probState[j] == 0) // || SScPair[j] == 8)
                    continue;
                if(emissionType[SScPair[i]][SScPair[j]]){
                        beta[0][windowedLength-1][i*numberOfStates + j] =
                            //probState[i] * probState[j];
                            probState_head[i] * probState_tail[j];
                }
            }
        }
    }
    else{
        for(ulint i=0; i < numberOfStates; i++)
         
            for(ulint j=0; j < numberOfStates; j++)
                if(emissionType[SScPair[i]][SScPair[j]] > 0)
                    beta[0][windowedLength-1][i*numberOfStates + j] = 1.0;
    }
    
    
    normalize(beta[0][windowedLength-1]);
    
    // 2.2 Induction
    
    for(int t=0; t < windowedLength-1; t++){
        
        for (int v=windowedLength-1; v >= t; v--){
 
            if(t==0 && v==windowedLength-1)
                continue;
            
            for(ulint prevLeft=0; prevLeft<numberOfStates; prevLeft++){
                 if(probState[prevLeft] == 0) continue;
                real probPrevLeft = (t == 0 ? probState_head[prevLeft] : probState[prevLeft]);
                real probPrevLeftTminus1 = (t == 1 ? probState_head[prevLeft] : probState[prevLeft]);
                // if(probPrevLeft == 0 || probPrevLeftTminus1 == 0) continue;
                
                for(ulint nextRight=0; nextRight<numberOfStates; nextRight++){
                     if(probState[nextRight] == 0) continue;
                    real probNextRight = (v == windowedLength - 1 ? probState_tail[nextRight] : probState[nextRight]);
                    real probNextRightVplus1 = (v == windowedLength - 2 ? probState_tail[nextRight] : probState[nextRight]);
                    // if(probNextRight == 0 || probNextRightTplus1 == 0) continue;
                    
                    real probPrevLeftNextRight = probPrevLeft * probNextRight;
                    
                    real emission = emissionType[SScPair[prevLeft]][SScPair[nextRight]];
                    if(emission == 0) continue;
                    
                    ulint parentDoubleIndex = prevLeft * numberOfStates + nextRight;
                
                    if(emission == 1 && v < windowedLength-1) { //} && v > t && probPrevLeft && probNextRightVplus1){
                                              //  && (v > t + 1 || nextLeft == prevRight)){ // && beta[t][v+1][parentDoubleIndex]){
                                
                    
                        for(uint tt=0; tt<numberOfSecondarySymbols; tt++){
                            ulint transRight = tt * numberOfStates + nextRight;
                            // if(probFragment[transRight] == 0) continue;
                            // if(probTrans[v][v+1][transRight] == 0) continue;
                            ulint prevRight = transRight /  numberOfSecondarySymbols;
                            if(probState[prevRight] == 0) continue;
                            ulint doubleIndexRight = prevRight * numberOfStates + nextRight;
                                
                            // if(probTrans[v][v+1][transRight] == 0) continue;
                             
                            ulint doubleIndex = prevLeft*numberOfStates + prevRight;
                                //if(probDoubleState[doubleIndex] == 0)
                                //    continue;
                                
                             
                                
                                
                            real deltaBeta = beta[t][v+1][parentDoubleIndex]
                                                // * alpha[v][v+1][doubleIndexRight]
                                                * probFragment[transRight] / probState[nextRight];
                            
                            if ( v > t + 1 ) {
                                
                                beta[t][v][doubleIndex] += deltaBeta * alpha[v][v+1][doubleIndexRight];
                                beta[v][v+1][doubleIndexRight] += deltaBeta * alpha[t][v][doubleIndex];
                            }
                            
                            else {          // i.e. if(v == t+1) {
                                
                                for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                                    ulint transLeft = prevLeft * numberOfSecondarySymbols + ss;
                                  
                                    ulint nextLeft = transLeft % numberOfStates;
                                    if(probState[nextLeft] == probState[prevRight]) {
                                        deltaBeta *= probFragment[transLeft] / probState[nextLeft];
                                        beta[t][v][doubleIndex] += deltaBeta * alpha[v][v+1][doubleIndexRight];
                                        beta[v][v+1][doubleIndexRight] += deltaBeta * alpha[t][v][doubleIndex];
                                    }
                                }
                               
                            
                            }
                                
                                /*
                                if(t == 1 && v == t+1){
                                    // deltaBetaTminus1 *= probState[nextLeft];
                                 deltaBeta *= probState[nextLeft];
                                }
                                */
                                // if(t != 1 || (t == 1 && v != 2))
                                
                                // if(v > t+1 || alpha[t][v][doubleIndex])
                            // beta[t][v][doubleIndex] += deltaBeta;  // * (t <= 2 && v <= 71
                                                                          //&& emissionType[SScPair[prevLe ft]][SScPair[prevRight]] == 2
                                                                          //   ? probFragment[transLeft] : 1);
                         
                                //this line corresponds to the termination step that was previously in separate loop
                         //   if(v > 0)
                         //       beta[v][v+1][doubleIndexRight] += beta[t][v+1][parentDoubleIndex]
                         //                                           * alpha[t][v][doubleIndex]
                         //                                           * probFragment[transRight] / probState[nextRight];
                                                               
                        }
                    }
                            
                            
                    else if(emission == 2 && t > 0 ) { // }&& tt==0) { //&& v > t && probPrevLeftTminus1 && probNextRight){
                                                  //  && (v > t + 1 || nextLeft == prevRight)){ // && beta[t-1][v][parentDoubleIndex]){
                    
                        for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                            ulint transLeft = prevLeft * numberOfSecondarySymbols + ss;
                            // if(probFragment[transLeft] == 0) continue;
                            // if(probTrans[t-1][t][transLeft] == 0) continue;
                            ulint nextLeft = transLeft % numberOfStates;
                            if(probState[nextLeft] == 0) continue;
                            ulint doubleIndexLeft = prevLeft * numberOfStates + nextLeft;
                            
                                // if(probTrans[t-1][t][transLeft] == 0) continue;
                                  
                                
                            ulint doubleIndex = nextLeft*numberOfStates + nextRight;
                                // if(probDoubleState[doubleIndex] == 0)
                                //    continue;
                                
                             

                            real deltaBeta = beta[t-1][v][parentDoubleIndex]
                                                    // * alpha[t-1][t][doubleIndexLeft]
                                                    * probFragment[transLeft] / probState[prevLeft];
                            
                            if ( v > t + 1 ) {
                                beta[t][v][doubleIndex] += deltaBeta * alpha[t-1][t][doubleIndexLeft];
                                beta[t-1][t][doubleIndexLeft] += deltaBeta * alpha[t][v][doubleIndex];
                            }
                                                   
                            else {                     // i.e. if ( v == t + 1) {
                            
                                for(uint tt=0; tt<numberOfSecondarySymbols; tt++){
                                                
                                    ulint transRight = tt * numberOfStates + nextRight;
                                    ulint prevRight = transRight / numberOfSecondarySymbols;
                                    if(prevRight == nextLeft) {
                                        deltaBeta *= probFragment[transRight] / probState[prevRight];
                                        beta[t][v][doubleIndex] += deltaBeta * alpha[t-1][t][doubleIndexLeft];
                                        beta[t-1][t][doubleIndexLeft] += deltaBeta * alpha[t][v][doubleIndex];
                                    }
                                }
                            }
                            /*
                                if(v == t+1) {
                                
                                    ulint transMiddle = nextLeft*numberOfSecondarySymbols +
                                                        nextRight%numberOfSecondarySymbols;
                                    deltaBeta *= probFragment[transMiddle]
                                                    / probState[nextLeft];
                                                    // / probState[nextRight];
                            
                                }
                                */
                                // deltaBeta /= (t == 1) ? probState_head[prevLeft] : probState[prevLeft];
                                // if(v > t+1)
                                    
                                // if(t == 1 && v == t+1 && probFragment[transRight])
                                //    deltaBeta *= probFragment[transRight];
                                
                                // if(v > t) // +1 || alpha[t][v][doubleIndex])
                                 //   beta[t][v][doubleIndex] += deltaBeta; // * (v >= 70
                                                                      //  && emissionType[SScPair[nextLeft]][SScPair[nextRight]] == 1 ? probState[prevRight] : 1);
                            
                                //this line corresponds to the termination step that was previously in separate loop
                                // real deltaBetaTminus1 = beta[t-1][v][parentDoubleIndex]
                                //                                * alpha[t][v][doubleIndex]
                                //                                * probFragment[transLeft] / probState[prevLeft];
                                                              
                                // if(t != 1)
                               // beta[t-1][t][doubleIndexLeft] += deltaBetaTminus1;
                                                                 // / probState[nextLeft]
                                                                 //    * probFragment[transLeft] / probState[prevLeft];
                                                                // * ( t >0 && t <= 3 && (prevLeft == 12 || prevLeft == 15 ) ? probState[nextLeft] : 1
                        }
                    }
                            
                            /*
                           if(0 < t && t==2 && v == t+1){
                              
                               beta[t-1][t][doubleIndexLeft] *= probState[nextLeft];
                           }
                            */
                            
                           
                    
                    else if(emission == 3  ) { //} && probState[prevLeft] && probState[nextRight] // ) {
                             // && t > 0 && v < windowedLength -1) {
                            
                                
                                // if(emissionType[SScPair[nextLeft]][SScPair[prevRight]] == 0)
                                //      continue;
                                
                        for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                            ulint transLeft = prevLeft * numberOfSecondarySymbols + ss;
                            // if(probFragment[transLeft] == 0) continue;
                            // if(probTrans[t-1][t][transLeft] == 0) continue;
                            ulint nextLeft = transLeft % numberOfStates;
                            if(probState[nextLeft] == 0) continue;
                            ulint doubleIndexLeft = prevLeft * numberOfStates + nextLeft;
                                    
                                     /*
                                    if(0 < t && t==2 && v == t+1){
                                       
                                        beta[t-1][t][doubleIndexLeft] *= probState[nextLeft];
                                    }
                                     */
                                
                            for(uint tt=0; tt<numberOfSecondarySymbols; tt++){
                                ulint transRight = tt * numberOfStates + nextRight;
                                // if(probFragment[transRight] == 0) continue;
                                // if(probTrans[v][v+1][transRight] == 0) continue;
                                ulint prevRight = transRight /  numberOfSecondarySymbols;
                                if(probState[prevRight] == 0 && emission != 2) continue;
                                ulint doubleIndexRight = prevRight * numberOfStates + nextRight;
                                        
                         //        if(v == t+1 && nextLeft != prevRight)
                         //           continue;
                                
                                ulint doubleIndex = nextLeft*numberOfStates + prevRight;
                                ulint doubleIndex1 = nextLeft*numberOfStates + nextRight;
                                ulint doubleIndex2 = prevLeft*numberOfStates + prevRight;
                                
                                
                                
                                // real deltaBeta = 0.0;
                                if(t > 0 && v < windowedLength - 1 // && v > t + 1
                                   && beta[t-1][v+1][parentDoubleIndex]){
                                    
                                    real transProb = beta[t-1][v+1][parentDoubleIndex]
                                                        * probFragment[transLeft]
                                                        * probFragment[transRight]
                                                        / probState[prevLeft]
                                                        / probState[nextRight];
                                    
                                    real deltaBeta  =   transProb // beta[t-1][v+1][parentDoubleIndex]
                                                        // * probTrans[t-1][t][transLeft]
                                                        // * probTrans[v][v+1][transRight]
                                                        * alpha[t-1][t][doubleIndexLeft]
                                                        * alpha[v][v+1][doubleIndexRight];
                                                        // * probFragment[transLeft]
                                                        // * probFragment[transRight]
                                                        // / probState[prevLeft]
                                                        // / probState[nextRight];
                              
                                    if( v == t && nextLeft == prevRight) {
                                        transProb /= probState[nextLeft];
                                        beta[t-1][t][doubleIndexLeft] += transProb
                                                                            * alpha[v][v+1][doubleIndexRight];
                                                                            // * alpha[t][v+1][doubleIndex1];
                                                                            
                                        beta[v][v+1][doubleIndexRight] += transProb
                                                                            * alpha[t-1][t][doubleIndexLeft];
                                                                            // * alpha[t-1][v][doubleIndex2];
                                                                
                                    }
                                    
                                    else if(v == t+1 && alpha[t][v][doubleIndex]) {
                                        
                                        // deltaBeta /= (probState[nextLeft] * probState[prevRight]);
                                         // /*
                                        ulint transMiddle = nextLeft*numberOfSecondarySymbols +
                                                            prevRight%numberOfSecondarySymbols;
                                        ulint transMiddle2 = nextLeft * numberOfSecondarySymbols
                                                                / numberOfStates * numberOfStates+ prevRight;
                                       // if(transMiddle == transMiddle2) {
                                            transProb *= probFragment[transMiddle]
                                                        / probState[nextLeft]
                                                        / probState[prevRight];
                                       // }
                                     
                                     beta[t][v][doubleIndex] += transProb
                                                                    * alpha[t-1][t][doubleIndexLeft]
                                                                    * alpha[v][v+1][doubleIndexRight];
                                    
                                    //these lines corresponds to the termination step
                                    
                                 
                                    // /*
                                        beta[t-1][t][doubleIndexLeft] += transProb // beta[t-1][v+1][parentDoubleIndex]
                                                                         // * alpha[t][v+1][doubleIndex1];
                                                                          * alpha[v][v+1][doubleIndexRight]
                                                                          * alpha[t][v][doubleIndex];
                                                                        // / probState[prevRight]
                                                                        // / probState[nextRight];
                                                                        // * probFragment[transRight] / probState[nextRight];
                                                                        // * probFragment[transLeft]
                                                                        // / probState[prevLeft];
                                        beta[v][v+1][doubleIndexRight] += transProb // beta[t-1][v+1][parentDoubleIndex]
                                                                        // * alpha[t-1][v][doubleIndex2];
                                                                        * alpha[t-1][t][doubleIndexLeft]
                                                                        * alpha[t][v][doubleIndex];
                                                                        // / probState[prevLeft]
                                                                        // / probState[nextLeft];
                                                                        // * probFragment[transRight]
                                                                        // / probState[nextRight];
                                    // */
                                
                                    
                                    }
                                    
                                    else if(v > t+1) {
                                        
                                        beta[t][v][doubleIndex] += deltaBeta;
                                    
                                    //these lines corresponds to the termination step
                                    
                                 
                                    // /*
                                        beta[t-1][t][doubleIndexLeft] += transProb // beta[t-1][v+1][parentDoubleIndex]
                                                                        // * alpha[t][v+1][doubleIndex1];
                                                                         * alpha[v][v+1][doubleIndexRight]
                                                                        * alpha[t][v][doubleIndex];
                                                                        // / probState[prevRight]
                                                                        // / probState[nextRight];
                                                                        // * probFragment[transRight] / probState[nextRight];
                                                                        // * probFragment[transLeft]
                                                                        // / probState[prevLeft];
                                        beta[v][v+1][doubleIndexRight] += transProb // beta[t-1][v+1][parentDoubleIndex]
                                                                        // * alpha[t-1][v][doubleIndex2]
                                                                         * alpha[t-1][t][doubleIndexLeft]
                                                                        * alpha[t][v][doubleIndex];
                                                                        // / probState[prevLeft]
                                                                        // / probState[nextLeft];
                                                                        // * probFragment[transRight]
                                                                        // / probState[nextRight];
                                    // */
                                    }
                                }
                                
                                
                            }
                        }
                     
                    }
                }
            }
            
//            if(v > t+1)
//                normalize(beta[t][v]);

        }
    }
    
//     for(int t = 0; t < windowedLength - 1; t++)
//        normalize(beta[t][t+1]);
  
    //''''''''''''''''''''''''''''''''''''''''''''''
    // 3. gamma
    //''''''''''''''''''''''''''''''''''''''''''''''
    

    // gamma
    gamma.clear();
    gamma.resize(windowedLength, vector<real>(numberOfStates,0.0));
    
    for(uint t=0; t < windowedLength-1; t++){
        
        for(ulint i=0; i<numberOfStates; i++){
            if(probState[i] == 0) continue;
            // real probI = (t == 0 ? probState_head[i] : probState[i]);
            // if(probI == 0)
            //    continue;
            real deltagamma = 0.0;
            for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                
                ulint trans = i*numberOfSecondarySymbols + ss;
                if(probFragment[trans] == 0)
                    continue;
                ulint next = trans % numberOfStates;
                if(probState[next] == 0) continue;
                // real probNext = (t == windowedLength - 2 ? probState_tail[next] : probState[next]);
                // if(probNext == 0)
                //    continue;
                ulint doubleindex = i*numberOfStates + next;

                deltagamma += beta[t][t+1][doubleindex] * alpha[t][t+1][doubleindex];
                                  // / probFragment[trans];
                                // /  probState[next];
                                 // / probNext;
                if(t == windowedLength - 2)
                    gamma[windowedLength - 1][next] += deltagamma; // / probState[i];
                
            }
            gamma[t][i] += deltagamma; // / probState[i];
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
                //byte sec = indexToSequence(windowSize-1, i)[0];
                //probSecondarySymbol[a][sec] += gamma[a][i];
                 byte sec = indexToSequence(windowSize-1, i)[halfWindow];
                 probSecondarySymbol[a+halfWindow][sec] += gamma[a][i];
            }
            // normalize(probSecondarySymbol[a]);
            normalize(probSecondarySymbol[a+halfWindow]);
        }

        /*
        for(uint a=windowedLength; a < fullLength; a++){
            for(ulint i=0; i < numberOfStates; i++){
                byte sec = indexToSequence(windowSize-1, i)[a - windowedLength + 1];
                probSecondarySymbol[a][sec] += gamma[windowedLength-1][i];
            }
            normalize(probSecondarySymbol[a]);
        }
         */
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
        
    // SScpair is intended to be the central pair of secondary structure symbols
    // at the two central positions of a State fragment.
    // In terms of the indexing used in this code they should
    // correspond to positions (halfWindow - 1) and halfWindow.
    // For left states, therefore, the residue identity of the first position has already been emitted
    // while for the second position (halfWindow) it is waiting to be emitted.
    // For right states the second position (halfWindow) has already been emitted while the first position
    // is waiting to be emitted.
    // Secondary structure symbols (0 for E, 1 for H, 2 for C)
    // are obtained from SecondarySymbol%3,
    // and therefore SecondarySymbols must be constructed accordingly (AFPA)
    
    uint SScPair[numberOfStates];
    
    
     ulint number1 = numberOfXtraHalfStates * numberOfSecondarySymbols;
     ulint number2 = numberOfXtraHalfStates;
     ulint number3 = numberOfXtraHalfStates / numberOfSecondarySymbols;
    

/*
    // no need to use the variable numberOfXtraHalfStates, which is only defined when -h is used in the command line

    ulint number3 = 1.0;
    for(uint i = 1; i < windowSize - halfWindow - 1; i++)
        number3 *= numberOfSecondarySymbols;
    ulint number2 = number3 * numberOfSecondarySymbols;
    ulint number1 = number2 * numberOfSecondarySymbols;
   */
    for(uint i = 0; i < numberOfStates; i++){
        uint SScPair1 = i % number1 / number2 % 3;
        uint SScPair2 = i % number2 / number3 % 3;
        SScPair[i] = 3 * SScPair1 + SScPair2;
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
    
    // exclusively double emission for XE-W-EX
    
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
    
    
    // No H anywhere !!! Double E and C and otherwise first LEFT and then RIGHT
    
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
    
    // No H anywhere !!! Double E and C EXCLUSEVELY !!!!
    
    /*
   byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 3, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 3, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 3}};
      */
    
    // No H anywhere !!! C from RIGHT and E from LEFT EXCLUSEVELY !!!!
    
    /*
   byte emissionType[9][9] ={{0, 0, 0, 0, 0, 0, 0, 0, 2},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 2},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 2},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 2}};
      */
    
     /*
    byte emissionType[9][9] ={{0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {1, 1, 1, 0, 0, 0, 1, 0, 1}};
       */
    
   
    // No H anywhere !!! Double E exclusively C first LEFT and then RIGHT !!!!
    
     /*
   byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 2, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 3, 0, 0, 0, 1, 0, 1},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 2, 0, 0}};
      */
    
  
    // No H anywhere !!! Double E exclusively C first RIGHT and then LEFT !!!!
    
    /*
   byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 2, 0, 0, 0, 1, 0, 1},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 3, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 2, 0, 0, 0, 0, 0, 0}};
      */
    
    // No H anywhere !!! Double E exclusively C LEFT exclusively !!!!
    
     /*
   byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 2, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 3, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 2, 0, 0}};
      */
    
  
    // No H anywhere !!! Double E exclusively C exclusively RIGHT !!!!
    
    /*
   byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 1, 0, 1},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 3, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0}};
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
    
     // /*
     byte emissionType[9][9] ={{0, 0, 0, 0, 0, 0, 0, 0, 2},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 2},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 2},
                               {0, 0, 0, 0, 0, 0, 0, 0, 2},
                               {0, 0, 0, 0, 0, 0, 0, 0, 2},
                               {0, 0, 0, 0, 0, 0, 0, 0, 2},
                               {0, 0, 0, 0, 0, 0, 0, 0, 2}};
      //   */
    

  
    // Alternating sides for C and nonC (E or H) fragments beginning at RIGHT, also equivalent to the HMM ?

    /*
   byte emissionType[9][9] ={{0, 0, 0, 0, 0, 0, 2, 2, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 2, 0, 0, 2, 1, 1, 1},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 2, 2, 0},
                             {0, 0, 2, 0, 0, 2, 1, 1, 1},
                             {1, 0, 1, 0, 1, 1, 2, 2, 0},
                             {1, 0, 1, 0, 1, 1, 2, 2, 0},
                             {0, 0, 2, 0, 0, 2, 0, 0, 0}};
     */
    
    // Alternating sides for C and nonC (E or H) fragments beginning at LEFT, also equivalent to the HMM ?

    /*
   byte emissionType[9][9] ={{0, 0, 2, 0, 0, 2, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {1, 0, 1, 0, 1, 1, 2, 2, 2},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 2, 0, 0, 2, 0, 0, 0},
                             {1, 0, 1, 0, 1, 1, 2, 2, 2},
                             {0, 0, 2, 0, 0, 2, 1, 1, 1},
                             {0, 0, 2, 0, 0, 2, 1, 1, 1},
                             {0, 0, 0, 0, 0, 0, 2, 2, 2}};
     */
    
    // Alternating sides for C and nonC (only E and no H) fragments beginning at RIGHT, also equivalent to the HMM ?

    /*
   byte emissionType[9][9] ={{0, 0, 0, 0, 0, 0, 2, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 2, 0, 0, 0, 1, 0, 1},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {1, 0, 1, 0, 0, 0, 2, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 2, 0, 0, 0, 0, 0, 0}};
     */
    
    // Alternating sides for C and nonC (only E and no H) fragments beginning at LEFT, also equivalent to the HMM ?

    /*
   byte emissionType[9][9] ={{0, 0, 2, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {1, 0, 1, 0, 0, 0, 2, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 2, 0, 0, 0, 1, 0, 1},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 2, 0, 0}};
    */
    
    
      // Everybody emits at right, also equivalent to the HMM
    
      /*
     byte emissionType[9][9] ={{0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {1, 0, 1, 0, 1, 1, 1, 1, 1}};
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
    
    
    
    //Everybody emits at left/right, again, also equivalent to the HMM
    
    /*
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
    
      /*
    //Everybody emits at left/right, again, but no HE or HE
    
    byte emissionType[9][9] ={{3, 0, 3, 0, 3, 3, 3, 3, 3},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {3, 0, 3, 0, 3, 3, 3, 3, 3},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {3, 0, 3, 0, 3, 3, 3, 3, 3},
                              {3, 0, 3, 0, 3, 3, 3, 3, 3},
                              {3, 0, 3, 0, 3, 3, 3, 3, 3},
                              {3, 0, 3, 0, 3, 3, 3, 3, 3},
                              {3, 0, 3, 0, 3, 3, 3, 3, 3}};
       */
    
   
    
    ulint numberOfDoubleStates = numberOfStates * numberOfStates;
    
    vector < vector < vector<real> > > alpha(windowedLength, vector < vector<real> > (windowedLength, vector <real> (numberOfDoubleStates, 0.0)));
    vector < vector < vector<real> > > beta(windowedLength, vector < vector<real> > (windowedLength, vector <real> (numberOfDoubleStates, 0.0)));
    vector < vector < vector<real> > > probTrans(windowedLength, vector < vector<real> > (windowedLength, vector <real> (numberOfFragments, 0.0)));
    
    vector < real > probDoubleState(numberOfDoubleStates, 0.0);
    for(ulint left = 0; left < numberOfStates; left ++){
        if(probState[left] == 0)
            continue;
        for(ulint right = 0; right < numberOfStates; right ++){
            if(probState[right] == 0)
                continue;
            byte emission = emissionType[SScPair[left]][SScPair[right]];
            if(emission == 0)
                continue;
            ulint doubleIndex = left * numberOfStates + right;
            probDoubleState[doubleIndex] = probState[left] * probState[right];
        }
    }
//    normalize(probDoubleState);
    
    
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
   /*
    for(int t=0; t < windowedLength-1; t++){
        for(ulint i=0; i<numberOfFragments; i++){
            ulint doubleIndex = i / numberOfSecondarySymbols * numberOfStates + i % numberOfStates;
            alpha[t][t+1][doubleIndex] = probFragment[i] *
                probFragmentEmitsPrimarySymbol[i][ seq[t+halfWindow] ];
        }
        normalize(alpha[t][t+1]);
    }
    */
    
    for(int t=0; t < windowedLength-1; t++){
        for(ulint trans=0; trans<numberOfFragments; trans++){
            
            ulint left = trans / numberOfSecondarySymbols;
            ulint right = trans % numberOfStates;
            if(probState[left] == 0 || probState[right] == 0)
                continue;
            ulint doubleIndex = left * numberOfStates + right;
     //       ulint doubleIndex = trans / numberOfSecondarySymbols * numberOfStates + trans % numberOfStates;
            probTrans[t][t+1][trans] = probFragment[trans] *
                probFragmentEmitsPrimarySymbol[trans][ seq[t+halfWindow] ];
           
            alpha[t][t+1][doubleIndex] = probFragmentEmitsPrimarySymbol[trans][ seq[t+halfWindow] ]
                                                 * probFragment[trans];
                                                 // / probState[left] / probState[right];
        }
        normalize(probTrans[t][t+1]);
        normalize(alpha[t][t+1]);
    }
    
    //
    // 1.2 Induction
    //
    for(int v=2; v < windowedLength; v++){
        
        for (int t=v-2; t >= 0; t--) {
        
            for(ulint nextLeft=0; nextLeft<numberOfStates; nextLeft++){
                if(probState[nextLeft] == 0) continue;
                // real probPrevLeft = (t > 0 ? probState[prevLeft] : probState_head[prevLeft]);
                
                for(ulint prevRight=0; prevRight<numberOfStates; prevRight++){
                    if(probState[prevRight] == 0) continue;
                    // real probNextRight = (t < windowedLength - 1 ? probState[nextRight] : probState_tail[nextRight]);
                    
                  //  byte emission = emissionType[SScPair[prevLeft]][SScPair[nextRight]];
                  //  if(emission == 0) continue;
                    
                    if(t == v-2 && nextLeft != prevRight)
                        continue;
                    ulint childDoubleIndex = nextLeft * numberOfStates + prevRight;
                        
                    for(uint tt=0; tt<numberOfSecondarySymbols; tt++){
                                
                        ulint transRight = prevRight * numberOfStates + tt;
                        ulint nextRight = transRight % numberOfStates;
                        // if(probState[nextRight] == 0) continue;
                        ulint doubleIndexRight = prevRight * numberOfStates + nextRight;
                            
                        if ( emissionType[SScPair[nextLeft]][SScPair[prevRight]] == 1 ){
                            // if ( emission == 1 && probNextRight ){
                            ulint doubleIndex = nextLeft*numberOfStates + nextRight;
                            
                            real deltaAlpha = alpha[t][v-1][childDoubleIndex]
                                                * probTrans[v-1][v][transRight]
                                                / probState[prevRight];
                                                // / probNextRight;
                            alpha[t][v][doubleIndex] += deltaAlpha;
                        }
                    }
                            
                    for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                                    
                        ulint transLeft = ss * numberOfSecondarySymbols + nextLeft;
                        ulint prevLeft = transLeft / numberOfSecondarySymbols;
                        // if(probState[prevLeft] == 0) continue;
                        ulint doubleIndexLeft = prevLeft*numberOfStates + nextLeft;
        
                        if ( emissionType[SScPair[nextLeft]][SScPair[prevRight]] == 2 ){
                            // else if ( emission == 2 && probPrevLeft ){
                            ulint doubleIndex = prevLeft*numberOfStates + prevRight;
                                
                            real deltaAlpha = alpha[t+1][v][childDoubleIndex]
                                                * probTrans[t][t+1][transLeft]
                                                / probState[nextLeft];
                                                // / probPrevLeft;
                            alpha[t][v][doubleIndex] += deltaAlpha;
                        }
                    }
                           
                    for(uint tt=0; tt<numberOfSecondarySymbols; tt++){
                                
                        ulint transRight = prevRight * numberOfStates + tt;
                        ulint nextRight = transRight % numberOfStates;
                        // if(probState[nextRight] == 0) continue;
                        ulint doubleIndexRight = prevRight * numberOfStates + nextRight;
                        
                        for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                                        
                            ulint transLeft = ss * numberOfSecondarySymbols + nextLeft;
                            ulint prevLeft = transLeft / numberOfSecondarySymbols;
                            // if(probState[prevLeft] == 0) continue;
                            ulint doubleIndexLeft = prevLeft*numberOfStates + nextLeft;
                            
                            if (emissionType[SScPair[prevLeft]][SScPair[nextRight]] == 3) {
                                ulint doubleIndex = prevLeft*numberOfStates + nextRight;
                                
                                real deltaAlpha = 0.0;
                                
                                if(t == v - 2 && nextLeft == prevRight) {
                                    deltaAlpha =  probTrans[t+1][v][transRight]
                                                    * probTrans[t][t+1][transLeft]
                                                    / probState[nextLeft]; // / probState[nextRight];
                                }
                                if( (v - t) > 2){
                                    //ulint childDoubleIndex = nextLeft*numberOfStates + prevRight;
                                    deltaAlpha = alpha[t+1][v-1][childDoubleIndex]
                                                    * probTrans[t][t+1][transLeft]
                                                    / probState[nextLeft]
                                                    * probTrans[v-1][v][transRight]
                                                    / probState[prevRight];
                                }
                                
                                 /*
                                else if(t == 0 && v == windowedLength - 1) {
                                
                                    deltaAlpha =  (alpha[t][v-1][childDoubleIndex]
                                                    * probTrans[v-1][v][transRight]
                                                    / probState[prevRight])
                                                    +
                                                    (alpha[t+1][v][childDoubleIndex]
                                                    * probTrans[t][t+1][transLeft]
                                                     / probState[nextLeft]); //; // / probState[nextRight];
                                   
                                }
                                 */
                                
                                alpha[t][v][doubleIndex] += deltaAlpha;

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
    // will have non zero probability at the beginning and end of the chain.
    //
    
    if(headAndTail){
        for(ulint i=0; i < numberOfStates; i++){
            if(probState[i] == 0) // || SScPair[i] != 2)
                continue;
            for(ulint j=0; j < numberOfStates; j++){
                if(probState[j] == 0) // || SScPair[j] != 6)
                    continue;
                if(emissionType[SScPair[i]][SScPair[j]]){
                        beta[0][windowedLength-1][i*numberOfStates + j] =
                            //probState[i] * probState[j];
                            probState_head[i] * probState_tail[j];
                }
            }
        }
    }
    else{
        for(ulint i=0; i < numberOfStates; i++)
         
            for(ulint j=0; j < numberOfStates; j++)
                if(emissionType[SScPair[i]][SScPair[j]] > 0)
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
                 
                    // if(probState[nextRight] == 0) continue;
                    
                    real emission = emissionType[SScPair[prevLeft]][SScPair[nextRight]];
                    if(emission == 0) continue;
                    
                    ulint parentDoubleIndex = prevLeft * numberOfStates + nextRight;
                    
                    for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                        ulint transLeft = prevLeft * numberOfSecondarySymbols + ss;
                        // if(probFragment[transLeft] == 0) continue;
                        // if(probTrans[t-1][t][transLeft] == 0) continue;
                        ulint nextLeft = transLeft % numberOfStates;
                        if(probState[nextLeft] == 0) continue;
                        ulint doubleIndexLeft = prevLeft * numberOfStates + nextLeft;
                    
                    
                        for(uint tt=0; tt<numberOfSecondarySymbols; tt++){
                            ulint transRight = tt * numberOfStates + nextRight;
                            // if(probFragment[transRight] == 0) continue;
                            // if(probTrans[v][v+1][transRight] == 0) continue;
                            ulint prevRight = transRight /  numberOfSecondarySymbols;
                            if(probState[prevRight] == 0) continue;
                            ulint doubleIndexRight = prevRight * numberOfStates + nextRight;
                            
                            if(v == t+1 && nextLeft != prevRight)
                                continue;
                    
                            if((emission == 1 && v < windowedLength - 1) ){
                               // || (t == 0 && v == windowedLength - 2) ){
                                              //  && (v > t + 1 || nextLeft == prevRight)){ // && beta[t][v+1][parentDoubleIndex]){
                                
                                if(probTrans[v][v+1][transRight] == 0) continue;
                             
                                ulint doubleIndex = prevLeft*numberOfStates + prevRight;
                                // if(probDoubleState[doubleIndex] == 0 && v > t+1)
                                //    continue;
                                
                                real probNextRight = (v == windowedLength - 2 ? probState_tail[nextRight] : probState[nextRight]);
                                
                                if(probNextRight){
                                    real deltaBeta = beta[t][v+1][parentDoubleIndex]
                                                    * probTrans[v][v+1][transRight]
                                                    // / probNextRight;
                                                     / probState[nextRight];
                                                    // / probState[prevLeft];
                                    beta[t][v][doubleIndex] += deltaBeta;
                                }
                         
                                //this line corresponds to the termination step that was previously in separate loop
                                // beta[v][v+1][doubleIndexRight] += beta[t][v+1][parentDoubleIndex]
                                //                                * alpha[t][v][doubleIndex]
                                //                                / probState[prevLeft];
                            }
                    
                            else if(emission == 2 && t > 0 ){
                                                  //  && (v > t + 1 || nextLeft == prevRight)){ // && beta[t-1][v][parentDoubleIndex]){
                    
                                if(probTrans[t-1][t][transLeft] == 0) continue;
                                  
                                ulint doubleIndex = nextLeft*numberOfStates + nextRight;
                                // if(probDoubleState[doubleIndex] == 0 && v > t+1)
                                //    continue;

                                real probPrevLeft = (t == 1 ? probState_head[prevLeft] : probState[prevLeft]);
                                if(probPrevLeft){
                                    real deltaBeta = beta[t-1][v][parentDoubleIndex]
                                                    * probTrans[t-1][t][transLeft]
                                                    // / probPrevLeft;
                                                    / probState[prevLeft];
                                                    // * probState[nextRight];
                                    beta[t][v][doubleIndex] += deltaBeta;
                                }
                            
                                //this line corresponds to the termination step that was previously in separate loop
                                // beta[t-1][t][doubleIndexLeft] += beta[t-1][v][parentDoubleIndex]
                                //                                * alpha[t][v][doubleIndex]
                                //                                / probState[nextRight];
                            }
                            
                    
                            else if(emission == 3) {
                                //   && t > 0 && v < windowedLength -1)
                            
                                
                                // if(emissionType[SScPair[nextLeft]][SScPair[prevRight]] == 0)
                                //      continue;
                                
                                ulint doubleIndex = nextLeft*numberOfStates + prevRight;
                                ulint doubleIndex1 = nextLeft*numberOfStates + nextRight;
                                ulint doubleIndex2 = prevLeft*numberOfStates + prevRight;
                                
                                // real deltaBeta = 0.0;
                                if(t > 0 && v < windowedLength - 1 // && v > t + 1
                                   && beta[t-1][v+1][parentDoubleIndex]){
                                    real deltaBeta  =   beta[t-1][v+1][parentDoubleIndex]
                                                        * probTrans[t-1][t][transLeft]
                                                        * probTrans[v][v+1][transRight]
                                                        / probState[prevLeft]
                                                        / probState[nextRight];
                                    beta[t][v][doubleIndex] += deltaBeta;
                                    
                                    //these lines corresponds to the termination step
                                    
                                    /*
                                    beta[t-1][t][doubleIndexLeft] += beta[t-1][v+1][parentDoubleIndex]
                                                                        * alpha[t-1][v+1][parentDoubleIndex]
                                                                        / probTrans[t-1][t][transLeft]
                                                                        * probState[nextLeft]
                                                                        / probState[prevRight]
                                                                        / probState[nextRight];
                                                                        // / probFragment[transRight];
                                    beta[v][v+1][doubleIndexRight] += beta[t-1][v+1][parentDoubleIndex]
                                                                        * alpha[t-1][v+1][parentDoubleIndex]
                                                                        / probTrans[v][v+1][transRight]
                                                                        * probState[prevRight]
                                                                        / probState[prevLeft]
                                                                        / probState[nextLeft];
                                                                        // / probFragment[transLeft];
                                    */
                                    
                                     /*
                                    beta[t-1][t][doubleIndexLeft] += beta[t-1][v+1][parentDoubleIndex]
                                                                        // * alpha[t][v+1][doubleIndex1]
                                                                        // * alpha[v][v+1][doubleIndexRight]
                                                                        * probTrans[v][v+1][transRight]
                                                                        * alpha[t][v][doubleIndex]
                                                                        / probState[prevRight]
                                                                        / probState[nextRight];
                                                                        // * probFragment[transLeft]
                                                                        // / probState[prevLeft];
                                    beta[v][v+1][doubleIndexRight] += beta[t-1][v+1][parentDoubleIndex]
                                                                        // * alpha[t-1][v][doubleIndex2]
                                                                        * probTrans[t-1][t][transLeft]
                                                                        * alpha[t][v][doubleIndex]
                                                                        / probState[prevLeft]
                                                                        / probState[nextLeft];
                                                                        // * probFragment[transRight]
                                                                        // / probState[nextRight];
                                     */
                                }
                                
                                /*
                                else if(t == 0 && v == windowedLength - 2){
                                    
                                    ulint parentDoubleIndex1 = nextLeft * numberOfStates + nextRight;
                                    real deltaBeta = beta[t][v+1][parentDoubleIndex1]
                                                        * probTrans[v][v+1][transRight]
                                                        / probState[nextRight];
                                    beta[t][v][doubleIndex] += deltaBeta;
                                    
                                    //this line corresponds to the termination step that was previously in separate loop
                                    beta[v][v+1][doubleIndexRight] += beta[t][v+1][parentDoubleIndex1]
                                                                        * alpha[t][v+1][parentDoubleIndex1]
                                                                        / probTrans[v][v+1][transRight]
                                                                        * probState[prevRight]
                                                                        / probState[nextLeft];
                                    
                                }
                                
                                else if(t == 1 && v == windowedLength - 1){
                                    
                                    ulint parentDoubleIndex2 = prevLeft * numberOfStates + prevRight;
                                    real deltaBeta = beta[t-1][v][parentDoubleIndex2]
                                                        * probTrans[t-1][t][transLeft]
                                                        / probState[prevLeft];
                                    beta[t][v][doubleIndex] += deltaBeta;
                                    
                                    //this line corresponds to the termination step that was previously in separate loop
                                    beta[t-1][t][doubleIndexLeft] += beta[t-1][v][parentDoubleIndex2]
                                                                        * alpha[t-1][v][parentDoubleIndex2]
                                                                        / probTrans[t-1][t][transLeft]
                                                                        * probState[nextLeft]
                                                                        / probState[prevRight];
                                    
                                }
                                 */
                                
                                /*
                                if(v == t+1 && nextLeft == prevRight){
                                    
                                    if(t > 0 && beta[t-1][v][parentDoubleIndex]){
                                        beta[t][t+1][doubleIndexRight] += beta[t-1][v][parentDoubleIndex]
                                                                        * alpha[t-1][t][doubleIndexLeft]
                                                                        / probState[prevLeft];
                                    
                                        //this line corresponds to the termination step
                                        beta[t-1][t][doubleIndexLeft] += beta[t-1][v][parentDoubleIndex]
                                                                            * alpha[t][v][doubleIndexRight]
                                                                            / probState[nextRight];
                                    }
                                    if(v < windowedLength-1 && beta[t][v+1][parentDoubleIndex]){
                                        beta[t][t+1][doubleIndexLeft] += beta[t][v+1][parentDoubleIndex]
                                                                        * alpha[v][v+1][doubleIndexRight]
                                                                        / probState[nextRight];
                                    
                                        //this line corresponds to the termination step
                                        beta[v][v+1][doubleIndexRight] += beta[t][v+1][parentDoubleIndex]
                                                                            * alpha[t][t+1][doubleIndexLeft]
                                                                            / probState[prevLeft];
                                    }
                                    
                                }
                                 */
                                
                                /*
                                if(v == t+1 && t > 0 && nextLeft == prevRight
                                                    && beta[t-1][v][parentDoubleIndex]){
                                    beta[t][v][doubleIndexRight] +=  beta[t-1][v][parentDoubleIndex]
                                                                * alpha[t-1][t][doubleIndexLeft]
                                                                / probState[prevLeft];
                                    
                                    //this line corresponds to the termination step
                                    beta[t-1][t][doubleIndexLeft] += beta[t-1][v][parentDoubleIndex]
                                                                        * alpha[t][v][doubleIndex]
                                                                        / probState[nextRight];
                                }
                                */
                                 
                         //       beta[t][v][doubleIndex] += deltaBeta;
                                
                            }
                        }
                    }
                }
            }
            
       //     if(v > t+1)
                normalize(beta[t][v]);

        }
    }
    
    for(int t = 0; t < windowedLength - 1; t++){
        
        for(ulint trans=0; trans < numberOfFragments; trans++){
            if(probFragment[trans] == 0)
                continue;
            ulint prev = trans / numberOfSecondarySymbols;
            ulint next = trans % numberOfStates;
            ulint doubleIndex0 = prev * numberOfStates + next;
    
            //
            // possible production with other possibly long fragments at left of t
            // current short fragment emits at right
            //
            
            for(ulint prevLeft = 0; prevLeft < numberOfStates; prevLeft++){
                if(probState[prevLeft] == 0)
                    continue;
                ulint doubleIndex2 = prevLeft*numberOfStates + next;
                ulint doubleIndex1 = prevLeft*numberOfStates + prev;
                    
                if(emissionType[SScPair[prevLeft]][SScPair[next]] == 1  ) {
                    // && emissionType[SScPair[prevLeft]][SScPair[prev]] == 1) {
                    
                  
             //       if(emissionType[SScPair[prevLeft]][SScPair[prev]] == 1 ) {
                        
                        real deltaBeta = 0.0;
                        for(int w = 1; w < t; w++){
                            deltaBeta += beta[w][t+1][doubleIndex2] *
                                            alpha[w][t][doubleIndex1];
                        }
                        
                        beta[t][t+1][doubleIndex0] += deltaBeta
                                                    // / probState[next];
                                                    * probFragment[trans] / probState[next]
                                                     / probState[prevLeft];
                  
                        if(t > 0 && probState_head[prevLeft]){
                            beta[t][t+1][doubleIndex0] += beta[0][t+1][doubleIndex2]
                                                        * alpha[0][t][doubleIndex1]
                                                        * probFragment[trans] / probState[next]
                                                        / probState_head[prevLeft];
                        }
                    
               //     }
                    
                }
                
                /*
                if(emissionType[SScPair[prevLeft]][SScPair[next]] == 2 && t > 1) {
                    real probLeft = (t == 1 ? probState_head[prevLeft] : probState[prevLeft]);
                    if(probLeft) {
                        beta[t][t+1][doubleIndex0] += beta[t-2][t+1][doubleIndex2]
                                                        * alpha[t-2][t][doubleIndex1]
                                                        / probLeft;
                    }
                }
                 */
                
                /*
                if(emissionType[SScPair[prevLeft]][SScPair[prev]] == 2 && t > 1) {
                    
                    real deltaBeta = 0.0;
                    for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                        ulint transLeft = prevLeft * numberOfSecondarySymbols + ss;
                        ulint nextLeft = transLeft % numberOfStates;
                        if(probState[nextLeft] == 0) continue;
                        ulint doubleIndex1 = nextLeft*numberOfStates + prev;
                        ulint doubleIndex3 = prevLeft*numberOfStates + nextLeft;
                        
                        deltaBeta += beta[t-2][t+1][doubleIndex2]
                                        * alpha[t-2][t-1][doubleIndex3]
                                        * alpha[t-1][t][doubleIndex1]
                                        / probState[nextLeft];
                    }
                        
                    beta[t][t+1][doubleIndex0] += deltaBeta / probState[prevLeft];
                }
                 */
                
                
                    
                if(emissionType[SScPair[prevLeft]][SScPair[next]] == 3){
                    
                    for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                        ulint transLeft = prevLeft * numberOfSecondarySymbols + ss;
                        ulint nextLeft = transLeft % numberOfStates;
                        if(probState[nextLeft] == 0) continue;
                        
                        ulint doubleIndex1 = nextLeft * numberOfStates + prev;
                        ulint doubleIndex3 = prevLeft * numberOfStates + nextLeft;
                    
                        real deltaBeta = 0.0;
                        for(int w = 0; w < t; w++){
                            deltaBeta += beta[w][t+1][doubleIndex2]
                                            * alpha[w][w+1][doubleIndex3]
                                            * alpha[w+1][t][doubleIndex1];
                        }
                        beta[t][t+1][doubleIndex0] += deltaBeta
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
                ulint doubleIndex2 = prev*numberOfStates + nextRight;
                ulint doubleIndex1 = next*numberOfStates + nextRight;
                        
                if(emissionType[SScPair[prev]][SScPair[nextRight]] == 2 // ) {
                     && emissionType[SScPair[next]][SScPair[nextRight]] == 2 ) {
                    
               //     if(emissionType[SScPair[next]][SScPair[nextRight]] == 2 ) {
                    if(t > 0) {
                        real deltaBeta = 0.0;
                        for(int w = windowedLength - 2; w > t+1; w--){
                            deltaBeta += beta[t][w][doubleIndex2]
                                            * alpha[t+1][w][doubleIndex1];
                        }
                        beta[t][t+1][doubleIndex0] += deltaBeta
                                                    // / probState[prev];
                                                    * probFragment[trans] / probState[prev]
                                                     / probState[nextRight];
                    // /*
                        if(t < windowedLength-2 && probState_tail[nextRight]){
                            beta[t][t+1][doubleIndex0] += beta[t][windowedLength-1][doubleIndex2]
                                                        * alpha[t+1][windowedLength-1][doubleIndex1]
                                                        * probFragment[trans] / probState[prev]
                                                        / probState_tail[nextRight];
                        }
                    // */
                    }
                    // /*
                    else if(probState_head[prev]){
                        real deltaBeta = 0.0;
                        for(int w = windowedLength - 2; w > t+1; w--){
                            deltaBeta += beta[t][w][doubleIndex2]
                                            * alpha[t+1][w][doubleIndex1];
                        }
                        beta[t][t+1][doubleIndex0] += deltaBeta
                                                    // / probState[prev];
                                                     * probFragment[trans] / probState_head[prev]
                                                     / probState[nextRight];
                
                        if(t < windowedLength-2 && probState_tail[nextRight]){
                            beta[t][t+1][doubleIndex0] += beta[t][windowedLength-1][doubleIndex2]
                                                        * alpha[t+1][windowedLength-1][doubleIndex1]
                                                         * probFragment[trans] / probState_head[prev]
                                                        / probState_tail[nextRight];
                        }
                    //
                    }
                    // */
                        
                }
                 
                /*
                if(t < windowedLength - 3 && emissionType[SScPair[prev]][SScPair[nextRight]] == 1 ) {
                    real probRight = (t == windowedLength - 4 ? probState_tail[nextRight] : probState[nextRight]);
                    if(probRight) {
                        beta[t][t+1][doubleIndex0] += beta[t][t+3][doubleIndex2]
                                                        * alpha[t+1][t+3][doubleIndex1]
                                                        / probRight;
                    }
                }
                 */
                
                /*
                if(emissionType[SScPair[prev]][SScPair[nextRight]] == 1 && t < windowedLength - 2) {
                    
                    ulint doubleIndex1 = next*numberOfStates + nextRight;
                    beta[t][t+1][doubleIndex0] += beta[t][t+2][doubleIndex2] *
                                    alpha[t+1][t+2][doubleIndex1];
                }
                */
                 
                if(emissionType[SScPair[prev]][SScPair[nextRight]] == 3){
                        
                    for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                        ulint transRight = nextRight+(numberOfStates * ss);
                        ulint prevRight = transRight / numberOfSecondarySymbols;
                        if (probState[prevRight] == 0) continue;
                            
                        ulint doubleIndex1 = next * numberOfStates + prevRight;
                        ulint doubleIndex3 = prevRight * numberOfStates + nextRight;
                        
                        real deltaBeta = 0.0;
                        for(int w = windowedLength - 1; w > t+1; w--){
                            deltaBeta += beta[t][w][doubleIndex2]
                                                * alpha[w-1][w][doubleIndex3]
                                               // / probState[nextRight]
                                                * alpha[t+1][w-1][doubleIndex1];
                                                // / probState[prevRight];
                        }
                        beta[t][t+1][doubleIndex0] += deltaBeta
                                                    / probState[nextRight]
                                                    / probState[prevRight];
                        if((next == prevRight) && (t+1 < windowedLength - 1)){
                            beta[t][t+1][doubleIndex0] += beta[t][t+2][doubleIndex2]
                                                            * alpha[t+1][t+2][doubleIndex3]
                                                            / probState[nextRight];
                        }
                    }
                }
            }
        
        }
    
        normalize(beta[t][t+1]);
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
                if(probFragment[trans] == 0)
                    continue;
                ulint next = trans % numberOfStates;
                if(probState[next] == 0)
                    continue;
                ulint doubleindex = i*numberOfStates + next;

                 deltagamma += beta[t][t+1][doubleindex] * alpha[t][t+1][doubleindex]
                                  // * probFragment[trans]
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
                //byte sec = indexToSequence(windowSize-1, i)[0];
                //probSecondarySymbol[a][sec] += gamma[a][i];
                 byte sec = indexToSequence(windowSize-1, i)[halfWindow];
                 probSecondarySymbol[a+halfWindow][sec] += gamma[a][i];
            }
            // normalize(probSecondarySymbol[a]);
            normalize(probSecondarySymbol[a+halfWindow]);
        }

        /*
        for(uint a=windowedLength; a < fullLength; a++){
            for(ulint i=0; i < numberOfStates; i++){
                byte sec = indexToSequence(windowSize-1, i)[a - windowedLength + 1];
                probSecondarySymbol[a][sec] += gamma[windowedLength-1][i];
            }
            normalize(probSecondarySymbol[a]);
        }
         */
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
    
    // SScpair is intended to be the central pair of secondary structure symbols
    // at the two central positions of a State fragment.
    // In terms of the indexing used in this code they should
    // correspond to positions (halfWindow - 1) and halfWindow.
    // For left states, therefore, the residue identity of the first position has already been emitted
    // while for the second position (halfWindow) it is waiting to be emitted.
    // For right states the second position (halfWindow) has already been emitted while the first position
    // is waiting to be emitted.
    // Secondary structure symbols (0 for E, 1 for H, 2 for C)
    // are obtained from SecondarySymbol%3,
    // and therefore SecondarySymbols must be constructed accordingly (AFPA)
  
    
    uint SScPair[numberOfStates];
    
    ulint number1 = numberOfXtraHalfStates * numberOfSecondarySymbols;
    ulint number2 = numberOfXtraHalfStates;
    ulint number3 = numberOfXtraHalfStates / numberOfSecondarySymbols;
    
    // no need to use the variable numberOfXtraHalfStates, which is only defined when -h is used in the command line
 
   /*
    ulint number3 = 1.0;
    for(uint i = 1; i < windowSize - halfWindow - 1; i++)
        number3 *= numberOfSecondarySymbols;
    ulint number2 = number3 * numberOfSecondarySymbols;
    ulint number1 = number2 * numberOfSecondarySymbols;
    */
    
    for(uint i = 0; i < numberOfStates; i++){
        uint SScPair1 = i % number1 / number2 % 3;
        uint SScPair2 = i % number2 / number3 % 3;
        SScPair[i] = 3 * SScPair1 + SScPair2;
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
    
    
    // exclusively double emission for XE-W-EX
    
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
    
    //  No H anywhere !!! E from RIGHT and C from LEFT
    
    /*
   byte emissionType[9][9] ={{2, 0, 2, 0, 0, 0, 2, 0, 2},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {2, 0, 2, 0, 0, 0, 2, 0, 2},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {2, 0, 2, 0, 0, 0, 2, 0, 2},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {2, 0, 2, 0, 0, 0, 2, 0, 2}};
      */
    
     /*
    byte emissionType[9][9] ={{0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 2, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {1, 0, 1, 0, 0, 0, 2, 0, 1}};
       */
    
    
    //  No H anywhere !!! Double E and C EXCLUSEVELY !!!!
    
    // /*
   byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 3, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 3, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 3}};
     // */
    
    // No H anywhere !!! Double E exclusively C first LEFT and then RIGHT !!!!
    
    /*
   byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 2, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 3, 0, 0, 0, 1, 0, 1},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 2, 0, 0, 0, 2, 0, 2}};
      */
    
  
   //  No H anywhere !!! Double E exclusively C first RIGHT then LEFT !!!!
    
    /*
   byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 2, 0, 0, 0, 1, 0, 1},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 3, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 2, 0, 0, 0, 0, 0, 3}};
      */
    
    // No H anywhere !!! Double E exclusively C LEFT then RIGHT, except for double in the the ends !!!!
    
     /*
   byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 2, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 3, 0, 0, 0, 1, 0, 1},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 2, 0, 2}};
      */
    
  
    // No H anywhere !!! Double E exclusively C exclusively RIGHT !!!!
    
    /*
   byte emissionType[9][9] ={{3, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 3, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 3, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 3}};
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
    
    // Everybody emits at left, equivalent to the HMM

    /*
    byte emissionType[9][9] ={{0, 0, 0, 0, 0, 0, 0, 0, 2},
                              {0, 0, 0, 0, 0, 0, 0, 0, 2},
                              {0, 0, 0, 0, 0, 0, 0, 0, 2},
                              {0, 0, 0, 0, 0, 0, 0, 0, 2},
                              {0, 0, 0, 0, 0, 0, 0, 0, 2},
                              {0, 0, 0, 0, 0, 0, 0, 0, 2},
                              {0, 0, 0, 0, 0, 0, 0, 0, 2},
                              {0, 0, 0, 0, 0, 0, 0, 0, 2},
                              {0, 0, 0, 0, 0, 0, 0, 0, 2}};
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
    
    
    //
    //Removes probabilities for states containing any 'EH' or 'HE' pair of secondary structures,
    //which are inconsistent with the grammar but might appear rarely in the training data bank.
    //
    
    /*
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
     */
    
    ulint numberOfDoubleStates = numberOfStates * numberOfStates;
    
    unordered_map<ulint, ulint> betaIndex;
 //   unordered_map<ulint, ulint> alphaIndex;
    
    ulint maxBetaIndex = 0;
    for(ulint i = 0; i < numberOfStates; i++){
        if(probState[i] == 0) continue;
        for(ulint j = 0; j < numberOfStates; j++){
            if(probState[j] == 0) continue;
            if(emissionType[SScPair[i]][SScPair[j]]) {
                ulint doubleIndex = i*numberOfStates +j;
                betaIndex[doubleIndex] = maxBetaIndex;
                maxBetaIndex ++;
            }
        }
    }
    
    /*
    ulint maxAlphaIndex = 0;
    for(ulint i = 0; i < numberOfStates; i++){
        if(probState[i] == 0) continue;
        for(ulint j = 0; j < numberOfStates; j++){
            if(probState[j] == 0) continue;
            
            ulint doubleIndex = i*numberOfStates +j;
            alphaIndex[doubleIndex] = maxAlphaIndex;
            maxAlphaIndex ++;
            }
        }
    }
    
    */
    
    vector < vector < vector<real> > > alpha(windowedLength, vector < vector<real> > (windowedLength, vector <real> (numberOfDoubleStates, 0.0)));
  //  vector < vector < vector<real> > > beta(windowedLength, vector < vector<real> > (windowedLength+1, vector <real> (numberOfDoubleStates, 0.0)));
    
    
 //   vector < vector < vector<real> > > alpha(windowedLength, vector < vector<real> > (windowedLength+1, vector <real> (maxIndex, 0.0)));
    vector < vector < vector<real> > > beta(windowedLength, vector < vector<real> > (windowedLength, vector <real> (maxBetaIndex, 0.0)));
    
 //   vector < vector <real> > alphaTTplusOne (windowedLength,
 //                                            vector <real> (numberOfFragments, 0.0));
    vector < vector <real> > betaTTPlusOne (windowedLength,
                                            vector <real> (numberOfFragments, 0.0));
     

    
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
            
            for(ulint prevLeft=0; prevLeft<numberOfStates; prevLeft++){
                if(probState[prevLeft] == 0) continue;
                for(ulint nextRight=0; nextRight<numberOfStates; nextRight++){
                    if(probState[nextRight] == 0) continue;
                    
                    byte emission = emissionType[SScPair[prevLeft]][SScPair[nextRight]];
                    if(emission == 0) continue;
                    
                    ulint doubleIndex = prevLeft * numberOfStates + nextRight;
                
                    if( emission == 1){
                        for(uint tt=0; tt<numberOfSecondarySymbols; tt++){
                                
                            ulint transRight = tt * numberOfStates + nextRight;
                            ulint prevRight = transRight / numberOfSecondarySymbols;
                            if(probState[prevRight] == 0) continue;
                            ulint doubleIndexRight = prevRight * numberOfStates + nextRight;
                            ulint childDoubleIndex = prevLeft*numberOfStates + prevRight;
                            
                            real deltaAlpha = alpha[t][v-1][childDoubleIndex]
                                                * alpha[v-1][v][doubleIndexRight]
                                                / probState[prevRight];
                            alpha[t][v][doubleIndex] += deltaAlpha;
                        
                        }
                    }
            
                    else if(emission ==2 ) {

                        for(uint ss=0; ss<numberOfSecondarySymbols; ss++){
                             
                            ulint transLeft = prevLeft * numberOfSecondarySymbols + ss;
                            ulint nextLeft = transLeft % numberOfStates;
                            if(probState[nextLeft] == 0) continue;
                            ulint doubleIndexLeft = prevLeft*numberOfStates + nextLeft;
                            ulint childDoubleIndex = nextLeft*numberOfStates + nextRight;
                                
                            real deltaAlpha = alpha[t+1][v][childDoubleIndex]
                                                * alpha[t][t+1][doubleIndexLeft]
                                                / probState[nextLeft];
                            alpha[t][v][doubleIndex] += deltaAlpha;
                        
                        }
                    }
                    
                    else if (emission == 3) {
                            
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
                                
                                ulint childDoubleIndex = nextLeft*numberOfStates + prevRight;
                                
                                real deltaAlpha = 0.0;
                                if(t == v - 2 && nextLeft == prevRight) {
                                    deltaAlpha =  alpha[t+1][v][doubleIndexRight]
                                                    * alpha[t][t+1][doubleIndexLeft]
                                                    / probState[nextLeft];
                                }
                                else if(t < v - 2){
                                    deltaAlpha = alpha[t+1][v-1][childDoubleIndex]
                                                    * alpha[t][t+1][doubleIndexLeft]
                                                    / probState[nextLeft]
                                                    * alpha[v-1][v][doubleIndexRight]
                                                    / probState[prevRight];
                                }
                               
                                alpha[t][v][doubleIndex] += deltaAlpha;
                            
                              //  /*
                                if(t == 0 && v == windowedLength - 1 && (windowedLength % 2)){
                                    alpha[t][v][doubleIndex] += (alpha[t][v-1][childDoubleIndex]
                                                                    * alpha[v-1][v][doubleIndexRight]
                                                                    /probState[prevRight])
                                                                +
                                                                (alpha[t+1][v][childDoubleIndex]
                                                                    * alpha[t][t+1][doubleIndexLeft]
                                                                 / probState[nextLeft]);
                                }
                              //   */
                                
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
            if(probState[i] == 0 || SScPair[i] != 8)
                continue;
            for(ulint j=0; j < numberOfStates; j++){
                if(probState[j] == 0 || SScPair[j] != 8)
                    continue;
                if(emissionType[SScPair[i]][SScPair[j]]){
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
                if(emissionType[SScPair[i]][SScPair[j]]){
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
                    
                    real emission = emissionType[SScPair[prevLeft]][SScPair[nextRight]];
                    if(emission == 0) continue;
                    
                    ulint parentDoubleIndex = prevLeft * numberOfStates + nextRight;
                    ulint indexOfParentDoubleIndex = betaIndex.at(parentDoubleIndex);
                    
                    if(emission == 1 && v < windowedLength - 1 && beta[t][v+1][indexOfParentDoubleIndex])
                    {
                        for(uint tt=0; tt<numberOfSecondarySymbols; tt++)
                        {
                            ulint transRight = tt * numberOfStates + nextRight;
                            ulint prevRight = transRight /  numberOfSecondarySymbols;
                            if(probState[prevRight] == 0) continue;
                            ulint doubleIndexRight = prevRight * numberOfStates + nextRight;
                            
                            if(emissionType[SScPair[prevLeft]][SScPair[prevRight]] == 0)
                                continue;
                            
                            ulint doubleIndex = prevLeft*numberOfStates + prevRight;
                                
                            real deltaBeta = beta[t][v+1][indexOfParentDoubleIndex]
                                                * alpha[v][v+1][doubleIndexRight]
                                                / probState[nextRight];
                            beta[t][v][betaIndex.at(doubleIndex)] += deltaBeta;
                            
                            //this line corresponds to the termination step that was previously in separate loop
                            betaTTPlusOne[v][transRight] += beta[t][v+1][indexOfParentDoubleIndex]
                                                                    * alpha[t][v][doubleIndex]
                                                                    / probState[prevLeft];
                        }
                    }
                    
                    else if(emission == 2 && t > 0 && beta[t-1][v][indexOfParentDoubleIndex])
                    {
                        for(uint ss=0; ss<numberOfSecondarySymbols; ss++)
                        {
                            ulint transLeft = prevLeft * numberOfSecondarySymbols + ss;
                            ulint nextLeft = transLeft % numberOfStates;
                            if(probState[nextLeft] == 0) continue;
                            ulint doubleIndexLeft = prevLeft * numberOfStates + nextLeft;
                            
                            if(emissionType[SScPair[nextLeft]][SScPair[nextRight]] == 0)
                                continue;
                            
                            ulint doubleIndex = nextLeft*numberOfStates + nextRight;
                          
                            real deltaBeta = beta[t-1][v][indexOfParentDoubleIndex]
                                                * alpha[t-1][t][doubleIndexLeft]
                                                / probState[prevLeft];
                            beta[t][v][betaIndex.at(doubleIndex)] += deltaBeta;
                            
                            //this line corresponds to the termination step that was previously in separate loop
                            betaTTPlusOne[t-1][transLeft] += beta[t-1][v][indexOfParentDoubleIndex]
                                                                    * alpha[t][v][doubleIndex]
                                                                    / probState[nextRight];
                        }
                    }
                    
                    else if(emission == 3)
                    //   && t > 0 && v < windowedLength -1)
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
                                //ulint doubleIndexLongLeft = prevLeft * numberOfStates + prevRight;
                                //ulint doubleIndexLongRight = nextLeft * numberOfStates + nextRight;
                                
                                if(t > 0 && v < windowedLength - 1  // && v > t+1
                                                    && beta[t-1][v+1][indexOfParentDoubleIndex]
                                                    && emissionType[SScPair[nextLeft]][SScPair[prevRight]] ) {
                                    
                                    real deltaBeta =  beta[t-1][v+1][indexOfParentDoubleIndex]
                                                        * alpha[t-1][t][doubleIndexLeft]
                                                        * alpha[v][v+1][doubleIndexRight]
                                                        / probState[prevLeft]
                                                        / probState[nextRight];
                                    beta[t][v][betaIndex.at(doubleIndex)] += deltaBeta;
                                    
                                    //these lines corresponds to the termination step
                                    betaTTPlusOne[t-1][transLeft] += beta[t-1][v+1][indexOfParentDoubleIndex]
                                                                        // * alpha[t][v+1][doubleIndexLongRight]
                                                                        * alpha[v][v+1][doubleIndexRight]
                                                                        * alpha[t][v][doubleIndex]
                                                                        / probState[prevRight]
                                                                        / probState[nextRight];
                                    betaTTPlusOne[v][transRight] += beta[t-1][v+1][indexOfParentDoubleIndex]
                                                                        // * alpha[t-1][v][doubleIndexLongLeft]
                                                                        * alpha[t-1][t][doubleIndexLeft]
                                                                        * alpha[t][v][doubleIndex]
                                                                        / probState[nextLeft]
                                                                        / probState[prevLeft];
                                                                            
                                }
                                
                                
                                if(v == t+1 && nextLeft == prevRight){
                                    
                                    if(t > 0 && beta[t-1][v][indexOfParentDoubleIndex]){
                                        betaTTPlusOne[t][transRight] += beta[t-1][v][indexOfParentDoubleIndex]
                                                                        * alpha[t-1][t][doubleIndexLeft]
                                                                        / probState[prevLeft];
                                    
                                        //this line corresponds to the termination step
                                        betaTTPlusOne[t-1][transLeft] += beta[t-1][v][indexOfParentDoubleIndex]
                                                                            * alpha[t][v][doubleIndexRight]
                                                                            / probState[nextRight];
                                    }
                                    if(v < windowedLength-1 && beta[t][v+1][indexOfParentDoubleIndex]){
                                        betaTTPlusOne[t][transLeft] += beta[t][v+1][indexOfParentDoubleIndex]
                                                                        * alpha[v][v+1][doubleIndexRight]
                                                                        / probState[nextRight];
                                    
                                        //this line corresponds to the termination step
                                        betaTTPlusOne[v][transRight] += beta[t][v+1][indexOfParentDoubleIndex]
                                                                            * alpha[t][t+1][doubleIndexLeft]
                                                                            / probState[prevLeft];
                                    }
                                    
                                }
                                
                                
                            }
                        }
                    }
             
                    // this line prevents single probabilities to beadded more than once
                    // if(emission != 3) continue;
                }
            }
            
            if(v > t+1)
                normalize(beta[t][v]);

        }
    }
    
    for(int t = 0; t < windowedLength - 1; t++) {
        
        for(ulint trans=0; trans < numberOfFragments; trans++){
            if(probFragment[trans] == 0)
                continue;
            ulint prev = trans / numberOfSecondarySymbols;
            ulint next = trans % numberOfStates;
          
            ulint doubleIndex = prev * numberOfStates + next;
                
            if(emissionType[SScPair[prev]][SScPair[next]])
                betaTTPlusOne[t][trans] += beta[t][t+1][betaIndex.at(doubleIndex)];
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
        
     //   for(ulint i=0; i<numberOfStates; i++){
     //
     //       if(probState[i] == 0)
     //           continue;
            
     //   real deltagamma = 0.0;
        
        for(uint trans=0; trans<numberOfFragments; trans++){
                
            if(probFragment[trans] == 0)
                continue;
            ulint left = trans / numberOfSecondarySymbols;
            ulint right = trans % numberOfStates;
            if(probState[left] == 0 || probState[right] == 0)
                continue;
            ulint doubleIndex = left * numberOfStates + right;
        
           //     ulint next = trans % numberOfStates;
           //     if(probState[next] == 0)
           //         continue;
               // if(emissionType[SScPair[i]][SScPair[next]] == 0)
               //     continue;
             //   ulint doubleindex = i*numberOfStates + next;
                //ulint indexOfDoubleindex = betaINDEX.at(doubleindex);
    
            real deltagamma = betaTTPlusOne[t][trans]
                                    * alpha[t][t+1][doubleIndex]
                                    * probFragment[trans]
                                    / probState[right]; // / probState[right];
                                    /// probState[left];
            
           gamma[t][left] += deltagamma;
        //    gamma[t+1][right] += deltagamma;
                if(t == windowedLength - 2)
                    gamma[windowedLength - 1][right] += deltagamma; // / probState[i];
        
        }
        
       // gamma[t][left] += deltagamma;
      //  gamma[t+1][right] += deltagamma;
        
        //   gamma[t][i] += deltagamma; // / probState[i];
        //}
        
        //normalize(gamma[t]);
        //normalize(gamma[t+1]);
    }
    
    for(uint t=0; t<windowedLength; t++)
        normalize(gamma[t]);
    
    //normalize(gamma[windowedLength - 1]);
            
    
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
                // byte sec = indexToSequence(windowSize-1, i)[0];
                // probSecondarySymbol[a][sec] += gamma[a][i];
                byte sec = indexToSequence(windowSize-1, i)[halfWindow];
                probSecondarySymbol[a+halfWindow][sec] += gamma[a][i];
            }
            // normalize(probSecondarySymbol[a]);
            normalize(probSecondarySymbol[a+halfWindow]);
        }
        
        /*
        for(uint a=windowedLength; a < fullLength; a++){
            for(ulint i=0; i < numberOfStates; i++){
                byte sec = indexToSequence(windowSize-1, i)[a - windowedLength + 1];
                probSecondarySymbol[a][sec] += gamma[windowedLength-1][i];
            }
            normalize(probSecondarySymbol[a]);
        }
         */
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
