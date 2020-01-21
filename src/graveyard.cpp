		//			for(uint a=windowedLength; a < fullLength; a++){
		//				for(ulint partialIndex=0; partialIndex < numberOfPartials; partialIndex++){
		//
		//					uint offset = a - windowedLength + 1;
		//
		//					ulint rightSize = offset;
		//					ulint midSize = tam;
		//					ulint leftSize = windowSize - 1 - rightSize - midSize;
		//
		//					ulint numberOfRights = (ulint) pow((real)numberOfSecondarySymbols,(int)(rightSize));
		//					ulint numberOfMids = (ulint) pow((real)numberOfSecondarySymbols,(int)(midSize));
		//					ulint numberOfLefts = (ulint) pow((real)numberOfSecondarySymbols,(int)(leftSize));
		//
		//					ulint rightMultSize = numberOfMids * numberOfLefts;
		//					ulint midMultSize = numberOfLefts;
		//
		//					for(ulint r=0; r<numberOfRights; r++){
		//						ulint right = r * rightMultSize;
		//
		//						for(ulint m=0; m<numberOfMids; m++){
		//							ulint right_mid = right + m*midMultSize;
		//
		//							for(ulint l=0; l<numberOfLefts; l++){
		//								ulint fullIndex = right_mid + l;
		//								partial[a][m] += gamma[a][fullIndex];
		//
		//
		////								/*{ */
		////								vector<byte> vecR = indexToSequence(rightSize, r, numberOfSecondarySymbols);
		////								vector<byte> vecM = indexToSequence(midSize, m, numberOfSecondarySymbols);
		////								vector<byte> vecL = indexToSequence(leftSize, l, numberOfSecondarySymbols);
		////								vector<byte> vecFull = indexToSequence(windowSize-1, fullIndex, numberOfSecondarySymbols);
		////
		////								cout << offset << " % <";
		////								for(ulint z=0; z<rightSize; z++)
		////									cout << (int)vecR[z];
		////								cout << "[";
		////								for(ulint z=0; z<midSize; z++)
		////									cout << (int)vecM[z];
		////								cout << "]";
		////								for(ulint z=0; z<leftSize; z++)
		////									cout << (int)vecL[z];
		////								cout << ">   ";
		////								for(ulint z=0; z<windowSize-1; z++)
		////									cout << (int)vecFull[z];
		////								cout << endl;
		////								/*} */
		//							}
		//						}
		//					}
		//
		//				}
		//
		//				normalize(partial[a]);
		//			}//for a
