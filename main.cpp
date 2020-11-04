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

#include <string>
#include <sstream>
#include <iostream>
#include <cmath>
#include <limits>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

#include "mappings.h"
#include "fasta.h"
#include "Hmm.h"
// #include "Tests.h"

using namespace std;
namespace po = boost::program_options;

const uint PROGRESS_BAR_WIDTH = 70;
void printPBFills();
void printPBPercent(real p);

void usage();

inline bool fileExists(string filename){
	ifstream ifile(filename.c_str());
	return (bool) ifile;
}

int main(int argc, char **argv) {
	uint t0 = clock();

	//
	// 1. Initialize char <-> byte Mappings
	//
	initMappings();

// #define TESTING
	
#ifdef TESTING
	//========================================================//
	//                      Run Tests                         //
	//========================================================//
	Test::TextOutput output(Test::TextOutput::Verbose);
	Tests tests;
	return tests.run(output, false) ? EXIT_SUCCESS : EXIT_FAILURE;
	//========================================================//
#else
	//
	// 2. Parse arguments
	//
	string trainingFile;
	string evalFile;
	string resultsFile;
	int windowSize;
	string priAlphabet;
	string secAlphabet;
	int nLayers;
	bool probabilities = false;
	bool useDoubleFastaEval = false;
	string redux;
	bool stats = false;
	int seed;
	int nreplicas;
	bool noPartials = false;
	bool headAndTail = false;
	bool loglikes = false;
    
    //arguments added by AFPA
    int predGrammar;
    bool seqAlignment = false;
    bool burialConstraints = false;
    
	
	bool printTable = false;

	po::positional_options_description p;
	po::options_description desc("");
	desc.add_options()
		("training-data,d", 	po::value<string>(&trainingFile ) -> default_value(""))
		("input-eval,i", 	po::value<string>(&evalFile     ) -> default_value(""))
		("output-results,o", 	po::value<string>(&resultsFile  ) -> default_value(""))
		("window-size,w", 	po::value<int>   (&windowSize   ) -> default_value(7))
		("primary-alphabet,p",	po::value<string>(&priAlphabet  ) -> default_value("Letter"))
		("secondary-alphabet,s",po::value<string>(&secAlphabet  ) -> default_value("ck"))
		("layers,l",	 	po::value<int>   (&nLayers      ) -> default_value(2))
		("redux,r",		po::value<string>(&redux        ) -> default_value(""))
		("probabilities,b",     "")
		("stats,t",             "")
		("seed,e",		po::value<int>   (&seed         ) -> default_value(-1))
		("nreplicas,n", 	po::value<int>   (&nreplicas    ) -> default_value(50))
		("no-partials,x",	"")
		("head-and-tail,h",	"")
		("loglikes,k",		"")
        ("grammar,g", po::value<int> (&predGrammar ) -> default_value(0))
        ("burialconstraints,c",    "")
        ("seqalignment,a",         "")
		;

	p.add("training-data", 1);
	p.add("input-eval", 1);
	p.add("output-results", 1);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
	po::notify(vm);

	probabilities = vm.count("probabilities");
	stats = vm.count("stats");
	noPartials = vm.count("no-partials");
	headAndTail = vm.count("head-and-tail");
	loglikes = vm.count("loglikes");
	vector<string> errors;
    
    seqAlignment = vm.count("seqalignment");
    burialConstraints = vm.count("burialconstraints");
    

	// 2.1. Files
	if(trainingFile == "")
		errors.push_back("TRAINING_DATA_FILE not specified.");
	else if(!fileExists(trainingFile))
		errors.push_back("'" + trainingFile + "' not found.");

	if(evalFile == "")
		errors.push_back("INPUT_EVALUATION_FILE not specified.");
	else if(evalFile == "PrintTable")
		printTable = true;
	else if(!fileExists(evalFile))
		errors.push_back("'" + evalFile + "' not found.");
	else
		//Use double fasta if eval filename ends with '+'
		useDoubleFastaEval = stats || evalFile[evalFile.length()-1] == '+';


	if(resultsFile == "")
		errors.push_back("OUTPUT_RESULTS_FILE not specified.");

	// 2.2. Alphabet
	map<char,byte> priMapping;
	map<byte,char> priRevMapping;
	uint nPriSymbols = 0;
	string priString;

	if(priAlphabet == "HP"){
		priMapping = HP;		priRevMapping = REV_HP;			nPriSymbols = 3;
		priString = priAlphabet;
	}
	else if(priAlphabet == "HPN"){
		priMapping = HPN;		priRevMapping = REV_HPN;		nPriSymbols = 4;
		priString = priAlphabet;
	}
	else if(priAlphabet == "HPNAC"){
		priMapping = HPNAC;		priRevMapping = REV_HPNAC;		nPriSymbols = 6;
		priString = priAlphabet;
	}
	else if(priAlphabet == "HPNhp"){
		priMapping = HPNhp;		priRevMapping = REV_HPNhp;		nPriSymbols = 6;
		priString = priAlphabet;
	}
	else if(priAlphabet == "HPNhpn"){
		priMapping = HPNhpn;	priRevMapping = REV_HPNhpn;		nPriSymbols = 7;
		priString = priAlphabet;
	}
	else if(priAlphabet == "Letter"){
		priMapping = AMINO_ACIDS;	priRevMapping = REV_AMINO_ACIDS;	nPriSymbols = AMINO_ACIDS_N;
		priString = "None (20 letters + 2)";
	}
	else
		errors.push_back("Invalid primary alphabet: '" + priAlphabet + "'");

	map<char,byte> secMapping;
	map<byte,char> secRevMapping;
	uint nSecSymbols = 0;
	string secString;

	if(secAlphabet == "ehl" || secAlphabet == "EHL"){
		secMapping = SEC_STRUCT_EHL;
		secRevMapping = REV_SEC_STRUCT;
		nSecSymbols = SEC_STRUCT_N;
		secString = "CASP/EVA (EHL)";
	}
	else if(secAlphabet == "ck" || secAlphabet == "CK"){
		secMapping = SEC_STRUCT_CK;
		secRevMapping = REV_SEC_STRUCT;
		nSecSymbols = SEC_STRUCT_N;
		secString = "Frishman & Argos (CK)";
	}

	else if(secAlphabet == "layers"){
		for(byte i=0; i<nLayers; i++){
			secMapping[i + '0'] = (byte)i;
			secRevMapping[(byte)i] = i + '0';
		}
		secMapping['X'] = nLayers;
		secRevMapping[nLayers] = 'X';
		nSecSymbols = nLayers + 1;

		secString = int2string(nLayers) + " layers";
	}

	else if(secAlphabet == "chars"){
		for(byte i=0; i<nLayers; i++){
			secMapping[i + 'a'] = (byte)i;
			secRevMapping[(byte)i] = i + 'a';
		}
		secMapping['X'] = nLayers;
		secRevMapping[nLayers] = 'X';
		nSecSymbols = nLayers + 1;

		stringstream ss;
		ss << int2string(nLayers) << " characters (";
		for(uint i=0; i<nSecSymbols; i++)
			ss << secRevMapping[i];
		ss << ")";
		secString = ss.str();

	}

	else if(secAlphabet == "CAdir" || secAlphabet == "cadir"){
		secMapping = CAdir;
		secRevMapping = REV_CAdir;
		nSecSymbols = CAdir_N;
		secString = "2 layers with CA + direction of CB (0^v1X)";
	}

	else
		errors.push_back("Invalid secondary alphabet: '" + priAlphabet + "'");

	if(errors.size()){
		cout << "Error" << (errors.size()==1? "" : "s") << ":\n";
		for(uint i=0; i<errors.size(); i++)
			cout << "  * " << errors[i] << "\n";
		cout << "\n";
		usage();
		return 1;
	}

	//
	// 3. Run Predicition
	//	
	try{
		DoubleFasta training(trainingFile, priMapping, secMapping);
		FastaFile* eval;
		if(printTable){
			eval = NULL;
		}
		else{
			eval = useDoubleFastaEval?
				(FastaFile*) new DoubleFasta(evalFile, priMapping, secMapping)
				:
				(FastaFile*) new SingleFasta(evalFile, priMapping);
		}


		stringstream header;
		header  << "# HmmPred 1.2\n"
				<< "# Primary alphabet mapping: " << priString << "\n"
				<< "# Database:   " << trainingFile << "\n"
				<< "# Query file: " << evalFile << "\n"
				<< "# Secondary structure translation: " << secString << "\n"
				<< "# Fragment length: " << windowSize << "\n"
                << "# Grammar: "<< predGrammar << "\n";
		if(stats)
			header << "# Resampling sequences: " << nreplicas << " Replicas\n";
		else
			header << "#\n";


		uint runTimes = stats? nreplicas+1 : 1;
		//uint runTimes = stats? 46 : 1;

		ofstream results;
// 		results.precision(std::numeric_limits<real>::digits10 + 1); //maximum precision
		if(!printTable){
			results.open(resultsFile.c_str());
			results << header.str()
				<< "#\n"
				<< "\n";
		}

		Hmm hmm((seed==-1)? myclock() : seed, nreplicas );
		hmm.calculatePartials = !noPartials;
		hmm.headAndTail = headAndTail;
		hmm.loglikes = loglikes;
        
        hmm.predGrammar = predGrammar;
        hmm.seqAlignment = seqAlignment;
        hmm.burialConstraints = burialConstraints;

		//Statistics
		for(uint r=0; r<runTimes; r++){
			hmm.readTrainingData(training, windowSize, nPriSymbols, nSecSymbols, priMapping, priRevMapping, secRevMapping, redux, r!=0);
			
			if(printTable){
				
				vector<string> fne = extract_filename_extendsion(resultsFile);
				const string basename  = fne[0];
				const string extension = fne[1];					
				
				ofstream tableTransitions( (basename + "_transitions" + extension).c_str());
				ofstream tableEmissions  ( (basename + "_emissions"   + extension).c_str());
				
  				hmm.printTableTransitions(tableTransitions);
 				hmm.printTableEmissions(tableEmissions);
				
				tableTransitions.close();
				tableEmissions.close();

				
// 				ofstream resultsH((resultsFile + "H").c_str());
// 				ofstream resultsP((resultsFile + "P").c_str());
// 				ofstream resultsN((resultsFile + "N").c_str());
// 				
// 				hmm.printTableForPrimarySymbol(resultsH, 0);
// 				hmm.printTableForPrimarySymbol(resultsP, 1);
// 				hmm.printTableForPrimarySymbol(resultsN, 2);
// 				resultsH.close();
// 				resultsP.close();
// 				resultsN.close();				
/*				
				for(int i=0; i<20; i++)
					hmm.calculateEntropyForPrimarySymbol(i);*/
				
// 				hmm.calculateEntropyForPrimarySymbol(0);
// 				hmm.calculateEntropyForPrimarySymbol(1);
// 				hmm.calculateEntropyForPrimarySymbol(2);
				
// 				ofstream resultsV2((resultsFile + ".predinf").c_str());
// 				hmm.printTableV2(resultsV2);
// 				resultsV2.close();
				
				break;
			}

			if(stats){
				for(uint r=0; r<hmm.stats.size(); r++){
					hmm.stats[r].on = true;
					hmm.stats[r].hmm = &hmm;
				}
			}

			if( r==0 ){
				cout << header.str()
					 << "* Number of Fragments: " << hmm.numberOfFragments << " (" << hmm.numberOfSecondarySymbols << "^" << hmm.windowSize << ")" << endl
					 << "* Display probabilities: " << (probabilities? "true" : "false") << endl
					 << "* Calculate statistics: " << (stats? "true" : "false") << endl
					 << "* Using double fasta in eval file: " << (useDoubleFastaEval? "true" : "false") << endl
					 << "* Redux mapping: " << (redux == "" ? "none" : redux) << endl;
				if(seed != -1)
					cout
					 << "* Using random number seed " << seed << "\n";
				cout << endl;
				printPBFills();
			}

			list<Sequence*>::iterator evalIt = eval->sequences.begin();
			list<string>::iterator commIt = eval->comments.begin();

			list<Sequence*>::iterator secIt;
			if(useDoubleFastaEval)
				secIt = ((DoubleFasta*)eval)->secondarySequences.begin();

			real count = 0;
			real total = eval->sequences.size();

			if(r != 0)
				results << "\n\n# --- Replica " << r << " ---" << endl;

		    results	<< "\n# --- Predictions ---\n\n";

			while(evalIt != eval->sequences.end()){
				vector<byte> seq = (*evalIt)->intSeq;
				string comment = *commIt;

				if(stats)
					printPBPercent( (r*total + count) / (runTimes*total) );
				else
					printPBPercent(count / total);

				results << ">" << comment << "\n"
						<< vec2string(seq, priRevMapping) << "\n";

				if(probabilities){
					vector< vector<real> > probs;
					if(useDoubleFastaEval)
						hmm.predict(seq, (*secIt)->intSeq, r, &probs);
					else
						hmm.predict(seq, &probs);
					
					results << "> Secondary symbol probabilities\n";

					byte jmax = hmm.numberOfReducedSymbols? hmm.numberOfReducedSymbols : hmm.numberOfSecondarySymbols;

					for(uint i=0; i < probs.size(); i++){
						for(byte j=0; j < jmax-1; j++)
							results << probs[i][j] << "  ";
						results << probs[i][jmax-1] << "\n";
					}
					results << "\n";
				}
				else{
					vector<byte> predicted = useDoubleFastaEval?
							hmm.predict(seq, (*secIt)->intSeq, r)
							:
							hmm.predict(seq);
					string predstring = vec2string(predicted,
							hmm.numberOfReducedSymbols? hmm.rev_reducedAlphabet : secRevMapping
					);

					results << "> Local Prediction\n"
							<< predstring << "\n"
							<< "\n";
				}

				++evalIt; ++commIt; ++count;
				if(useDoubleFastaEval)
					++secIt;
			}

			if(stats)
				hmm.stats[r].print(results,hmm.calculatePartials);
		}//for r

		if(!printTable)
			printPBPercent(1);

		if(stats)
			hmm.printResampledStats(results);

		////
		if(!printTable){
			cout << endl;
			
			uint time = (clock() - t0) / CLOCKS_PER_SEC;
			uint hour = time/3600;
			time %= 3600;
			uint min = time/60;
			time %= 60;
			uint sec = time;
			cout    << "Run time: " << hour << " hours, " << min << " min, " << sec << " sec" << endl;
			results << "Run time: " << hour << " hours, " << min << " min, " << sec << " sec" << endl;
		}
		
	}
	catch(const char* msg){
		cerr << "Error:\n  * " << msg << endl;
	}

}

/**
 * Progress bar functions.
 */
void printPBFills(){
	cout << "      ";
	for(uint i=0; i<PROGRESS_BAR_WIDTH; i++)
		cout << "-";
	for(uint i=0; i<PROGRESS_BAR_WIDTH; i++)
		cout << "\b";
	cout << "\b\b\b\b\b\b" << flush;
}

real lastP = -1;
void printPBPercent(real p){
	assert(p <= 1.0);

	if(p!=1.0 && abs(p-lastP) < 0.01)
		return;

	uint fills = (uint) (p * PROGRESS_BAR_WIDTH);
	for(uint i=0; i<fills; i++)
		cout << "\b";

	cout << "\b\b\b\b\b\b";
	cout << " " << setprecision(0) << setw( 3 ) << setfill(' ') << fixed << p*100.0 << "% ";

	for(uint i=0; i<fills; i++)
		cout << "#";

	cout << flush;

	lastP = p;
#endif
}

/**
 * Program usage.
 */
void usage(){
	cout
	<< "---------------------------------------------------------------------------------------------------------\n"
	<< " HmmPred v.1.2\n"
	<< " Copyright (C) 2014, Marx Gomes van der Linden.\n"
	<< "---------------------------------------------------------------------------------------------------------\n"
	<< "\n"
	<< "Usage: HmmPred [-d] TRAINING_DATA_FILE [-i] INPUT_EVALUATION_FILE [-o] OUTPUT_RESULTS_FILE\n"
	<< "\n"
	<< "Options:\n"
	<< "        -w  --window              Fragment window size (default: 7)\n"
	<< "\n"
	<< "        -p  --primary-alphabet    Alphabet for the primary structure. Should be one of:\n"
	<< "                                  * HP\n"
	<< "                                  * HPN \n"
	<< "                                  * HPNAC\n"
	<< "                                  * HPNhp\n"
	<< "                                  * HPNhpn\n"
	<< "                                  * Letter (Default)\n"
	<< "\n"
	<< "        -s  --secondary-alphabet  Alphabet for whatever is to be predicted according to the\n"
	<< "                                  primary structure. Should be one of:\n"
	<< "                                  * ehl      EHL mapping for secondary structure\n"
	<< "                                  * ck       CK  mapping for secondary structure   (Default)\n"
	<< "                                  * layers   Burial layers (0123..., can be used with the -l option)\n"
	<< "                                  * chars    Character alphabet (abcd..., can be used with the -l option)\n"
	<< "                                  * cadir    2 burial layers with CA->CB direction (0^v1)\n"
	<< "\n"
	<< "        -l  --layers=N            Number of symbols of the secondary alphabet, if --secondary-alphabet is\n"
	<< "                                  set to 'layers' or 'chars' \n"
	<< "                                  Default: 2\n"
	<< "\n"
	<< "        -r  --redux=???...        If secondary alphabet is layers or chars, reduces output to given\n"
	<< "                                  mapping.\n"
	<< "\n"
	<< "        -b  --probabilities       Output symbols with their probabilities.\n"
	<< "\n"
	<< "        -t  --stats               Calculate statistics about the results.\n"
	<< "                                  Requires the evaluation file to be in double-fasta format.\n"
	<< "\n"
	<< "        -n  --nreplicas=N         Set the number of bootstrapping replicas when using --stats.\n"
	<< "                                  Default: 50\n"	
	<< "\n"
	<< "        -e  --seed=N              Use a fixed seed number for the random number generator.\n"
	<< "\n"
	<< "        -x  --no-partials         Do not calculate partial fragment log likelihoods.\n"
	<< "                                  Use this option to save some RAM when you don't need these data.\n"
	<< "\n"
	<< "        -h  --head-and-tail       Treat probabilities at the beginning and end of sequences separately.\n"
	<< "\n"
	<< "        -k  --loglikes            Write out log likelihoods per residue in the output file.\n"
	<< "\n"	
	<< "HmmPred was compiled on " << __DATE__ << ", " << __TIME__
#ifdef __INTEL_COMPILER
	<< " using the Intel C++ Compiler version " << 	__INTEL_COMPILER << "." << __INTEL_COMPILER_BUILD_DATE
#else
	#ifdef __GNUC__
		<< " using the GNU C++ compiler version " << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__
	#endif
#endif
#ifdef __LP64__
	<< " (64-bit)"
#else
	<< " (32-bit)"
#endif
	".\n";
}
