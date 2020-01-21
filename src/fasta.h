/*
 * Copyright 2012, Marx Gomes van der Linden.
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

#ifndef FASTA_H_
#define FASTA_H_

#include <string>
#include <list>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include "misc.h"
#include "Sequence.h"

using namespace std;

/**
 * An input file with sequences in Fasta format.
 */
class FastaFile {
public:
	int size;
	list<string> comments;
	list<Sequence*> sequences;
};

/**
 * A fasta file with single pairs of comment-sequences.
 */
class SingleFasta : public FastaFile{
public:
	SingleFasta(){}

	SingleFasta(string filename, map<char,byte>& mapping){
		ifstream fin( filename.c_str() );

		string comment, seq, tmp;
		while( getline(fin,comment) ) {
			getline(fin,seq);
			getline(fin, tmp);
			try{
				comment = comment.substr(1);
			}
			catch(std::out_of_range){
				stringstream ss;
				ss << "Invalid input file \"" << filename << "\". (Should be in single fasta format)\n";
				throw ss.str().c_str();
			}
			comments.push_back(trim(comment));
			sequences.push_back(new Sequence(seq, mapping));
		}

		size = comments.size();
	}
};

/**
 * A fasta file with double pairs of comment-sequences.
 */
class DoubleFasta  : public FastaFile {
public:
	list<Sequence*> secondarySequences;

	DoubleFasta(){}

	DoubleFasta(string filename, map<char,byte>& priMapping, map<char,byte>& secMapping){
		ifstream fin( filename.c_str() );

		string comment1, priSeq, tmp, secSeq;
		while( getline(fin,comment1) ) {
			if(!getline(fin,priSeq)) break;
			if(!getline(fin, tmp)  ) break;
			if(!getline(fin,secSeq)) break;
			getline(fin, tmp);

			try{
				comment1 = comment1.substr(1);
			}
			catch(std::out_of_range){
				stringstream ss;
				ss << "Invalid input file \"" << filename << "\". (Should be in double fasta format)\n";
				throw ss.str().c_str();
			}
			comments.push_back(trim(comment1));
			sequences.push_back(new Sequence(priSeq, priMapping));
			secondarySequences.push_back(new Sequence(secSeq, secMapping));
		}

		size = comments.size();
	}
};

#endif /* FASTA_H_ */
