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

#ifndef TESTS_H_
#define TESTS_H_

#include <cpptest.h>
#include <sys/stat.h>
#include "misc.h"

class Tests : public Test::Suite {
public:
	Tests(){

		results[0] = "HHHHHHLLHHHHHHHHHHHHHHHHHLLEEEEEELLLLLLLLLLHHHHHHHHHHLLLLLLLLLLHHHHHHHLLLLLLLEEEEEEELLLLLLLHHHHHHHHHHHHHHHHHHLLLEEEE";
		results[1] = "ELLLLLLLHHHHHHHHHHHHHHHHHHHHHHHHLLLLLLLLLLLLLLLHHHHHHHHHHHHHHHHHHHHHHLHHHHHHHHHHLLLLLHHLLHHHHHHHHHHHHHHHHHHHHHHHHHLLLLLLLLLLLLLLHHHHHHHHHHHHHLLLLLLLHHHHHHHHLLEEEE";
		results[2] = "EELLLLLLLLLLLLLLLEELLLLLLLLLLLLLLLLEEEEEEELLLLLLEEEEELLLLLHHHHHHHHHHHHHLLLLHHHHHHHHHHLLLLLLLLLLEEEEELLLLLEELLLLLEEEEEEL";
		results[3] = "LLHHHHHHHLHHHHHHHHHHHLLLLLLLLHHHHHHHHLLLLLHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH";
		results[4] = "HHHLLLLHHHHHHHHHHHHHLLLEEEELLLLLLEEEEELLLLLLLLLLLLLLLLLLLLLLLHHHHHHHHHHHLLLHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHLLLLLEEHHHHHHHLLLLLLLLLLLLLLLLLLLLLLLLLLLL";

		TEST_ADD(Tests::compareWithSecondHMM_1);
		TEST_ADD(Tests::compareWithSecondHMM_2);
	}
private:

	string results[5];

	//Tests
	void compareWithSecondHMM_1();
	void compareWithSecondHMM_2();

	//Auxiliary
	bool fileExists(string strFilename) {
		struct stat stFileInfo;
		return !stat(strFilename.c_str(),&stFileInfo);
	}


};

#endif /* TESTS_H_ */
