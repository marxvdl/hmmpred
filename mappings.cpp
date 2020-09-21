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

#include "mappings.h"

map<char,byte> SEC_STRUCT_EHL;
map<char,byte> SEC_STRUCT_CK;
map<byte,char> REV_SEC_STRUCT;

map<char,byte> AMINO_ACIDS;
map<byte,char> REV_AMINO_ACIDS;

map<char,byte> DEBUG_01;
map<byte,char> REV_DEBUG_01;
map<byte,char> REV_DEBUG_000111;
map<byte,char> REV_DEBUG_ABC;

map<byte,char> NULL_REV_MAP;

map<char,byte> CAdir;
map<byte,char> REV_CAdir;

map<char, byte> HP;
map<char, byte> HPN;
map<char, byte> HPNAC;
map<char, byte> HPNhp;
map<char, byte> HPNhpn;

map<byte,char> REV_HP;
map<byte,char> REV_HPN;
map<byte,char> REV_HPNAC;
map<byte,char> REV_HPNhp;
map<byte,char> REV_HPNhpn;

map<byte,char>& reverseOf(map<char,byte>& map){
	if(&map == &AMINO_ACIDS)
		return REV_AMINO_ACIDS;
	else if(&map == &SEC_STRUCT_EHL)
		return REV_SEC_STRUCT;

	return NULL_REV_MAP;
}

bool initMappingsAlready = false;

void initMappings(){
	if(initMappingsAlready)
		return;

	SEC_STRUCT_EHL['E'] = SEC_STRUCT_EHL['B'] = SEC_STRUCT_EHL['b'] = 0;
	SEC_STRUCT_EHL['H'] = SEC_STRUCT_EHL['G'] = SEC_STRUCT_EHL['I'] = 1;
	SEC_STRUCT_EHL['T'] = SEC_STRUCT_EHL['S'] = SEC_STRUCT_EHL['C'] =
	SEC_STRUCT_EHL[' '] = SEC_STRUCT_EHL['_'] = SEC_STRUCT_EHL['-'] =
	                                            SEC_STRUCT_EHL['L'] = 2;
	SEC_STRUCT_EHL['X'] = SEC_STRUCT_EHL['?']                       = 3;

	SEC_STRUCT_CK['E']                                              = 0;
	SEC_STRUCT_CK['H']                                              = 1;
	  SEC_STRUCT_CK['G'] = SEC_STRUCT_CK['I']
	= SEC_STRUCT_CK['B'] = SEC_STRUCT_CK['b'] = SEC_STRUCT_CK['T']
	= SEC_STRUCT_CK['S'] = SEC_STRUCT_CK['C'] = SEC_STRUCT_CK[' ']
	= SEC_STRUCT_CK['_'] = SEC_STRUCT_CK['-'] = SEC_STRUCT_CK['L']  = 2;
	SEC_STRUCT_CK['X'] = SEC_STRUCT_CK['?']                         = 3;

	REV_SEC_STRUCT[0] = 'E';
	REV_SEC_STRUCT[1] = 'H';
	REV_SEC_STRUCT[2] = 'L';
	REV_SEC_STRUCT[3] = 'X';

	DEBUG_01['0'] = 0;
	DEBUG_01['1'] = 1;
	DEBUG_01['X'] = 2;

	REV_DEBUG_01[0] = '0';
	REV_DEBUG_01[1] = '1';
	REV_DEBUG_01[2] = 'X';
	
	REV_DEBUG_ABC[0] = 'a';
	REV_DEBUG_ABC[1] = 'b';
	REV_DEBUG_ABC[2] = 'c';
	REV_DEBUG_ABC[3] = 'd';
	REV_DEBUG_ABC[4] = 'e';
	REV_DEBUG_ABC[5] = 'f';
	REV_DEBUG_ABC[6] = 'g';
	REV_DEBUG_ABC[7] = 'h';
	REV_DEBUG_ABC[8] = 'X';
	
	REV_DEBUG_000111[0] = '0';
	REV_DEBUG_000111[1] = '0';
	REV_DEBUG_000111[2] = '0';
	REV_DEBUG_000111[3] = '1';
	REV_DEBUG_000111[4] = '1';
	REV_DEBUG_000111[5] = '1';
	REV_DEBUG_000111[6] = 'X';

	AMINO_ACIDS['A'] = 0;
	AMINO_ACIDS['C'] = 1;
	AMINO_ACIDS['D'] = 2;
	AMINO_ACIDS['E'] = 3;
	AMINO_ACIDS['F'] = 4;
	AMINO_ACIDS['G'] = 5;
	AMINO_ACIDS['H'] = 6;
	AMINO_ACIDS['I'] = 7;
	AMINO_ACIDS['K'] = 8;
	AMINO_ACIDS['L'] = 9;
	AMINO_ACIDS['M'] = 10;
	AMINO_ACIDS['N'] = 11;
	AMINO_ACIDS['P'] = 12;
	AMINO_ACIDS['Q'] = 13;
	AMINO_ACIDS['R'] = 14;
	AMINO_ACIDS['S'] = 15;
	AMINO_ACIDS['T'] = 16;
	AMINO_ACIDS['V'] = 17;
	AMINO_ACIDS['W'] = 18;
	AMINO_ACIDS['Y'] = 19;
	AMINO_ACIDS['X'] = 20;

	REV_AMINO_ACIDS[0]  = 'A';
	REV_AMINO_ACIDS[1]  = 'C';
	REV_AMINO_ACIDS[2]  = 'D';
	REV_AMINO_ACIDS[3]  = 'E';
	REV_AMINO_ACIDS[4]  = 'F';
	REV_AMINO_ACIDS[5]  = 'G';
	REV_AMINO_ACIDS[6]  = 'H';
	REV_AMINO_ACIDS[7]  = 'I';
	REV_AMINO_ACIDS[8]  = 'K';
	REV_AMINO_ACIDS[9]  = 'L';
	REV_AMINO_ACIDS[10] = 'M';
	REV_AMINO_ACIDS[11] = 'N';
	REV_AMINO_ACIDS[12] = 'P';
	REV_AMINO_ACIDS[13] = 'Q';
	REV_AMINO_ACIDS[14] = 'R';
	REV_AMINO_ACIDS[15] = 'S';
	REV_AMINO_ACIDS[16] = 'T';
	REV_AMINO_ACIDS[17] = 'V';
	REV_AMINO_ACIDS[18] = 'W';
	REV_AMINO_ACIDS[19] = 'Y';
	REV_AMINO_ACIDS[20] = 'X';

	CAdir['0'] = 0;
	CAdir['^'] = 1;
	CAdir['v'] = 2;
	CAdir['1'] = 3;
	CAdir['X'] = 4;

	REV_CAdir[0] = '0';
	REV_CAdir[1] = '^';
	REV_CAdir[2] = 'v';
	REV_CAdir[3] = '1';
	REV_CAdir[4] = 'X';

	HP['A'] = HP['C'] = HP['F'] = HP['G'] = HP['I'] = HP['L'] = HP['M'] = HP['V'] = HP['W'] = HP['Y'] = 0;
	HP['D'] = HP['E'] = HP['H'] = HP['K'] = HP['N'] = HP['P'] = HP['Q'] = HP['R'] = HP['S'] = HP['T'] = 1;
	HP['X']                                                                                           = 2;

	HPN['C'] = HPN['F'] = HPN['I'] = HPN['L'] = HPN['M'] = HPN['V'] = HPN['W'] = HPN['Y'] = 0;
	HPN['D'] = HPN['E'] = HPN['K'] = HPN['N'] = HPN['P'] = HPN['Q'] = HPN['R']            = 1;
	HPN['A'] = HPN['G'] = HPN['H'] = HPN['S'] = HPN['T']                                  = 2;
	HPN['X']                                                                              = 3;

	HPNAC['C'] = HPNAC['I'] = HPNAC['L'] = HPNAC['M'] = HPNAC['V']      = 0;
	HPNAC['N'] = HPNAC['P'] = HPNAC['Q']                                = 1;
	HPNAC['A'] = HPNAC['G'] = HPNAC['H'] = HPNAC['S'] = HPNAC['T']      = 2;
	HPNAC['F'] = HPNAC['W'] = HPNAC['Y']                                = 3;
	HPNAC['D'] = HPNAC['E'] = HPNAC['K'] = HPNAC['R']                   = 4;
	HPNAC['X']                                                          = 5;

	HPNhp['F'] = HPNhp['I'] = HPNhp['L'] = HPNhp['V'] = HPNhp['W']      = 0;
	HPNhp['D'] = HPNhp['E'] = HPNhp['K']                                = 1;
	HPNhp['A'] = HPNhp['G'] = HPNhp['H'] = HPNhp['S'] = HPNhp['T']      = 2;
	HPNhp['C'] = HPNhp['M'] = HPNhp['Y']                                = 3;
	HPNhp['N'] = HPNhp['P'] = HPNhp['Q'] = HPNhp['R']                   = 4;
	HPNAC['X']                                                          = 5;

	HPNhpn['F'] = HPNhpn['I'] = HPNhpn['L'] = HPNhpn['V'] = HPNhpn['W'] = 0;
	HPNhpn['D'] = HPNhpn['E'] = HPNhpn['K']                             = 1;
	HPNhpn['A'] = HPNhpn['H']                                           = 2;
	HPNhpn['G'] = HPNhpn['S'] = HPNhpn['T']                             = 3;
	HPNhpn['C'] = HPNhpn['M'] = HPNhpn['Y']                             = 4;
	HPNhpn['N'] = HPNhpn['P'] = HPNhpn['Q'] = HPNhpn['R']               = 5;
	HPNhpn['X']                                                         = 6;

	REV_HP[0] = 'H'; REV_HP[1] = 'P'; REV_HP[2] = 'X';

	REV_HPN[0] = 'H'; REV_HPN[1] = 'P'; REV_HPN[2] = 'N'; REV_HPN[3] = 'X';

	REV_HPNAC[0] = 'H'; REV_HPNAC[1] = 'P'; REV_HPNAC[2] = 'N';
	REV_HPNAC[3] = 'A'; REV_HPNAC[4] = 'C'; REV_HPNAC[5] = 'X';

	REV_HPNhp[0] = 'H'; REV_HPNhp[1] = 'P'; REV_HPNhp[2] = 'N';
	REV_HPNhp[3] = 'h'; REV_HPNhp[4] = 'p'; REV_HPNhp[5] = 'X';

	REV_HPNhpn[0] = 'H'; REV_HPNhpn[1] = 'P'; REV_HPNhpn[2] = 'N';
	REV_HPNhpn[3] = 'h'; REV_HPNhpn[4] = 'p'; REV_HPNhpn[4] = 'n';
	REV_HPNhpn[5] = 'X';

	initMappingsAlready = true;
}
