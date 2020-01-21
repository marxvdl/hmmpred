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

#ifndef MAPPINGS_H_
#define MAPPINGS_H_

#include <map>
#include <vector>

#include "misc.h"

using namespace std;

const int AMINO_ACIDS_N = 21;
const int SEC_STRUCT_N = 4;

extern map<char,byte> SEC_STRUCT_EHL;
extern map<char,byte> SEC_STRUCT_CK;
extern map<byte,char> REV_SEC_STRUCT;

extern map<char,byte> DEBUG_01;
extern map<byte,char> REV_DEBUG_01;
extern map<byte,char> REV_DEBUG_000111;
extern map<byte,char> REV_DEBUG_ABC;

extern map<char,byte> AMINO_ACIDS;
extern map<byte,char> REV_AMINO_ACIDS;

extern map<byte,char> NULL_REV_MAP;

extern map<char,byte> CAdir;
extern map<byte,char> REV_CAdir;
const int CAdir_N = 5;

extern vector< map<char,byte> > BURIALS;
extern vector< map<byte,char> > REV_BURIALS;

extern map<char, byte> HP;
extern map<char, byte> HPN;
extern map<char, byte> HPNAC;
extern map<char, byte> HPNhp;
extern map<char, byte> HPNhpn;

extern map<byte,char> REV_HP;
extern map<byte,char> REV_HPN;
extern map<byte,char> REV_HPNAC;
extern map<byte,char> REV_HPNhp;
extern map<byte,char> REV_HPNhpn;

extern map<byte,char> REV_ABC;

extern map<string, uint> ALPHABET_SIZE;

extern bool initMappingsAlready;

map<byte,char>& reverseOf(map<char,byte>& map);

void initMappings();

#endif /* MAPPINGS_H_ */
