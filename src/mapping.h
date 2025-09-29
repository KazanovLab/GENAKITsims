//
//  mapping.h
//  MutaGena
//
//  Created by Gennady on 8/9/25.
//  Copyright © 2025 Gennady. All rights reserved.
//

#ifndef mapping_h
#define mapping_h

#include "vector"
//#include "string"
#include <unordered_map>

#include "xrosoma.h"
#include "genbase.h"

using namespace std;

//#define MOTKEY_SIZE 8
#define KEY_SET_RESERV 2305
typedef unordered_map<string, int> KEY_MAP;

#define _iKEY_GENE      3
#define _iKEY_CODING    4
#define _iKEY_RT        5
#define _iKEY_STRAND    6


void setKEY_GENES(char *key, int g );
void setKEY_CODING(char *key, int c );
void setKEY_RT(char *key, int t );
void setKEY_STRAND(char *key, const char *pTag, int invert);

void make_keySet();
void reservIndexSpace();
int checkKeySet(vector <INDEX_BAS> &vIndexG);

void rewindMotif(const char *pKey, char *pCmpl);
int xtract_Motif(const char *pXr, char *key);
int formKEY (const char *motif, char *pTag, char *Key );     // returns invert
//void setKEYbyTag (char *key, char *pTag, unsigned char invertM);     // returns invert

#endif /* mapping_h */
