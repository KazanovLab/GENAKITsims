//
//  genbase.h
//  MutaGena
//
//  Created by Gennady on 9/1/25.
//  Copyright © 2025 Gennady. All rights reserved.
//

#ifndef genbase_h
#define genbase_h

#include <stdio.h>
#include "vector"
#include "string"

//#include "mapping.h"
using namespace std;

#define MOTKEY_SIZE 8
typedef pair< char, unsigned int>  INV_POS;

struct MOTIF_vPOS {
    char motKey__ [MOTKEY_SIZE];                 // non used nowhere 
    vector < INV_POS > vPosM; // < invert, posM >
    MOTIF_vPOS( ) { strcpy(motKey__, "???"); };
    MOTIF_vPOS( char *pk) { strncpy(motKey__, pk, MOTKEY_SIZE-1); };
};
bool lesser_POS ( const INV_POS &x1, const INV_POS &x2 );

struct INDEX_BAS {
    char mKey [MOTKEY_SIZE];
    int amtPos;
    long offsetF;
    INDEX_BAS() { memset(mKey, '\0', MOTKEY_SIZE); amtPos=0;  offsetF=-1L;   };
    INDEX_BAS(char *mK) { strncpy(mKey, mK, MOTKEY_SIZE); };
    INDEX_BAS( char *pK, int cnt )  { strncpy(mKey, pK, MOTKEY_SIZE); amtPos=cnt; offsetF=-1L; };
};

int buildGenBase( );
int readGenBase( );
int readData(INDEX_BAS &Motif, INV_POS *bufData, FILE *binFile);
//int readMotif (INDEX_BAS &Motif, MOTIF_vPOS &dataPos, FILE *binFile );
int readIndex ( vector <INDEX_BAS> &vIndex, FILE *binFile );

#endif /* genbase_h */
//================================================================================
//================================================================================
//     Структура базы позиций Генома :
/*
 --------------------   HEADER  ------------------------------------
 offset <индекс_хромосом>            [long]
 offset <индекс_Ключей>              [long]
 offset <индекс_Footer>              [long]
 
 --------------------   INDEX ------------------------------------
 
 ---> <индекс_хромосом>:
 число хромосом                  [int]
 {   XroID                       [XRO_ID_SIZE]       char[]
 #1_поз. в Геноме            [unsigned int]
 #посл_позиции в Геноме      [unsigned int]
 Xsize                       [unsigned int]
 } [число хромосом]
 
 ---> <индекс_Ключей>:
 keySet.size()                   [int]
 vIndexG                         [keySet.size()*sizeof(INDEX_BAS)]
 {  KEY_ID  [MOTKEY_SIZE]
    amt позиций    [int]
    offset позиций [long]
 } [число Ключей]
 
 --------------------     DATA ------------------------------------
 
 ---> <offset позиции> Ключа_1:
 vPosM                           [amt1 * sizeof(INV_POS)]
 { invert [char]
   posM   [unsigned int]
 }

 ---> <offset позиции> Ключа_2:
 vPosM                           [amt2 * sizeof(INV_POS)]
 
 ---> ........................
 
 --------------------   FOOTER  ------------------------------------
 длина имени Геном.файла            [int]
 имя Геном.файла                    []                   char[]
 
 struct {
 long  ixXRO;
 long  ixMOTIF;
 long  ixGen_ID;
 
 
 }
 
 
 */
