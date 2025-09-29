//
//  nucles.cpp
//  mClust_ins
//
//  Created by Gennady on 7/23/25.
//  Copyright © 2025 Gennady. All rights reserved.
//

#include <stdio.h>
#include "xrosoma.h"

char Nucleos[6] = "ACTGN";
int _A=0, _C=1, _T=2, _G=3, _N=4;

////////////////////////////////////////////////////////////////////////////////////

char getNuc ( int NucID )
{
    return Nucleos[NucID];
}
////////////////////////////////////////////////////////////////////////////////////

int getNucID ( const char Nuc )
{
    switch (Nuc) {
        case 'A':
            return _A;
        case 'C':
            return _C;
        case 'T':
            return _T;
        case 'G':
            return _G;
            
        default:
            break;
    }
    return _N;
}

////////////////////////////////////////////////////////////////////////////////////

char getCmpl_Nuc( const char Nuc)
{
    switch (Nuc) {
        case 'A':
            return 'T';
        case 'C':
            return 'G';
        case 'T':
            return 'A';
        case 'G':
            return 'C';
            
        default:
            break;
    }
    return 'N';
}
////////////////////////////////////////////////////////////////////////////////////

char getCmpl_Nuc( int NucID )
{
    const char Nuc = getNuc ( NucID );
    
    return getCmpl_Nuc( Nuc);
}
////////////////////////////////////////////////////////////////////////////////////

int getCmpl_NucId ( const char Nuc )
{
    switch (Nuc) {
        case 'A':
            return _T;
        case 'C':
            return _G;
        case 'T':
            return _A;
        case 'G':
            return _C;
            
        default:
            break;
    }
    return _N;
}
////////////////////////////////////////////////////////////////////////////////////

int getCmpl_NucId ( int NucID )
{
    const char Nuc = getNuc ( NucID );
    
    return getCmpl_NucId (Nuc);
}
////////////////////////////////////////////////////////////////////////////////////

void testNuc()
{                   //"ACTGN"
    int idN;
    char Nuc;
    
    printf( "Nuc\tId\tcN>N\tcN>i\tci>N\tci>i\n");
    for (int n=0; n<=4; n++ )   {
        Nuc = Nucleos[n];
        idN = getNucID(Nuc);
        
        printf( "%c\t%d\t%c\t%d\t%c\t%d\n", Nuc, idN, getCmpl_Nuc(Nuc), getCmpl_NucId(Nuc), getCmpl_Nuc(idN), getCmpl_NucId(Nuc));
        
    }
}
////////////////////////////////////////////////////////////////////////////////////



