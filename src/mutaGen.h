//
//  mutaGen.h
//  MutaGena
//
//  Created by Gennady on 9/5/25.
//  Copyright © 2025 Gennady. All rights reserved.
//

#ifndef mutaGen_h
#define mutaGen_h

#define MUTGEN_ID_SIZE 16
typedef pair< char[MUTGEN_ID_SIZE], float>  MUTGEN_FRACT;
//------------------------------------------------------
struct MUTSIG {     // contains all types of munations (96) from inp.file
    char muttype[8];        // "A[T>G]A"
    char motif[4];          // "ATA"
    char cREF;              // 'T'      ==motif[1]
    char cALT;              // 'G'
//    int rnd_cnt;
    MUTSIG() { memset(muttype,'\0', sizeof(muttype)); memset(motif,'\0', sizeof(motif)); cREF='\0'; cALT='\0'; };
    MUTSIG(char *t, char *m, char r, char a) { strncpy(muttype, t, 7); strncpy(motif, m, 3); cREF=r; cALT=a; };
};

//------------------------------------------------------------------

#define RT_TAB_SIZE 8
#define DOUBLE_ONE 0.99999
struct MUTAGEN {
    char MGname[MUTGEN_ID_SIZE];
    int amtMut;         // %
    vector < double > fract_mut;    //=  MUTSIG.size()
    vector < double > fract_rt;     //  = 0:7  (RT0 ... RT7)
    double fract_cod[2];            //  [0]=noncoding; [1]=coding
    double fract_gen[2];            //  [0]=intergenes; [1]=genes
    double fract_strand[3];         //  [0]=undefined; [1]=leading; [2]=lagging
    MUTAGEN( char *n ) { strncpy(MGname, n, 15); amtMut=0;
                        fract_cod[0]=0; fract_cod[1]=0; fract_gen[0]=0; fract_gen[1]=0;
                        fract_strand[0]=0; fract_strand[1]=0; fract_strand[2]=0;
    };

    int splitMut2Key( );
};

struct FILE_RECRD {
    int mPos;
    string Recrd;
    FILE_RECRD() { mPos=0; };
};

int makeMutTab( );
int openXtempFiles( );
int mergeXtempFiles( );

#endif /* mutaGen_h */
