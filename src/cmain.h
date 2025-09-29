
#ifndef _CMAIN_H__
#define _CMAIN_H__

#include <dirent.h>
#include "vector"
#include "string"
#include <unordered_map>

using namespace std;

#include "mutaGen.h"

//#define DEBUG_TRACE_TSLD

#define SRAND_VALUE 199
#define RANDOM_CYCLES 1000

//==========================================================
// arguments processing

#define _ARG_I 0x01     // run with <genome.bin>
#define _ARG_G 0x02     // build motif_base using ref to <genome.fa>
//#define _ARG_D 0x04
//#define _ARG_M 0x08
//==========================================================
struct PROGARGS {
    string HUGpath;     //genome
    string GENpath;   //genes
    string REPTIpath;
    string BASEpath;    //genome.bin
    string OUTdir;
    string INdir;
    FILE *fGENES;
    FILE *fREPTI;
    unsigned char argTAG;
    int amtMut;
    vector < MUTGEN_FRACT > vMGfract;
    
    PROGARGS() { argTAG='\0'; amtMut=0; };
    int procArg( int argc, char* argv[] );
    int openOutFiles ( );//const char *vcf_Fname );
    void closeOutFiles (  );
    bool isArg_I() { return ( (argTAG & _ARG_I) > 0); };
    bool isArg_G() { return ( (argTAG & _ARG_G) > 0); };
   
};
//==========================================================

void  xtrFileName ( const char  *Fpath, string &sName);
bool is_dir(const char  *path );
bool is_file(const char  *path );

int fgets_ShortRec( char *shortRec, int sizeRec, FILE *f_in );
int readCTGR ( );
int testHeader(char *hdr, char *fName);
int loadMutRanges(FILE *f_CAT, char *fName);
int loadRtRanges(FILE *f_CAT, char *fName);
int loadCodRanges(FILE *f_CAT, char *fName);
int loadGenRanges(FILE *f_CAT, char *fName);
int loadStrandRanges(FILE *f_CAT, char *fName);

void print_N_zone( );


#endif

