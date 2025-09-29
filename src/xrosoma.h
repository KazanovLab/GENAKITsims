#ifndef _XROSOMA_H__
#define _XROSOMA_H__


#include "vector"
#include "string"

#include "genbase.h"
#include "mutaGen.h"
using namespace std;

/////////////////////////////////////////////////////////////////

#define CODING_TAG 0x40
#define GEN_TAG 0x20

#define LEADING_TAG 0x10
#define LAGGING_TAG 0x08
#define RT_TAG  0x07


#define SET_GEN_TAG(p)    *(p) |= GEN_TAG
#define GET_GEN_TAG(p)  (((*(p) & GEN_TAG)==0) ? 0 : 1 )

#define SET_CODING_TAG(p) *(p) |= CODING_TAG
#define GET_CODING_TAG(p) (((*(p) & CODING_TAG)==0) ? 0 : 1 )

#define SET_LEADING_TAG(p) *(p) |= LEADING_TAG
#define GET_LEADING_TAG(p) (((*(p) & LEADING_TAG)==0) ? 0 : 1 )

#define SET_LAGGING_TAG(p) *(p) |= LAGGING_TAG
#define GET_LAGGING_TAG(p) (((*(p) & LAGGING_TAG)==0) ? 0 : 1 )

#define SET_RT_VAL(p, v)  *(p) |= (RT_TAG & v)
#define GET_RT_VAL(p)    (*(p) & RT_TAG)

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

#define XRO_ID_SIZE 64

class  XROSOMA {
private:
//		vector <HUGEN> :: iterator begSearch;
public:
    static int  chrIDmode;
    static int  bufDsizeX;
    static long maxXsize;

//
//  ===== real body of XRO will be read from files =============================
//
    string XroID;
    unsigned int XstartPos;  // Absolut #First Nuc_pos across GENOME
    unsigned int XstopPos;   // Absolut #Last Nuc_pos across GENOME
    unsigned int Xsize;
    char *Xbody;
    char *Xtag;
    vector <INDEX_BAS> vIndex_X;    // copying from HUGEN.vIndexG at processing Xro
    static vector <MOTIF_vPOS> vX_Data;
    
    long f_GENstart;  // start pos in GENESfile
    long f_RTstart;   // start pos in RTfile
    
    string XtempFilePath;
    FILE *XtempFile;
    int cntOutMut;  // real count mutations at my VCF file
    
    
// testing
    int cntGcod;
    int cntGnoc;
    int cntRT[RT_TAB_SIZE];
    
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
///////////
    XROSOMA() { XroID="?"; Xsize=0L; XstartPos=0; XstopPos=0; Xbody=NULL; Xtag=NULL; XtempFile=NULL; f_GENstart=-1L; f_RTstart=-1L; cntOutMut=0;};
    XROSOMA(char *xName) { XroID=xName; Xsize=0L; XstartPos=0; XstopPos=0; Xbody=NULL; Xtag=NULL; XtempFile=NULL; f_GENstart=-1L; f_RTstart=-1L;cntOutMut=0;};
    XROSOMA(unsigned int s) { XstopPos = s; };  //for finding 
    
    int testValidDNK( int Pos, char chREF );
    int APOtest( long Pos, char chREF, char chALT );
    int setGENzone( );
    int setRTzone( );
    int motivator( );
    int writeBas_X ( );
//    int readMotif_X ( char *pKey, vPOS_MOTIF &motSet );
//---------------------------------------
//  testing
    void test_xtract( );
    int test_readBasX( );
};
//--------------------------------------------
bool lesser_XRO ( const XROSOMA &x1, const XROSOMA &x2 );

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
///////////////////////////////////////////////////////////////
//                              for 2 stage: contain counts rnd_mut for each position of HUGEN.vIndexG
struct ALT_MUT_CNTS    {    //               MOIF >A, >C, >T, >G
     unsigned int cntP[4];
    ALT_MUT_CNTS() { cntP[0]=0; cntP[1]=0; cntP[2]=0; cntP[3]=0; };
};
//-------------
class  HUGEN {
public:
    string GenName;    // without extension
    unsigned int Gsize;
    vector <INDEX_BAS> vIndexG;     // etalon for any
    vector < ALT_MUT_CNTS > vAltCnt; // size == vIndexG.size()
    
//    MOTIF_vPOS vData;               // only for 1 motif. Not used at 1-st stage ( only writeBas()  )
    INV_POS *BufDataG; 
    unsigned int sizeBufDataG;
    
    FILE *GenBasFile;
    int writeGenBase ( );
    int readIndXro ( vector <XROSOMA> &vecX );
    HUGEN() { Gsize=0; BufDataG=NULL; sizeBufDataG=0; };
    
    void initALT_MUT_CNTS ( );
    void addALT_MUT_CNTS(int iKey, char cALT);
    int rndMutation( const char *mgName );
//---------------------------------------
//  testing
    int testRNDmut( );
};
///////////////////////////////////////////////////////////////

int LoadXromoSet( const char *fPath );
int loadGENEs (  );
int loadRT (  );
int findXroByID( char *xID, int say=1 );
XROSOMA *findXroByPOS( INV_POS posG, unsigned int &posX );

///////////////////////////////////////////////

//char * skipN( char *pXro );
int setNucleID(char Nuc);
char getNuc ( int NucID );
int getNucID(const char Nuc);
char getCmpl_Nuc( int NucID );
char getCmpl_Nuc( const char Nuc);
int getCmpl_NucId(const char Nuc);
int getCmpl_NucId (int NucID);
void testNuc();

int selectXid ( char *buff );       // see it at 'mClust' prog

#endif











