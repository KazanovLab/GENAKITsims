


#include <string.h>
#include <time.h>
#include <ctype.h>
#include <algorithm>
#include <functional>      // For greater<int>( )

#include "cmain.h"
#include "xrosoma.h"

extern PROGARGS ArgKit;
HUGEN GENOM;

vector < XROSOMA > vecDNK;
int XROSOMA::chrIDmode =-1;
int XROSOMA::bufDsizeX = 0;
long XROSOMA:: maxXsize = 0;

extern FILE *Ftrace;

///////////////////////////////////////////////////////

bool lesser_XRO ( const XROSOMA &x1, const XROSOMA &x2 )
{
    return x1.XstopPos < x2.XstopPos;
}
///////////////////////////////////////////////////////

void xtrctXID( char *pBuff)
{
    char xID[XRO_ID_SIZE];
    char *pX = xID ;
    char *pB = pBuff ;
    
    pB++;   // >
    while ( *pB && *pB==' ') pB++;
    strncpy(xID, pB, XRO_ID_SIZE-1);
    
    while ( *pX && *pX > ' ') pX++;
    *pX = '\0';
    strcpy(pBuff,  xID);
    
    return;
}
///////////////////////////////////////////////////////

int LoadXromoSet( const char *fPath )
{
    FILE *fHuGen=NULL;
    char fBuff[1024];
    char *pb;
    int skipLoad=0;
    string sGen;
//    char ValXro[128] = { ">chr20&>chr21&>chr22& "};    //short variant
    printf("\nLoading genome from '%s'.....\n", fPath);
    
    if ( ( fHuGen=fopen(fPath, "r"))==NULL )    {
        printf("Cannt OPEN reference genome: '%s'\n",fPath);
        return -1;
    }
    fprintf(Ftrace,"Loading genome from '%s'.....\n", fPath);
    
    int cntXro = 0;
    XROSOMA *pCurXro = NULL;

    while ( 1 )    {
        if( fgets(fBuff, sizeof(fBuff)-1, fHuGen) == NULL )
            break;
        if ( (pb = strchr(fBuff, '\r')) )   //for Windows
            *pb = '\0';
        else
            if ( (pb = strchr(fBuff, '\n')) )
                *pb = '\0';
        if ( fBuff[0] != '>' )    {
            if ( skipLoad )
                continue;
            pb = fBuff;
            while ( *pb )   {
                *pb = toupper(*pb);
                pb++;
            }
            sGen += fBuff;
            continue;
        }
// new XRO
        if ( pCurXro  )  {
            pCurXro->Xsize = (unsigned int) sGen.size() ;
            pCurXro->Xbody =  new char[pCurXro->Xsize+1];
            strcpy(pCurXro->Xbody, sGen.c_str() );
            pCurXro->Xtag =  new char[pCurXro->Xsize+1];
            memset(pCurXro->Xtag, '\0', pCurXro->Xsize+1 );
            if ( pCurXro->Xsize > pCurXro->maxXsize )
                pCurXro->maxXsize = pCurXro->Xsize;
            cntXro++;
            printf("%s  ", pCurXro->XroID.c_str() );
            fprintf(Ftrace, "%s  ", pCurXro->XroID.c_str() );
//            fprintf("%s=[%ld]\tSum=%u\n", pCurXro->XroID.c_str(), pCurXro->Xsize, cntNucALL );
            pCurXro = NULL;
        }
        
//        char tstXro[32];                              // short variant
//        sprintf(tstXro,"%s&",fBuff);                  // short variant
//        if ( ! strstr(ValXro, tstXro) ) {             // short variant
        if ( strlen(fBuff) > 6 )    {
            if ( ! skipLoad )
                fprintf(Ftrace, "\n" );
            fprintf(Ftrace, "%s SKIPPED\n", fBuff);
            skipLoad = 1;
            continue;
        }
        else
            skipLoad = 0;
            
        
        xtrctXID( fBuff );
        if ( findXroByID(fBuff, 0) >= 0 ) {
            printf("\nRedefinition chroID=%s", fBuff );
            fclose (fHuGen );
            return -1;
        }
        
//        if ( (vecDNK.size() == 0) || ( vecDNK.back().chrNum < Xro_ID )  ) {
//            n = (int)vecDNK.size() - 1;
//        } else {
//            for ( n=0; n<vecDNK.size(); n++ )
//                if ( vecDNK[n].chrNum > Xro_ID )
//                    break;
//            vecDNK.insert(vecDNK.begin()+n, XROSOMA(Xro_ID));
//        }
        
        vecDNK.push_back(XROSOMA(fBuff));
        pCurXro = &vecDNK[vecDNK.size()-1];
        if ( pCurXro->chrIDmode < 0 )      // first
            pCurXro->chrIDmode =  (strstr(fBuff, "chr") == fBuff) ? 1 : 0;

        sGen.clear();
    }
//
    if ( pCurXro  )  {
        pCurXro->Xsize = (unsigned int) sGen.size() ;
        pCurXro->Xbody =  new char[pCurXro->Xsize+1];//+3];
        strcpy(pCurXro->Xbody, sGen.c_str() );
        pCurXro->Xtag =  new char[pCurXro->Xsize+1];
        memset(pCurXro->Xtag, '\0', pCurXro->Xsize+1 );
        if ( pCurXro->Xsize > pCurXro->maxXsize )
            pCurXro->maxXsize = pCurXro->Xsize;
        cntXro++;
        printf("%s\n", pCurXro->XroID.c_str() );
        fprintf(Ftrace, "%s\n", pCurXro->XroID.c_str());
//        fprintf(Ftrace, "%s=[%ld]\tSum=%u\n", pCurXro->XroID.c_str(), pCurXro->Xsize, cntNucALL );
    }
    
    unsigned int cntNucGen = 0;
    for ( int n=0; n<vecDNK.size(); n++ )   {
        vecDNK[n].XstartPos = cntNucGen; 
        cntNucGen += vecDNK[n].Xsize;
        vecDNK[n].XstopPos = cntNucGen - 1;      // start from 0
    }
    printf("=== Loaded %d chromos. cntNuc=%u\n", cntXro, cntNucGen );
    fprintf(Ftrace,"=== Loaded %d chromos. cntNuc=%u\n", cntXro, cntNucGen );
    
    fclose (fHuGen );
    
    return cntXro;
}
///////////////////////////////////////////////////////////////////////////////////

int loadGENEs (  )
{
    char fBuff[1024];
    char typRec[512];
    int startPos, stopPos;
    char curXroId[XRO_ID_SIZE] = "\0\0";
    char coding;
    char readXid[XRO_ID_SIZE];
    
    ArgKit.fGENES = fopen(ArgKit.GENpath.c_str(), "r");
    
    printf("\nLoading indGENEs:\t... \n" );
    memset(fBuff, '\0', sizeof(fBuff));

    XROSOMA *pCurXro = NULL;
    int cntGEN=0;
    while ( 1 ) {
        if ( ! fgets_ShortRec( fBuff, sizeof(fBuff)-2, ArgKit.fGENES) )
            break;
        if ( fBuff[0] == '#' )
            continue;
        if ( sscanf(fBuff, "%s\t%*s\t%s\t%d\t%d\t%*s\t%c", readXid, typRec, &startPos, &stopPos, &coding) != 5 ) {
            printf("loadGENEs():: Err.sscanf('%s')!=5 \n", fBuff );
            return -1;
        }
        if ( strcmp(typRec, "gene") )
            continue;
        if ( coding == '+' || coding == '-' )   {}
        else {
            printf("setGENzone():: CodingID=%c '%s' \n", coding, fBuff);
            return -1;
        }
        
        if ( ! pCurXro || strcmp(readXid, curXroId) != 0 )  {
            int iXr = findXroByID(readXid, 0);
            if ( iXr < 0 )
                continue;
            if ( pCurXro )
                printf("\t%s: +cnt=%d -cnt=%d\n", pCurXro->XroID.c_str(), pCurXro->cntGcod, pCurXro->cntGnoc );
            pCurXro = &vecDNK[iXr];
            pCurXro->cntGcod =0;
            pCurXro->cntGnoc =0;
            strcpy(curXroId, readXid);
        }
        cntGEN++;
        char *pT = pCurXro->Xtag + startPos-1;
        for ( int n=startPos-1; n<=stopPos-1; n++, pT++ )   {
            if ( GET_GEN_TAG(pT) )      // cross GENEs
                continue;
            SET_GEN_TAG(pT);
            if (coding == '+')
                SET_CODING_TAG(pT); 
            if ( coding == '+' )
                pCurXro->cntGcod++;
            else
                pCurXro->cntGnoc++;
        }
    }
    if ( pCurXro )
        printf("\t%s: +cnt=%d -cnt=%d\n", pCurXro->XroID.c_str(), pCurXro->cntGcod, pCurXro->cntGnoc );
    
    fclose(ArgKit.fGENES);
    printf("\t=====cntGENEs = %d\n",  cntGEN);
    
    return cntGEN;
}
///////////////////////////////////////////////////////////////////////////////////

int loadRT (  )
// RT && STRAND (Leading Lagging)
{
    char fBuff[1024];
    int startPos, stopPos;
    float rV;
    char strand[16];
    int leading, lagging;
    char curXroId[XRO_ID_SIZE] = "\0\0";
    unsigned char rtN;
    char readXid[XRO_ID_SIZE];
    struct RT_ID    {
        int idRT;
        float minV;
        float maxV;
    } RTzone[RT_TAB_SIZE-1] = {
        { 1,    -100,       -1.18464 },
        { 2,    -1.18464,   -0.606658 },
        { 3,    -0.606658,  -0.152577 },
        { 4,    -0.152577,  0.259648 },
        { 5,    0.259648,   0.641621 },
        { 6,    0.641621,   1.09141 },
        { 7,    1.09141,    100 }
    };
    
    ArgKit.fREPTI = fopen(ArgKit.REPTIpath.c_str(), "r");
    
    printf("\nLoading REP_TIME:\t...\n" );
    memset(fBuff, '\0', sizeof(fBuff));
    
    XROSOMA *pCurXro = NULL;
    int cntRT=0;
    
    fgets_ShortRec( fBuff, sizeof(fBuff)-2, ArgKit.fGENES);  // header
    while ( 1 ) {
        if ( ! fgets_ShortRec( fBuff, sizeof(fBuff)-2, ArgKit.fGENES) )
            break;
        if ( fBuff[0] == '#' )
            continue;
        if ( sscanf(fBuff, "%s\t%d\t%d\t%*s\t%f\t%s\n", readXid, &startPos, &stopPos, &rV, strand) != 5 ) {
            printf("loadRT():: Err.sscanf('%s')!=6 \n", fBuff );
            return -1;
        }
        
        if ( ! pCurXro || strcmp(readXid, curXroId) != 0 )  {
            int iXr = findXroByID(readXid,0);
            if ( iXr < 0 )
                continue;
            if ( pCurXro )
                printf("\t%s: RT1=%d RT2=%d RT3=%d RT4=%d RT5=%d RT6=%d RT7=%d\n", pCurXro->XroID.c_str(),
                       pCurXro->cntRT[1], pCurXro->cntRT[2], pCurXro->cntRT[3], pCurXro->cntRT[4],
                       pCurXro->cntRT[5], pCurXro->cntRT[6], pCurXro->cntRT[7] );
            pCurXro = &vecDNK[iXr];
            for ( int t=0; t<RT_TAB_SIZE; t++ )
                pCurXro->cntRT[t] = 0;
            strcpy(curXroId, readXid);
        }
        rtN = RT_TAB_SIZE;
        for ( int n=0; n<RT_TAB_SIZE; n++)
            if ( rV >= RTzone[n].minV && rV <= RTzone[n].maxV ) {
                rtN = RTzone[n].idRT;
                break;
            }
        if ( rtN >= (RT_TAB_SIZE) )  {
            printf("loadRT(%s):: Err.read RTvalue=%f\n", pCurXro->XroID.c_str(), rV );
            return -1;
        }
        leading = ! strcmp(strand, "leading");
        lagging = ! strcmp(strand, "lagging");
        cntRT++;
        char *pT = pCurXro->Xtag + startPos-1;

        for ( int n=startPos-1; n<stopPos-1; n++, pT++ )    {
            SET_RT_VAL(pT, rtN);
            if ( leading )
                SET_LEADING_TAG(pT);
            if ( lagging )
                SET_LAGGING_TAG(pT);
            
            pCurXro->cntRT[rtN] +=1;
        }

    }
    if ( pCurXro )
        printf("\t%s: RT1=%d RT2=%d RT3=%d RT4=%d RT5=%d RT6=%d RT7=%d\n", pCurXro->XroID.c_str(),
               pCurXro->cntRT[1], pCurXro->cntRT[2], pCurXro->cntRT[3], pCurXro->cntRT[4],
               pCurXro->cntRT[5], pCurXro->cntRT[6], pCurXro->cntRT[7] );
    
    fclose(ArgKit.fREPTI);
    printf("\t=====cntRTsegm cnt=%d\n",  cntRT);
    
    return cntRT;
}
///////////////////////////////////////////////////////////////////////////////////

int XROSOMA::APOtest( long Pos, char chREF, char chALT )
{
    // APOBEC:
    //    a) TCx ---> 'C' > ['G' | 'T' ]
    //    b) xGA ---> 'G' > ['C' | 'A' ]
    // Returns: 1=APOBEC; 2=APOBEC in ZONE; 0= any other
    //
    
    char *pB = Xbody+Pos-1;
    int retC=0;
    
    if (  *pB != chREF ) {
        printf( "\n !!! NUCLEO_mismatch:: CHR_%s[%ld]='%c';  MUTATION: '%c' > '%c'\n",
               XroID.c_str(), Pos, *pB, chREF, chALT);
        return 0;
    }
    
    switch ( chREF ) {
        case 'C':
            if ( ! (chALT == 'G' || chALT == 'T') )
                return 0;
            if ( *(pB-1) != 'T' )
                return 0;
            retC = 1;
            break;
        case 'G':
            if ( ! (chALT == 'C' || chALT == 'A') )
                return 0;
            if ( *(pB+1) != 'A' )
                return 0;
            retC = 1;
            break;
        default:
            return 0;
    }
    
    return retC;
}
////////////////////////////////////////////////////////////////////////////

int XROSOMA::testValidDNK( int Pos, char chREF ) 
{
    //    char *pB = Xbody+Pos-1;
    
    if (  *(Xbody+Pos-1) != chREF ) {
        printf( "\nReference base mismatch at %s:%d. Reference genome has '%c', but VCF REF = '%c'.\n",
               XroID.c_str(), Pos,*(Xbody+Pos-1), chREF);
        return 0;
    }
    return 1;
}
////////////////////////////////////////////////////////////////////////////////////

int findXro( char *xID )    {
    for ( int indX=0; indX<vecDNK.size(); indX++ )  {
        if ( vecDNK[indX].XroID == xID )
            return indX;
    }
    return -1;
}
////////////////////////////////////////////////////////////////////////////////////

int  findXroByID( char *xID, int say )
{
    int indX;
    int stat;
    char my_xID[XRO_ID_SIZE];
    
    if  ( strstr(xID,"chr")==xID )
        stat = ( vecDNK[0].chrIDmode ) ? 1 : 2;
    else
        stat = ( vecDNK[0].chrIDmode ) ? 3 : 4;
//     vcf     genome
// 1=(chrID) : (chrID)
// 2=(chrID) : (   ID)
// 3=(   ID) : (chrID)
// 1=(   ID) : (   ID)
    
    strncpy(my_xID, xID, XRO_ID_SIZE-1);
    switch (stat) {
            
        case 2:         // (chrID) : (   ID)
            if ( (indX = findXro((xID+3)) ) >= 0 )    // (chrID->ID) : ( ID )
                break;
            
            if ( strcmp(xID, "chrM")==0 )           // (chrM->MT) : ( ID )
                strcpy(my_xID, "MT");
            if ( (indX = findXro(my_xID)) >= 0 )
                break;
            
            if ( (indX = findXro(xID)) >= 0 )      // (chrID) : ( ID )   for finding in patches
                break;
            
            indX =-1;
            break;
            
        case 3:         // (   ID) : (chrID)
            strcpy(my_xID,"chr");
            strcat(my_xID, xID);
            if ( (indX = findXro(my_xID)) >= 0 )      // (ID->chrID) : (chrID)
                break;
            
            if ( strcmp(xID, "MT")==0 )
                strcpy(my_xID, "chrM");
            if ( (indX = findXro(my_xID)) >= 0 )      // (MT->chrM) : ( chrID )
                break;
            
            if ( (indX = findXro(xID)) >= 0 )      // (ID) : (chrID)   for finding in patches
                break;
            
            indX =-1;
            break;

        default:        // case1:: (chrID) : (chrID)
                        // case4:: (   ID) : (   ID)
            indX = findXro(xID);
            break;

    }
    
    if (say && indX < 0 )
        printf("ERR: XroID='%s' NOT FOUND\n", xID);
    
    return indX;
}
////////////////////////////////////////////////////////////////////////////////////

XROSOMA *findXroByPOS( INV_POS posG, unsigned int &posX )
{
    vector < XROSOMA > ::iterator  iterXr;
    int indXro;
    XROSOMA *pXro;
    
    iterXr = lower_bound( vecDNK.begin(), vecDNK.end(), posG.second, lesser_XRO);
    if ( iterXr == vecDNK.end() ) {
        printf("NOT found Xromo for Pos.%d\n", posG.second );
        return NULL;
    }
    indXro = (int)(iterXr - vecDNK.begin() );
    pXro = &vecDNK[indXro];
    posX = posG.second - pXro->XstartPos;
    
    return pXro;
}
////////////////////////////////////////////////////////////////////////////////////

void HUGEN:: initALT_MUT_CNTS ( )
{
     if ( vAltCnt.empty() )  {
         vAltCnt.reserve( vIndexG.size() );
         for ( int n=0; n<vIndexG.size(); n++ )
             vAltCnt.push_back(ALT_MUT_CNTS());
         return;
     }
    for ( int n=0; n<vAltCnt.size(); n++ )
        for ( int c=0; c<4; c++ )
            vAltCnt[n].cntP[c] = 0;
 return;
 }
 //----------------------------------------------
 
 void HUGEN:: addALT_MUT_CNTS(int iKey, char cALT)
 {
     int iNuc = getNucID ( cALT );
     vAltCnt[iKey].cntP[iNuc] += 1;
     
     return;
 }
//////////////////////////////////////////////////////////////////////////


