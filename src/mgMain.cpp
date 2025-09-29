//
//  muMain.cpp
//  MutaGena
//
//  Created by Gennady on 8/8/25.
//  Copyright © 2025 Gennady. All rights reserved.
//

#include <stdio.h>
#include <time.h>
#include <dirent.h>
#include <sys/stat.h>

#include "cmain.h"
#include "xrosoma.h"
#include "mapping.h"

using namespace std;

PROGARGS ArgKit;
extern HUGEN GENOM;
extern vector < XROSOMA > vecDNK;
extern vector < MUTAGEN > MutaGen;
vector <MOTIF_vPOS> XROSOMA:: vX_Data;

FILE *Ftrace=NULL;

void printAllInd( );
int test_readBasGen( );
void testTagGEN_RT ( );
void testFileBin( );
//void compareKeys( ); //testing
/////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
    char Fpath[1024];
    
    double duration;
    clock_t startMain = clock();
    
    if ( ArgKit.procArg( argc, argv) < 0 )
        return -1;
    
    sprintf(Fpath, "%smutrace.txt", ArgKit.OUTdir.c_str() );
    Ftrace = fopen(Fpath, "w");
    
    make_keySet();
    
    if ( ArgKit.isArg_G() ) {
        if ( buildGenBase( ) < 0 )
            return -1;
        return 0;
    }
//--------------===========================================

    if ( makeMutTab( ) < 0 )
        return -1;
    if ( readGenBase( ) < 0 )
        return  -1;
    if ( readCTGR ( ) < 0 )
        return -1;
    
    if ( openXtempFiles()<0 )
        return -1;
    
//    GENOM.testRNDmut( );
//    return 0;
    
    for ( int n=0; n<MutaGen.size(); n++ )  {
        clock_t startMG = clock();
        printf("Processing Mutagen %s (muts=%d).... \n", MutaGen[n].MGname, MutaGen[n].amtMut );
        int rc = MutaGen[n].splitMut2Key();
        if (  rc < 0 )
            return -1;
        int rc1 = GENOM.rndMutation( MutaGen[n].MGname);
        if ( rc1 < 0 )
            return -1;
        
        clock_t stoptMG = clock();
        duration = (double)(stoptMG - startMG) / CLOCKS_PER_SEC;
        printf( "---> (%d,%d) finish dT==%5.2f\n", rc, rc1, duration );
    }
    
    mergeXtempFiles( );
    
/*
// =================================================
   for ( int indX=0; indX<vecDNK.size(); indX++ ) {
        vecDNK[indX].motivator( );
//       vecDNK[indX].test_xtract( );
        if (  vecDNK[indX].writeBas_X() < 0 )
            return -1;
//        vecDNK[indX].test_readBasX( );
    }
// ====================================================
    printAllInd( );
    
    if ( GENOM.writeGenBase() < 0 )
        return -1;
//    test_readBasGen( );
    
//  delete temp_base  !!!!!!!!!!!!!!!!!
    for ( int nXr=0; nXr<vecDNK.size(); nXr++ ) {
        if ( vecDNK[nXr].XtempFile )
            fclose(vecDNK[nXr].XtempFile);
        remove(vecDNK[nXr].XtempFilePath.c_str());
    }
  */
    
    clock_t stopMain = clock();
    duration = (double)(stopMain - startMain) / CLOCKS_PER_SEC;
    printf( "\nTotal time  dT==%5.2f\n", duration );
    
    return 0;
}
/////////////////////////////////////////////////////////////////////////

int buildGenBase( )
{
    xtrFileName ( ArgKit.HUGpath.c_str(), GENOM.GenName);
    if ( LoadXromoSet( ArgKit.HUGpath.c_str() ) <= 0 )
        return -1;
    if ( loadGENEs(  ) < 0 )
        return -1;
    if ( loadRT (  ) < 0 )
        return -1;
    //      testTagGEN_RT( );  return -1;
    
    reservIndexSpace();
    // =================================================
    for ( int indX=0; indX<vecDNK.size(); indX++ ) {
        vecDNK[indX].motivator( );
//       vecDNK[indX].test_xtract( );
        if (  vecDNK[indX].writeBas_X() < 0 )
            return -1;
//        vecDNK[indX].test_readBasX( );
    }
    GENOM.sizeBufDataG = vecDNK[0].bufDsizeX;
    GENOM.BufDataG =new INV_POS [GENOM.sizeBufDataG];
// ====================================================
//    printAllInd( );
    
    if ( GENOM.writeGenBase() < 0 )
        return -1;
    
//    test_readBasGen( );
    
/*  delete temp_base  !!!!!!!!!!!!!!!!!
    for ( int nXr=0; nXr<vecDNK.size(); nXr++ ) {
        if ( vecDNK[nXr].XtempFile )
            fclose(vecDNK[nXr].XtempFile);
        remove(vecDNK[nXr].XtempFilePath.c_str());
    }
*/
    return 0;
}
/////////////////////////////////////////////////////////////////////////

int PROGARGS:: procArg( int argc, char* argv[] )
{
    
    if ( argc < 1 ) {
        printf ("\nUsage: MutaGen -o <output_directory> -g <reference_genome> \n\n");
        printf ("Options:\n");
        printf ("  -o <output_directory>\tDirectory to save output files (required).\n");
        printf ("  -g <reference_genome>\tPath to the reference genome (required).\n");
        printf ("  -b build motif base for Genome.\n");
        return -1;
    }
    
    char what[16];
    char parList[] = "-b@0 -g@1 -gen@2 -rt@3 -o@4 -i@5 -mut@6 -ap@7 -uv@8 -sm@9 -cl@10 ";
    int parmN;
    MUTGEN_FRACT m_frakt;
    
    int nP=1;
    while ( nP+1 < argc ) {
        snprintf(what, sizeof(what)-4, "%s@", argv[nP] );
        char *pL = strstr(parList, what);
        if ( ! pL ) {
            parmN = -1;
        }
        else {
            while ( *pL != '@' ) pL++;
            pL++;
            parmN = atoi(pL);
        }
        switch ( parmN )    {
            case 0:     // -i BASEpath :   processing with using BASE.bin
                argTAG |= _ARG_I;
                if ( argv[nP+1][0] == '-' ) // skipped <output directory>
                    break;
                nP++;
                BASEpath = argv[nP];
                if ( ! is_file(argv[nP]) )    {
                    printf ("<reference_genome_base.bin> does not exist: '%s'\n", BASEpath.c_str());
                    return -1;
                }
                break;
            case 1:     // -g build BASE using  <gen_name>.fa
                argTAG |= _ARG_G;
                if ( argv[nP+1][0] == '-' ) // skipped <reference_genome>
                    break;
                nP++;
                HUGpath = argv[nP];
                if ( ! is_file(argv[nP]) )    {
                    printf ("<reference_genome> does not exist: '%s'\n", HUGpath.c_str());
                    return -1;
                }
                break;

            case 2:     //genas
                if ( argv[nP+1][0] == '-' ) // skipped <output directory>
                    break;
                nP++;
                GENpath = argv[nP];
                if ( ! is_file(argv[nP]) )    {
                    printf ("<reference_genes> does not exist: '%s'\n", GENpath.c_str());
                    return -1;
                }
                break;
            case 3:     //RT
                if ( argv[nP+1][0] == '-' ) // skipped <output directory>
                    break;
                nP++;
                REPTIpath = argv[nP];
                if ( ! is_file(argv[nP]) )    {
                    printf ("<reference_RT> does not exist: '%s'\n", REPTIpath.c_str());
                    return -1;
                }
                break;
            case 4:          // -o <output directory>
                if ( argv[nP+1][0] == '-' ) // skipped <output directory>
                    break;
                nP++;
                OUTdir = argv[nP];
                if (OUTdir.back() != '/' )
                    OUTdir += '/';
                if ( ! is_dir(argv[nP]) )    {
                    printf ("<output directory> does not exist: '%s'\n", OUTdir.c_str());
                    return -1;
                }
                break;
            case 5:     // -i dir for input files
                if ( argv[nP+1][0] == '-' ) // skipped <output directory>
                    break;
                nP++;
                INdir = argv[nP];
                if (INdir.back() != '/' )
                    INdir += '/';
                if ( ! is_dir(argv[nP]) )    {
                    printf ("<output directory> does not exist: '%s'\n", INdir.c_str());
                    return -1;
                }
                break;
            case 6:     // -mut  amount mutation
                if ( argv[nP+1][0] == '-' ) { // Mutations
                    printf ("'%s' missed Mutations value\n", argv[nP]);
                    break;
                }
                nP++;
                amtMut = atoi(argv[nP]);
                break;
            case 7:     // -ap APOBEC % of -mut
            case 8:     // -uv ultraviolet % of -mut
            case 9:     // -sm SMOKING % of -mut
            case 10:     // -cl CLOCK % of -mut
                if ( argv[nP+1][0] == '-' ) { // fraction
                    printf ("'%s' missed fraction value\n", argv[nP]);
                    return -1;
                }
                strncpy(m_frakt.first, argv[nP], MUTGEN_ID_SIZE-2);
                m_frakt.second = (float)(atof(argv[nP+1]) / 100);
                if ( m_frakt.second <= 0 )  { // fraction
                    printf ("'%s' inv. fraction value\n", argv[nP]);
                    return -1;
                }
                vMGfract.push_back(m_frakt);
                nP++;
                break;
                
            default:    //-mut@6 -ap@7 -uv@8 -sm@9 -cl@10
                printf ("Invalid option: '%s'\n", argv[nP]);
                break;
        }
        nP++;
    }
    
    if ( isArg_I() && isArg_G() )   {
        printf ("May be only -i  or -g \n");
        return -1;
    }
    
    if ( isArg_G() ) {
        if  ( HUGpath.empty() ) {
            printf ("Not defined:  -g <reference_genome> \n");
            return -1;
        }
        if ( GENpath.empty() )  {
            printf ("Not defined:  -gen <reference_genes> \n");
            return -1;
        }
        if ( REPTIpath.empty() )  {
            printf ("Not defined:  -rt <reference_RT> \n");
            return -1;
        }
    }
    if ( isArg_I() && BASEpath.empty() )    {
        printf ("Not defined:  -i <reference_genome_base.bin> \n");
        return -1;
    }
    if ( OUTdir.empty() )    {
        printf ("Not defined:  -o <output_directory> \n");
        return -1;
    }
    if ( amtMut <= 0 )  {
        printf ("Not defined:  -mut value \n");
        return -1;
    }

    return nP;
}
/////////////////////////////////////////////////////////////////////////

int readGenBase( )
{
    clock_t Tstart = clock();
    printf( "readGenBase( ) ... \n");
    
    GENOM.GenBasFile = fopen(ArgKit.BASEpath.c_str(), "rb");
    if ( ! GENOM.GenBasFile )   {
        printf ( "Cannt OPEN  '%s'\n", ArgKit.BASEpath.c_str() );
        return -1;
    }
// index_XRO read  ----------------
    GENOM.readIndXro ( vecDNK );
    
// index read  --------------------
    if ( readIndex (GENOM.vIndexG, GENOM.GenBasFile ) < 0 )
        return -1;

//    readFooter( );
    
    if ( ! checkKeySet( GENOM.vIndexG ) )
        return -1;
    
    clock_t Tstop = clock();
    double  duration = (double)(Tstop - Tstart) / CLOCKS_PER_SEC;
    printf( "------> end load dT==%5.2f\n", duration );

    return 0;
}
/////////////////////////////////////////////////////////////////////////

bool is_dir(const char  *path ) {
    struct stat statbuf;
    if (stat(path, &statbuf) != 0) {
        return false; // Ошибка при получении информации о файле
    }
    return S_ISDIR(statbuf.st_mode);
}
//////////////////////////////////////

bool is_file(const char *path ) {
    struct stat statbuf;
    if (stat(path, &statbuf) != 0) {
        return false; // Ошибка при получении информации о файле
    }
    return S_ISREG(statbuf.st_mode);
}
/////////////////////////////////////////////////////////////////////////

int fgets_ShortRec( char *shortRec, int sizeRec, FILE *f_in )
{
    char buff[4096];
    int lng;
    
    if ( fgets(shortRec, sizeRec, f_in)==NULL )
        return 0;
    lng = (int)strlen(shortRec);
    if ( *(shortRec+lng-1) != '\n' )
        while ( 1 ) {
            if ( fgets(buff, sizeof(buff), f_in)==NULL )
                break;
            if ( buff[strlen(buff)-1] == '\n' )
                break;
        }
    
    return lng;
}

////////////////////////////////////////////////////// //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void  xtrFileName ( const char  *Fpath, string &sName)
{
    const char *pFp = Fpath + strlen(Fpath) - 1;
    
    while ( *pFp != '/' ) {
        if ( pFp == Fpath )
            break;
        pFp--;
    }
    if ( *pFp == '/' )
        pFp++;
    
    //    const char *pExt = strstr(pFp, ".vcf");
    const char *pExt = Fpath + strlen(Fpath) - 1;
    while ( *pExt != '.' ) {
        if ( pExt == Fpath )
            break;
        pExt--;
    }
    if ( *pExt == '.' )
        pExt--;
    
    sName.clear();
    for ( const char *p=pFp; p<=pExt;  p++ )
        sName += *p;
    
    return;
}
////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
