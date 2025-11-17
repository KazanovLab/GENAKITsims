//
//  genom.cpp
//  MutaGena
//
//  Created by Gennady on 9/3/25.
//  Copyright © 2025 Gennady. All rights reserved.
//

#include <stdio.h>
#include <vector>
#include <string>
#include <time.h>
#include <random>

#include "cmain.h"
#include "xrosoma.h"
#include "mapping.h"
#include "mutaGen.h"

extern PROGARGS ArgKit;
extern HUGEN GENOM;
extern vector < XROSOMA > vecDNK;
extern KEY_MAP keySet;
extern FILE *Ftrace;

vector < MUTAGEN > MutGenTab;
vector < MUTSIG >  vMutType;

/////////////////////////////////////////////////////////////////////////

bool lesser_RECRD (  const FILE_RECRD &x1, const FILE_RECRD &x2 )
{
    return x1.mPos < x2.mPos;
}
/////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

int openXtempFiles( )
{
    XROSOMA *pXro;
    for ( int nX=0; nX<vecDNK.size(); nX++ )  {
        pXro = &vecDNK[nX];
        pXro->XtempFilePath = ArgKit.OUTdir + "temp_" + pXro->XroID + ".txt";
        pXro->XtempFile = fopen(pXro->XtempFilePath.c_str(), "w+");
        if ( ! pXro->XtempFile )   {
            printf ( "Cannt OPEN  '%s'\n", pXro->XtempFilePath.c_str() );
            return -1;
        }
    }
    return 0;
}
//////////////////////////////////////////////////////////////////////////

int mergeXtempFiles( )
{    
    int nXro;
    int maxFileSize=0;
    FILE * fOut = NULL;
    char Buffr[1024];
    
    printf ("mergeXtempFiles( )\n");
    sprintf(Buffr,"%smyVCF_%d.txt", ArgKit.OUTdir.c_str(), ArgKit.amtMut);
    fOut = fopen(Buffr, "w");
    if ( ! fOut )   {
        printf("Inv.OPENop '%s'\n", Buffr);
        return -1;
    }
    fprintf(fOut, "CHR\tPOS\tREF\tALT\tMOTIF\tCoding\tGenes\tRT\tSTRAND\tMutaGen\n");
    for ( nXro=0; nXro<vecDNK.size(); nXro++ )    {
        if ( maxFileSize < vecDNK[nXro].cntOutMut )
            maxFileSize = vecDNK[nXro].cntOutMut;
    }
    
    FILE_RECRD curRec;
    vector <FILE_RECRD> vecRec;
    vecRec.reserve(maxFileSize);
    vector <FILE_RECRD> ::iterator itR;
    
    int cntR;
    int pos;
    for ( nXro=0; nXro<vecDNK.size(); nXro++ )    {
        cntR = 0;
        vecRec.clear();
        if ( fseek( vecDNK[nXro].XtempFile, 0, SEEK_SET) ) {
            printf ("mergeXtempFiles( ): Fseek failed for Xromo.%s tempFile\n", vecDNK[nXro].XroID.c_str() );
            return -1;
        }
        while ( fgets(Buffr, sizeof(Buffr)-1, vecDNK[nXro].XtempFile) ) {
            pos = 0;
            
            sscanf(Buffr, "%*s\t%d\t", &pos);
            if ( pos <= 0 ) {
                printf ("mergeXtempFiles( ): POS=%d for Xromo.%s tempFile '%s'\n", pos,
                        vecDNK[nXro].XroID.c_str(), Buffr );
                return -1;
            }
            curRec.mPos = pos;
            curRec.Recrd = Buffr;
            vecRec.push_back(curRec);
            cntR++;
        }
        sort(vecRec.begin(), vecRec.end(), lesser_RECRD );
        for ( itR = vecRec.begin(); itR != vecRec.end(); itR++ )
            fprintf(fOut, "%s", itR->Recrd.c_str() );
    }
    fclose(fOut);
    
    return 0;
}
//////////////////////////////////////////////////////////////////////////

int makeMutTab( )
{
    int n, v;
    struct INIT_MUT {
        char shortID[8];
        char fullID [MUTGEN_ID_SIZE];
    }
    iniMut[4] = { {"-ap", "APOBEC"}, {"-uv", "UV"},  {"-sm", "SMOKING"},  {"-cl", "CLOCK"} };
    
    MutGenTab.clear();
    for ( n=0; n<4; n++)
        MutGenTab.push_back(MUTAGEN(iniMut[n].fullID));
    
    for ( n=0; n<ArgKit.vMGfract.size(); n++ )  {
        for ( v=0; v<MutGenTab.size(); v++ )    {
            if ( strcmp(ArgKit.vMGfract[n].first, iniMut[v].shortID) == 0 )
                break;
        }
        if ( v >= 4 )   {
            printf("UnKnown full name for '%s'\n", ArgKit.vMGfract[n].first );
            return -1;
        }
//        MutGenTab.push_back(MUTAGEN(iniMut[v].fullID));
        MutGenTab[v].amtMut = ArgKit.amtMut * ArgKit.vMGfract[n].second;
    }
    
    return 0;
}
//////////////////////////////////////////////////////////////////////////

int MUTAGEN:: splitMut2Key( )
{
    vector < double > ::iterator ITfract;
    KEY_MAP ::iterator itKey;
    int indM, indRT, coding, gen, strnd;
    double fract;
    char Tag, invert ='\0';
    char Key[MOTKEY_SIZE];
    char Motif[MOTKEY_SIZE];
    int nMut;
    
    GENOM.initALT_MUT_CNTS();
    
    for ( nMut=0; nMut<amtMut; nMut++ ) {   // amtMut=%
        Tag = '\0';
// ----- mutation
        fract = (double)rand()/(double)RAND_MAX;
        ITfract = lower_bound(fract_mut.begin(), fract_mut.end(), fract);
        indM = (int)(ITfract - fract_mut.begin());
            strcpy(Motif, vMutType[indM].motif);
// ----- Gene
        fract = (double)rand()/(double)RAND_MAX;
        gen = (fract <= fract_gen[0]) ? 0 : 1;
        if ( gen )
            SET_GEN_TAG(&Tag);
// ----- RT
        fract = (double)rand()/(double)RAND_MAX;
        if ( gen )  {
            ITfract = lower_bound(fract_rt_gene.begin(), fract_rt_gene.end(), fract);
            indRT = (int)(ITfract-fract_rt_gene.begin());
        }
        else    {
            ITfract = lower_bound(fract_rt_intg.begin(), fract_rt_intg.end(), fract);
            indRT = (int)(ITfract-fract_rt_intg.begin());
        }
            SET_RT_VAL(&Tag, indRT);
// ----- Coding
        fract = (double)rand()/(double)RAND_MAX;
        coding = (fract <= fract_cod[0]) ? 0 : 1;
            if ( coding )
                SET_CODING_TAG(&Tag);
// ----- Strand
        if ( indRT > 0 )    {
            fract = (double)rand()/(double)RAND_MAX;
            strnd = (fract<=fract_strand[0]) ? 0 : ( (fract<=fract_strand[1]) ? 1 : 2 );
            if ( strnd==1)
                SET_LEADING_TAG(&Tag);
            if ( strnd==2)
                SET_LAGGING_TAG(&Tag);
        }

//        setKEYbyTag (Key, &Tag, invert);
        invert = formKEY (Motif, &Tag, Key);        // = 0 alwey
        itKey = keySet.find(Key);
        if ( itKey == keySet.end() )  {
            printf( "splitMut2Key( ): KEY [%s] not found at MAP\n", Key);
            return -1;
        }
        GENOM.addALT_MUT_CNTS(itKey->second, vMutType[indM].cALT);
    }
    
    return nMut;
}
//////////////////////////////////////////////////////////////////////////
// int comparePos_( int nMot, INV_POS &PosG);   //  mark='comparePos_Xro ()'

int HUGEN:: rndMutation( const char *mgName ) 
{
    unsigned int  Rnd_pos;
    INV_POS PosG;
    char cREF, cALT;
    char cod, gen, rt, strnd;
    char Coding[12], Genes[12], Strand[12];
    char motif[4];
    unsigned int PosXmotif;
    XROSOMA *pXro;
    int cntMut=0;
    int RetC;
    
    memset(motif, '\0', sizeof(motif));
    for ( int nK=0; nK<vIndexG.size(); nK++)    {
        if ( vIndexG[nK].amtPos <= 0 )
            continue;
        RetC = vAltCnt[nK].cntP[0] + vAltCnt[nK].cntP[1] + vAltCnt[nK].cntP[2] + vAltCnt[nK].cntP[3];
        if ( RetC==0 )
            continue;
        
        RetC = readData (vIndexG[nK], BufDataG, GenBasFile );
        
        cod = vIndexG[nK].mKey[_iKEY_CODING];
        switch (cod) {
            case 'c':
                strcpy(Coding, "Coding");
                break;
            case 'n':
                strcpy(Coding, "NonCoding");
                break;
            default:
                strcpy(Coding, "-");
                break;
        }
        gen    = vIndexG[nK].mKey[_iKEY_GENE];
        switch (gen) {
            case 'g':
                strcpy(Genes, "Genes");
                break;
            case 'i':
                strcpy(Genes, "InterGenes");
                break;
            default:
                strcpy(Genes, "-");
                break;
        }
        strnd = vIndexG[nK].mKey[_iKEY_STRAND];
        switch (strnd) {
            case '>':
                strcpy(Strand, "leading");
                break;
            case '<':
                strcpy(Strand, "lagging");
                break;
            default:
                strcpy(Strand, "-");
                break;
        }
        rt     = vIndexG[nK].mKey[_iKEY_RT];
        
        for ( int nALT=0; nALT<4; nALT++ )  {
            if ( vAltCnt[nK].cntP[nALT]==0 )
                continue;
            cALT = getNuc (nALT);
            for ( int cnt=0; cnt<vAltCnt[nK].cntP[nALT]; cnt++ )    {
                Rnd_pos = ( (double)rand() / (double)RAND_MAX ) * (double)(vIndexG[nK].amtPos-1);
                PosG = BufDataG[Rnd_pos];
                pXro = findXroByPOS( PosG, PosXmotif );
                if ( ! pXro )   {
                    printf("rndMutation(%s):: Pos=%d, lastPosGENOME=%d\n", mgName,
                           PosG.second, vecDNK.back().XstopPos);
                    return -1;
                }
//                if ( comparePos_( nK, PosG) < 0 )     //  mark='comparePos_Xro ()'
//                    continue;                         //  mark='comparePos_Xro ()'
                
                if ( PosG.first )
                    rewindMotif(GENOM.vIndexG[nK].mKey, motif);
                else
                    strncpy(motif, GENOM.vIndexG[nK].mKey, MOTKEY_SIZE);
                cREF = motif[1];
                cALT = ( PosG.first ) ? getCmpl_Nuc(getNuc (nALT)) : getNuc (nALT);
//  mark='comparePos_Xro ()'
//                if ( pXro->testValidDNK( (int)PosXmotif+1, cREF ) == 0 )   {
//                    printf( "Key[%s] Ref'%c'\n",vIndexG[nK].mKey, cREF);
//                    continue;
//                }   //  mark='comparePos_Xro ()'

               motif[3] = '\0';
//              CHR POS REF ALT MOTIF CODING GENES RT STRAND MUTAGEN
                fprintf(pXro->XtempFile, "%s\t%u\t%c\t%c\t%s\t%s\t%s\t%c\t%s\t%s\n",
                        pXro->XroID.c_str(), PosXmotif+1, cREF,cALT, motif,
                        Coding, Genes, rt, Strand, mgName );
                pXro->cntOutMut++;
                cntMut++;
            }
        }
        
    }
    
    return cntMut;
}
//////////////////////////////////////////////////////////////////////////

char *scanNextMG( char *pB, char *mutID )
{
    *mutID = '\0';
    if ( ! *pB  || *pB=='\n' || *pB=='\r' )
        return NULL;
        
    while ( *pB > ' ' ) *mutID++ = *pB++;
    *mutID = '\0';
    if ( *pB == '\t' )
        pB++;
    
    return pB;
}
//////////////////////////////////////////////////////////////////////////

int testHeader(char *hdr, int indMG[], char *fName)
{
    char *pB = hdr;
    char MutaID[64];
    int nMG;
    int nCol=0;

    for ( int n=0; n<MutGenTab.size(); n++ )
        indMG[n] = -1;
        
    while ( pB )    {
        pB = scanNextMG( pB, MutaID );
        if ( ! MutaID[0] )
            continue;
//        if ( nMG >= MutGenTab.size()  )   {
//            nMG++;
//            continue;
//        }
        for ( nMG=0; nMG<MutGenTab.size(); nMG++ ) {
            if ( strcmp(MutaID, MutGenTab[nMG].MGname) == 0 )
                break;
        }
        if ( nMG >= MutGenTab.size() )    {
            printf("%s :: Unknown  Mutagen_Name '%s'\n", fName, MutaID);
            continue;
        }
        indMG[nCol] = nMG;
        nCol++;
    }
    if ( nCol==0 )  {
        printf("%s :: Empty table header in file '%s'\n", fName, MutaID);
        return -1;
    }
/*    if ( nMG != MutGenTab.size() )    {
        printf("%s :: Mismatch of Mutagen count=%d : must be = %d\n",
               fName, nMG, (int)MutGenTab.size() );
        return -1;
    }
 */
    return 0;
}
/////////////////////////////////////////////////////////////////////////

int parsRecrd(char *Buff, char *Rtyp, double Rvals[], int refMG[] )
{
    int nMG;
    char *pB=Buff;
    int ref;
    
    
    if ( ! Buff ) {
        printf ("Not enough lines at file  ... ");
        return 0;
    }
    for ( int n=0; n<MAX_MUTAGEN; n++ )
        Rvals[n] = 0;
    
    while ( *pB == ' ' ) pB++;
    if ( *pB < ' ' )    {
        printf ("Not defined RecType '%s' at file ...", Buff );
        return 0;
    }
    
    while ( *pB != '\t')
        *Rtyp++ = *pB++;
    *Rtyp = '\0';
    pB++;
    for ( nMG=0; nMG<MutGenTab.size(); nMG++ )    {
        if ( (ref = refMG[nMG]) < 0 )
            break;
        if ( !*pB ||  *pB=='\n' )
            break;
        
        while ( *pB==' ') pB++;
        if (  *pB != '\t' )
            sscanf(pB, "%lf", &Rvals[ref]);
        while ( *pB>=' ') pB++;
        pB++;
    }
    if ( nMG==0 )
        printf ("Not Empty Rec '%s' at file ...", Buff );
    
    return nMG;
}
//////////////////////////////////////////////////////////////////////////

int readCTGR( )
{
    FILE * f_CAT;
    string fPath;
    char Buff[512];
    char f_motif[]  = "mutsig.txt";
    char f_coding[] = "TranscriptionStrand.txt";
    char f_gen[]    = "genes.txt";
    char f_rt_gene[]     = "ReplicationTiming_gene.txt";
    char f_rt_intg[]     = "ReplicationTiming_inter.txt";
    char f_strand[] = "ReplicationStrand.txt";
    int refMG[MAX_MUTAGEN];

// ----- mutation -------------------------------------------------
    fPath = ArgKit.INdir + f_motif;
    if ( ! (f_CAT = fopen(fPath.c_str(), "r")) ) {
        printf ( "Cannt OPEN  '%s'\n", fPath.c_str() );
        return -1;
    }
    fgets(Buff, sizeof(Buff)-2, f_CAT);
    if ( testHeader(Buff, refMG, f_motif) < 0 )
        return -1;
    if ( loadMutRanges(f_CAT, refMG, f_motif) < 0 )
        return -1;
    fclose(f_CAT);
// ----- RT_gene -------------------------------------------------
    fPath = ArgKit.INdir + f_rt_gene;
    if ( ! (f_CAT = fopen(fPath.c_str(), "r")) ) {
        printf ( "Cannt OPEN  '%s'\n", fPath.c_str() );
        return -1;
    }
    fgets(Buff, sizeof(Buff)-2, f_CAT);
    if ( testHeader(Buff, refMG, f_rt_gene) < 0 )
        return -1;
    if ( loadRtRanges(f_CAT, refMG, f_rt_gene) < 0 )
        return -1;
    fclose(f_CAT);
// ----- RT_interGene ----------------------------------------------
    fPath = ArgKit.INdir + f_rt_intg;
    if ( ! (f_CAT = fopen(fPath.c_str(), "r")) ) {
        printf ( "Cannt OPEN  '%s'\n", fPath.c_str() );
        return -1;
    }
    fgets(Buff, sizeof(Buff)-2, f_CAT);
    if ( testHeader(Buff, refMG, f_rt_intg) < 0 )
        return -1;
    if ( loadRtRanges(f_CAT, refMG, f_rt_intg) < 0 )
        return -1;
    fclose(f_CAT);
    
// ----- Coding ----------------------------------------------
    fPath = ArgKit.INdir + f_coding;
    if ( ! (f_CAT = fopen(fPath.c_str(), "r")) ) {
        printf ( "Cannt OPEN  '%s'\n", fPath.c_str() );
        return -1;
    }
    fgets(Buff, sizeof(Buff)-2, f_CAT);
    if ( testHeader(Buff, refMG, f_coding) < 0 )
        return -1;
    if ( loadCodRanges(f_CAT, refMG, f_coding) < 0 )
        return -1;
    fclose(f_CAT);
// ----- Genes ----------------------------------------------
    fPath = ArgKit.INdir + f_gen;
    if ( ! (f_CAT = fopen(fPath.c_str(), "r")) ) {
        printf ( "Cannt OPEN  '%s'\n", fPath.c_str() );
        return -1;
    }
    fgets(Buff, sizeof(Buff)-2, f_CAT);
    if ( testHeader(Buff, refMG, f_gen) < 0 )
        return -1;
    if ( loadGenRanges(f_CAT, refMG, f_gen) < 0 )
        return -1;
    fclose(f_CAT);
// ----- Strand ----------------------------------------------
    fPath = ArgKit.INdir + f_strand;
    if ( ! (f_CAT = fopen(fPath.c_str(), "r")) ) {
        printf ( "Cannt OPEN  '%s'\n", fPath.c_str() );
        return -1;
    }
    fgets(Buff, sizeof(Buff)-2, f_CAT);
    if ( testHeader(Buff, refMG, f_gen) < 0 )
        return -1;
    if ( loadStrandRanges(f_CAT, refMG, f_gen) < 0 )
        return -1;
    fclose(f_CAT);
// ===========================================================
    
    for ( int nMu=0; nMu<MutGenTab.size(); nMu++ )  {
        if ( MutGenTab[nMu].amtMut==0 )
            continue;
        if ( MutGenTab[nMu].beta_dirichlet( ) < 0 )
            return -1;
    }

    return 0;
}
//////////////////////////////////////////////////////////////////////////

int loadMutRanges(FILE *f_CAT, int refMG[], char *fName)
{
    char Buff[512];
    char *pB;
    MUTSIG recrd;
    char MutType[64];
    double portion;
    int nMG, ref;
    
    while ( 1 )     {
        if ( ! ( pB=fgets(Buff, sizeof(Buff)-2, f_CAT) ) )
            break;
// ------ create table of mutation types -------------
        if ( sscanf(pB, "%s\t", MutType ) != 1 )    {
            printf ( "%s :: Inv.format motif in rec: '%s'\n", fName, Buff);
            return -1;
        }
        if ( strlen(MutType) != 7 )  {
            printf ( "%s :: Inv.format motif='%s' in rec: '%s'\n", fName, MutType, Buff);
            return -1;
        }
        sscanf ( MutType,"%c[%c>%c]%c", &recrd.motif[0], &recrd.motif[1], &recrd.cALT, &recrd.motif[2] );
        recrd.cREF = recrd.motif[1];
        strcpy(recrd.muttype, MutType);
        vMutType.push_back(recrd);
        while ( *pB != '\t' ) pB++;
        
//-------- fill field fract_mut at Mutagen
        for ( nMG=0; nMG<MutGenTab.size(); nMG++ )  {
            if ( (ref = refMG[nMG]) < 0 )
                break;
            if ( ! *pB )
                break;
            portion = 0;
            if ( (sscanf(pB, "\t%lf", &portion)==0) &&
                  MutGenTab[ref].amtMut > 0  )    {
                printf ( "%s :: '%s' empty rangeValue for %s : '%s'\n",
                        fName, MutType, MutGenTab[ref].MGname, Buff);
                return -1;
            }
            if ( MutGenTab[ref].amtMut == 0 )
                continue;
            
            MutGenTab[ref].fract_mut.push_back (portion);
            
            pB++;   //// '\t'
            while ( *pB && (*pB != '\t') )
                pB++;
        }

    }
    
    for ( nMG=0; nMG<MutGenTab.size(); nMG++ )  {
        if ( (ref = refMG[nMG]) < 0 )
            break;
        if ( MutGenTab[ref].amtMut == 0 )
            continue;
        for ( int p=1; p<MutGenTab[nMG].fract_mut.size(); p++ )
            MutGenTab[nMG].fract_mut[p] += MutGenTab[nMG].fract_mut[p-1];
    }
    for ( nMG=0; nMG<MutGenTab.size(); nMG++ )  {
        if ( (ref = refMG[nMG]) < 0 )
            break;
        if ( MutGenTab[ref].amtMut == 0 )
            continue;
        if ( MutGenTab[nMG].fract_mut.back() < DOUBLE_ONE )
            printf ( "%s :: Check distribution of shares for %s (sum<1)",
                    fName, MutGenTab[nMG].MGname );
        MutGenTab[nMG].fract_mut.back() = (double) 1;
    }

    return 0;
}
//////////////////////////////////////////////////////////////////////////

int loadRtRanges(FILE *f_CAT, int refMG[], char *fName)
{
    char RecTypes[] =  "1RT1_2RT2_3RT3_4RT4_5RT5_6RT6_7RT7_0undefRT_";
    char Buff[512];
    int recTyp[RT_TAB_SIZE];
    
    double rVals[MAX_MUTAGEN];
    char scanTyp[32], *pTyp;
    char *pB;
    int RT_gene;

    RT_gene =  ( strstr(fName, "inter") ) ? 0 : 1;
    for ( int n=0; n<MutGenTab.size(); n++ )
        for ( int p=0; p<RT_TAB_SIZE; p++ ) {
            if ( RT_gene )
                MutGenTab[n].fract_rt_gene.push_back (0);
            else
                MutGenTab[n].fract_rt_intg.push_back (0);
        }
    for ( int Line=0; Line<RT_TAB_SIZE; Line++)   {
        pB = fgets(Buff, sizeof(Buff)-2, f_CAT);
        if ( parsRecrd(pB, scanTyp, rVals, refMG) == 0 )  {
            printf ("%s\n", fName);
            return -1;
        }
        if ( ! ( pTyp=strstr(RecTypes, scanTyp)) )    {
            printf ("Inv. RT_Type '%s' at file  '%s'\n", scanTyp, fName );
            return -1;
        }
        recTyp[Line] = atoi(pTyp-1);
        for ( int n=0; n<MutGenTab.size(); n++ )    {
            if ( RT_gene )
                MutGenTab[n].fract_rt_gene[recTyp[Line]] = rVals[n];
            else
                MutGenTab[n].fract_rt_intg[recTyp[Line]] = rVals[n];
        }
    }

    for ( int nMu=0; nMu<MutGenTab.size(); nMu++ )  {
        if ( MutGenTab[nMu].amtMut==0 )
            continue;
        for ( int n=0; n<RT_TAB_SIZE-1; n++ ) {
            rVals[0] = ( RT_gene )  ? MutGenTab[nMu].fract_rt_gene[n]
                                    : MutGenTab[nMu].fract_rt_intg[n];
            if ( rVals[0]==0 )   {
                printf ("%s :: MutGen='%s' : all RTvalues must be > 0 \n",
                        fName, MutGenTab[nMu].MGname);
                return -1;
            }
        }
    }
    
    return 0;
}
/////////////////////////////////////////////////////////////////////////

int loadCodRanges(FILE *f_CAT, int refMG[], char *fName)
{
    char RecTypes[] =  "1Coding_0Noncoding_";
    char Buff[512];
    int recTyp[2];
    double rVals[MAX_MUTAGEN];
    char scanTyp[32], *pTyp;
    char *pB;

    for ( int Line=0; Line<2; Line++)   {
        pB = fgets(Buff, sizeof(Buff)-2, f_CAT);
        if ( parsRecrd(pB, scanTyp, rVals, refMG) == 0 )  {
            printf ("%s\n", fName);
            return -1;
        }
        if ( ! ( pTyp=strstr(RecTypes, scanTyp)) )    {
            printf ("Inv. CodingType '%s' at file  '%s'\n", scanTyp, fName );
            return -1;
        }
        recTyp[Line] = atoi(pTyp-1);
        for ( int n=0; n<MutGenTab.size(); n++ )
            MutGenTab[n].fract_cod[recTyp[Line]] = rVals[n];
    }
    if ( recTyp[0]==recTyp[1])  {
        printf ("%s :: redefined CodingType : '%s'\n", fName, Buff);
        return -1;
    }
    
    for ( int n=0; n<MutGenTab.size(); n++ )  {
        if ( MutGenTab[n].amtMut==0 )
            continue;
        if ( MutGenTab[n].fract_cod[0]==0 || MutGenTab[n].fract_cod[1]==0 ) {
            printf ("%s :: MutGen='%s' : Coding and non-Coding must be > 0 : %lf %lf\n",
                    fName, MutGenTab[n].MGname, MutGenTab[n].fract_cod[1], MutGenTab[n].fract_cod[0]);
            return -1;
        }
    }
    
    return 0;
}
/////////////////////////////////////////////////////////////////////////

int loadGenRanges(FILE *f_CAT, int refMG[], char *fName)
{
    char RecTypes[] =  "1Genes_0Intergenes_";
    char Buff[512];
    int recTyp[2];
    double rVals[MAX_MUTAGEN];
    char scanTyp[32], *pTyp;
    char *pB;
    
    for ( int Line=0; Line<2; Line++)   {
        pB = fgets(Buff, sizeof(Buff)-2, f_CAT);
        if ( parsRecrd(pB, scanTyp, rVals, refMG) == 0 )  {
            printf ("%s\n", fName);
            return -1;
        }
        if ( ! ( pTyp=strstr(RecTypes, scanTyp)) )    {
            printf ("Inv. GeneType '%s' at file  '%s'\n", scanTyp, fName );
            return -1;
        }
        recTyp[Line] = atoi(pTyp-1);
        for ( int n=0; n<MutGenTab.size(); n++ )
            MutGenTab[n].fract_gen[recTyp[Line]] = rVals[n];
    }
    if ( recTyp[0]==recTyp[1])  {
        printf ("%s :: redefined GeneType : '%s'\n", fName, Buff);
        return -1;
    }
    
    for ( int n=0; n<MutGenTab.size(); n++ )  {
        if ( MutGenTab[n].amtMut==0 )
            continue;
        if ( MutGenTab[n].fract_gen[0]==0 || MutGenTab[n].fract_gen[1]==0 ) {
            printf ("%s :: MutGen='%s' : Gene and interGene must be > 0 : %lf %lf\n",
                    fName, MutGenTab[n].MGname, MutGenTab[n].fract_gen[1], MutGenTab[n].fract_gen[0]);
            return -1;
        }
    }
    
    return 0;
}
/////////////////////////////////////////////////////////////////////////

int loadStrandRanges(FILE *f_CAT, int refMG[], char *fName)
{
    char RecTypes[] =  "0undefStrend_1Leading_2Lagging_";
    char Buff[512];
    int recTyp[3];
    double rVals[MAX_MUTAGEN];
    char scanTyp[32], *pTyp;
    char *pB;
    
    for ( int n=0; n<MutGenTab.size(); n++ )   {
        for ( int p=0; p<3; p++ )
            MutGenTab[n].fract_strand.push_back (0);
    }
    for ( int Line=0; Line<3; Line++)   {
        pB = fgets(Buff, sizeof(Buff)-2, f_CAT);
        if ( parsRecrd(pB, scanTyp, rVals, refMG) == 0 )  {
            printf ("%s\n", fName);
            return -1;
        }
        if ( ! ( pTyp=strstr(RecTypes, scanTyp)) )    {
            printf ("Inv. StrandType '%s' at file  '%s'\n", scanTyp, fName );
            return -1;
        }
        recTyp[Line] = atoi(pTyp-1);
        for ( int n=0; n<MutGenTab.size(); n++ )
            MutGenTab[n].fract_strand[recTyp[Line]] = rVals[n];
    }
    if ( recTyp[0]==recTyp[1] || recTyp[0]==recTyp[2] || recTyp[1]==recTyp[2])  {
        printf ("%s :: redefined StrandType \n", fName);
        return -1;
    }
    
    for ( int n=0; n<MutGenTab.size(); n++ )  {
        if ( MutGenTab[n].amtMut==0 )
            continue;
        if ( MutGenTab[n].fract_strand[0]==0 || MutGenTab[n].fract_strand[1]==0 ||
             MutGenTab[n].fract_strand[2]==0 ) {
            printf ("%s :: MutGen='%s' : Leading, Lagging, undefStrend must be > 0 : %lf %lf %lf\n",
                    fName, MutGenTab[n].MGname,
                    MutGenTab[n].fract_strand[1], MutGenTab[n].fract_strand[2], MutGenTab[n].fract_strand[0]);
            return -1;
        }
    }
    
    return 0;
}
/////////////////////////////////////////////////////////////////////////

double rbeta(double alpha, double beta, mt19937_64& rng)
{
    if (alpha <= 0.0 || beta <= 0.0)    {
        printf ("rbeta:: alpha and beta must be > 0");
        return (double)-1;
    }
    
    gamma_distribution<double> gA(alpha, 1.0);
    gamma_distribution<double> gB(beta, 1.0);
    
    double x = gA(rng);
    double y = gB(rng);
    return x / (x + y);  // Beta(alpha, beta)
}
/////////////////////////////////////////////////////////////////////////

vector<double> rdirichlet(const vector<double>& alpha, mt19937_64& rng)
{
    vector<double> x(alpha.size());
    double sum = 0.0;
    
// генерируем гаммы для каждой компоненты
    for (size_t i = 0; i < alpha.size(); ++i) {
        if (alpha[i] <= 0.0)    {
            printf ("rdirichlet:: all alpha[i] must be > 0 ");
            x[0] = -1;
            return x;
        }
        gamma_distribution<double> g(alpha[i], 1.0);
        x[i] = g(rng);
        sum += x[i];
    }
    
// нормализация
    for (size_t i = 0; i < alpha.size(); ++i)
        x[i] /= sum;
    
    return x; // вектор вероятностей (сумма = 1)
}
//////////////////////////////////////////////////////////////////////////

int MUTAGEN:: beta_dirichlet( )
{
    random_device rd;
    mt19937_64 rng(rd());
    vector < double > pv;
/*
 vector < double > fract_rt;     //  = 0:7  (RT0 ... RT7)
 double fract_cod[2];            //  [0]=noncoding; [1]=coding
 double fract_gen[2];            //  [0]=intergenes; [1]=genes
 vector < double > fract_strand; //  [0]=undefined; [1]=leading; [2]=lagging
*/

// Gene/intergene
    if (  (fract_gen[0]=rbeta(fract_gen[0], fract_gen[1], rng) ) < 0 )    {
        printf ( " for Gene/intergene\n");
        return -1;
    }
    fract_gen[1] = (double) 1;
// Coding/non-coding
    if (  (fract_cod[0]=rbeta(fract_cod[0], fract_cod[1], rng) ) < 0)   {
        printf ( " for Coding/non-coding\n");
        return -1;
    }
    fract_cod[1] = (double) 1;
    
// Leading/lagging
    pv = rdirichlet(fract_strand, rng);
    if ( pv[0] < 0.0 )  {
        printf ( " for Leading/lagging\n");
        return -1;
    }
    fract_strand[0] = pv[0];
    fract_strand[1] = fract_strand[0] + pv[1];
    fract_strand[2] = fract_strand[1] + pv[2];
    
// RT_gene
    pv.clear();
    pv = rdirichlet(fract_rt_gene, rng);
    if ( pv[0] < 0.0 )  {
        printf ( " for RT_gene\n");
        return -1;
    }
    fract_rt_gene = pv;
    for ( int n=1; n<RT_TAB_SIZE; n++ )
        fract_rt_gene[n] += fract_rt_gene[n-1];
    
// RT_intergene
    pv.clear();
    pv = rdirichlet(fract_rt_intg, rng);
    if ( pv[0] < 0.0 )  {
        printf ( " for RT_interGene\n");
        return -1;
    }
    fract_rt_intg = pv;
    for ( int n=1; n<RT_TAB_SIZE; n++ )
        fract_rt_intg[n] += fract_rt_intg[n-1];

    return 0;
}
//////////////////////////////////////////////////////////////////////////











