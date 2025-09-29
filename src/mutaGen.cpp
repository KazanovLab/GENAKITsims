//
//  genom.cpp
//  MutaGena
//
//  Created by Gennady on 9/3/25.
//  Copyright © 2025 Gennady. All rights reserved.
//

#include <stdio.h>
#include "vector"
#include "string"
#include <time.h>

#include "cmain.h"
#include "xrosoma.h"
#include "mapping.h"
#include "mutaGen.h"

extern PROGARGS ArgKit;
extern HUGEN GENOM;
extern vector < XROSOMA > vecDNK;
extern KEY_MAP keySet;

vector < MUTAGEN > MutaGen;
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
    
    for ( n=0; n<ArgKit.vMGfract.size(); n++ )  {
        for ( v=0; v<4; v++ )
            if ( strcmp(ArgKit.vMGfract[n].first, iniMut[v].shortID) == 0 )
                break;
        if ( v >= 4 )   {
            printf("UnKnown full name for '%s'\n", ArgKit.vMGfract[n].first );
            return -1;
        }
        MutaGen.push_back(MUTAGEN(iniMut[v].fullID));
        MutaGen.back().amtMut = ArgKit.amtMut * ArgKit.vMGfract[n].second;
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
        fract = (double)rand()/(double)RAND_MAX;
        ITfract = lower_bound(fract_mut.begin(), fract_mut.end(), fract);
        indM = (int)(ITfract - fract_mut.begin());
            strcpy(Motif, vMutType[indM].motif);

        fract = (double)rand()/(double)RAND_MAX;     
        ITfract = lower_bound(fract_rt.begin(), fract_rt.end(), fract);
        indRT = (int)(ITfract-fract_rt.begin());
            SET_RT_VAL(&Tag, indRT);
        
        fract = (double)rand()/(double)RAND_MAX;
        gen = (fract <= fract_gen[0]) ? 0 : 1;
            if ( gen )
                SET_GEN_TAG(&Tag);
        
        fract = (double)rand()/(double)RAND_MAX;
        coding = (fract <= fract_cod[0]) ? 0 : 1;
            if ( coding )
                SET_CODING_TAG(&Tag);
        
        fract = (double)rand()/(double)RAND_MAX;
        strnd = (fract<=fract_strand[0]) ? 0 : ( (fract<=fract_strand[1]) ? 1 : 2 );
        if ( strnd==1)
            SET_LEADING_TAG(&Tag);
        if ( strnd==2)
            SET_LAGGING_TAG(&Tag);

//        setKEYbyTag (Key, &Tag, invert);
        invert = formKEY (Motif, &Tag, Key);
        itKey = keySet.find(Key);
        if ( itKey == keySet.end() )  {
            printf( "splitMut2Key( ): KEY [%s] not found at MAP\n", Key);
            return -1;
        }
//      unsigned int  Rnd_pos = ( (double)rand() / (double)RAND_MAX ) * (double)(GENOM.vIndexG[itKM->second].amtPos-1);
        GENOM.addALT_MUT_CNTS(itKey->second, vMutType[indM].cALT);
    }
    
    return nMut;
}
//////////////////////////////////////////////////////////////////////////

int HUGEN:: rndMutation( const char *mgName ) 
{
    unsigned int  Rnd_pos;
    INV_POS PosG;
    char cREF, cALT;
    char cod, gen, rt, strnd;
    char Coding[12], Genes[12], Strand[12];
    char motif[4];
    unsigned int PosX;
    XROSOMA *pXro;
    int cntMut=0;
    int RetC;
    
    memset(motif, '\0', sizeof(motif));
    for ( int nK=0; nK<vIndexG.size(); nK++)    {
        if ( vIndexG[nK].amtPos <= 0 )
            continue;
//        vData.vPosM.clear();   
        RetC = readData (vIndexG[nK], BufDataG, GenBasFile );
        if ( RetC != vIndexG[nK].amtPos )   {
            printf( "rndMutation(%s)::readData(%s)=%d != %d\n", mgName,
                   vIndexG[nK].mKey, RetC, vIndexG[nK].amtPos);
            return -1;
        }
        cod = vIndexG[nK].mKey[_iKEY_CODING];
        switch (cod) {
            case 'c':
                strcpy(Coding, "Coding");
                break;
            case 'n':
                strcpy(Coding, "NonCoding");
                break;
            default:
                Coding[0] = '\0';
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
                Genes[0] = '\0';
                break;
        }
        strnd = vIndexG[nK].mKey[_iKEY_STRAND];
        switch (strnd) {
            case '<':
                strcpy(Strand, "leading");
                break;
            case '>':
                strcpy(Strand, "lagging");
                break;
            default:
                Strand[0] = '\0';
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
                pXro = findXroByPOS( PosG, PosX );
                if ( ! pXro )
                    return -1;
                cREF = ( PosG.first ) ? getCmpl_Nuc( vIndexG[nK].mKey[1] ) : vIndexG[nK].mKey[1];
                cALT = ( PosG.first ) ? getCmpl_Nuc(getNuc (nALT)) : getNuc (nALT);
                strncpy(motif, vIndexG[nK].mKey, 3);
 //               motif[3] = '\0';
//              CHR POS REF ALT MOTIF CODING GENES RT STRAND MUTAGEN
                fprintf(pXro->XtempFile, "%s\t%u\t%c\t%c\t%s\t%s\t%s\t%c\t%s\t%s\n",
                        pXro->XroID.c_str(), PosX, cREF,cALT, motif, //vIndexG[nK].mKey,
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

int testHeader(char *hdr, char *fName)
{
    char *pB = hdr;
    char MutaID[64];
    int nMG=0;

    while ( pB )    {
        pB = scanNextMG( pB, MutaID );
        if ( ! MutaID[0] )
            continue;
        if ( nMG >= MutaGen.size()  )   {
            nMG++;
            continue;
        }
        if ( strcmp(MutaID, MutaGen[nMG].MGname) != 0 )    {
            printf("%s :: Mismatch of order at Mutagen List '%s'\n", fName, MutaID);
            return -1;
        }
        nMG++;
    }
    if ( nMG != MutaGen.size() )    {
        printf("%s :: Mismatch of Mutagen count=%d : must be = %d\n",
               fName, nMG, (int)MutaGen.size() );
        return -1;
    }
    return 0;
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
    char f_rt[]     = "ReplicationTiming.txt";
    char f_strand[] = "ReplicationStrand.txt";
//    int nMG=0;
//    MUTSIG recrd;
//    char MutType[64];
//    double portion;

    fPath = ArgKit.INdir + f_motif;
    if ( ! (f_CAT = fopen(fPath.c_str(), "r")) ) {
        printf ( "Cannt OPEN  '%s'\n", fPath.c_str() );
        return -1;
    }
    fgets(Buff, sizeof(Buff)-2, f_CAT);
    if ( testHeader(Buff,f_rt) < 0 )
        return -1;
/*
    char *pB = Buff;
//  --------------load header  (mutagen IDs) -----------------
    while ( pB )    {
        pB = scanNextMG( pB, MutType );
        if ( ! MutType[0] )
            continue;
        MutaGen.push_back(MUTAGEN(MutType));
    }
*/
    if ( loadMutRanges(f_CAT, f_motif) > 0 )
        return -1;
    fclose(f_CAT);
// ===========================================================
    fPath = ArgKit.INdir + f_rt;
    if ( ! (f_CAT = fopen(fPath.c_str(), "r")) ) {
        printf ( "Cannt OPEN  '%s'\n", fPath.c_str() );
        return -1;
    }
    fgets(Buff, sizeof(Buff)-2, f_CAT);
    if ( testHeader(Buff,f_rt) < 0 )
        return -1;
    if ( loadRtRanges(f_CAT,f_rt) < 0 )
        return -1;
    fclose(f_CAT);
// ===========================================================
    fPath = ArgKit.INdir + f_coding;
    if ( ! (f_CAT = fopen(fPath.c_str(), "r")) ) {
        printf ( "Cannt OPEN  '%s'\n", fPath.c_str() );
        return -1;
    }
    fgets(Buff, sizeof(Buff)-2, f_CAT);
    if ( testHeader(Buff,f_coding) < 0 )
        return -1;
    if ( loadCodRanges(f_CAT,f_coding) < 0 )
        return -1;
    fclose(f_CAT);
// ===========================================================
    fPath = ArgKit.INdir + f_gen;
    if ( ! (f_CAT = fopen(fPath.c_str(), "r")) ) {
        printf ( "Cannt OPEN  '%s'\n", fPath.c_str() );
        return -1;
    }
    fgets(Buff, sizeof(Buff)-2, f_CAT);
    if ( testHeader(Buff,f_gen) < 0 )
        return -1;
    if ( loadGenRanges(f_CAT,f_gen) < 0 )
        return -1;
    fclose(f_CAT);
// ===========================================================
    fPath = ArgKit.INdir + f_strand;
    if ( ! (f_CAT = fopen(fPath.c_str(), "r")) ) {
        printf ( "Cannt OPEN  '%s'\n", fPath.c_str() );
        return -1;
    }
    fgets(Buff, sizeof(Buff)-2, f_CAT);
    if ( testHeader(Buff,f_gen) < 0 )
        return -1;
    if ( loadStrandRanges(f_CAT,f_gen) < 0 ) 
        return -1;
    fclose(f_CAT);

    return 0;
}
//////////////////////////////////////////////////////////////////////////

int loadMutRanges(FILE *f_CAT, char *fName)
{
    char Buff[512];
    char *pB;
    MUTSIG recrd;
    char MutType[64];
    double portion;
    int nMG;
    
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
        for ( nMG=0; nMG<MutaGen.size(); nMG++ )  {
            if ( ! *pB )
                break;
            if ( sscanf(pB, "\t%lf", &portion)==0 )    {
                printf ( "%s :: '%s' empty rangeValue for %s : '%s'\n",
                        fName, MutType, MutaGen[nMG].MGname, Buff);
                return -1;
            }
            pB++;   //// '\t'
            MutaGen[nMG].fract_mut.push_back (portion);
            while ( *pB && (*pB != '\t') )
                pB++;
        }
        if ( nMG != MutaGen.size() )    {
            printf("%s :: Mismatch of Mutagen count=%d : must be = %d '%s'\n",
                   fName, nMG, (int)MutaGen.size(), Buff);
            return -1;
        }
    }
    
    for ( nMG=0; nMG<MutaGen.size(); nMG++ )  {
        for ( int p=1; p<MutaGen[nMG].fract_mut.size(); p++ )
            MutaGen[nMG].fract_mut[p] += MutaGen[nMG].fract_mut[p-1];
    }
    for ( nMG=0; nMG<MutaGen.size(); nMG++ )  {
        if ( MutaGen[nMG].fract_mut.back() < DOUBLE_ONE )
            printf ( "%s :: Check distribution of shares for %s (sum<1)",
                    fName, MutaGen[nMG].MGname );
        MutaGen[nMG].fract_mut.back() = (double) 1;
    }
/*
    vector < double > ::iterator itD;

     for ( int p=0; p<MutaGen[0].fract_mut.size(); p++ )  {
         portion = MutaGen[0].fract_mut[p];
         itD = lower_bound(MutaGen[0].fract_mut.begin(), MutaGen[0].fract_mut.end(), portion);
         printf ( "%d. --> %d\n", p, (int)(itD-MutaGen[0].fract_mut.begin()));
     }

    for ( int n=0; n<100; n++ )  {
        portion = (double)rand()/(double)RAND_MAX;
        itD = lower_bound(MutaGen[0].fract_mut.begin(), MutaGen[0].fract_mut.end(), portion);
        printf ( "%d. rnd=%lf ---> %d. %lf\n", n, portion, (int)(itD-MutaGen[0].fract_mut.begin()), *itD);
    }
*/
    return 0;
}
//////////////////////////////////////////////////////////////////////////

int loadRtRanges(FILE *f_CAT, char *fName)
{
    char Buff[512];
    char *pB;
    int RTn;
    double portion;
    int nMG;
    
    for ( int n=0; n<MutaGen.size(); n++ )   {
        for ( int p=0; p<RT_TAB_SIZE; p++ )
            MutaGen[n].fract_rt.push_back (0);
    }
    
    while ( 1 )     {
        if ( ! ( pB=fgets(Buff, sizeof(Buff)-2, f_CAT) ) )
            break;
        if ( sscanf(pB, "RT%d\t", &RTn ) != 1 )    {
            printf ( "%s :: Need RTn in rec: '%s'\n", fName, Buff);
            return -1;
        }
        if ( RTn >= RT_TAB_SIZE )  {
            printf ( "%s :: ignored RTn=%d  [0:7] in rec: '%s'\n", fName, RTn, Buff);
            continue;
        }
        while (*pB != '\t' ) pB++;
        
        for ( nMG=0; nMG<MutaGen.size(); nMG++ )  {
            if ( ! *pB )
                break;
            if ( sscanf(pB, "\t%lf", &portion)==0 )    {
                printf ( "%s :: 'RT%d' empty rangeValue for %s : '%s'\n",
                        fName, RTn, MutaGen[nMG].MGname, Buff);
                return -1;
            }
            pB++;   //// '\t'
            MutaGen[nMG].fract_rt[RTn] = portion;
            while ( *pB && *pB != '\t' ) pB++;
        }
        if ( nMG != MutaGen.size() )    {
            printf("%s :: Mismatch of Mutagen count=%d : must be = %d '%s'\n",
                   fName, nMG, (int)MutaGen.size(), Buff);
            return -1;
        }
    }
    for ( int n=0; n<MutaGen.size(); n++ )  {
        for ( int p=1; p<MutaGen[n].fract_rt.size(); p++ )
            MutaGen[n].fract_rt[p] += MutaGen[n].fract_rt[p-1];
        
    }
    for ( nMG=0; nMG<MutaGen.size(); nMG++ )  {
        if ( MutaGen[nMG].fract_rt.back() < DOUBLE_ONE )
            printf ( "%s :: Check distribution of shares for %s (sum<1)\n",
                    fName, MutaGen[nMG].MGname );
        MutaGen[nMG].fract_rt.back() = (double) 1;
    }
    return 0;
}
/////////////////////////////////////////////////////////////////////////

int parsRec(char *Buff, const char *RecTypes, vector < double > &vRanges)
{
    const char *pB;
    int nMG, typ;
    double portion;
    char rType[64];
    
    if ( ! Buff ) {
        printf ("Not enough lines at file  ... ");
        return -1;
    }
    vRanges.clear();
    sscanf(Buff, "%s\t", rType );
    if ( ! ( pB=strstr(RecTypes, rType)) )    {
        printf ("Inv. Type line '%s' at file  ... ", rType );
        return -1;
    }
    typ = atoi(pB-1);
    pB = Buff;
    while (*pB != '\t' ) pB++;
    for ( nMG=0; nMG<MutaGen.size(); nMG++ )  {
        if ( ! *pB )
            break;
        if ( sscanf(pB, "\t%lf", &portion)==0 )    {
            printf ( "Empty rangeValue for %s '%s' ... ",
                    MutaGen[nMG].MGname, Buff);
            return -1;
        }
        pB++;   //// '\t'
        vRanges.push_back(portion);
        while ( *pB && *pB != '\t' ) pB++;
    }
    if ( vRanges.size() != MutaGen.size() )    {
        printf("Mismatch of Mutagen count=%d : must be = %d '%s' ... ",
               nMG, (int)MutaGen.size(), Buff);
        return -1;
    }
    
    return typ;
}
/////////////////////////////////////////////////////////////////////////

int loadCodRanges(FILE *f_CAT, char *fName)
{
    char RecTypes[] =  "1Coding_0Noncoding_";
    char Buff[512];
    int rTyp1, rTyp2;
    vector < double > vRanges;
    
    char *pB = fgets(Buff, sizeof(Buff)-2, f_CAT);
    if ( (rTyp1=parsRec(pB, RecTypes, vRanges)) < 0 )   {
        printf ("%s\n", fName);
        return -1;
    }
    for ( int n=0; n<MutaGen.size(); n++ )
        MutaGen[n].fract_cod[rTyp1] = vRanges[n];
//-------
    pB = fgets(Buff, sizeof(Buff)-2, f_CAT);
    if ( (rTyp2=parsRec(pB, RecTypes, vRanges)) < 0 )   {
        printf ("%s\n", fName);
        return -1;
    }
    if ( rTyp2==rTyp1)  {
        printf ("%s :: redefined line : '%s'\n", fName, Buff);
        return -1;
    }
    for ( int n=0; n<MutaGen.size(); n++ )  {
        MutaGen[n].fract_cod[rTyp2] = vRanges[n];
        MutaGen[n].fract_cod[1] += MutaGen[n].fract_cod[0];
    }
    for ( int n=0; n<MutaGen.size(); n++ )  {
        if ( MutaGen[n].fract_cod[1] < DOUBLE_ONE )
            printf ( "%s :: Check distribution of shares for %s (sum<1)\n",
                    fName, MutaGen[n].MGname );
        MutaGen[n].fract_cod[1] = (double) 1;
    }
    
    return 0;
}
/////////////////////////////////////////////////////////////////////////

int loadGenRanges(FILE *f_CAT, char *fName)
{
    char RecTypes[] =  "1Genes_0Intergenes_";
    char Buff[512];
    int rTyp1, rTyp2;
    vector < double > vRanges;
    
    char *pB = fgets(Buff, sizeof(Buff)-2, f_CAT);
    if ( (rTyp1=parsRec(pB, RecTypes, vRanges)) < 0 )   {
        printf ("%s\n", fName);
        return -1;
    }
    for ( int n=0; n<MutaGen.size(); n++ )
        MutaGen[n].fract_gen[rTyp1] = vRanges[n];
//-------
    pB = fgets(Buff, sizeof(Buff)-2, f_CAT);
    if ( (rTyp2=parsRec(pB, RecTypes, vRanges)) < 0 )   {
        printf ("%s\n", fName);
        return -1;
    }
    if ( rTyp2==rTyp1)  {
        printf ("%s :: redefined line : '%s'\n", fName, Buff);
        return -1;
    }
    for ( int n=0; n<MutaGen.size(); n++ )  {
        MutaGen[n].fract_gen[rTyp2] = vRanges[n];
        MutaGen[n].fract_gen[1] += MutaGen[n].fract_gen[0];
    }
    for ( int n=0; n<MutaGen.size(); n++ )  {
        if ( MutaGen[n].fract_gen[1] < DOUBLE_ONE )
            printf ( "%s :: Check distribution of shares for %s (sum<1)\n",
                    fName, MutaGen[n].MGname );
        MutaGen[n].fract_gen[1] = (double) 1;
    }
    
    return 0;
}
/////////////////////////////////////////////////////////////////////////

int loadStrandRanges(FILE *f_CAT, char *fName)
{
    char RecTypes[] =  "1Leading_2Lagging_";
    char Buff[512];
    int rTyp1, rTyp2;
    vector < double > vRanges;
    
    char *pB = fgets(Buff, sizeof(Buff)-2, f_CAT);
    if ( (rTyp1=parsRec(pB, RecTypes, vRanges)) < 0 )   {
        printf ("%s\n", fName);
        return -1;
    }
    for ( int n=0; n<MutaGen.size(); n++ )
        MutaGen[n].fract_strand[rTyp1] = vRanges[n];
    //-------
    pB = fgets(Buff, sizeof(Buff)-2, f_CAT);
    if ( (rTyp2=parsRec(pB, RecTypes, vRanges)) < 0 )   {
        printf ("%s\n", fName);
        return -1;
    }
    if ( rTyp2==rTyp1)  {
        printf ("%s :: redefined line : '%s'\n", fName, Buff);
        return -1;
    }
    for ( int n=0; n<MutaGen.size(); n++ )  {
        MutaGen[n].fract_strand[rTyp2] = vRanges[n];
        MutaGen[n].fract_strand[1] += MutaGen[n].fract_strand[0];
        MutaGen[n].fract_strand[2] += MutaGen[n].fract_strand[1];
    }
    for ( int n=0; n<MutaGen.size(); n++ )  {
        if ( MutaGen[n].fract_strand[2] < DOUBLE_ONE )
            printf ( "%s :: Check distribution of shares for %s (sum<1)\n",
                    fName, MutaGen[n].MGname );
        MutaGen[n].fract_strand[2] = (double) 1;
    }
    
    return 0;
}
/////////////////////////////////////////////////////////////////////////
/*
int loadCodRanges(FILE *f_CAT, char *fName)
{
    char Buff[512];
    char *pB;
    char rType[64];
    char RecTypes[2][16] = { "Coding", "Noncoding"};
    int nMG, typ;
    double portion;
    
    if ( ! ( pB=fgets(Buff, sizeof(Buff)-2, f_CAT) ) ) {
        printf ("%s :: must be '%s' & '%s' lines at file\n", fName, RecTypes[0], RecTypes[1]);
        return -1;
    }
    sscanf(pB, "%s\t", rType );
    if ( strcmp(rType, RecTypes[0])==0 )
        typ = 0;
    else
        if ( strcmp(rType,  RecTypes[0])==0 )
            typ = 1;
        else {
            printf ("%s :: must be '%s' & '%s' lines at file : '%s'\n", fName, RecTypes[0], RecTypes[1], Buff);
            return -1;
        }
    while (*pB != '\t' ) pB++;
    
    for ( nMG=0; nMG<MutaGen.size(); nMG++ )  {
        if ( ! *pB )
            break;
        if ( sscanf(pB, "\t%lf", &portion)==0 )    {
            printf ( "%s ::  empty rangeValue for %s : '%s'\n",
                    fName, MutaGen[nMG].MGname, Buff);
            return -1;
        }
        MutaGen[nMG].fract_cod[typ] = portion;
        while ( *pB && *pB != '\t' ) pB++;
    }
    if ( nMG != MutaGen.size() )    {
        printf("%s :: Mismatch of Mutagen count=%d : must be = %d '%s'\n",
               fName, nMG, (int)MutaGen.size(), Buff);
        return -1;
    }
//----------------------------------------------------------
    if ( ! ( pB=fgets(Buff, sizeof(Buff)-2, f_CAT) ) ) {
        printf ("%s :: must be '%s' & '%s' lines at file\n", fName, RecTypes[0], RecTypes[1]);
        return -1;
    }
    sscanf(pB, "%s\t", rType);

    int typ1;
    if ( strcmp(rType, RecTypes[0])==0 )
        typ1 = 0;
    else
        if ( strcmp(rType,  RecTypes[1])==0 )
            typ1 = 1;
        else {
            printf ("%s :: must be '%s' & '%s' lines at file : '%s'\n", fName, RecTypes[0], RecTypes[1], Buff);
            return -1;
        }
    if ( typ1 == typ )  {
        printf ("%s :: redefined line '%s' in file : '%s'\n", fName, rType, Buff);
        return -1;
    }
    while ( *pB != '\t' ) pB++;
    
    for ( nMG=0; nMG<MutaGen.size(); nMG++ )  {
        if ( ! *pB )
            break;
        if ( sscanf(pB, "\t%lf", &portion)==0 )    {
            printf ( "%s ::  empty rangeValue for %s : '%s'\n",
                    fName, MutaGen[nMG].MGname, Buff);
            return -1;
        }
        MutaGen[nMG].fract_cod[typ1] = portion;
        while ( *pB && *pB != '\t' ) pB++;
    }
    if ( nMG != MutaGen.size() )    {
        printf("%s :: Mismatch of Mutagen count=%d : must be = %d '%s'\n",
               fName, nMG, (int)MutaGen.size(), Buff);
        return -1;
    }
    
    for ( int n=0; n<MutaGen.size(); n++ )
        MutaGen[n].fract_cod[1] += MutaGen[n].fract_cod[0];
        
    
    return 0;
}
//////////////////////////////////////////////////////////////////////////

int loadGenRanges(FILE *f_CAT, char *fName)
{
    char Buff[512];
    char *pB;
    char rType[64];
    char RecTypes[2][16] = { "Genes", "Intergenes"};
    int nMG, typ;
    double portion;
    
    if ( ! ( pB=fgets(Buff, sizeof(Buff)-2, f_CAT) ) ) {
        printf ("%s :: must be '%s' & '%s' lines at file\n", fName, RecTypes[0], RecTypes[1]);
        return -1;
    }
    sscanf(pB, "%s\t", rType );
    if ( strcmp(rType, RecTypes[0])==0 )
        typ = 0;
    else
        if ( strcmp(rType,  RecTypes[0])==0 )
            typ = 1;
        else {
            printf ("%s :: must be '%s' & '%s' lines at file : '%s'\n", fName, RecTypes[0], RecTypes[1], Buff);
            return -1;
        }
    while ( *pB != '\t' ) pB++;
    
    for ( nMG=0; nMG<MutaGen.size(); nMG++ )  {
        if ( ! *pB )
            break;
        if ( sscanf(pB, "\t%lf", &portion)==0 )    {
            printf ( "%s ::  empty rangeValue for %s : '%s'\n",
                    fName, MutaGen[nMG].MGname, Buff);
            return -1;
        }
        MutaGen[nMG].fract_gen[typ] = portion;
        while ( *pB && *pB != '\t' ) pB++;
    }
    //----------------------------------------------------------
    if ( ! ( pB=fgets(Buff, sizeof(Buff)-2, f_CAT) ) ) {
        printf ("%s :: must be '%s' & '%s' lines at file\n", fName, RecTypes[0], RecTypes[1]);
        return -1;
    }
    sscanf(pB, "%s\t", rType);
    int typ1;
    if ( strcmp(rType, RecTypes[0])==0 )
        typ1 = 0;
    else
        if ( strcmp(rType,  RecTypes[1])==0 )
            typ1 = 1;
        else {
            printf ("%s :: must be '%s' & '%s' lines at file : '%s'\n", fName, RecTypes[0], RecTypes[1], Buff);
            return -1;
        }
    if ( typ1 == typ )  {
        printf ("%s :: redefined line '%s' in file : '%s'\n", fName, rType, Buff);
        return -1;
    }
    while ( *pB != '\t' ) pB++;
    
    for ( nMG=0; nMG<MutaGen.size(); nMG++ )  {
        if ( ! *pB )
            break;
        if ( sscanf(pB, "\t%lf", &portion)==0 )    {
            printf ( "%s ::  empty rangeValue for %s : '%s'\n",
                    fName, MutaGen[nMG].MGname, Buff);
            return -1;
        }
        MutaGen[nMG].fract_gen[typ1] = portion;
        while ( *pB && *pB != '\t' ) pB++;
    }
    for ( int n=0; n<MutaGen.size(); n++ )
        MutaGen[n].fract_gen[1] += MutaGen[n].fract_gen[0];

    return 0;
}
*/
//////////////////////////////////////////////////////////////////////////

