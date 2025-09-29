//
//  genbase.cpp
//  MutaGena
//
//  Created by Gennady on 9/1/25.
//  Copyright © 2025 Gennady. All rights reserved.
//

#include "cmain.h"
#include "mapping.h"
#include "xrosoma.h"
#include "genbase.h"

extern PROGARGS ArgKit;
extern KEY_MAP keySet;
extern HUGEN GENOM;
extern vector < XROSOMA > vecDNK;

extern FILE *Ftrace;

/////////////////////////////////////////////////////////////////////////

bool lesser_POS (  const INV_POS &x1, const INV_POS &x2 )
{
    return x1.second < x2.second;
}
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

int XROSOMA:: writeBas_X ( )
{
    fprintf(Ftrace, "\t%s.Make tempBase\n", XroID.c_str());
    
    XtempFilePath = ArgKit.OUTdir + "temp_" + XroID + ".bin";
    XtempFile = fopen(XtempFilePath.c_str(), "w+b");
    if ( ! XtempFile )   {
        printf ( "Cannt OPEN  '%s'\n", XtempFilePath.c_str() );
        return -1;
    }
    
    //    --------------------   HEADER  ------------------------------------
    //    offset <индекс_хромосом>            [long]
    //    offset <индекс_Motif>               [long]
    //    offset <индекс_Footer>              [long]
    long offs_indX=0L;
    long offs_indM=0L;
    long offs_Foot=0L;
    
    fwrite(&offs_indX, sizeof(long), 1, XtempFile);
    fwrite(&offs_indM, sizeof(long), 1, XtempFile);
    fwrite(&offs_Foot, sizeof(long), 1, XtempFile);
//  ----------------------------------------------------------------------
//     ---> <индекс_хромосом>:
//     SKIPPED for Bas_X
//  ----------------------------------------------------------------------
//    ---> <индекс_Ключей>:
//    keySet.size()                   [int]
//   {   Motif_ID                     [MOTKEY_SIZE]      char[]
//        число позиций               [int]
//        offset позиций              [long]
//    } [число Motif_ов]
    
    offs_indM = ftell(XtempFile);
    int cntM = (int)keySet.size();
    
    fwrite(&cntM, sizeof(int), 1, XtempFile);
    
    INDEX_BAS *pIndex = vIndex_X.data();
    size_t rc = fwrite(pIndex, sizeof(INDEX_BAS), cntM, XtempFile); //// vIndex_X[n].offsetF don't set now!!!!!!!!!!!!
    if ( rc != cntM )   {
        printf("writeBas_X ( ):: IndexWrite ERROR rc=%ld, size=%d", rc, cntM);
    }

//    --------------------     DATA ------------------------------------
//    ---> <позиции для Motif_1>:
//         #поз. в Геноме для Motif_1      [unsigned int] * число позиций
    INV_POS *pIpos;
    for ( int n=0; n<keySet.size(); n++ ) {
        vIndex_X[n].offsetF = ftell(XtempFile);
        
        for ( int p=0; p<vX_Data[n].vPosM.size(); p++ )
            vX_Data[n].vPosM[p].second += XstartPos;
        pIpos = vX_Data[n].vPosM.data();
        fwrite(pIpos, sizeof(INV_POS), vX_Data[n].vPosM.size(), XtempFile);
    }
//  --------------------   FOOTER  ------------------------------------
    offs_Foot = ftell(XtempFile);  //==EOF
//  SKIPPED for Bas_X
//  ----------------------------------------------------------------------
//====================================================================================
//
// ---> rewrite  <индекс_Ключей>:
    if ( fseek( XtempFile, offs_indM+sizeof(int), SEEK_SET) ) {
        printf ( "Fseek failed\n" );
        return -1;
    };
    pIndex = vIndex_X.data();
    rc = fwrite(pIndex, sizeof(INDEX_BAS), cntM, XtempFile);
    if ( rc != cntM )   {
        printf("writeBas_X ( ):: Index_RE_Write ERROR rc=%ld, size=%d", rc, cntM);
    }
//  ---> rewrite   HEADER  -------------------------------------------------
    if ( fseek( XtempFile, sizeof(long), SEEK_SET) ) {
        printf ( "Fseek failed\n" );
        return -1;
    };
    fwrite(&offs_indM, sizeof(long), 1, XtempFile);
    
    fclose(XtempFile);
    XtempFile = fopen(XtempFilePath.c_str(), "rb");
    if ( ! XtempFile )   {
        printf ( "Cannt reOPEN  '%s'\n", XtempFilePath.c_str() );
        return -1;
    }
// -----------test writing -------------------                                  // tttttttttttttttttttttttt
/*
    long offsM;
    int cntInd;
    fseek( XtempFile, sizeof(long), SEEK_SET);
    fread(&offsM, sizeof(long), 1, XtempFile);
    fseek( XtempFile, offsM, SEEK_SET);
    fread(&cntInd, sizeof(int), 1, XtempFile);
    
    INDEX_BAS Indx_1;
    for ( int n=0; n<10; n++ )  {
        fread(&Indx_1, sizeof(INDEX_BAS),  1, XtempFile);
        printf ( "%d. '%s' %d %ld\n", n, Indx_1.mKey, Indx_1.amtPos, Indx_1.offsetF);
    }
*/
    return 0;
}
/////////////////////////////////////////////////////////////////////////

int HUGEN:: writeGenBase ( )
{
    string BasPath;
    char XroId[XRO_ID_SIZE];
    int RetC;
    
    printf("\nCreating GENOME Base for '%s' ... \n", GenName.c_str());
        clock_t start = clock();
    
    BasPath = ArgKit.OUTdir + GenName +"_Base.bin";
    fprintf(Ftrace, "Create GENOME Base for '%s'\n", GenName.c_str());
    GenBasFile = fopen(BasPath.c_str(), "w+b");
    if ( ! GenBasFile )   {
        printf ( "Cannt OPEN  '%s'\n", BasPath.c_str() );
        return -1;
    }
//    --------------------   HEADER  ------------------------------------
    long offs_indX=3*sizeof(long);
    long offs_indM=0L;
    long offs_Foot=0L;
    
    fwrite(&offs_indX, sizeof(long), 1, GenBasFile);
    fwrite(&offs_indM, sizeof(long), 1, GenBasFile);
    fwrite(&offs_Foot, sizeof(long), 1, GenBasFile);
//  ----------------------------------------------------------------------
//     ---> <индекс_хромосом>:
    offs_indX = ftell(GenBasFile);
    int nXr = (int)vecDNK.size();
    fwrite(&nXr, sizeof(int), 1, GenBasFile);
    for ( nXr=0; nXr<vecDNK.size(); nXr++ ) {
        strcpy(XroId, vecDNK[nXr].XroID.c_str());
        fwrite(XroId, 1, XRO_ID_SIZE, GenBasFile);
        fwrite(&vecDNK[nXr].XstartPos, sizeof(int), 1, GenBasFile);
        fwrite(&vecDNK[nXr].XstopPos, sizeof(int), 1, GenBasFile);
        fwrite(&vecDNK[nXr].Xsize, sizeof(int), 1, GenBasFile);
    }
    //  ----------------------------------------------------------------------
    // ---> <индекс_Ключей>:
    offs_indM = ftell(GenBasFile);
    int cntM = (int)keySet.size();
    fwrite(&cntM, sizeof(int), 1, GenBasFile);
    INDEX_BAS *pIndex = vIndexG.data();
    fwrite(pIndex, sizeof(INDEX_BAS), cntM, GenBasFile);
//
//  --------------------     DATA ------------------------------------
    for ( int n=0; n<keySet.size(); n++ ) {
        vIndexG[n].offsetF = ftell(GenBasFile); //  ---> <позиции для Ключа_1>:
//        vData.vPosM.clear();
        int cntPos = 0;
        for ( nXr=0; nXr<vecDNK.size(); nXr++ ) {
            RetC = readData(vecDNK[nXr].vIndex_X[n],  BufDataG,  vecDNK[nXr].XtempFile);
            if ( RetC < 0 )
                return -1;
            cntPos += vecDNK[nXr].vIndex_X[n].amtPos;
            RetC = (int)fwrite(BufDataG, sizeof(INV_POS), vecDNK[nXr].vIndex_X[n].amtPos, GenBasFile);
            if ( RetC != vecDNK[nXr].vIndex_X[n].amtPos ) {
                printf("writeGenBase():: %s.fwrite(data)=%d : must be %d\n",
                       vecDNK[nXr].XroID.c_str(), RetC, vecDNK[nXr].vIndex_X[n].amtPos);
            }
                
        }
        vIndexG[n].amtPos = cntPos;
    }
    // --------------------   FOOTER  ------------------------------------------
    offs_Foot = ftell(GenBasFile);
    int sz = (int)GenName.size()+1;
    fwrite(&sz, sizeof(int), 1, GenBasFile);
    fwrite(GenName.c_str(), 1, sz, GenBasFile);
    
    //=======================================================================================
    //
    // ---> rewrite  <индекс_Ключей>:
    if ( fseek( GenBasFile, offs_indM+sizeof(int), SEEK_SET) ) {
        printf ( "Fseek failed\n" );
        return -1;
    };
    pIndex = vIndexG.data();
    RetC = (int)fwrite(pIndex, sizeof(INDEX_BAS), cntM, GenBasFile);
    if ( RetC != cntM )   {
        printf("writeBas_X ( ):: Index_RE_Write ERROR rc=%d, size=%d", RetC, cntM);
    }
    //  ---> rewrite   HEADER  -------------------------------------------------
    if ( fseek( GenBasFile, sizeof(long), SEEK_SET) ) {
        printf ( "Fseek failed\n" );
        return -1;
    };
    fwrite(&offs_indM, sizeof(long), 1, GenBasFile);
    fwrite(&offs_Foot, sizeof(long), 1, GenBasFile);
    
    fclose(GenBasFile);
    GenBasFile = NULL;
    
    clock_t finish = clock();
    float duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "\t\tCreated '%s'  dT=%5.2fsec\n", GenName.c_str(), duration );
    
    return 0;
}
//===============================================================================
//===============================================================================

int HUGEN:: readIndXro ( vector <XROSOMA> &vecX )
{
    int cntX;
    long offs_indX;
    //    INDEX_BAS Indx1;
    
    vecX.clear();
    
    fseek( GenBasFile, 0L, SEEK_SET);
    fread(&offs_indX, sizeof(long), 1, GenBasFile);            // <индекс_хромосом>
    if ( fseek( GenBasFile, offs_indX, SEEK_SET) ) {
        printf ("readIndXro():: INVop.fseek(offs_indX=%ld)\n", offs_indX );
        return -1;
    }
    fread(&cntX, sizeof(int), 1, GenBasFile);
    
    // index read  -----------------------------------------------------------
    char chrName[XRO_ID_SIZE];
    Gsize = 0;
    memset(chrName, '\0', XRO_ID_SIZE);
    for ( int nX=0; nX<cntX; nX++ )    {
        fread(chrName, 1, XRO_ID_SIZE, GenBasFile);
        vecX.push_back(XROSOMA(chrName));
        fread(&vecX[nX].XstartPos, sizeof(int), 1, GenBasFile);
        fread(&vecX[nX].XstopPos, sizeof(int), 1, GenBasFile);
        fread(&vecX[nX].Xsize, sizeof(int), 1, GenBasFile);
        Gsize += vecX[nX].Xsize;
    }
    Gsize -= 1;     // from 0
    return 0;
}
/////////////////////////////////////////////////////////////////////////

int readIndex ( vector <INDEX_BAS> &vIndex, FILE *binFile )
// not needed for XRO, becouse it's already exist           !!!!!!!================================================
// may be for testing only
{
    int cntM;
    long offs_indM;
    INDEX_BAS Indx_1;
    
    fseek( binFile, sizeof(long), SEEK_SET);        // <индекс_хромосом>
    fread(&offs_indM, sizeof(long), 1, binFile);
    if ( fseek( binFile, offs_indM, SEEK_SET) ) {
        printf ("readIndex():: INVop.fseek(offs_indM=%ld)\n", offs_indM );
        return -1;
    }
    fread(&cntM, sizeof(int), 1, binFile);
    if ( cntM != keySet.size() )    {
        printf("readIndex():: ErrREADing(keySet.size): %d  [%d]\n", cntM, (int)keySet.size() );
        return -1;
    }
    vIndex.clear();
    if ( vIndex.capacity() < cntM )
        vIndex.reserve(cntM);
    
    // index read  -----------------------------------------------------------
    int maxAmt=0;
    for ( int n=0; n<cntM; n++ )    {
        fread(&Indx_1, sizeof(INDEX_BAS),  1, binFile);
        vIndex.push_back(Indx_1);
        if ( maxAmt < Indx_1.amtPos )
            maxAmt = Indx_1.amtPos;
    }
    if (  maxAmt > GENOM.sizeBufDataG )    {
        if ( GENOM.sizeBufDataG > 0 )
            delete [] GENOM.BufDataG;
        GENOM.sizeBufDataG = maxAmt;
        GENOM.BufDataG = new INV_POS [maxAmt];
    }
    
    return 0;
}
/////////////////////////////////////////////////////////////////////////

int readData ( INDEX_BAS &KeyIndx, INV_POS *dataP, FILE *binFile )
// read all positions for 'Motif.pKey' into array dataP
{
//    long offs;
//    unsigned int posMotif;
//    char invert;
    
    if ( dataP==NULL || KeyIndx.amtPos > GENOM.sizeBufDataG )   {   // sizeof dataP
        printf ( "readData(%s) nP=%d:: LOW Buffer Size=%d\n", KeyIndx.mKey, KeyIndx.amtPos, GENOM.sizeBufDataG );
        return -1;
    }
//    offs = Motif.offsetF;
    if ( fseek( binFile, KeyIndx.offsetF, SEEK_SET) ) {
        printf ( "readData(%s):: :: INVop.fseek(vIndex.offsetF=%ld) \n", KeyIndx.mKey, KeyIndx.offsetF );
        return -1;
    }
    int rc = (int)fread(dataP, sizeof(INV_POS), KeyIndx.amtPos, binFile);
    if ( rc != KeyIndx.amtPos )   {
        printf ( "readData(%s) nP=%d:: READED=%d\n", KeyIndx.mKey, KeyIndx.amtPos, rc );
        return -1;
    }

    return rc;
}
/////////////////////////////////////////////////////////////////////////
/*
int readMotif ( INDEX_BAS &Motif, MOTIF_vPOS &dataPos, FILE *binFile )
// read all positions for 'Motif.pKey'
// and add it to dataPos.vPosM
{
    long offs;
    unsigned int posMotif;
    char invert;
    
    offs = Motif.offsetF;
    if ( fseek( binFile, offs, SEEK_SET) ) {
        printf ( "readMotif(%s):: :: INVop.fseek(vIndex.offsetF=%ld) \n", Motif.mKey, offs );
        return -1;
    }
    for ( int p=0; p<Motif.amtPos; p++ )  {
//    for ( int p=0; p<10; p++ )  {
        fread(&invert, 1, 1, binFile);
        fread(&posMotif, sizeof(unsigned int), 1, binFile);
        dataPos.vPosM.push_back( INV_POS(invert,posMotif) );
    }

    return 1;
}*/
/////////////////////////////////////////////////////////////////////////

int  getMotif_Pos( const char *Key, const INV_POS PosM,  char *Motif )
{
// not used now
    int indXro=-1;
    unsigned int pos;
    vector < XROSOMA > ::iterator  iterXr;
    
    strcpy(Motif, "???");
    iterXr = lower_bound( vecDNK.begin(), vecDNK.end(), XROSOMA(PosM.second), lesser_XRO);
    if ( iterXr == vecDNK.end() ) {
        printf("getMotif_X(%s, {%d,%u} : not found Xromo\n", Key, PosM.first, PosM.second);
        return -1;
    }
    indXro = (int)(iterXr - vecDNK.begin() );
    pos = PosM.second - vecDNK[indXro].XstartPos;
    char *pXr = vecDNK[indXro].Xbody + pos -1;       // pos==(middl char)
    xtract_Motif(pXr, Motif);
    
    return indXro;
    
}
/////////////////////////////////////////////////////////////////////////

/*
 int XROSOMA:: readIndex_X ( vector <INDEX_BAS> &vIndex )
 // not needed for XRO, becouse it's already exist           !!!!!!!================================================
 // may be for testing only
 {
 int cntM;
 int rc;
 long offs_indM;
 INDEX_BAS Indx;
 //    offset <индекс_хромосом>            [long]
 //    offset <индекс_Motif>               [long]
 
 fseek( XtempFile, sizeof(long), SEEK_SET);
 rc = (int)fread(&offs_indM, sizeof(long), 1, XtempFile);
 if ( fseek( XtempFile, offs_indM, SEEK_SET) ) {
 printf ("readBegIndex_X(%s): Fseek %ld failed\n", XroID.c_str(), offs_indM );
 return -1;
 }
 rc = (int)fread(&cntM, sizeof(int), 1, XtempFile);
 if ( cntM != keySet.size() )    {
 printf("readIndex_X(%s): Error fread(&cntM)=%d  [%d]\n", XroID.c_str(), cntM, (int)keySet.size() );
 return -1;
 }
 vIndex.clear();
 
 // index read  -----------------------------------------------------------
 for ( int n=0; n<cntM; n++ )    {
 fread(Indx.mKey, 1,  MOTKEY_SIZE, XtempFile);
 fread(&Indx.amtPos, sizeof(int), 1, XtempFile);
 fread(&Indx.offsetF, sizeof(long), 1, XtempFile);
 vIndex.push_back(Indx);
 }
 
 return 0;
 }
 /////////////////////////////////////////////////////////////////////////
 */


