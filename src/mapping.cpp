//
//  mapping.cpp
//  MutaGena
//
//  Created by Gennady on 8/9/25.
//  Copyright © 2025 Gennady. All rights reserved.
//
#include "string"

#include "cmain.h"
#include "mapping.h"
#include "xrosoma.h"
#include "genbase.h"

extern PROGARGS ArgKit;
extern HUGEN GENOM;
extern vector < XROSOMA > vecDNK;

extern FILE *Ftrace;

KEY_MAP keySet;

/////////////////////////////////////////////////////////////////////////

//-------------------------------------
void setKEY_GENES(char *key, int g )
{
    *(key + _iKEY_GENE)  =  (g) ? 'g' : 'i';
    if ( ! g )
        *(key + _iKEY_CODING) = '-';
    return;
}
//---------------------------------
void setKEY_CODING(char *key, int c )
{
    if ( *(key + _iKEY_GENE) == 'g' )
        *(key + _iKEY_CODING)  =  (c) ? 'c' : 'n';
    else
        *(key + _iKEY_CODING) = '-';
    return;
}
//---------------------------------
void setKEY_RT(char *key, int r )
{
    if ( r < 1 || r > 7 )  {
        *(key + _iKEY_RT) = '0';
        return;
    }
    char rt = r + '0';
    *(key + _iKEY_RT) = rt;
    return;
}
//---------------------------------
void setKEY_STRAND(char *key, char *pTag, int invert)
{
    *(key + _iKEY_STRAND) = '-';
    
    if ( GET_LEADING_TAG(pTag) )    {
        *(key + _iKEY_STRAND) = ( invert ) ? '<' : '>';
    }
    
    if ( GET_LAGGING_TAG(pTag) )    {
        *(key + _iKEY_STRAND) = ( invert ) ? '>' : '<';
    }
    return;
}
//---------------------------------
void make_keySet()
{
    char key[MOTKEY_SIZE];
    char Motif[MOTKEY_SIZE];
    char tag;
//    char Nucs[5]  = "ACGT";
//    unsigned char invert;
    int indx;
    KEY_MAP::iterator itK;
 
//--------------------------------------------
//    keySet.clear();
    keySet.reserve(KEY_SET_RESERV);
    indx=0;
    for (int m=0; m<4; m++ )    {
        Motif[1] = getNuc(m); //Nucs[m];
        for (int l=0; l<4; l++ )    {
            Motif[0] = getNuc(l); // Nucs[l];
            for (int r=0; r<4; r++ )    {
                Motif[2] = getNuc(r); //Nucs[r];
                for (int g=1; g>=0; g-- )    {          //"gi"
                    for (int c=1; c>=0; c-- )    {      //"cn-"
                        for ( int s=2; s>=0; s--)   {   // strand = "<>-" (leading, lagging, undef)
                            for (int t=0; t<=RT_TAB_SIZE-1; t++ )    {     //  0 : 7
                                tag = '\0';
                                if ( g )
                                    SET_GEN_TAG(&tag);
                                if ( c )
                                    SET_CODING_TAG(&tag);
                                if ( s==1 )
                                    SET_LEADING_TAG(&tag);
                                if ( s==2 )
                                    SET_LAGGING_TAG(&tag);
                                SET_RT_VAL(&tag, t);
                                formKEY (Motif, &tag, key);
//                                setKEYbyTag (key, &tag, invert);
                                if ( indx > 0 ) {
                                    itK = keySet.find(key);
                                    if ( itK != keySet.end() )
                                        continue;
                                }
                                keySet[key]=indx++;
                            }
                        }
                    }
                }
            }
        }
    }
//-----------------------------------------------
/*    int n=0;
    for ( auto it = keySet.begin(); it != keySet.end(); ++it, n++ )
        fprintf(Ftrace,"%d\t%s\t%d\n", n, it->first.c_str(), it->second);
*/
    return;
}
////////////////////////////////////////////////////////////////////////////////////

void reservIndexSpace()
{
    char key[MOTKEY_SIZE] = "???";
    KEY_MAP::iterator itK;
    
    GENOM.vIndexG.reserve(keySet.size());
    vecDNK[0].vX_Data.reserve(keySet.size());
    
    for (int t=0; t<keySet.size(); t++ )    {
        GENOM.vIndexG.push_back(INDEX_BAS(key));
        vecDNK[0].vX_Data.push_back(MOTIF_vPOS(key));        //prepare dataSet for each index
    }
    
    for ( itK = keySet.begin(); itK != keySet.end(); itK++ )
        strcpy ( GENOM.vIndexG[itK->second].mKey, itK->first.c_str() );
    
//    GENOM.sizeBufDataG = vecDNK[0].bufDsizeX;        ==0 now
//    GENOM.BufDataG =new INV_POS [GENOM.sizeBufDataG];
    
    return;
}
/////////////////////////////////////////////////////////////////////////

void rewindMotif(const char *pKey, char *pCmpl)
{
    //    strncpy(pCmpl, pKey, 3);
    *pCmpl     = getCmpl_Nuc( *(pKey+2) );
    *(pCmpl+1) = getCmpl_Nuc( *(pKey+1) );
    *(pCmpl+2) = getCmpl_Nuc( *pKey );
    *(pCmpl+3) = '\0';
    
    return;
}
/////////////////////////////////////////////////////////////////////////

int xtract_Motif(const char *pXr, char *key)
{
// RC:  < 0 invalid key; ==0 orig; ==1 inverted
    int rc=7;
    
    strncpy(key, pXr, 3);
    if ( *pXr=='N' || *(pXr+1)=='N' || *(pXr+2)=='N' )
        return rc;
    

    switch ( *(pXr+1) )  {
        case 'G':
        case 'A':
            rewindMotif(pXr,key);
//            *key     = getCmpl_Nuc( *(pXr+2) );
//            *(key+1) = getCmpl_Nuc( *(pXr+1) );
//            *(key+2) = getCmpl_Nuc( *pXr );
            rc = 1;
            break;
        default: //case 'C': case 'T':
 //           strncpy(key, pXr, 3);
//            *key = *pXr;
//            *(key+1) = *(pXr+1);
//            *(key+2) = *(pXr+2);
            rc = 0;
            break;
    }
    
    return rc;
}
/////////////////////////////////////////////////////////////////////////

int formKEY (const char *motif, char *pTag, char *Key ) 
{
    int coding;
    int invert = xtract_Motif(motif, Key);
    
    for ( int n = 3; n<MOTKEY_SIZE; n++)
        *(Key+n) = '\0';
    
    setKEY_GENES(Key, GET_GEN_TAG(pTag) );
    if ( GET_GEN_TAG(pTag) )    {
        coding =    ( ! invert &&   GET_CODING_TAG(pTag) ) ||
                    (   invert && ! GET_CODING_TAG(pTag      ) );
        setKEY_CODING(Key, coding );
    }
    setKEY_RT(Key, GET_RT_VAL(pTag) );
    setKEY_STRAND(Key, pTag, invert);
    
    return invert;
}
/////////////////////////////////////////////////////////////////////////
/*
void setKEYbyTag (char *key, char *pTag, unsigned char invertM)
{
    int coding;
    
    for ( int n = 3; n<MOTKEY_SIZE; n++)
        *(key+n) = '\0';
    setKEY_GENES(key, GET_GEN_TAG(pTag) );
    if ( GET_GEN_TAG(pTag) )    {
        coding =    ( ! invertM &&   GET_CODING_TAG(pTag) ) ||
                    (   invertM && ! GET_CODING_TAG(pTag) );
        setKEY_CODING(key, coding );
    }
    setKEY_RT(key, GET_RT_VAL(pTag) );
    return;
} */
/////////////////////////////////////////////////////////////////////////
 
int XROSOMA:: motivator( )
{
    unsigned int posMot;
    unsigned char  invert;
    KEY_MAP::iterator itM ;
    char *pTag;
    
    char *pXr = Xbody;
    int cntM = 0;
    char Key[MOTKEY_SIZE];
    char Motif[MOTKEY_SIZE];
//    memset(key, '\0', MOTKEY_SIZE);
    
    printf( "%s.Motivator", XroID.c_str());
    clock_t start = clock();
    
    for ( int n=0; n<GENOM.vIndexG.size(); n++ )   {
        GENOM.vIndexG[n].amtPos = 0;
        vX_Data[n].vPosM.clear();
    }
    while ( *(pXr+2) )  {
        while ( *pXr=='N' ) pXr++;
        if ( ! *pXr )
            break;
        strncpy(Motif, pXr, 3); Motif[3] = '\0';
        if ( strchr(Motif, 'N') ) {
            pXr++;
            continue;
        }

        posMot = (unsigned int)( pXr+1 - Xbody );       // +1 ---> middl_NUC
        pTag = Xtag + posMot;
        invert = formKEY (Motif, pTag, Key);
        
        itM = keySet.find(Key);
        if ( itM != keySet.end() )  {
            int sec = itM->second;              /// ttttttteeeeeeeeesssssssstttttttt
            string  fst = itM->first;              /// ttttttteeeeeeeeesssssssstttttttt
            vX_Data[itM->second].vPosM.push_back(INV_POS (invert,posMot));
            cntM++;
//            pXr++;
        }
        else    {
            printf( "KEY [%s] not found at MAP\n", Key);
            return -1;
        }
        pXr++;
    }
//  -------------------------
//  building Index_X now
    
    vIndex_X = GENOM.vIndexG;
    for ( int n=0; n<GENOM.vIndexG.size(); n++ )    {
        vIndex_X[n].amtPos = (int)vX_Data[n].vPosM.size();
        if ( bufDsizeX < vIndex_X[n].amtPos )
            bufDsizeX = vIndex_X[n].amtPos;
    }
    
    clock_t finish = clock();
    float duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "\tXsize=%u selected %d motif.Pos dT=%5.2fsec\n", Xsize, cntM, duration );
    
    return cntM;
}
/////////////////////////////////////////////////////////////////////////

int checkKeySet(vector <INDEX_BAS> &vIndx)
{
    printf( "checkKeySet()\n");
    if ( vIndx.size() != keySet.size() )    {
        printf ( "checkKeySet() : Difference keySet.size=%ld and Index.size=%ld\n",
                keySet.size(), vIndx.size() );
        return 0;
    }

    KEY_MAP ::iterator iterK;
    for ( int n=0; n<keySet.size(); n++ )   {
//        char tKey[MOTKEY_SIZE];                         ////tttttttttt
//        strcpy(tKey, vIndx[n].mKey);                    ////ttttttttttmmm
        iterK = keySet.find( vIndx[n].mKey );
        if ( iterK == keySet.end() )    {
            printf ( "checkKeySet() : Motif '%s' NOT found at keySet\n", vIndx[n].mKey );
            return 0;
        }
        if ( iterK->second != n )   {
            printf ( "checkKeySet() : Mism.Motif '%s' keySet->%d Index=%d\n", vIndx[n].mKey,
                    iterK->second, n);
            return 0;
        }
    }

    return 1;
}
/////////////////////////////////////////////////////////////////////////



