/*******************************************************************
 * Copyright (C) 2003 University at Buffalo
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Author:
 * Description:
 *
 *******************************************************************/
#ifndef SFC_H_
#define SFC_H_

#include "constant.h"
#include <ostream>

/*typedef uint64_t sfc_key;

constexpr uint64_t key_zero=0;
constexpr uint64_t oldkey_1mask=~0U;

#define SET_OLDKEY(oldkey, newkey) {oldkey[1]=(oldkey_1mask&newkey);oldkey[0]=(newkey>>32)}

#define GET_OLDKEY_PART(newkey,part) (part==0)?(newkey>>32):(oldkey_1mask&newkey)

#define SET_NEWKEY(oldkey) ((((uint64_t) oldkey[0]) << 32) | oldkey[1])*/

#if USE_ARRAY_SFC_KEY
class SFC_Key
{
public:
    SFC_Key(){
        //for(int i=0;i<KEYLENGTH;++i)key[i]=0;
    }
    //!copy constructor at this point implicit copy constructor should work fine
    SFC_Key(const SFC_Key& other ){
        for(int i=0;i<KEYLENGTH;++i)key[i]=other.key[i];
    }
    SFC_Key(const int & ikey)
    {
        key[0]=ikey;
        for(int i=1;i<KEYLENGTH;++i)
            key[i]=0;
    }


    SFC_Key& operator=(const int & ikey)
    {
        key[0]=ikey;
        for(int i=1;i<KEYLENGTH;++i)
            key[i]=0;
        return *this;
    }
    unsigned key[KEYLENGTH];

    //uint64_t key conversion
    //uint64_t get_ukey(){return (((uint64_t) key[0]) << 32) | key[1];}
};
inline bool operator >(const SFC_Key &L, const SFC_Key &R)
{
    for(int i=0;i<KEYLENGTH;++i){
        if(L.key[i]<R.key[i])return false;
        if(L.key[i]>R.key[i])return true;
    }
    return false;
}
inline bool operator <(const SFC_Key &L, const SFC_Key &R)
{
    for(int i=0;i<KEYLENGTH;++i){
        if(L.key[i]>R.key[i])return false;
        if(L.key[i]<R.key[i])return true;
    }
    return false;
}
inline bool operator ==(const SFC_Key &L, const SFC_Key &R)
{
    for(int i=0;i<KEYLENGTH;++i){
        if(L.key[i]!=R.key[i])return false;
    }
    return true;
}
inline bool operator !=(const SFC_Key &L, const SFC_Key &R)
{
    for(int i=0;i<KEYLENGTH;++i){
        if(L.key[i]!=R.key[i])return true;
    }
    return false;
}
inline std::ostream& operator<< (std::ostream& os, const SFC_Key& obj)
{
    for(int i=0;i<KEYLENGTH;++i){
        os<<obj.key[i];
        if(i<KEYLENGTH-1)os<<" ";
    }
    return os;
}
inline void fprintf_sfc_key(FILE *fout, const SFC_Key& obj)
{
    for(int i=0;i<KEYLENGTH;++i){
        fprintf(fout,"%u",obj.key[i]);
        if(i<KEYLENGTH-1)fprintf(fout," ");
    }

}

//constexpr int sfc_key_zero=0;
//#define sfc_key_null SFC_Key(0)
extern SFC_Key sfc_key_zero;
extern SFC_Key sfc_key_null;

inline uint64_t get_ukey_from_sfc_key(const SFC_Key& key){return (((uint64_t) key.key[0]) << 32) | key.key[1];};


#define SET_OLDKEY(oldkey, newkey) {oldkey[0]=newkey.key[0];oldkey[1]=newkey.key[1];}
#define SET_NEWKEY(newkey, oldkey) {newkey.key[0]=oldkey[0];newkey.key[1]=oldkey[1];}

inline SFC_Key sfc_key_from_oldkey(const unsigned ikey[] ){
    SFC_Key key;
    for(int i=0;i<KEYLENGTH;++i)key.key[i]=ikey[i];
    return key;
}

//#define GET_OLDKEY_PART(newkey,part) (part==0)?(newkey>>32):(oldkey_1mask&newkey)

//back compatibility with elements serialization as soon as it will be improved it should be gone
inline void sfc_key_write_to_space(const SFC_Key &key, unsigned writespace[],int & pos)
{
    for(int i=0;i<KEYLENGTH;++i){
        writespace[pos++] = key.key[i];
    }
}
//back compatibility with elements serialization as soon as it will be improved it should be gone
inline void sfc_key_read_from_space(SFC_Key &key, const unsigned readspace[],int & pos)
{
    for(int i=0;i<KEYLENGTH;++i){
        key.key[i]=readspace[pos++];
    }
}
inline SFC_Key sfc_key_read_from_space(const unsigned readspace[],int & pos)
{
    SFC_Key tmp;
    SET_NEWKEY(tmp,(readspace+pos));
    pos+=KEYLENGTH;
    return tmp;
}
#else
#include <stdio.h>
#include <stdint.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
typedef uint64_t SFC_Key;

inline void fprintf_sfc_key(FILE *fout, const SFC_Key& obj)
{
    fprintf(fout,"%" PRIu64 "\n",obj);
}

extern SFC_Key sfc_key_zero;
extern SFC_Key sfc_key_null;

inline uint64_t get_ukey_from_sfc_key(const SFC_Key& key)
{
    return key;
};

constexpr uint64_t oldkey_1mask=~0U;

#define SET_OLDKEY(oldkey, newkey) {oldkey[1]=(oldkey_1mask&newkey);oldkey[0]=(newkey>>32);}
#define SET_NEWKEY(newkey,oldkey) {newkey=(((uint64_t) oldkey[0]) << 32) | oldkey[1];}

inline SFC_Key sfc_key_from_oldkey(const unsigned ikey[] ){
    SFC_Key key;
    SET_NEWKEY(key,ikey)
    return key;
}

//back compatibility with elements serialization as soon as it will be improved it should be gone
inline void sfc_key_write_to_space(const SFC_Key &key, unsigned writespace[],int & pos)
{
    unsigned oldkey[KEYLENGTH];
    SET_OLDKEY(oldkey,key);
    for(int i=0;i<KEYLENGTH;++i){
        writespace[pos++] = oldkey[i];
    }
    
}
//back compatibility with elements serialization as soon as it will be improved it should be gone
inline void sfc_key_read_from_space(SFC_Key &key, const unsigned readspace[],int & pos)
{
    unsigned oldkey[KEYLENGTH];
    for(int i=0;i<KEYLENGTH;++i){
        oldkey[i]=readspace[pos++];
    }
    SET_NEWKEY(key,oldkey);
}
inline SFC_Key sfc_key_read_from_space(const unsigned readspace[],int & pos)
{
    SFC_Key tmp;
    SET_NEWKEY(tmp,(readspace+pos));
    pos+=KEYLENGTH;
    return tmp;
}

#endif
/*class SFC_Key_Hash
{
public:
    SFC_Key_Hash(){};
    ~SFC_Key_Hash(){};

};*/
#endif /* SFC_H_ */
