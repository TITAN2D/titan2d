/* 
 * File:   tivector.h
 * Author: mikola
 *
 * Created on August 30, 2015, 7:46 PM
 */

#ifndef TIVECTOR_H
#define	TIVECTOR_H

#include <string.h>
#include <assert.h>

typedef int tisize_t;
typedef int ti_ndx_t;

constexpr ti_ndx_t ti_ndx_unknown=-1;
constexpr ti_ndx_t ti_ndx_doesnt_exist=-2;
inline bool ti_ndx_not_negative(const ti_ndx_t &ndx){return (ndx>=0);}
inline bool ti_ndx_negative(const ti_ndx_t &ndx){return (ndx<0);}

//#define LOW_MEMORY

class tivector_interface
{
public:
    
};

#ifdef LOW_MEMORY
//low memory footprint implementation
//used to save memory and identify deallocated pointers usage

template<typename T>
class tivector
{
protected:
public:
    tisize_t size_;
    tisize_t reserved_size_;
    T *array_;
    tisize_t size_old_;
    tisize_t reserved_size_old_;
    T *array_old_;
public:
    tivector(){
        size_=0;
        size_old_=0;
        reserved_size_=1;
        reserved_size_old_=reserved_size_;
        
        array_=new T[reserved_size_];
        array_old_=nullptr;
        //for(size_t i=0;i<reserved_size_;i++)array_[i]=0;
    }
    ~tivector(){
        if(array_!=nullptr)delete [] array_;
        if(array_old_!=nullptr)delete [] array_old_;
    }
    const tisize_t& size() const {return size_;}

public:
    void resize(size_t new_size,bool movecontent=true)
    {
        if(new_size>reserved_size_){
            //move to array_old_
            //delete [] array_old_;
            size_old_=size_;
            reserved_size_old_=reserved_size_;
            array_old_=array_;
            //allocate new
            reserved_size_=new_size;//2*(new_size/reserved_size_)*reserved_size_;
            array_=new T[reserved_size_];
            //for(tisize_t i=0;i<size_old_;++i)
            //    array_[i]=array_old_[i];
            if(movecontent)
                memcpy(array_, array_old_, sizeof(T)*size_old_);

            delete [] array_old_;
            array_old_=nullptr;
        }
        size_=new_size;
    }
    void insert(ti_ndx_t pos){
        assert(pos<=size_);

        resize(size_+1);

        for(tisize_t i=size_-1;i>pos;--i)
            array_[i]=array_[i-1];
        //array_[pos]=NULL;
    }
    void remove(ti_ndx_t pos){
        assert(pos<size_);

        for(tisize_t i=pos;i<size_-1;++i)
            array_[i]=array_[i+1];
        resize(size_-1);
    }
    ti_ndx_t push_back()
    {
        resize(size_+1);
        return size_-1;
    }
    
    void reorder(ti_ndx_t *new_order, tisize_t new_size)
    {
    	array_old_=new T[size_];
    	            //for(tisize_t i=0;i<size_old_;++i)
    	            //    array_[i]=array_old_[i];
		memcpy(array_old_, array_, sizeof(T)*size_);


        size_=new_size;
        
        for(ti_ndx_t i=0;i<new_size;i++)
        {
            array_[i]=array_old_[new_order[i]];
        }

        delete [] array_old_;
		array_old_=nullptr;
    }
    void set(const T& value)
    {
        for(ti_ndx_t i=0;i<size_;i++)
        {
            array_[i]=value;
        }
    }
    
    T& operator[](ti_ndx_t i){ return *(array_ + i);}
    const T& operator[](ti_ndx_t i) const { return *(array_ + i);}
};

#else
template<typename T>
class tivector
{
protected:
public:
    tisize_t size_;
    tisize_t reserved_size_;
    T *array_;
    tisize_t size_old_;
    tisize_t reserved_size_old_;
    T *array_old_;
public:
    tivector(tisize_t reserved_size=1000000){
        size_=0;
        size_old_=0;
        reserved_size_=reserved_size;
        reserved_size_old_=reserved_size_;
        
        array_=new T[reserved_size_];
        array_old_=new T[reserved_size_old_];
        //for(size_t i=0;i<reserved_size_;i++)array_[i]=0;
    }
    ~tivector(){
        if(array_!=nullptr)delete [] array_;
        if(array_old_!=nullptr)delete [] array_old_;
    }
    const tisize_t& size() const {return size_;}
public:
    void swap_arrays(){
        tisize_t tmp_size_=size_;
        tisize_t tmp_reserved_size_=reserved_size_;
        T *tmp_array_=array_;
        
        array_=array_old_;
        size_=size_old_;
        reserved_size_=reserved_size_old_;
        
        array_old_=tmp_array_;
        size_old_=tmp_size_;
        reserved_size_old_=tmp_reserved_size_;
    }
public:
    void resize(size_t new_size,bool movecontent=true)
    {
        if(new_size>reserved_size_){
            //move to array_old_
            delete [] array_old_;
            size_old_=size_;
            reserved_size_old_=reserved_size_;
            array_old_=array_;
            //allocate new
            reserved_size_=2*(new_size/reserved_size_)*reserved_size_;
            array_=new T[reserved_size_];
            //for(tisize_t i=0;i<size_old_;++i)
            //    array_[i]=array_old_[i];
            if(movecontent)
                memcpy(array_, array_old_, sizeof(T)*size_old_);
        }
        size_=new_size;
    }
    void insert(ti_ndx_t pos){
        assert(pos<=size_);
        swap_arrays();
        //reallocate array_ if needed
        if(size_old_+1>reserved_size_){
            //move to array_old_
            delete [] array_;
            //allocate new
            size_=size_old_+1;
            reserved_size_=2*(size_/reserved_size_)*reserved_size_;
            array_=new T[reserved_size_];
        }
        //size_=size_+1;
        //for(tisize_t i=0;i<pos;++i)
        //    array_[i]=array_old_[i];
        memcpy(array_, array_old_, sizeof(T)*pos);
        //for(tisize_t i=pos+1;i<size_;++i)
        //    array_[i]=array_old_[i-1];
        memcpy(array_+pos+1, array_old_+pos, sizeof(T)*(size_old_-pos));
    }
    void remove(ti_ndx_t pos){
        assert(pos<size_);
        swap_arrays();
        //reallocate array_ if needed
        if(size_old_-1>reserved_size_){
            //move to array_old_
            delete [] array_;
            //allocate new
            size_=size_old_-1;
            reserved_size_=2*(size_/reserved_size_)*reserved_size_;
            array_=new T[reserved_size_];
        }
        size_=size_old_-1;
        //size_=size_+1;
        //for(tisize_t i=0;i<pos;++i)
        //    array_[i]=array_old_[i];
        memcpy(array_, array_old_, sizeof(T)*pos);
        //for(tisize_t i=pos+1;i<size_;++i)
        //    array_[i]=array_old_[i-1];
        memcpy(array_+pos, array_old_+pos+1, sizeof(T)*(size_-pos));
    }
    ti_ndx_t push_back()
    {
        resize(size_+1);
        return size_-1;
    }
    
    void reorder(ti_ndx_t *new_order, tisize_t new_size)
    {
        swap_arrays();
        if(new_size>reserved_size_){
            //move to array_old_
            delete [] array_;
            //allocate new
            reserved_size_=2*(new_size/reserved_size_)*reserved_size_;
            array_=new T[reserved_size_];
        }
        size_=new_size;
        
        for(ti_ndx_t i=0;i<new_size;i++)
        {
            array_[i]=array_old_[new_order[i]];
        }
    }
    void set(const T& value)
    {
        for(ti_ndx_t i=0;i<size_;i++)
        {
            array_[i]=value;
        }
    }
    
    T& operator[](ti_ndx_t i){ return *(array_ + i);}
    const T& operator[](ti_ndx_t i) const { return *(array_ + i);}
};
#endif
#endif	/* TIVECTOR_H */

