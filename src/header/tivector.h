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
        reserved_size_=10;
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
protected:
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
    void resize(size_t new_size)
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
            memcpy(array_, array_old_, sizeof(T)*size_old_);
        }
        size_=new_size;
    }
    void insert(size_t pos){
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
    
    T& operator[](tisize_t i){ return *(array_ + i);}
    const T& operator[](tisize_t i) const { return *(array_ + i);}
};

#endif	/* TIVECTOR_H */

