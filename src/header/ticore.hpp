/*******************************************************************
 * Copyright (C) 2015 University at Buffalo
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
 *******************************************************************
 */

#ifndef SRC_HEADER_TICORE_HPP_
#define SRC_HEADER_TICORE_HPP_

#include <cstddef>

#if defined(__INTEL_COMPILER)
#include <malloc.h>
#else
#include <mm_malloc.h>
#define NO_MM_MALLOC
#endif // defined(__GNUC__)

//!Memory management
#ifdef NO_MM_MALLOC
//this one c++ way but no alignment
#define TI_ALLOC(element_type, number_of_elements) new element_type[number_of_elements]
#define TI_FREE(pointer) if(pointer!=nullptr){delete [] pointer;pointer=nullptr;}
#define TI_ASSUME_ALIGNED(pointer)

#define AlignmentAllocator allocator

#define FREE_VAR_IF_NOT_NULLPTR(variable) if(variable!=nullptr){delete variable;variable=nullptr;}

#else
//this one is with alignment
#define TI_ALIGNMENT 64
#define TI_ALLOC(element_type, number_of_elements) (element_type*)_mm_malloc(sizeof(element_type)*number_of_elements,TI_ALIGNMENT)
#define TI_FREE(pointer) _mm_free(pointer)
#define TI_ASSUME_ALIGNED(pointer) __assume_aligned(pointer,TI_ALIGNMENT)

#define FREE_VAR_IF_NOT_NULLPTR(variable) if(variable!=nullptr){delete variable;variable=nullptr;}



template<typename T, std::size_t N = TI_ALIGNMENT>
class AlignmentAllocator
{
public:
    typedef T value_type;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;

    typedef T * pointer;
    typedef const T * const_pointer;

    typedef T & reference;
    typedef const T & const_reference;

public:
    inline AlignmentAllocator() throw ()
    {
    }

    template<typename T2>
    inline AlignmentAllocator(const AlignmentAllocator<T2, N> &) throw ()
    {
    }

    inline ~AlignmentAllocator() throw ()
    {
    }

    inline pointer adress(reference r)
    {
        return &r;
    }

    inline const_pointer adress(const_reference r) const
    {
        return &r;
    }

    inline pointer allocate(size_type n)
    {
        return (pointer) _mm_malloc(n * sizeof(value_type), N);
    }

    inline void deallocate(pointer p, size_type)
    {
        _mm_free(p);
    }

    inline void construct(pointer p, const value_type & wert)
    {
        new (p) value_type(wert);
    }

    inline void destroy(pointer p)
    {
        p->~value_type();
    }

    inline size_type max_size() const throw ()
    {
        return size_type(-1) / sizeof(value_type);
    }

    template<typename T2>
    struct rebind
    {
        typedef AlignmentAllocator<T2, N> other;
    };

    bool operator!=(const AlignmentAllocator<T, N>& other) const
    {
        return !(*this == other);
    }

    // Returns true if and only if storage allocated from *this
    // can be deallocated from other, and vice versa.
    // Always returns true for stateless allocators.
    bool operator==(const AlignmentAllocator<T, N>& other) const
    {
        return true;
    }
};
#endif

#endif /* SRC_HEADER_TICORE_HPP_ */
