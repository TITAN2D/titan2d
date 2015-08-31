/* 
 * File:   tivector_test.C
 * Author: mikola
 *
 * Created on August 30, 2015, 7:50 PM
 */

#include <stdio.h>
#include <stdlib.h>

#include <iostream>
using namespace std;

#include "../header/tivector.h"

/*
 * 
 */
int main(int argc, char** argv) {

    cout<<"tivector test\n";
    tivector<double> v;
    cout<<"resize and init\n";
    v.resize(10);
    for(int i=0;i<v.size();i++)
        v[i]=i;
    
    for(int i=0;i<v.size();i++)
        cout<<i<<" "<<v[i]<<"\n";
    
    cout<<"insert\n";
    v.insert(5);
    v[5]=55;
    v.insert(11);
    v[11]=11;
    //for(int i=9;i<v.size();i++)
    //    v[i]=i*10;
    for(int i=0;i<v.size();i++)
        cout<<i<<" "<<v[i]<<"\n";
    cout<<"remove\n";
    v.remove(5);
    for(int i=0;i<v.size();i++)
        cout<<i<<" "<<v[i]<<"\n";
    
    v.remove(10);
    for(int i=0;i<v.size();i++)
        cout<<i<<" "<<v[i]<<"\n";
    return (EXIT_SUCCESS);
}

