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
 *******************************************************************
 * $Id: boundary.C 206 2009-01-26 17:32:10Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <iostream>
#include <fstream>
using namespace std;

#include <stdio.h>
#include "boundary.h"
#include "useful_lib.h"
#include "../header/FileFormat.h"


Boundary::Boundary(){}


void Boundary::setparameters(Node* n, double xv, double yv, int t)
{

  node=n;
  x_value=xv;
  y_value=yv;
  type=t;

}

void Boundary::write_b_data_bin(FILE *fp){
  fwriteI(fp,type);
  fwriteI(fp,node->key[0]);
  fwriteI(fp,node->key[1]);
#ifdef WRITEDOUBLEASFLOAT
  fwriteF(fp,x_value);
  fwriteF(fp,y_value);
#else
  fwriteD(fp,x_value);
  fwriteD(fp,y_value);
#endif
  return;
}

void Boundary::write_b_data(ofstream* out)
{
  *out<<type<<' '<<node->key[0]<<' '<<node->key[1]<<' '<<x_value<<' '<<y_value<<' '<<'\n';

}
