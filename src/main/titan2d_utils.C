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
 */


#include "../header/hd5calls.h"
#include "../header/titan2d_utils.h"

H5::EnumType datatypeElementType;

vector<string> coord_names;
vector<string> SinglePhaseVarNames;
vector<string> TwoPhasesVarNames;

void init_TiH5()
{
    //init enums
    datatypeElementType=H5::EnumType(sizeof(ElementType));
    ElementType val;
    val=ElementType::UnknownElementType;
    datatypeElementType.insert("UnknownElementType",&val);
    val=ElementType::SinglePhase;
    datatypeElementType.insert("SinglePhase",&val);
    val=ElementType::TwoPhases;
    datatypeElementType.insert("TwoPhases",&val);

    //init indexes to char* conversions
    coord_names.push_back("X");
    coord_names.push_back("Y");
    coord_names.push_back("Z");

    SinglePhaseVarNames.push_back("h");
    SinglePhaseVarNames.push_back("hVx");
    SinglePhaseVarNames.push_back("hVy");

    TwoPhasesVarNames.push_back("h");
    TwoPhasesVarNames.push_back("h_liq");
    TwoPhasesVarNames.push_back("hVx_sol");
    TwoPhasesVarNames.push_back("hVy_sol");
    TwoPhasesVarNames.push_back("hVx_liq");
    TwoPhasesVarNames.push_back("hVy_liq");
}



