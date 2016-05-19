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
#include "../header/properties.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>


H5::EnumType datatypeElementType;
H5::EnumType datatypeInterface_Capturing_Type;
H5::EnumType datatypePileType;
H5::EnumType datatypeOutLineInitSize;

vector<string> coord_names;
//vector<string> SinglePhaseVarNames;
vector<string> SinglePhaseHeuristicVarNames;
vector<string> SinglePhaseLevelSetVarNames;
vector<string> SinglePhasePhaseFieldVarNames;
vector<string> TwoPhasesVarNames;

void init_TiH5()
{
    //init enums
    ElementType val;
    Interface_Capturing_Type val_;

    datatypeElementType=H5::EnumType(sizeof(ElementType));
    val=ElementType::UnknownElementType;
    datatypeElementType.insert("UnknownElementType",&val);
    val=ElementType::SinglePhase;
    datatypeElementType.insert("SinglePhase",&val);
    val=ElementType::TwoPhases;
    datatypeElementType.insert("TwoPhases",&val);

    datatypeInterface_Capturing_Type=H5::EnumType(sizeof(Interface_Capturing_Type));
    val_=Interface_Capturing_Type::Heuristic;
    datatypeInterface_Capturing_Type.insert("Heuristic",&val_);
    val_=Interface_Capturing_Type::LevelSet;
    datatypeInterface_Capturing_Type.insert("LevelSet",&val_);
    val_=Interface_Capturing_Type::PhaseField;
    datatypeInterface_Capturing_Type.insert("PhaseField",&val_);

    datatypePileType=H5::EnumType(sizeof(PileProps::PileType));
    PileProps::PileType pile_type;
    pile_type=PileProps::PARABALOID;
    datatypePileType.insert("PARABALOID",&pile_type);
    pile_type=PileProps::CYLINDER;
    datatypePileType.insert("CYLINDER",&pile_type);
    pile_type=PileProps::PLANE;
    datatypePileType.insert("PLANE",&pile_type);
    pile_type=PileProps::CASITA;
    datatypePileType.insert("CASITA",&pile_type);
    pile_type=PileProps::POPO;
    datatypePileType.insert("POPO",&pile_type);
    pile_type=PileProps::ID1;
    datatypePileType.insert("ID1",&pile_type);
    pile_type=PileProps::ID2;
    datatypePileType.insert("ID2",&pile_type);

    datatypeOutLineInitSize=H5::EnumType(sizeof(OutLine::OutLineInitSize));
    OutLine::OutLineInitSize init_size;
    init_size=OutLine::AMR;
    datatypeOutLineInitSize.insert("AMR",&init_size);
    init_size=OutLine::DEM;
    datatypeOutLineInitSize.insert("DEM",&init_size);

    //init indexes to char* conversions
    coord_names.push_back("X");
    coord_names.push_back("Y");
    coord_names.push_back("Z");

    SinglePhaseHeuristicVarNames.push_back("h");
    SinglePhaseHeuristicVarNames.push_back("hVx");
    SinglePhaseHeuristicVarNames.push_back("hVy");


    SinglePhaseLevelSetVarNames.push_back("h");
    SinglePhaseLevelSetVarNames.push_back("hVx");
    SinglePhaseLevelSetVarNames.push_back("hVy");
    SinglePhaseLevelSetVarNames.push_back("phi");
    SinglePhaseLevelSetVarNames.push_back("junk1");
    SinglePhaseLevelSetVarNames.push_back("junk2");

    SinglePhasePhaseFieldVarNames.push_back("h");
    SinglePhasePhaseFieldVarNames.push_back("hVx");
    SinglePhasePhaseFieldVarNames.push_back("hVy");
    SinglePhasePhaseFieldVarNames.push_back("phi");
    SinglePhasePhaseFieldVarNames.push_back("junk1");
    SinglePhasePhaseFieldVarNames.push_back("junk2");


    TwoPhasesVarNames.push_back("h");
    TwoPhasesVarNames.push_back("h_liq");
    TwoPhasesVarNames.push_back("hVx_sol");
    TwoPhasesVarNames.push_back("hVy_sol");
    TwoPhasesVarNames.push_back("hVx_liq");
    TwoPhasesVarNames.push_back("hVy_liq");
}

bool ti_dir_exists(const char *filename)
{
    struct stat info;

    if(stat( filename, &info ) != 0)
        return false;
    else if(info.st_mode & S_IFDIR)
        return true;
    else
        return false;
}
void ti_mkdir(const char *filename)
{
    mkdir(filename, S_IRWXU | S_IRGRP|S_IXGRP | S_IROTH | S_IXOTH);
}
