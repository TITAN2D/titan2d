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

#ifndef OUTLINE_H
#define	OUTLINE_H

#if 0
//! this function updates the maximum throughout time pileheight in every cell covered by an arbitrary element
inline void OutLine::update(ElementsHashTable* ElemTable, NodeHashTable* NodeTable)
{
    int Nelms = ElemTable->size();
    
    ContentStatus *status_=&(ElemTable->status_[0]);
    int *adapted_=&(ElemTable->adapted_[0]);
    
    //int Nelms = ElemTable->getNumberOfLocalElements();
    //Element** Elms = (Element**) ElemTable->getLocalElementsValues();

    for(int i = 0; i < Nelms; i++)
    {
        if(status_[i]<0||adapted_[i]<=0)
            continue;
        
        Element* Curr_El = &(ElemTable->elenode_[i]);
        
        double edx = Curr_El->dx(0);
        double edy = Curr_El->dx(1);
        
        double xstart=Curr_El->coord(0) - 0.5 * edx;
        double xstop=Curr_El->coord(0) + 0.5 * edx;
        double ystart=Curr_El->coord(1) - 0.5 * edy;
        double ystop=Curr_El->coord(1) + 0.5 * edy;


        int ixstart = (int) ((xstart - xminmax[0]) / dx + 0.5);
        int ixstop = (int) ((xstop - xminmax[0]) / dx + 0.5);
        int iystart = (int) ((ystart - yminmax[0]) / dy + 0.5);
        int iystop = (int) ((ystop - yminmax[0]) / dy + 0.5);

        ixstart=max(ixstart, 0);
        
        if(ixstop == ixstart)
        {
            ixstart = (int) ((xstart - xminmax[0]) / dx);
            ixstop = ixstart + 1;
        }
        ixstop=min(Nx,ixstop);

        iystart=max(iystart, 0);
        if(iystop == iystart)
        {
            iystart = (int) ((ystart - yminmax[0]) / dy);
            iystop = iystart + 1;
        }
        iystop=min(Ny,iystop);

        double ke = 0.0;
        if(Curr_El->state_vars(0) > 1.0E-04)
        {
            ke = 0.5 * (Curr_El->state_vars(1) * Curr_El->state_vars(1) + Curr_El->state_vars(2) * Curr_El->state_vars(2)) / Curr_El->state_vars(0);
        }


        for(int iy = iystart; iy < iystop; iy++)
        {
            for(int ix = ixstart; ix < ixstop; ix++)
            {
                cum_kinergy[iy][ix] += ke;
                if(Curr_El->state_vars(0) > pileheight[iy][ix])
                    pileheight[iy][ix] = Curr_El->state_vars(0);
                if(ke > max_kinergy[iy][ix])
                    max_kinergy[iy][ix] = ke;
            }
        }
    }
}
#endif
#if 0
#include "tivector.h"

//! this function updates the maximum throughout time pileheight in every cell covered by an arbitrary element
inline void OutLine::update(ElementsHashTable* ElemTable, NodeHashTable* NodeTable)
{
    int Nelms = ElemTable->getNumberOfLocalElements();
    Element** Elms = (Element**) ElemTable->getLocalElementsValues();

    for(int i = 0; i < Nelms; i++)
    {
        Element* Curr_El = Elms[i];
        double edx = Curr_El->dx(0);
        double edy = Curr_El->dx(1);
        
        double xstart=Curr_El->coord(0) - 0.5 * edx;
        double xstop=Curr_El->coord(0) + 0.5 * edx;
        double ystart=Curr_El->coord(1) - 0.5 * edy;
        double ystop=Curr_El->coord(1) + 0.5 * edy;


        int ixstart = (int) ((xstart - xminmax[0]) / dx + 0.5);
        int ixstop = (int) ((xstop - xminmax[0]) / dx + 0.5);
        int iystart = (int) ((ystart - yminmax[0]) / dy + 0.5);
        int iystop = (int) ((ystop - yminmax[0]) / dy + 0.5);

        ixstart=max(ixstart, 0);
        
        if(ixstop == ixstart)
        {
            ixstart = (int) ((xstart - xminmax[0]) / dx);
            ixstop = ixstart + 1;
        }
        ixstop=min(Nx,ixstop);

        iystart=max(iystart, 0);
        if(iystop == iystart)
        {
            iystart = (int) ((ystart - yminmax[0]) / dy);
            iystop = iystart + 1;
        }
        iystop=min(Ny,iystop);

        double ke = 0.0;
        if(Curr_El->state_vars(0) > 1.0E-04)
        {
            ke = 0.5 * (Curr_El->state_vars(1) * Curr_El->state_vars(1) + Curr_El->state_vars(2) * Curr_El->state_vars(2)) / Curr_El->state_vars(0);
        }


        for(int iy = iystart; iy < iystop; iy++)
        {
            for(int ix = ixstart; ix < ixstop; ix++)
            {
                cum_kinergy[iy][ix] += ke;
                if(Curr_El->state_vars(0) > pileheight[iy][ix])
                    pileheight[iy][ix] = Curr_El->state_vars(0);
                if(ke > max_kinergy[iy][ix])
                    max_kinergy[iy][ix] = ke;
            }
        }
    }
}
#endif
#if 0
//! this function updates the maximum throughout time pileheight in every cell covered by an arbitrary element
inline void OutLine::update(ElementsHashTable* ElemTable, NodeHashTable* NodeTable)
{
    int Nelms = ElemTable->size();
    
    ContentStatus *status_=&(ElemTable->status_[0]);
    int *adapted_=&(ElemTable->adapted_[0]);
    double *dx_=&((ElemTable->dx_[0])[0]);
    double *dy_=&((ElemTable->dx_[1])[0]);
    
    double *coord_x=&((ElemTable->coord_[0])[0]);
    double *coord_y=&((ElemTable->coord_[1])[0]);
    
    double *h   = &((ElemTable->state_vars_[0])[0]);
    double *hVx = &((ElemTable->state_vars_[1])[0]);
    double *hVy = &((ElemTable->state_vars_[2])[0]);
    
    
    for(int i = 0; i < Nelms; i++)
    {
        if(status_[i]<0||adapted_[i]<=0)
            continue;
        
        double xstart=coord_x[i] - 0.5 * dx_[i];
        double xstop=coord_x[i] + 0.5 * dx_[i];
        double ystart=coord_y[i] - 0.5 * dy_[i];
        double ystop=coord_y[i] + 0.5 * dy_[i];


        int ixstart = (int) ((xstart - xminmax[0]) / dx + 0.5);
        int ixstop = (int) ((xstop - xminmax[0]) / dx + 0.5);
        int iystart = (int) ((ystart - yminmax[0]) / dy + 0.5);
        int iystop = (int) ((ystop - yminmax[0]) / dy + 0.5);

        ixstart=max(ixstart, 0);
        
        if(ixstop == ixstart)
        {
            ixstart = (int) ((xstart - xminmax[0]) / dx);
            ixstop = ixstart + 1;
        }
        ixstop=min(Nx,ixstop);

        iystart=max(iystart, 0);
        if(iystop == iystart)
        {
            iystart = (int) ((ystart - yminmax[0]) / dy);
            iystop = iystart + 1;
        }
        iystop=min(Ny,iystop);

        double ke = 0.0;
        if(h[i] > 1.0E-04)
        {
            ke = 0.5 * (hVx[i] * hVx[i] + hVy[i] * hVy[i]) / h[i];
        }


        for(int iy = iystart; iy < iystop; iy++)
        {
            for(int ix = ixstart; ix < ixstop; ix++)
            {
                cum_kinergy[iy][ix] += ke;
                pileheight[iy][ix] = max(pileheight[iy][ix],h[i]);
                max_kinergy[iy][ix] = max(max_kinergy[iy][ix],ke);
            }
        }
    }
}
#endif
inline void OutLine::update(ElementsHashTable* ElemTable, NodeHashTable* NodeTable)
{
    int Nelms = ElemTable->getNumberOfLocalElements();
    Element** Elms = (Element**) ElemTable->getLocalElementsValues();

    for(int i = 0; i < Nelms; i++)
    {
        Element* Curr_El = Elms[i];
        double edx = Curr_El->dx(0);
        double edy = Curr_El->dx(1);
        //update the record of maximum pileheight in the area covered by this element
        double height = Curr_El->state_vars(0);
        //if (hheight > 0 && hheight < 0)
        //	;
        double xstart=Curr_El->coord(0) - 0.5 * edx;
        double xstop=Curr_El->coord(0) + 0.5 * edx;
        double ystart=Curr_El->coord(1) - 0.5 * edy;
        double ystop=Curr_El->coord(1) + 0.5 * edy;

        double hv[4];
        if(elementType == ElementType::TwoPhases)
        {
            hv[0]=Curr_El->state_vars(2);
            hv[1]=Curr_El->state_vars(3);
            hv[2]=Curr_El->state_vars(4);
            hv[3]=Curr_El->state_vars(5);
        }
        if(elementType == ElementType::SinglePhase)
        {
            hv[0]=Curr_El->state_vars(1);
            hv[1]=Curr_El->state_vars(2);
        }
        //update(xstart, xstop, ystart, ystop, height, hv);

        int ixstart = (int) ((xstart - xminmax[0]) / dx + 0.5);
        int ixstop = (int) ((xstop - xminmax[0]) / dx + 0.5);
        int iystart = (int) ((ystart - yminmax[0]) / dy + 0.5);
        int iystop = (int) ((ystop - yminmax[0]) / dy + 0.5);

        if(ixstart < 0)
            ixstart = 0;
        if(ixstop == ixstart)
        {
            ixstart = (int) ((xstart - xminmax[0]) / dx);
            ixstop = ixstart + 1;
        }
        if(ixstop > Nx)
            ixstop = Nx;

        if(iystart < 0)
            iystart = 0;
        if(iystop == iystart)
        {
            iystart = (int) ((ystart - yminmax[0]) / dy);
            iystop = iystart + 1;
        }
        if(iystop > Ny)
            iystop = Ny;

        double ke = 0.;
        if(height > 1.0E-04)
        {
            if(elementType == ElementType::TwoPhases)
            {
                //@TODO: Check the correctness of TwoPhases ke calculation
                ke = 0.5 * (hv[0] * hv[0] + hv[1] * hv[1] + hv[2] * hv[2] + hv[3] * hv[3]) / height;
            }
            if(elementType == ElementType::SinglePhase)
            {
                ke = 0.5 * (hv[0] * hv[0] + hv[1] * hv[1]) / height;
            }
        }


        for(int iy = iystart; iy < iystop; iy++)
        {
            for(int ix = ixstart; ix < ixstop; ix++)
            {
                cum_kinergy[iy][ix] += ke;
                if(height > pileheight[iy][ix])
                    pileheight[iy][ix] = height;
                if(ke > max_kinergy[iy][ix])
                    max_kinergy[iy][ix] = ke;
            }
        }
    }
}
#endif	/* OUTLINE_H */

