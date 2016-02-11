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
#ifndef ELEMENTS_PROPERTIES_INLINE_HPP
#define	ELEMENTS_PROPERTIES_INLINE_HPP

#include "hashtab.h"
#include "titan2d_utils.h"


inline void ElementsProperties::calc_flux(ti_ndx_t ndx,FluxProps *fluxprops, TimeProps *timeprops)
{
    int no_of_sources = fluxprops->no_of_sources;
    double temp_coef;
    double major, minor, dswap, xcoord, ycoord;
    int inode, isrc;
    double curr_time = timeprops->cur_time;
    double start_time, end_time;
    double m_Influx_[3]={0.0,0.0,0.0};

    for(inode = 0; inode < 9; inode++)
    {
        double temp_coef2 = temp_coef = 0.0;
        double sum_flux_xmom_ymom[3] ={ 0.0, 0.0, 0.0 };
        ti_ndx_t nd_ndx;
        if(inode==8)
            nd_ndx=node_bubble_ndx_[ndx];
        else
            nd_ndx=node_key_ndx_[inode][ndx];
        ASSERT3(ti_ndx_not_negative(nd_ndx));

        xcoord = node_coord_[0][nd_ndx];
        ycoord = node_coord_[1][nd_ndx];

        for(int isrc = 0; isrc < no_of_sources; isrc++)
        {
            // normalize start and end time
            start_time = fluxprops->start_time[isrc];
            end_time = fluxprops->end_time[isrc];

            if((fluxprops->start_time[isrc] <= curr_time) && (fluxprops->end_time[isrc] >= curr_time))
            {

                /*
                 if(timeprops->timesec()>30.0){
                 printf("flux not stopped after 30 seconds\n");
                 assert(0);
                 }
                 */

                // "undo" elliptical pile rotation ... from (x,y)->(major,minor)
                major = xcoord - fluxprops->xCen[isrc];
                minor = ycoord - fluxprops->yCen[isrc];
                dswap = (major * fluxprops->cosrot[isrc] + minor * fluxprops->sinrot[isrc]) / fluxprops->majorrad[isrc];
                minor = (-major * fluxprops->sinrot[isrc] + minor * fluxprops->cosrot[isrc]) / fluxprops->minorrad[isrc];
                major = dswap;
                if(major * major + minor * minor < 1.0)
                {
                    if(inode < 4)
                        temp_coef = 1.0 / 16.0;
                    if(inode >= 4 && inode < 8)
                        temp_coef = 1.0 / 8.0;
                    if(inode == 8)
                    {
                        temp_coef = 1.0 / 4.0;
                        //printf("xcoord=%g ycoord=%g isrc=%d major=%g minor=%g\n",
                        //   xcoord,ycoord,isrc,major,minor);
                    }
                }
                /*
                 * distance from cell-center is used regardless of cell-center
                 * being inside or outside the source region
                 */

                //hydrograph flux starts at max rate decays linearly to zero at end of duration
                double tempflux = (fluxprops->influx[isrc])
                        * (1.0 - (curr_time - fluxprops->start_time[isrc])
                                / (fluxprops->end_time[isrc] - fluxprops->start_time[isrc]));
                temp_coef *= (1.0 - major * major - minor * minor) * tempflux;
                //temp_coef*=(1.0-major*major-minor*minor)*(fluxprops->influx[isrc]);

                if(temp_coef > 0.0)
                {
                    sum_flux_xmom_ymom[0] += temp_coef;
                    sum_flux_xmom_ymom[1] += temp_coef * fluxprops->xVel[isrc];
                    sum_flux_xmom_ymom[2] += temp_coef * fluxprops->yVel[isrc];
                }

                if(temp_coef > temp_coef2)
                {
                    temp_coef2 = temp_coef;

                    //printf("xcoord=%g ycoord=%g isrc=%d major=%g minor=%g\n",
                    //   xcoord,ycoord,isrc,major,minor);
                }

            }

        }
        if(temp_coef2 > 0.0)
        {
            m_Influx_[0] += temp_coef2;

            //if(sum_flux_xmom_ymom[0]>0.0) {
            m_Influx_[1] += sum_flux_xmom_ymom[1] * temp_coef2 / sum_flux_xmom_ymom[0];
            m_Influx_[2] += sum_flux_xmom_ymom[2] * temp_coef2 / sum_flux_xmom_ymom[0];
            //}
        }

    }

    Influx_[0][ndx]=m_Influx_[0];
    Influx_[1][ndx]=m_Influx_[1];
    Influx_[2][ndx]=m_Influx_[2];

    ASSERT2(Influx_[0][ndx] >= 0.0);//"error in Influx[0] assignment\n"
    Awet_[ndx]=(Influx_[0][ndx] > 0.0) ? 1.0 : 0.0;

    return;
}



#endif	/* ELEMENTS_PROPERTIES_INLINE_HPP */

