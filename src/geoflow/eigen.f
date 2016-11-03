C*******************************************************************
C* Copyright (C) 2003 University at Buffalo
C*
C* This software can be redistributed free of charge.  See COPYING
C* file in the top distribution directory for more details.
C*
C* This software is distributed in the hope that it will be useful,
C* but WITHOUT ANY WARRANTY; without even the implied warranty of
C* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
C*
C* Author: 
C* Description: 
C*
C*******************************************************************
C* $Id: eigen.f 128 2007-06-07 19:51:52Z dkumar $ 
C*

C***********************************************************************
      subroutine eigen(Uvec, eigenvxmax, eigenvymax, evalue, tiny,
     1     kactxy, gravity, v_solid, v_fluid, eps, flowtype)
C***********************************************************************

      include 'rnr.h'
      integer flowtype
      double precision eigenvxmax, eigenvymax, v_solid(2), v_fluid(2)
      double precision evalue, tiny, Uvec(6), kactxy(2), gravity(3)
      double precision eps
      double precision sound_speed

      if (Uvec(1) .gt. tiny) then
c     iverson and denlinger
         if(kactxy(1) .lt. 0.d0) then
            kactxy(1) = -kactxy(1)
         endif

         if (flowtype.eq.1) then
           sound_speed = dsqrt(uvec(1)*kactxy(1)*gravity(3))
         else if (flowtype.eq.2) then
           sound_speed = dsqrt(uvec(1)*gravity(3))
         else
           sound_speed = dsqrt(uvec(2)*gravity(3)*kactxy(1)+
     $                         (uvec(1)-uvec(2))*gravity(3))
         endif

!        x-direction
         eigenvxmax=dmax1(dabs(v_solid(1)+sound_speed), 
     $                    dabs(v_fluid(1)+sound_speed))

!        y-direction
         eigenvymax=dmax1(dabs(v_solid(2)+sound_speed),
     $                    dabs(v_fluid(2)+sound_speed))
      else
         eigenvxmax=tiny
         eigenvymax=tiny
      endif
      evalue=dmax1(eigenvxmax,eigenvymax)

      if ( evalue .lt. 0. ) then
         write (*,*) "ERROR: Eigen-values -ve"
         write (*,*) sound_speed, v_solid, v_fluid
      endif

      return
      end
