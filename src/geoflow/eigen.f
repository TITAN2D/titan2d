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
      subroutine eigen(Uvec,eigenvxmax,eigenvymax,evalue,tiny,
     1     kactxy,gravity,VxVy)
C***********************************************************************

      include 'rnr.h'
      double precision eigenvxmax,eigenvymax, VxVy(2)
      double precision evalue, tiny, Uvec(3), kactxy(2),gravity(3)
         
      if (Uvec(1) .gt. tiny) then
c     iverson and denlinger
         if(kactxy(1) .lt. 0.d0) then
c	      write(*,*) 'negative kactxy'
  	      kactxy(1) = -kactxy(1)
	 endif
         eigenvxmax=dabs(VxVy(1)) +
     +        dsqrt(kactxy(1)*gravity(3)*Uvec(1))
         eigenvymax=dabs(VxVy(2)) +
     +        dsqrt(kactxy(1)*gravity(3)*Uvec(1))

      else
         eigenvxmax=tiny
         eigenvymax=tiny
      endif
         
      evalue=dmax1(eigenvxmax,eigenvymax)
         
      return
      end
