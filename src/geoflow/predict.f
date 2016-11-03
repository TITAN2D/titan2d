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
C* $Id: predict.f 143 2007-06-25 17:58:08Z dkumar $ 
C*

C***********************************************************************
      subroutine predict(Uvec,dUdx,dUdy,Uprev,tiny,kactxy,dt2, g, 
     1     curv, bedfrictang, intfrictang,
     2     dgdx, frict_tiny, order_flag, VxVy, 
     3     IF_STOPPED,fluxsrc)
C***********************************************************************
c     this routine only calculates Uvec for a half step

***********************************************************************
******NOTE:  d(g(3)*Uvec(1))/dx is approximated by g(3)*dUvec(1)/dx !!!
***********************************************************************
      integer order_flag, IF_STOPPED
      double precision dUdx(6),dUdy(6),kactxy,Uvec(6)
      double precision tiny,dt2,Uprev(6),dgdx(2)
      double precision c_sq, h_inv, g(3),curv(2)
      double precision intfrictang, bedfrictang, sgn, frict_tiny
      double precision sgn_dudy, sgn_dvdx, tmp

!     local variables
      double precision speed,unitvx,unitvy,VxVy(2)
      double precision forcegrav, forceintx, forceinty
      double precision forcebedmax, forcebedequil
      double precision forcebedx, forcebedy
      double precision tanbed,fluxsrc(3)
      external sgn


      end subroutine predict
