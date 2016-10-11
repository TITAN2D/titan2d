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
C* $Id: getcoef.f 128 2007-06-07 19:51:52Z dkumar $ 
C*

C***********************************************************************
      subroutine gmfggetcoef(Uvec, dUdx, dUdy, dx, bedfrictang,
     1     intfrictang, Kactx, Kacty, tiny, epsilon)
C***********************************************************************
        
      include 'rnr.h'
      double precision dUdx(3), dUdy(3), Uvec(3), coefABCD(4), dx(2)
      double precision Kactx, Kacty, zeta, epsilonx, tiny
      double precision vel,cosphi,tandel,w,sgn,ka, epsilon
      
      double precision bedfrictang,intfrictang, tttest
      integer i,j 
c     vel is used to determine if velocity gradients are converging or diverging
      external sgn
      
c     COEFFICIENTS
      cosphi=dcos(intfrictang)
      tandel=dtan(bedfrictang)
      if(uvec(1) .gt. tiny) then
c         tttest = sgn(vel, tiny)
         vel=(dUdx(2)/uvec(1) - uvec(2)*dudx(1)/uvec(1)**2+
     1        dUdy(3)/uvec(1) - uvec(3)*dudy(1)/uvec(1)**2)
         Kactx=(2.d0/cosphi**2)*(1.d0-sgn(vel,tiny)*
     1        dsqrt(dabs(1.d0-(1.d0+tandel**2)*cosphi**2) )) -1.d0
         Kacty=(2.d0/cosphi**2)*(1.d0-sgn(vel,tiny)*
     1        dsqrt(dabs(1.d0-(1.d0+tandel**2)*cosphi**2) )) -1.d0

c     if there is no yielding...
         if(dabs(uvec(2)/uvec(1)) .lt. tiny .and.
     1        dabs(uvec(3)/uvec(1)) .lt. tiny) then
            kactx = 1.d0
            kacty = 1.d0
         endif
      else 
         vel = 0.d0
         kactx = 1.d0
         kacty = 1.d0
      endif
      
      kactx = epsilon * kactx
      kacty = epsilon * kacty

      return
      end
