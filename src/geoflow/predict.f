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
     2     dgdx, frict_tiny, order_flag, VxVyB, 
     3     IF_STOPPED,fluxsrc)
C***********************************************************************
c     this routine only calculates Uvec for a half step

***********************************************************************
******NOTE:  d(g(3)*Uvec(1))/dx is approximated by g(3)*dUvec(1)/dx !!!
***********************************************************************
      include 'rnr.h'
      integer order_flag, IF_STOPPED
      double precision dUdx(3),dUdy(3),kactxy,Uvec(3)
      double precision tiny,dt2, Uprev(3), dgdx(2)
      double precision c_sq, h_inv, g(3),curv(2)
      double precision intfrictang, bedfrictang, sgn, frict_tiny
      double precision sgn_dudy, sgn_dvdx, tmp
c      double precision dnorm
      double precision speed,unitvx,unitvy,VxVy(2),VxVyB(2),VxVyS(2)
      double precision forcegrav, forceintx, forceinty
      double precision forcebedmax, forcebedequil
      double precision forcebedx, forcebedy
      double precision tanbed,fluxsrc(3)
      external sgn


c     curv := inverse of radius of curvature = second derivative of 
c     position normal to tangent with respect to distance along tangent,
c     if dz/dx=0 curve=d2z/dx2, otherwise rotate coordinate system so 
c     dz/dx=0, that is mathematical definition of curvature I believe
c     laercio returns d2z/dx2 whether or not dz/dx=0 in his GIS functions

c TEST********** order_flag = 1 

c      write(*,*) "order_flag in predict=", order_flag 
      
      uprev(1) = uvec(1)
      uprev(2) = uvec(2)
      uprev(3) = uvec(3)

      if(IF_STOPPED.eq.2) then
         VxVy(1)=0.0
         VxVy(2)=0.0
         VxVyS(1)=0.0
         VxVyS(2)=0.0
      else
         VxVy(1)=VxVyB(1)
c     uvec(2)/uvec(1)
         VxVy(2)=VxVyB(2)
c     uvec(3)/uvec(1)
         VxVyS(1)=VxVyB(1)
         VxVyS(2)=VxVyB(2)
      endif

      if (order_flag .eq. 2) then
         c_sq = kactxy*g(3)*Uvec(1)
c     h_inv := 1/Uvec(1)                         
         
         Uvec(1)=Uvec(1) - dt2*(dUdx(2)+dUdy(3)+fluxsrc(1))
         Uvec(1)=dmax1(Uvec(1),0.d0)
         
c     dF/dU, dG/dU and S terms if Uvec(1) > TINY !
         if(Uprev(1) .gt. TINY) then
            h_inv = 1.d0/Uprev(1)
            tanbed=dtan(bedfrictang)

c     here speed is speed squared
            speed=VxVy(1)**2+VxVy(2)**2
            if(speed.GT.0.0) then
c     here speed is speed
               speed=dsqrt(speed)
               unitvx=VxVy(1)/speed
               unitvy=VxVy(2)/speed
            else
               unitvx=0.0
               unitvy=0.0
            endif

c            dnorm=dsqrt(uprev(2)**2+uprev(3)**2+tiny**2)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ****** X-dir ******
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     dF/dU and dG/dU terms
            Uvec(2) = Uvec(2)- 
     1           dt2*((c_sq-VxVy(1)**2)*dUdx(1)+ 
     2           2.d0*VxVy(1)*dUdx(2)-
     3           VxVy(1)*VxVy(2)*dUdy(1)+
     4           VxVy(2)*dUdy(2)+ 
     5           VxVy(1)*dUdy(3)+
     6           fluxsrc(2)) 
            
c     x direction source terms

c     the gravity force in the x direction
            forcegrav=g(1)*Uvec(1) 

c     the internal friction force
            tmp = h_inv*(dudy(2)-VxVyS(1)*dudy(1))
            sgn_dudy = sgn(tmp, frict_tiny)
            forceintx=sgn_dudy*Uvec(1)*kactxy*(g(3)*dUdy(1)
     $           +dgdx(2)*Uvec(1))*dsin(intfrictang) 
            
c     the bed friction force for fast moving flow 
            forcebedx=unitvx*
     $           dmax1(g(3)*Uvec(1)+VxVyS(1)*Uvec(2)*curv(1),0.0d0)
     $           *tanbed

            if(IF_STOPPED.eq.2.and.1.eq.0) then
c     the bed friction force for stopped or nearly stopped flow

c     the static friction force is LESS THAN or equal to the friction
c     coefficient times the normal force but it can NEVER exceed the 
c     NET force it is opposing

c     maximum friction force the bed friction can support
               forcebedmax=
     $              dmax1(g(3)*Uvec(1)+VxVyS(1)*Uvec(2)*curv(1),0.0d0)
     $              *tanbed

c     the NET force the bed friction force is opposing 
               forcebedequil=forcegrav
c     $              -kactxy*g(3)*Uvec(1)*dUdx(1)
     $              -forceintx

c     the "correct" stopped or nearly stopped flow bed friction force 
c     (this force is not entirely "correct" it will leave a "negligible"
c     (determined by stopping criteria) amount of momentum in the cell
               forcebedx=sgn(forcebedequil,dmin1(forcebedmax,
     $              dabs(forcebedx)+dabs(forcebedequil)))
c               forcebedx=
c     $              sgn(forcebed2,dmin1(forcebed1,dabs(forcebed2)))
c            else
            endif

c     all the x source terms
            Uvec(2) = Uvec(2) + dt2*(forcegrav -forcebedx -forceintx)
c            write(*,*) 'int', forceintx, 'bed', forcebedx

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ****** Y-dir ******
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     dF/dU and dG/dU terms
            Uvec(3) = Uvec(3) - 
     1           dt2*((c_sq-VxVy(2)**2)*dUdy(1) +
     2           2.d0*VxVy(2)*dUdy(3) -
     3           VxVy(1)*VxVy(2)*dUdx(1)+
     4           VxVy(2)*dUdx(2)+
     5           VxVy(1)*dUdx(3)+
     6           fluxsrc(3)) 

c     the gravity force in the y direction
            forcegrav=g(2)*Uvec(1) 

c     the internal friction force
            tmp = h_inv*(dudx(3)-VxVyS(2)*dudx(1))
            sgn_dvdx = sgn(tmp, frict_tiny)
            forceinty=sgn_dvdx*Uvec(1)*kactxy*(g(3)*
     $           dudx(1)+dgdx(1)*Uvec(1))*dsin(intfrictang)

c     the bed friction force for fast moving flow 
               forcebedy=unitvy*
     $           dmax1(g(3)*Uvec(1)+VxVyS(2)*Uvec(3)*curv(2),0.0d0)
     $           *dtan(bedfrictang)

            if(If_STOPPED.eq.2.and.1.eq.0) then
c     the bed friction force for stopped or nearly stopped flow

               forcebedmax=
     $              dmax1(g(3)*Uvec(1)+VxVyS(2)*Uvec(3)*curv(2),0.0d0)
     $              *tanbed

c     the NET force the bed friction force is opposing 
               forcebedequil=forcegrav
c     $              -kactxy*g(3)*Uvec(1)*dUdy(1)
     $              -forceinty

c     the "correct" stopped or nearly stopped flow bed friction force 
c     (this force is not entirely "correct" it will leave a "negligible"
c     (determined by stopping criteria) amount of momentum in the cell
               forcebedy=sgn(forcebedequil,dmin1(forcebedmax,
     $              dabs(forcebedy)+dabs(forcebedequil)))

c               forcebedy=sgn(forcebed2,dmin1(forcebed1,dabs(forcebed2)))
c            else
            endif

c     all the y source terms
            Uvec(3) = Uvec(3) + dt2*(forcegrav -forcebedy -forceinty)

         endif
      endif
      
      return
      end
