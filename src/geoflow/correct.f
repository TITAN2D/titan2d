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
C* $Id: correct.f 143 2007-06-25 17:58:08Z dkumar $ 
C*


C***********************************************************************
      subroutine correct2ph(uvec, uprev, fluxxp, fluxyp, fluxxm, fluxym,
     1     tiny, dtdx, dtdy, dt, dUdx, dUdy, xslope, yslope,
     2     curv, intfrictang, bedfrictang, g, kactxy, frict_tiny,
     3     forceint, forcebed, DO_EROSION, eroded, v_solid, v_fluid,
     4     den_solid, den_fluid, terminal_vel, eps, IF_STOPPED, fluxsrc)
C***********************************************************************

      implicit none
      double precision forceint, forcebed, eroded, speed
      double precision forceintx, forceinty
      double precision forcebedx, forcebedy
      double precision forcebedmax, forcebedequil, forcegrav
      double precision unitvx, unitvy, v_solid(2), v_fluid(2)
      double precision den_frac, den_solid, den_fluid
      double precision alphaxx, alphayy, alphaxy, alphaxz, alphayz
      double precision tanbed, terminal_vel

      double precision fluxxp(6),fluxyp(6),tiny, uprev(6), ustore(6)
      double precision fluxxm(6), fluxym(6)
      double precision uvec(6), dUdx(6), dUdy(6)
      double precision h_inv, hphi_inv, curv(2), frict_tiny
      double precision intfrictang, bedfrictang, kactxy, dgdx(2)
      double precision dtdx, dtdy, dt, g(3), sgn_dudy, sgn_dvdx, tmp
      double precision dnorm, fluxsrc(6)
      double precision xslope,yslope,slope
      double precision t1, t2, t3, t4, t5
      double precision erosion_rate,threshold,es,totalShear
      double precision eps, drag(4)

!     function calls
      double precision sgn

      integer i
      integer DO_EROSION, IF_STOPPED
      parameter(threshold=1.0D-02,erosion_rate=0.1)

c     initialize to zero
      forceintx=0.0
      forcebedx=0.0
      forceinty=0.0
      forcebedy=0.0
      unitvx=0.0
      unitvy=0.0
      eroded=0.0

      slope=dsqrt(xslope*xslope+yslope*yslope)
      den_frac = den_fluid/den_solid
      do 10 i = 1,6
10    ustore(i)=uprev(i)+dt*fluxsrc(i)
     1     -dtdx*(fluxxp(i)-fluxxm(i))
     2     -dtdy*(fluxyp(i)-fluxym(i))

      if(ustore(1).gt.tiny) then
c     Source terms ...
c     here speed is speed squared
         speed=v_solid(1)**2+v_solid(2)**2
         if(speed.gt.0.0) then
c     here speed is speed
            speed=dsqrt(speed)
            unitvx=v_solid(1)/speed
            unitvy=v_solid(2)/speed
         else
            unitvx=0.0
            unitvy=0.0
         endif
         tanbed=dtan(bedfrictang)
         h_inv = 1.d0/uvec(1)
         hphi_inv = 1.0/uvec(2)
         alphaxx = kactxy
         alphayy = kactxy
         den_frac = den_fluid/den_solid
         call calc_drag_force (uvec, v_solid, v_fluid, den_solid, 
     &                         den_fluid, terminal_vel, drag, tiny)
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        solid fraction x-direction source terms
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        alphaxy -- see pitman-le (2005)
         tmp = hphi_inv*(dudy(3)-v_solid(1)*dudy(2))
         sgn_dudy = sgn(tmp, frict_tiny)
         alphaxy = sgn_dudy*dsin(intfrictang)*kactxy

c        alphaxz (includes centrifugal effects)
         alphaxz = -unitvx*tanbed*(1.0+v_solid(1)**2*curv(1)/g(3))

c        evaluate t1 
         t1 = (1.0-den_frac)*(-alphaxx*xslope
     $         -alphaxy*yslope + alphaxz)*uvec(2)*g(3)
c        evaluate t2
         t2 = eps*den_frac*uvec(2)*g(3)*dudx(1)
c        evaluate t3
         t3 = eps*den_frac*uvec(2)*g(3)*xslope
c        evaluate t4
         t4 = uvec(2)*g(1)
c        evaluate drag
         t5 = drag(1)

c        update ustore
         ustore(3) = ustore(3) + dt*(t1-t2-t3+t4+t5)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     solid fraction y-direction source terms
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        alphaxy -- see pitman-le (2005) for definitions
         tmp = hphi_inv*(dudx(4)-v_solid(2)*dudx(2))
         sgn_dvdx = sgn(tmp, frict_tiny)
         alphaxy = sgn_dvdx*dsin(intfrictang)*kactxy

c        alphayz
         alphayz = -unitvy*tanbed*(1.0+v_solid(2)**2*curv(2)/g(3))

c        evaluate t1
         t1 = (1.0-den_frac)*(-alphaxy*xslope
     $         -alphayy*yslope + alphayz)*uvec(2)*g(3)
c        evaluate t2
         t2 = eps*den_frac*uvec(2)*dudy(1)
c        evaluate t3
         t3 = eps*den_frac*uvec(2)*g(3)*yslope
c        evaluate t4 ( gravity along y-dir )
         t4 = uvec(2)*g(2)
c        drag term
         t5 = drag(2)
         ustore(4) = ustore(4) + dt*(t1-t2-t3+t4+t5)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    fluid fraction x-direction source terms
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        gravity on fluid
         t4 = uvec(1)*g(1)
c        drag force on fluid
         t5 = drag(3)
         ustore(5) = ustore(5) + dt*(t4 - t5)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    fluid fraction y-direction source terms
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        gravity on fluid
         t4 = uvec(1)*g(2)
c        drag force on fluid
         t5 = drag(4)
         ustore(6) = ustore(6) + dt*(t4 - t5)
      endif

c     computation of magnitude of friction forces for statistics
      forceint=unitvx*forceintx+unitvy*forceinty
      forcebed=unitvx*forcebedx+unitvy*forcebedy

c     update the state variables
      do 20 i=1,6
20       uvec(i) = ustore(i)

      return
      end


C***********************************************************************
      subroutine correct(Uvec,Uprev,fluxxp,fluxyp, fluxxm, fluxym,
     1     tiny,dtdx,dtdy,dt,dUdx,dUdy,xslope,yslope,
     2     curv, intfrictang, bedfrictang,g,kactxy, dgdx,
     3     frict_tiny,forceint,forcebed,DO_EROSION,eroded, VxVy,
     4     IF_STOPPED,fluxsrc)
C***********************************************************************

      include 'rnr.h'
      double precision forceint, forcebed, eroded, speed
      double precision forceintx, forceinty
      double precision forcebedx, forcebedy
      double precision forcebedmax, forcebedequil, forcegrav
      double precision unitvx, unitvy, VxVy(2)
      double precision tanbed

      double precision fluxxp(3),fluxyp(3),tiny, Uprev(3),Ustore(3)
      double precision fluxxm(3), fluxym(3)
      double precision uvec(3), dUdx(3), dUdy(3)
      double precision sgn, h_inv, curv(2), frict_tiny
      double precision intfrictang, bedfrictang, kactxy, dgdx(2)
      double precision dtdx, dtdy, dt, g(3), sgn_dudy, sgn_dvdx, tmp
      double precision dnorm,fluxsrc(3)
      double precision xslope,yslope,slope
      double precision erosion_rate,threshold,es,totalShear
      integer DO_EROSION, IF_STOPPED
      parameter(threshold=1.0D-02,erosion_rate=0.1)

c      parameter(erosion_rate=0.015)(data for this is in ...DATA dir)
c      parameter(erosion_rate=0.125).....2.983063e+05.........5.718950e+05.....dec_2_erosion_4 c      parameter(erosion_rate=0.250) ....2.983063e+05.........2.582245e+06.....dec_2_erosion_5
c      parameter(erosion_rate=0.19)......2.983063e+05........ 1.183863e+06.....dec_10_erosion_2.
c      parameter(erosion_rate=0.17)......2.983063e+05........ 9.570268e+05.....dec_10_erosion_3.
c      parameter(erosion_rate=0.0) ......2.983063e+05.........2.983062e+05.....dec_11_erosion1.
c      parameter(erosion_rate=0.16)......2.983063e+05.........9.306538e+05.....dec_11_erosion2.
c      parameter(erosion_rate=0.14)......2.983063e+05.........6.282101e+05.....dec_11_erosion3.
c       parameter(erosion_rate=0.15)
c      ......2.983063e+05.........7.039735e+05.....dec_11_erosion4.
c      parameter(erosion_rate=0.155).....2.983063e+05.........7.390704e+05.....dec_11_erosion5.
c       parameter(erosion_rate=0.156)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      external sgn            

      slope=dsqrt(xslope*xslope+yslope*yslope)
c
c      threshold=1.278820338*dcos(slope)*
C     &     dmax1((1-dtan(slope)/dtan(intfrictang)),0.0d0)**2

      Ustore(1)=Uprev(1)
     1     -dtdx*(fluxxp(1)-fluxxm(1))
     2     -dtdy*(fluxyp(1)-fluxym(1))
     3     +dt*fluxsrc(1)
      Ustore(1)=dmax1(Ustore(1),0.0d0)
      
      Ustore(2)=Uprev(2)
     1     -dtdx*(fluxxp(2)-fluxxm(2))
     2     -dtdy*(fluxyp(2)-fluxym(2))
     3     +dt*fluxsrc(2)
      
      Ustore(3)=Uprev(3)
     1     -dtdx*(fluxxp(3)-fluxxm(3))
     2     -dtdy*(fluxyp(3)-fluxym(3))
     3     +dt*fluxsrc(3)

c     initialize to zero
      forceintx=0.0
      forcebedx=0.0
      forceinty=0.0
      forcebedy=0.0
      unitvx=0.0
      unitvy=0.0
      eroded=0.0

      if(Uvec(1). gt. tiny) then
c     S terms
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
         tanbed=dtan(bedfrictang)
         h_inv = 1.d0/Uvec(1)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     x direction source terms
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     the gravity force in the x direction
         forcegrav=g(1)*Uvec(1) 

c     the internal friction force
         tmp = h_inv*(dudy(2)-VxVy(1)*dudy(1))
         sgn_dudy = sgn(tmp, frict_tiny)
         forceintx=sgn_dudy*Uvec(1)*kactxy*(g(3)*dUdy(1)
     $        +dgdx(2)*Uvec(1))*dsin(intfrictang) 
     
c     the bed friction force for fast moving flow 
         forcebedx=unitvx*
     $        dmax1(g(3)*Uvec(1)+VxVy(1)*Uvec(2)*curv(1),0.0d0)
     $        *tanbed

         if(IF_STOPPED.eq.2.and.1.eq.0) then
c     the bed friction force for stopped or nearly stopped flow

c     the static friction force is LESS THAN or equal to the friction
c     coefficient times the normal force but it can NEVER exceed the 
c     NET force it is opposing

c     maximum friction force the bed friction can support
            forcebedmax=g(3)*Uvec(1)*tanbed

c     the NET force the bed friction force is opposing 
            forcebedequil=forcegrav
c     $           -kactxy*g(3)*Uvec(1)*dUdx(1)
     $           -forceintx

c     the "correct" stopped or nearly stopped flow bed friction force 
c     (this force is not entirely "correct" it will leave a "negligible"
c     (determined by stopping criteria) amount of momentum in the cell
            forcebedx=sgn(forcebedequil,dmin1(forcebedmax,
     $           dabs(forcebedx)+dabs(forcebedequil)))
c     forcebedx=sgn(forcebed2,dmin1(forcebed1,dabs(forcebed2)))

c     not really 1 but this makes friction statistics accurate
            unitvx=1.0
c     else
            
         endif

         Ustore(2) = Ustore(2) + dt*(forcegrav -forcebedx -forceintx)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     y direction source terms
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     the gravity force in the y direction
         forcegrav=g(2)*Uvec(1) 

c     the internal friction force
         tmp = h_inv*(dudx(3)-VxVy(2)*dudx(1))
         sgn_dvdx = sgn(tmp, frict_tiny)
         forceinty=sgn_dvdx*Uvec(1)*kactxy*(g(3)*
     $        dudx(1)+dgdx(1)*Uvec(1))*dsin(intfrictang)

c     the bed friction force for fast moving flow 
         forcebedy=unitvy
     $        *dmax1(g(3)*Uvec(1)+VxVy(2)*Uvec(3)*curv(2),0.0d0)
     $        *tanbed
         if(IF_STOPPED.eq.2.and.1.eq.0) then
c     the bed friction force for stopped or nearly stopped flow

c     the NET force the bed friction force is opposing 
            forcebedequil=forcegrav
c     $           -kactxy*g(3)*Uvec(1)*dUdy(1)
     $           -forceinty

c     the "correct" stopped or nearly stopped flow bed friction force 
c     (this force is not entirely "correct" it will leave a "negligible"
c     (determined by stopping criteria) amount of momentum in the cell
            forcebedy=sgn(forcebedequil,dmin1(forcebedmax,
     $           dabs(forcebedy)+dabs(forcebedequil)))

c     not really 1 but this makes friction statistics accurate
            unitvy=1.0
c         else
         endif

         Ustore(3) = Ustore(3) + dt*(forcegrav -forcebedy -forceinty)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c (erosion terms) this is Camil's logic, Keith changed some variable 
c names for clarity
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         if((.false.).and.(DO_EROSION.ne.0).and.(IF_STOPPED.eq.0)) then
            totalShear=sqrt(forcebedx**2+forcebedy**2)
            if ((totalShear.gt.threshold).and.
     &           (Uvec(1).gt.0.004)) then
            
               es = erosion_rate*dsqrt(dabs(totalShear-threshold))
               eroded=dt*es
               Ustore(1)= Ustore(1) + eroded
               Ustore(2)= Ustore(2) + eroded*VxVy(1)
               Ustore(3)= Ustore(3) + eroded*VxVy(2)
               !write (*,*) 'Doing Keith Erosion Model'
            endif
         endif
         if ((DO_EROSION.ne.0).and.(uvec(1).gt.threshold)) then
            es = erosion_rate*dsqrt(uvec(2)**2+uvec(3)**2)/uvec(1)
            ustore(1) = ustore(1) + dt*es
            ustore(2) = ustore(2) + dt*es*ustore(2)
            ustore(3) = ustore(3) + dt*es*ustore(3)
            !write (*,*) 'Doing Camil Erosion Model'
         endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      endif

c     computation of magnitude of friction forces for statistics
      forceint=unitvx*forceintx+unitvy*forceinty
      forcebed=unitvx*forcebedx+unitvy*forcebedy

c     update the state variables
      Uvec(1) = Ustore(1)
      Uvec(2) = Ustore(2)
      Uvec(3) = Ustore(3)
      
      return
      end subroutine correct

