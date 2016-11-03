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
      subroutine correct(uvec, uprev, fluxxp, fluxyp, fluxxm, fluxym,
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
