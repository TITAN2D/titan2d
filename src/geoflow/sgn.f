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
C* $Id: sgn.f 2 2003-08-13 19:26:11Z sorokine $ 
C*

C***********************************************************************
      double precision function sgn(zz, tiny)
C***********************************************************************

      double precision zz, tiny
      
      if (dabs(zz) .lt. tiny) then
         sgn=0.d0
      else if (zz .ge. tiny) then
         sgn=1.d0
      else
         sgn=-1.d0
      endif
      
      return
      end
