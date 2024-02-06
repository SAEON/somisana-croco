      function o2sato(T,S)
      implicit none
      real    o2sato
      real   T
      real   S
      real A0, A1, A2, A3, A4, A5, B0, B1, B2, B3, C0
      parameter (A0 = 2.00907D0, A1 = 3.22014D0, A2 = 4.05010D0,
     &           A3 = 4.94457D0, A4 = -2.56847D-1, A5 = 3.88767D0,
     &           B0=-6.24523D-3, B1=-7.37614D-3, B2=-1.03410D-2,
     &           B3=-8.17083D-3, C0=-4.88682D-7)
      real    TT
      real    TK
      real    TS, TS2, TS3, TS4, TS5
      real    CO
      TT  = 298.15D0-T
      TK  = 273.15D0+T
      TS  = LOG(TT/TK)
      TS2 = TS**2
      TS3 = TS**3
      TS4 = TS**4
      TS5 = TS**5
      CO  = A0 + A1*TS + A2*TS2 + A3*TS3 + A4*TS4 + A5*TS5
     $     + S*(B0 + B1*TS + B2*TS2 + B3*TS3)
     $     + C0*(S*S)
      o2sato = EXP(CO)
      o2sato = o2sato/22.3916D0*1000.0D0
      return
      end
