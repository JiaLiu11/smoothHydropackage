     
c Setting of global paramters
c
      subroutine params
      implicit none
      real*8 A0, chi
      include 'coms.f'

c     gw = 0.25 fm^-2 width of the gaussian

      logSky = .true.
      logYuk = .true.
      logCb  = .true.
      logPau = .false.

      gw     = 0.25
      sgw    = sqrt(gw)
      Cb0    = 1.44
      Yuk0   = 0.0 !-85.0
      gamYuk = 1.4
      drPau  = 9.0
      dpPau  = 0.0144
      Pau0   = 0.0 !99.5*(hqc/sqrt(drPau*dpPau))**3

C
C hard Skyrme EOS (usual stuff)
C
c      Sky30  = 70.5 
c      gamSky = 2.0 
c      A0     = -124.2 * 0.5
C
C hard Skyrme EOS (JK parametrisation corrected for Gausswidth gw=0.25)
C
       Sky30  = 125.93
       gamSky = 1.676
       A0     = -87.67
       chi    = 0.93 

c      Sky30  = 303.0 !188.18
c      gamSky = 7.0/6.0 !1.457
c      A0 = -356.0 * 0.5

c      print *, 'Chi?'
c      read(5,*) chi

      Sky20 = chi*2.0*A0
      Yuk0 = (1.0-chi)/(2.0*pi*gamYuk**2)*A0

claw
c      print *, Sky20, Yuk0,
c     +         0.5*Sky20 + 2.0*pi*gamYuk**2*Yuk0

c      Sky20 = 0.0d0
c      Sky30 = 0.0d0
c      Yuk0  = 0.0d0
c      Cb0   = 0.0d0

      delr = 0.2
      fdel = delr*delr/6.0
      da = -(1.0/delr)
      db = -da

      cutdww = 20.0
      cutPau = 20.0
      cutYuk = 20.0
      cutCb  = 20.0

      dtimestep=0.2	
      dt  = 0.02
c      dt2 = 0.5*dt
c      dt6 = dt/6.0

      return
      end


c Reset of all indexed variables
c
      subroutine set0
      implicit none
      integer i, j, k
      include 'coms.f'

      do 10 i=1,nspl
        spPauy(i) = 0.0 
        outPau(i) = 0.0 
        spCby(i)  = 0.0  
        outCb(i)  = 0.0
        spYuky(i) = 0.0 
        outYuk(i) = 0.0
        spSkyy(i) = 0.0 
        outSky(i) = 0.0
        spdwwy(i) = 0.0 
        outdww(i) = 0.0
  10  continue

      do 20 j=1,nmax
        spin(j)  = 0
        iso3(j)   = 0
        ncoll(j) = 0
        rx(j)    = 0.0
        ry(j)    = 0.0
        rz(j)    = 0.0
        px(j)    = 0.0
        py(j)    = 0.0
        pz(j)    = 0.0
        fmass(j) = 0.0
        do 30 k=1,nmax
c          iPau(j,k) = .false.
c          rjk(j,k)  = 0.0
c          pjk(j,k)  = 0.0
c          rxjk(j,k)  = 0.0
c          pxjk(j,k)  = 0.0
c          ryjk(j,k)  = 0.0
c          pyjk(j,k)  = 0.0
c          rzjk(j,k)  = 0.0
c          pzjk(j,k)  = 0.0
  30    continue
  20  continue
      return
      end


      subroutine shift
      implicit none
      integer j
      real*8 x_cm, y_cm, z_cm
      real*8 px_cm, py_cm, pz_cm
      include 'coms.f'

c move center of mass to (0, 0, 0)
      x_cm = 0.0
      y_cm = 0.0
      z_cm = 0.0
      px_cm = 0.0
      py_cm = 0.0
      pz_cm = 0.0
      do 10 j=1,npart
        x_cm = x_cm + rx(j)
        y_cm = y_cm + ry(j)
        z_cm = z_cm + rz(j)
        px_cm = px_cm + px(j)
        py_cm = py_cm + py(j)
        pz_cm = pz_cm + pz(j)
   10 continue
      do 20 j=1,npart
        rx(j) = rx(j)-x_cm/dble(npart)
        ry(j) = ry(j)-y_cm/dble(npart)
        rz(j) = rz(j)-z_cm/dble(npart)
        px(j) = px(j)-px_cm/dble(npart) 
        py(j) = py(j)-py_cm/dble(npart)
        pz(j) = pz(j)-pz_cm/dble(npart)
   20 continue
      end

      subroutine angmom
      implicit none
      integer j
      real*8 angx, angy, angz, angtot
      real*8 Jxx, Jxy, Jxz, Jyy, Jyz, Jzz
      include 'coms.f'

      angx = 0.0
      angy = 0.0
      angz = 0.0
      Jxx = 0.0
      Jxy = 0.0
      Jxz = 0.0
      Jyy = 0.0
      Jyz = 0.0
      Jzz = 0.0
      do 10 j=1,npart
        angx = angx + ry(j)*pz(j)-rz(j)*py(j)
        angy = angy + rz(j)*px(j)-rx(j)*pz(j)
        angz = angz + rx(j)*py(j)-ry(j)*px(j)
        Jxx = Jxx + fmass(j)*(ry(j)*ry(j)+rz(j)*rz(j))
        Jxy = Jxy + fmass(j)*rx(j)*ry(j)
        Jxz = Jxz + fmass(j)*rx(j)*rz(j)
        Jyy = Jyy + fmass(j)*(rx(j)*rx(j)+rz(j)*rz(j))
        Jyz = Jyz + fmass(j)*ry(j)*rz(j)
        Jzz = Jzz + fmass(j)*(rx(j)*rx(j)+ry(j)*ry(j))
   10 continue
      angtot = sqrt(angx**2+angy**2+angz**2)/hqc
      print 11, angx/hqc, angy/hqc, angz/hqc, angtot
   11 format('angular momentum', 4F8.3)
c      print *, Jxx, Jxy, Jxz, Jyy, Jyz, Jzz
      end



      subroutine denscheck
      implicit none
      integer j, k, l
      real*8 radius, padius, d_t2, p_t2, rr2
      real*8 x_t(6,0:15), y_t(6,0:15), z_t(6,0:15), rho_r(0:15,6)
      real*8 px_t(6,0:15), py_t(6,0:15), pz_t(6,0:15), rho_p(0:15)
      real*8 nro(0:15), nucrad
      include 'coms.f'

c move center of mass to (0, 0, 0)
      call shift
      
      radius = (0.75*dble(npart)/rho0/pi)**0.33333333
      padius = 0.3/radius

      rr2 = 0.0
      do 80 l=1,npart
        rr2 = rr2+rx(l)**2+ry(l)**2+rz(l)**2
 80   continue
      rr2 = sqrt(rr2/dble(npart)+0.75/gw)
      print *, 'RMS-Radius:', rr2
      print *, 'nucrad:    ', nucrad(npart)    
      do 10 j=0,15
        nro(j) = 0.0
        x_t(1,j) = dble(j)*0.1*radius      
        y_t(1,j) = 0.0
        z_t(1,j) = 0.0
        x_t(2,j) = -(dble(j)*0.1*radius)     
        y_t(2,j) = 0.0
        z_t(2,j) = 0.0
        x_t(3,j) = 0.0
        y_t(3,j) = dble(j)*0.1*radius      
        z_t(3,j) = 0.0
        x_t(4,j) = 0.0
        y_t(4,j) = -(dble(j)*0.1*radius)    
        z_t(4,j) = 0.0
        x_t(5,j) = 0.0
        y_t(5,j) = 0.0
        z_t(5,j) = dble(j)*0.1*radius      
        x_t(6,j) = 0.0
        y_t(6,j) = 0.0
        z_t(6,j) = -(dble(j)*0.1*radius)
        do 20 k=1,6
          px_t(k,j) = x_t(k,j)*padius
          py_t(k,j) = y_t(k,j)*padius
          pz_t(k,j) = z_t(k,j)*padius
   20   continue
   10 continue
      do 30 k=0,15
        do 33 l=1,6
          rho_r(k,l) = 0.0
   33   continue
        rho_p(k) = 0.0
        do 40 l=1,npart
          if (rx(l)**2+ry(l)**2+rz(l)**2.le.x_t(1,k)**2) 
     +       nro(k) = nro(k)+1.0
          do 50 j=1,6
            d_t2 = (rx(l)-x_t(j,k))**2 + (ry(l)-y_t(j,k))**2 +
     +             (rz(l)-z_t(j,k))**2
            rho_r(k,j) = rho_r(k,j) + exp(-(2.0*gw*d_t2))
            p_t2 = (px(l)-px_t(j,k))**2 + (py(l)-py_t(j,k))**2 +
     +             (pz(l)-pz_t(j,k))**2
            rho_p(k) = rho_p(k) + exp(-(0.5/gw/hqc/hqc*p_t2))
   50     continue
   40   continue
        do 44 l=1,6
          rho_r(k,l) = rho_r(k,l)*(2.0*gw/pi)**1.5 ! fm^-3
   44   continue
        rho_p(k) = rho_p(k)/6.0 * (0.5/gw/pi)**1.5 ! fm^3
        nro(k) = nro(k)*0.75/pi/x_t(1,k)**3
        print 15, dble(k)*0.1*radius, (rho_r(k,l), l=1,6)
        write(88,15) dble(k)*0.1*radius, (rho_r(k,l), l=1,6)
c        print 14, k, x_t(1,k), rho_r(k), nro(k), px_t(1,k), rho_p(k)
c        write(99,14) k, x_t(1,k), rho_r(k), nro(k), px_t(1,k), rho_p(k)
c   14   format(i4,5f10.5)
   15   format(7f10.5)
   30 continue
      return
      end



      subroutine derivs(row)
      implicit none
      integer j, k, index, row
      real*8 spu, spo, outu, outo, tmp, a, b, dy, dp, drdp, dpj
      real*8 rxjku, ryjku, rzjku, rjku, pxjku, pyjku, pzjku, pjku
      logical iPau
      include 'coms.f'
  
C aopx(j,?) = -dH/dx_j
C aorx(j,?) = dH/dpx_j

      do 10 j=1,nbar
        rww(j) = 0.0
   10 continue  
      if(logSky) then 
      do 20 j=1,nbar
        do 30 k=j+1,nbar
          rxjku = (airx(j)-airx(k))
          ryjku = (airy(j)-airy(k))
          rzjku = (airz(j)-airz(k))
          rjku = sqrt(rxjku**2+ryjku**2+rzjku**2)
          if(rjku.lt.cutdww) then
            index = int(rjku/delr)+1
            a = dble(index) - rjku/delr
            b = 1.0 - a 
            tmp = a*spdwwy(index)+b*spdwwy(index+1)+
     +       ((a**3-a)*outdww(index)+(b**3-b)*outdww(index+1))*fdel
            rww(j) = rww(j) + tmp
            rww(k) = rww(k) + tmp
          end if
  30    continue
  20  continue
      end if

      do 40 j=1,nbar
        aopx(j,row) = 0.0
        aopy(j,row) = 0.0
        aopz(j,row) = 0.0
        dpj = 1.0/sqrt(aipx(j)*aipx(j)+aipy(j)*aipy(j)+aipz(j)*aipz(j)+
     +             fmass(j)*fmass(j))
        aorx(j,row) = aipx(j)*dpj
        aory(j,row) = aipy(j)*dpj
        aorz(j,row) = aipz(j)*dpj
  40  continue

      do 50 j=1,nbar
        do 60 k=j+1,nbar
          rxjku = (airx(j)-airx(k))
          ryjku = (airy(j)-airy(k))
          rzjku = (airz(j)-airz(k))
          rjku = sqrt(rxjku**2+ryjku**2+rzjku**2)
          if (rjku.ge.1.0E-8) then
            rxjku = rxjku/rjku
            ryjku = ryjku/rjku
            rzjku = rzjku/rjku
          else
            rxjku = 0.0
            ryjku = 0.0
            rzjku = 0.0
          end if
          spu  = 0.0
          spo  = 0.0
          outu = 0.0
          outo = 0.0
          dy = 0.0
          index = int(rjku/delr)+1
          a = dble(index)-rjku/delr
          b = 1.0-a
          if(logYuk.and.rjku.lt.cutYuk) then
            spu  = spu  + spYuky(index)
            spo  = spo  + spYuky(index+1)
            outu = outu + outYuk(index)
            outo = outo + outYuk(index+1)
          end if 
          if(logSky.and.rjku.lt.cutdww) then
            tmp = Sky20 + Sky30*gamSky/(gamSky+1.0)*
     *            (rww(j)**(gamSky-1.0)+rww(k)**(gamSky-1.0))
            spu  = spu  + spdwwy(index)*tmp
            spo  = spo  + spdwwy(index+1)*tmp
            outu = outu + outdww(index)*tmp
            outo = outo + outdww(index+1)*tmp
          end if    
          dy = da*spu+db*spo+
     +      ((3.0*a**2-1.0)*da*outu+(3.0*b**2-1.0)*db*outo)*fdel
          if(logCb) then
          if(rjku.lt.cutCb) then
            dy = dy + (da*spCby(index)+db*spCby(index+1)+
     +                ((3.0*a**2-1.0)*da*outCb(index)+
     +                (3.0*b**2-1.0)*db*outCb(index+1))*fdel)*
     *                dble(charge(j)*charge(k))
          else 
            dy = dy - Cb0/rjku/rjku*dble(charge(j)*charge(k))
          end if
          end if
          if(logPau.and.iPau(j,k)) then
            pxjku = (aipx(j)-aipx(k))
            pyjku = (aipy(j)-aipy(k))
            pzjku = (aipz(j)-aipz(k))
            pjku = sqrt(pxjku**2+pyjku**2+pzjku**2)
            if (pjku.ge.1.0E-8) then
              pxjku = pxjku/pjku
              pyjku = pyjku/pjku
              pzjku = pzjku/pjku
            else
              pxjku = 0.0
              pyjku = 0.0
              pzjku = 0.0
            end if
            drdp = 0.5*(pjku*pjku/dpPau+rjku*rjku/drPau)
            if(drdp.lt.cutPau) then
              index = int(drdp/delr)+1   
              a = dble(index)-drdp/delr
              b = 1.0-a
              tmp = da*spPauy(index)+db*spPauy(index+1)+
     +              ((3.0*a**2-1.0)*da*outPau(index)+
     +               (3.0*b**2-1.0)*db*outPau(index+1))*fdel
              dy = dy+tmp*rjku/drPau
              dp = tmp*pjku/dpPau*0.001
              aorx(j,row) = aorx(j,row)+dp*pxjku
              aory(j,row) = aory(j,row)+dp*pyjku
              aorz(j,row) = aorz(j,row)+dp*pzjku
              aorx(k,row) = aorx(k,row)-dp*pxjku
              aory(k,row) = aory(k,row)-dp*pyjku
              aorz(k,row) = aorz(k,row)-dp*pzjku
            end if
          end if 
          dy = -(0.001*dy)
          aopx(j,row) = aopx(j,row)+dy*rxjku
          aopy(j,row) = aopy(j,row)+dy*ryjku
          aopz(j,row) = aopz(j,row)+dy*rzjku
          aopx(k,row) = aopx(k,row)-dy*rxjku
          aopy(k,row) = aopy(k,row)-dy*ryjku
          aopz(k,row) = aopz(k,row)-dy*rzjku
  60    continue
  50  continue
      return
      end

      real*8 function Ekintot() 
      implicit none
      integer j
      real*8 Ekin
      include 'coms.f'

      Ekintot = 0.0
      do 3 j=1,npart
         Ekintot= Ekintot+Ekin(j)
 3    continue
      return
      end

      real*8 function EtotJK() 
      implicit none
      real*8 Etot
      integer j
      include 'coms.f'

      EtotJK = Etot()
      do 3 j=1,npart
         EtotJK= EtotJK-fmass(j)
 3    continue
      EtotJK = EtotJK/npart
      return
      end


      real*8 function Etot() 
      implicit none
      integer j, k, index
claw
c      real*8 Etot, a, b, y, drdp, tp, tr, tmp, Ekintot
      real*8 a, b, y, drdp, tp, tr, tmp, Ekintot
      real*8 Ekinbar, Ekinmes, ESky2, ESky3, EYuk, ECb, EPau, Ekin
      real*8 rxjku, ryjku, rzjku, rjku, pxjku, pyjku, pzjku, pjku
      logical iPau
      include 'coms.f'
      common /energies/ Ekinbar, Ekinmes, ESky2, ESky3, EYuk, ECb, EPau

      Etot = 0.0
      Ekinbar = 0.0
      Ekinmes = 0.0
      ESky2 = 0.0
      ESky3 = 0.0
      EYuk = 0.0
      ECb = 0.0
      EPau = 0.0

c2811      call dists(rx,ry,rz,px,py,pz)
      if(EoS.eq.0) then
csab
c CASCADE mode
         Etot=Ekintot()
         return
      else
c with potentials
c kinetic energies of mesons first
         do 4 j=nbar+1,npart
            Etot= Etot+Ekin(j)
            Ekinmes = Ekinmes+Ekin(j)
 4       continue

      do 10 j=1,nbar
        rww(j) = 0.0
   10 continue  
      if(logSky) then 
      do 20 j=1,nbar
        do 30 k=j+1,nbar
          rxjku = (rx(j)-rx(k))
          ryjku = (ry(j)-ry(k))
          rzjku = (rz(j)-rz(k))
          rjku = sqrt(rxjku**2+ryjku**2+rzjku**2)
          if(rjku.lt.cutdww) then
            index = int(rjku/delr)+1
            a = dble(index) - rjku/delr
            b = 1.0 - a 
            tmp = a*spdwwy(index)+b*spdwwy(index+1)+
     +       ((a**3-a)*outdww(index)+(b**3-b)*outdww(index+1))*fdel
            rww(j) = rww(j) + tmp
            rww(k) = rww(k) + tmp
          end if
  30    continue
  20  continue
      end if

      do 40 j=1,nbar
        Etot = Etot + Ekin(j) + 0.0005*Sky20*rww(j) +
     +         0.001*Sky30/(gamSky+1.0)*rww(j)**gamSky
        Ekinbar = Ekinbar + Ekin(j)
        ESky2 = ESky2 + 0.0005*Sky20*rww(j)
        ESky3 = ESky3 + 0.001*Sky30/(gamSky+1.0)*rww(j)**gamSky
        do 50 k=j+1,nbar
          rxjku = (rx(j)-rx(k))
          ryjku = (ry(j)-ry(k))
          rzjku = (rz(j)-rz(k))
          rjku = sqrt(rxjku**2+ryjku**2+rzjku**2)
          index = int(rjku/delr)+1
          a = dble(index)-rjku/delr
          b = 1.0-a
          if(logYuk.and.rjku.lt.cutYuk) then
            y = a*spYuky(index)+b*spYuky(index+1)+
     +          ((a**3-a)*outYuk(index)+(b**3-b)*outYuk(index+1))*fdel
            Etot = Etot + 0.001*y
            EYuk = EYuk + 0.001*y
          end if
          if(logCb) then
          if(rjku.lt.cutCb) then
            y = (a*spCby(index)+b*spCby(index+1)+
     +          ((a**3-a)*outCb(index)+(b**3-b)*outCb(index+1))*fdel)*
     *          dble(charge(j)*charge(k))
          else
            y = Cb0/rjku*dble(charge(j)*charge(k))
          end if
          Etot = Etot + 0.001*y
          ECb = ECb + 0.001*y
          end if
          if(logPau.and.iPau(j,k)) then
            pxjku = (px(j)-px(k))
            pyjku = (py(j)-py(k))
            pzjku = (pz(j)-pz(k))
            pjku = sqrt(pxjku**2+pyjku**2+pzjku**2)
            tp = pjku
            tr = rjku
            drdp = 0.5*(pjku*pjku/dpPau+rjku*rjku/drPau)
            if(drdp.lt.cutPau) then
              index = int(drdp/delr)+1
              a = dble(index)-drdp/delr
              b = 1.0-a
              y = a*spPauy(index)+b*spPauy(index+1)+
     +            ((a**3-a)*outPau(index)+(b**3-b)*outPau(index+1))*fdel
              Etot = Etot + 0.001*y
              EPau = EPau + 0.001*y
            end if
          end if
  50    continue
  40  continue
cJFK      Etot = Etot/dble(nbar)
      Ekinbar = Ekinbar/dble(nbar)
      Ekinmes = Ekinmes/dble(max(1,npart-nbar))
      ESky2 = ESky2/dble(nbar)
      ESky3 = ESky3/dble(nbar)
      EYuk = EYuk/dble(nbar)
      ECb = ECb/dble(nbar)
      EPau = EPau/dble(nbar)
      end if
      return
      end
          

      subroutine cascstep(tim,dtime)
      implicit none
      real*8 tim,dtime,energ
c     real*8 swapi
      integer j
      include 'coms.f'
      include 'boxinc.f'
      include 'options.f'

cdebug
c      write(6,'(a8,2f10.3)') '    casc', tim, dtime
c	write(6,*) 'tim,dtime = ',tim,dtime,npart

      do 1 j=1,npart
         energ = sqrt(px(j)**2+py(j)**2+pz(j)**2+fmass(j)**2)
         r0(j) = r0(j) + dtime
         rx(j) = rx(j) + px(j)/energ*dtime
         rz(j) = rz(j) + pz(j)/energ*dtime
         ry(j) = ry(j) + py(j)/energ*dtime
1     continue          
      return
      end
 
      subroutine proprk(tim,dtime)
c Modifikationen von M. Hofmann: CASCADE und variable Propagationszeit
      implicit none
      real*8 tim,dtime,energ,dt2, dt6
c     real*8 swapi
      integer j
      include 'coms.f'
      include 'boxinc.f'
      include 'options.f'

c      write(6,'(a8,2f10.3)') 'prop    ', tim, dtime

      if (EoS.eq.0) then
c  cascade mode
c	write(6,*) 'tim,dtime = ',tim,dtime,npart
         do 1 j=1,npart
               energ = p0(j)    ! sqrt(px(j)**2+py(j)**2+pz(j)**2+fmass(j)**2)
               r0(j) = r0(j) + dtime
               rx(j) = rx(j) + px(j)/energ*dtime
               ry(j) = ry(j) + py(j)/energ*dtime
               rz(j) = rz(j) + pz(j)/energ*dtime
1        continue          
         return
      else
c propagation with potentials
c propagate mesons on straight lines
         do 2 j=nbar+1,npart
            energ = p0(j) ! sqrt(px(j)**2+py(j)**2+pz(j)**2+fmass(j)**2)
            r0(j) = r0(j) + dtime
            rx(j) = rx(j) + px(j)/energ*dtime
            ry(j) = ry(j) + py(j)/energ*dtime
            rz(j) = rz(j) + pz(j)/energ*dtime
2        continue          

c propagate baryons
c adjust time-step parameters
      dt = dtime
      dt2 = dtime/2.0d0
      dt6 = dtime/6.0d0

      do 10 j=1,nbar
        airx(j) = rx(j) 
        airy(j) = ry(j) 
        airz(j) = rz(j) 
        aipx(j) = px(j)
        aipy(j) = py(j)
        aipz(j) = pz(j)
  10  continue
        
c2811      call dists(rx,ry,rz,px,py,pz)
      call derivs(1)
      do 20 j=1,nbar
        airx(j) = rx(j) + dt2*aorx(j,1) 
        airy(j) = ry(j) + dt2*aory(j,1) 
        airz(j) = rz(j) + dt2*aorz(j,1) 
        aipx(j) = px(j) + dt2*aopx(j,1)
        aipy(j) = py(j) + dt2*aopy(j,1)
        aipz(j) = pz(j) + dt2*aopz(j,1)
  20  continue

c2811      call dists(airx,airy,airz,aipx,aipy,aipz)
      call derivs(2)
        
      do 30 j=1,nbar
        airx(j) = rx(j) + dt2*aorx(j,2) 
        airy(j) = ry(j) + dt2*aory(j,2) 
        airz(j) = rz(j) + dt2*aorz(j,2) 
        aipx(j) = px(j) + dt2*aopx(j,2)
        aipy(j) = py(j) + dt2*aopy(j,2)
        aipz(j) = pz(j) + dt2*aopz(j,2)
  30  continue
      
c2811      call dists(airx,airy,airz,aipx,aipy,aipz)
      call derivs(3)

      do 40 j=1,nbar
        airx(j) = rx(j) + dt*aorx(j,3)
        airy(j) = ry(j) + dt*aory(j,3)
        airz(j) = rz(j) + dt*aorz(j,3)
        aipx(j) = px(j) + dt*aopx(j,3)
        aipy(j) = py(j) + dt*aopy(j,3)
        aipz(j) = pz(j) + dt*aopz(j,3)
  40  continue
   
c2811      call dists(airx,airy,airz,aipx,aipy,aipz)
      call derivs(4)

      do 50 j=1,nbar
         r0(j) = r0(j) + dtime
        rx(j)=rx(j)+dt6*(aorx(j,1)+2.0*(aorx(j,2)+aorx(j,3))+aorx(j,4))
        ry(j)=ry(j)+dt6*(aory(j,1)+2.0*(aory(j,2)+aory(j,3))+aory(j,4))
        rz(j)=rz(j)+dt6*(aorz(j,1)+2.0*(aorz(j,2)+aorz(j,3))+aorz(j,4))
        px(j)=px(j)+dt6*(aopx(j,1)+2.0*(aopx(j,2)+aopx(j,3))+aopx(j,4))
        py(j)=py(j)+dt6*(aopy(j,1)+2.0*(aopy(j,2)+aopy(j,3))+aopy(j,4))
        pz(j)=pz(j)+dt6*(aopz(j,1)+2.0*(aopz(j,2)+aopz(j,3))+aopz(j,4))
        p0(j)=sqrt(px(j)**2+py(j)**2+pz(j)**2+fmass(j)**2)
 50   continue
      end if
c MB
c            if ((mbox.gt.0).and.(mbflag.eq.1)) then
c              write(*,*) 'Hallo, Welt3 !'
c		rx(j)=swapi(rx(j),lboxhalbe)
c		ry(j)=swapi(ry(j),lboxhalbe)
c		rz(j)=swapi(rz(j),lboxhalbe)
c            endif
cc 
      return
      end


      subroutine potPau
      implicit none
      integer i, ncut, index
      real*8 Ecut, dr, abl0, abln, a, b, y, dy, Pau
      include 'coms.f'

      rx(1) = 0.0d0
      ry(1) = 0.0d0
      rz(1) = 0.0d0
      ry(2) = 0.0d0
      rz(2) = 0.0d0
      px(1) = 0.0d0
      py(1) = 0.0d0
      pz(1) = 0.0d0
      px(1) = 0.0d0
      py(2) = 0.0d0
      pz(2) = 0.0d0
      Ecut = 1.0E-5
      i = 0
  99  i = i+1    
      dr = delr*dble(i-1)
      rx(2) = sqrt(2.0*dr*drPau)
      spx(i) = dr
      spPauy(i) = Pau(1,2)
      if(spPauy(i).lt.Ecut) then
        spPauy(i) = 0.0
        cutPau = dr
        abl0 = -Pau0 
        abln = 0.0
        ncut = i
      else
        goto 99
      end if
      call spline(spx,spPauy,ncut,abl0,abln,outPau)

      write(6,'(''Pauli-Potential    '',e10.3,i5,f7.1)') 
     +      Ecut, ncut, cutPau

      do 10 i=0,20
        dr = 0.323*dble(i)
        if(dr.gt.cutPau) then
          y = 0.0
          dy = 0.0
        else
        rx(2) = dr
        dr = 0.5*dr*dr/drPau
        index = int(dr/delr)+1
        a = dble(index) - dr/delr
        b = 1.0 - a
        y = a*spPauy(index)+b*spPauy(index+1)+
     +     ((a**3-a)*outPau(index)+
     +      (b**3-b)*outPau(index+1))*fdel
        dy = da*spPauy(index)+db*spPauy(index+1)+
     +        ((3.0*a**2-1.0)*da*outPau(index)+
     +         (3.0*b**2-1.0)*db*outPau(index+1))*fdel
        dy = dy*sqrt(2.0*dr*drPau)/drPau
        end if
c        write(6,'(1p,4e15.5)') dr, dPaur(1,1), dy, dPaur(1,1)-dy
c        write(6,'(1p,4e15.5)') dr, Pau(1,1), y, Pau(1,1)-y
  10  continue
      return
      end

      subroutine potCb
      implicit none
      integer i, ncut, index
      real*8 Ecut, dr, abl0, abln, a, b, y, dy, dCb, Cb
      include 'coms.f'

      rx(1) = 0.0d0
      ry(1) = 0.0d0
      rz(1) = 0.0d0
      ry(2) = 0.0d0
      rz(2) = 0.0d0
      Ecut = 1.0E-5
      iso3(1) = 1
      iso3(2) = 1
      i = 0
  99  i = i+1    
      dr = delr*dble(i-1)
      rx(2) = dr
      spx(i) = dr
      spCby(i) = Cb(1,2)
      if(abs(spCby(i)*dr-Cb0)/max(dr,1.0d-5).lt.Ecut) then
        spCby(i) = Cb0/dr
        cutCb  = dr
        abln = dCb(1,2)
        abl0 = 0.0 
        ncut = i
      else
        goto 99
      end if
      call spline(spx,spCby,ncut,abl0,abln,outCb)

      write(6,'(''Coulomb-Potential  '',e10.3,i5,f7.1)')
     +      Ecut, ncut, cutCb

      do 10 i=0,20
        dr = 0.2*dble(i)+0.01212
        rx(2) = dr
        if(dr.ge.cutCb) then
          y = Cb0/dr
          dy = -(Cb0/dr/dr)
        else
          index = int(dr/delr)+1
          a = dble(index) - dr/delr
          b = 1.0 - a
          y = a*spCby(index)+b*spCby(index+1)+
     +       ((a**3-a)*outCb(index)+(b**3-b)*outCb(index+1))*fdel
          dy = da*spCby(index)+db*spCby(index+1)+
     +         ((3.0*a**2-1.0)*da*outCb(index)+
     +           (3.0*b**2-1.0)*db*outCb(index+1))*fdel
        end if
c        write(6,'(1p,4e15.5)') dr, dCb(1,2), dy, dCb(1,2)-dy
c        write(6,'(1p,4e15.5)') dr, Cb(1,2), y, Cb(1,2)-y
  10  continue
      return
      end

      subroutine potYuk
      implicit none
      integer i, ncut, index
      real*8 Ecut, dr, abl0, abln, a, b, y, dy
      real*8 Yuk
      include 'coms.f'

      rx(1) = 0.0d0
      ry(1) = 0.0d0
      rz(1) = 0.0d0
      ry(2) = 0.0d0
      rz(2) = 0.0d0
      Ecut = 1.0E-5
      i = 0
  99  i = i+1    
      dr = delr*dble(i-1)
      rx(2) = dr 
      spx(i) = dr
      spYuky(i) = Yuk(1,2)
      if(abs(spYuky(i)).lt.Ecut) then
        spYuky(i) = 0.0
        cutYuk = dr
        abl0 = 0.0 
        abln = 0.0
        ncut = i
      else
        goto 99
      end if
      call spline(spx,spYuky,ncut,abl0,abln,outYuk)

      write(6,'(''Yukawa-Potential   '',e10.3,i5,f7.1)') 
     +      Ecut, ncut, cutYuk

      do 10 i=0,40
        dr = 0.2*dble(i)
        rx(2) = dr
        if(dr.gt.cutYuk) then
          y = 0.0
          dy = 0.0
        else
        index = int(dr/delr)+1
        a = dble(index) - dr/delr
        b = 1.0 - a
        y = a*spYuky(index)+b*spYuky(index+1)+
     +     ((a**3-a)*outYuk(index)+(b**3-b)*outYuk(index+1))*fdel
        dy = da*spYuky(index)+db*spYuky(index+1)+
     +        ((3.0*a**2-1.0)*da*outYuk(index)+
     +         (3.0*b**2-1.0)*db*outYuk(index+1))*fdel
        end if
c        write(6,'(1p,4e15.5)') dr, dYuk(1,2), dy, dYuk(1,2)-dy
c        write(6,'(1p,4e15.5)') dr, Yuk(1,2), y, Yuk(1,2)-y
  10  continue
      return
      end

c      subroutine potSky
c      implicit none
c      integer i, ncut, index
c      real*8 Ecut, dr, abl0, abln, a, b, y, dy
c      include 'coms.f'
c
c      Ecut = 1.0E-5
c      i = 0
c  99  i = i+1    
c      dr = delr*dble(i-1)
c      rjk(1,1) = dr 
c      spx(i) = dr
c      spSkyy(i) = Sky(1,1)
c      if(abs(spSkyy(i)).lt.Ecut) then
c        spSkyy(i) = 0.0
c        cutSky = dr
c        abl0 = 0.0 
c        abln = 0.0
c        ncut = i
c      else
c        goto 99
c      end if
c      call spline(spx,spSkyy,ncut,abl0,abln,outSky)
c
c      write(6,*) Ecut, ncut, cutSky
c
c      do 10 i=0,20
c        dr = 0.295*dble(i)
c        rjk(1,1) = dr
c        if(dr.gt.cutSky) then
c          y = 0.0
c          dy = 0.0
c        else
c        index = int(dr/delr)+1
c        a = dble(index) - dr/delr
c        b = 1.0 - a
c        y = a*spSkyy(index)+b*spSkyy(index+1)+
c     +     ((a**3-a)*outSky(index)+(b**3-b)*outSky(index+1))*fdel
c        dy = da*spSkyy(index)+db*spSkyy(index+1)+
c     +        ((3.0*a**2-1.0)*da*outSky(index)+
c     +         (3.0*b**2-1.0)*db*outSky(index+1))*fdel
c        end if
cc        write(6,'(1p,4e15.5)') dr, dSky(1,1), dy, dSky(1,1)-dy
c        write(6,'(1p,4e15.5)') dr, Sky(1,1), y, Sky(1,1)-y
c  10  continue
c      return
c      end


      subroutine potdww
      implicit none
      integer i, ncut, index
      real*8 Ecut, dr, abl0, abln, a, b, y, dy
      real*8 dww
      include 'coms.f'

      rx(1) = 0.0d0
      ry(1) = 0.0d0
      rz(1) = 0.0d0
      ry(2) = 0.0d0
      rz(2) = 0.0d0
      Ecut = 1.0E-8
      i = 0
  99  i = i+1    
      dr = delr*dble(i-1)
      rx(2) = dr 
      spx(i) = dr
      spdwwy(i) = dww(1,2)
      if(abs(spdwwy(i)).lt.Ecut) then
        spdwwy(i) = 0.0
        cutdww = dr
        abl0 = 0.0 
        abln = 0.0
        ncut = i
      else
        goto 99
      end if
      call spline(spx,spdwwy,ncut,abl0,abln,outdww)

      write(6,'(''Interaction-Density'',e10.3,i5,f7.1)') 
     +      Ecut, ncut, cutdww

      do 10 i=0,20
        dr = 0.295*dble(i)
        rx(2) = dr
        if(dr.gt.cutdww) then
          y = 0.0
          dy = 0.0
        else
        index = int(dr/delr)+1
        a = dble(index) - dr/delr
        b = 1.0 - a
        y = a*spdwwy(index)+b*spdwwy(index+1)+
     +     ((a**3-a)*outdww(index)+(b**3-b)*outdww(index+1))*fdel
        dy = da*spdwwy(index)+db*spdwwy(index+1)+
     +        ((3.0*a**2-1.0)*da*outdww(index)+
     +         (3.0*b**2-1.0)*db*outdww(index+1))*fdel
        end if
c        write(6,'(1p,4e15.5)') dr, dSky(1,1), dy, dSky(1,1)-dy
c        write(6,'(1p,4e15.5)') dr, dww(1,1), y, dww(1,1)-y
  10  continue
      return
      end

c Evaluation of rjk, pjk and rww
c
c      subroutine dists(aarx,aary,aarz,aapx,aapy,aapz)
c      implicit none
c      include 'coms.f'      
c      integer j, k, index
c      real*8 aarx(nmax),aary(nmax),aarz(nmax),
c     +     aapx(nmax),aapy(nmax),aapz(nmax),
c     +     a, b, tmp, dpj
c
c      do 10 j=1,nbar
c        do 20 k=j+1,nbar
c          rjk(j,k) = sqrt((aarx(j)-aarx(k))*(aarx(j)-aarx(k))+
c     +                    (aary(j)-aary(k))*(aary(j)-aary(k))+
c     +                    (aarz(j)-aarz(k))*(aarz(j)-aarz(k)))
c          pjk(j,k) = sqrt((aapx(j)-aapx(k))*(aapx(j)-aapx(k))+
c     +                    (aapy(j)-aapy(k))*(aapy(j)-aapy(k))+
c     +                    (aapz(j)-aapz(k))*(aapz(j)-aapz(k)))
c          rjk(k,j) = rjk(j,k)
c          pjk(k,j) = pjk(j,k)
c          if (rjk(j,k).ge.1.0E-8) then
c            rxjk(j,k) = (aarx(j)-aarx(k))/rjk(j,k)
c            ryjk(j,k) = (aary(j)-aary(k))/rjk(j,k)
c            rzjk(j,k) = (aarz(j)-aarz(k))/rjk(j,k)
c          else
c            rxjk(j,k) = 0.0
c            ryjk(j,k) = 0.0
c            rzjk(j,k) = 0.0
c          end if
c          if (pjk(j,k).ge.1.0E-8) then
c            pxjk(j,k) = (aapx(j)-aapx(k))/pjk(j,k)
c            pyjk(j,k) = (aapy(j)-aapy(k))/pjk(j,k)
c            pzjk(j,k) = (aapz(j)-aapz(k))/pjk(j,k)
c          else
c            pxjk(j,k) = 0.0
c            pyjk(j,k) = 0.0
c            pzjk(j,k) = 0.0
c          end if
c          rxjk(k,j) = -rxjk(j,k)
c          ryjk(k,j) = -ryjk(j,k)
c          rzjk(k,j) = -rzjk(j,k)
c          pxjk(k,j) = -pxjk(j,k)
c          pyjk(k,j) = -pyjk(j,k)
c          pzjk(k,j) = -pzjk(j,k)
c  20    continue
c  10  continue
c      do 30 j=1,nbar
c        dpj = 1.0/sqrt(aapx(j)*aapx(j)+aapy(j)*aapy(j)+aapz(j)*aapz(j)+
c     +             fmass(j)*fmass(j)) 
c        pxjk(j,j) = aapx(j)*dpj
c        pyjk(j,j) = aapy(j)*dpj
c        pzjk(j,j) = aapz(j)*dpj
c        rww(j) = 0.0
c        rtest(j) = 1.0
c  30  continue
c      do 40 j=1,nbar
c        do 50 k=j+1,nbar
c          rtest(j) = rtest(j) + exp(-2.0*gw*rjk(j,k)**2)
c          rtest(k) = rtest(k) + exp(-2.0*gw*rjk(j,k)**2)
c          if(rjk(j,k).lt.cutdww) then
c            index = int(rjk(j,k)/delr)+1
c            a = dble(index) - rjk(j,k)/delr
c            b = 1.0 - a 
c            tmp = a*spdwwy(index)+b*spdwwy(index+1)+
c     +       ((a**3-a)*outdww(index)+(b**3-b)*outdww(index+1))*fdel
c            rww(j) = rww(j) + tmp
c            rww(k) = rww(k) + tmp
c          end if
c  50    continue
c  40  continue
c      do 60 j=1,nbar
c        rtest(j) = rtest(j)*(2.0*gw/pi)**1.5
c  60  continue  
c      return
c      end


c Kinetic Energy
c
      function Ekin(j)
      implicit none
      integer j
      real*8 Ekin
      include 'coms.f'
 
      Ekin = sqrt((px(j)+ffermpx(j))*(px(j)+ffermpx(j))+
     +            (py(j)+ffermpy(j))*(py(j)+ffermpy(j))+
     +            (pz(j)+ffermpz(j))*(pz(j)+ffermpz(j))+
     +            fmass(j)*fmass(j)) 

      return
      end

c Derivative for Kinetic Energy
c
      function dEkin(j)
      implicit none
      integer j
      real*8 dEkin
      include 'coms.f'
      
      dEkin = 1.0/sqrt(px(j)*px(j)+py(j)*py(j)+pz(j)*pz(j)+
     +                 fmass(j)*fmass(j)) 
      return
      end

c Skyrme Potential (3-body) rwwterm
c
      function dww(j,k)
      implicit none
      integer j, k
      real*8 dww, rjk
      include 'coms.f'
      
      dww = gw/pi*sqrt(gw/pi)*exp(-(gw*rjk(j,k)*rjk(j,k)))/
     /      rho0
c     /      (rho0-(gw/pi)**1.5)
      return
      end

cc Skyrme Potential
cc
      function Sky(j,k)
      implicit none
      integer j, k
      real*8 Sky, rjk
      include 'coms.f'
      
      Sky = Sky20*gw/pi*sqrt(gw/pi)*exp(-(gw*rjk(j,k)*rjk(j,k))) 
      return
      end
c
cc Derivative for Skyrme Potential
cc
c      function dSky(j,k)
c      implicit none
c      integer j, k
c      include 'coms.f'
c      
c      dSky = -Sky20*gw/pi*sqrt(gw/pi)*exp(-gw*rjk(j,k)*rjk(j,k))*
c     *       2.0*gw*rjk(j,k)
c      return
c      end

c Coulomb Potential
c
      function Cb(j,k)
      implicit none
      integer j, k
      real*8 Cb, rjk
      real*8 erf
      include 'coms.f'
      
      if (iso3(j).eq.1.and.iso3(k).eq.1) then
        if (rjk(j,k).lt.eps) then
          Cb = Cb0*er0*sgw
        else
          Cb = Cb0/rjk(j,k)*erf(sgw*rjk(j,k))
        end if
      else 
        Cb = 0.0
      end if 
      return
      end

c Derivative for Coulomb Potential
c
      function dCb(j,k)
      implicit none
      integer j, k
      real*8 dCb, rjk
      real*8 erf
      include 'coms.f'
      
      if (iso3(j).eq.1.and.iso3(k).eq.1) then
        if (rjk(j,k).lt.eps) then
          dCb = 0.0
        else
          dCb = Cb0*(er0*exp(-(gw*rjk(j,k)*rjk(j,k)))*sgw*rjk(j,k)-
     +               erf(sgw*rjk(j,k)))/rjk(j,k)/rjk(j,k)
        end if
      else 
        dCb = 0.0
      end if 
      return
      end

c Yukawa Potential
c
      function Yuk(j,k)
      implicit none
      integer j, k
      real*8 Yuk, rjk
      real*8 erf
      include 'coms.f'
      
      if(rjk(j,k).lt.eps) then
        Yuk = Yuk0*(er0*sgw-exp(0.25/gamYuk/gamYuk/gw)/gamYuk*
     *              (1.0-erf(0.5/gamYuk/sgw)))
      else
        Yuk = Yuk0*0.5/rjk(j,k)*exp(0.25/gamYuk/gamYuk/gw)*
     *           (exp(-(rjk(j,k)/gamYuk))*
     +            (1.0-erf(0.5/gamYuk/sgw-sgw*rjk(j,k)))-
     -            exp(rjk(j,k)/gamYuk)*
     +            (1.0-erf(0.5/gamYuk/sgw+sgw*rjk(j,k))))
      end if
      return
      end

c Derivative for Yukawa Potential
c
      function dYuk(j,k)
      implicit none
      integer j, k
      real*8 dYuk, rjk
      real*8 erf
      include 'coms.f'
      
      if(rjk(j,k).lt.eps) then
        dYuk = 0.0
      else
        dYuk = 0.5*Yuk0/rjk(j,k)*( exp(0.25/gamYuk/gamYuk/gw)*(
     *          (-(1.0/rjk(j,k))-1.0/gamYuk)*exp(-(rjk(j,k)/gamYuk))*
     *             (1.0-erf(0.5/gamYuk/sgw-sgw*rjk(j,k))) +
     *          (1.0/rjk(j,k)-1.0/gamYuk)*exp(rjk(j,k)/gamYuk)*
     *             (1.0-erf(0.5/gamYuk/sgw+sgw*rjk(j,k))) ) +
     +          sgw*er0*2.0*exp(-(gw*rjk(j,k)*rjk(j,k))) )
      end if
      return
      end


c Pauli Potential
c
      function Pau(j,k)
      implicit none
      integer j, k
      real*8 Pau, pjk, rjk
      include 'coms.f'
      
      Pau = Pau0*exp(-(0.5*rjk(j,k)*rjk(j,k)/drPau))*
     *           exp(-(0.5*pjk(j,k)*pjk(j,k)/dpPau))
      return
      end

c Derivative (p) for Pauli Potential
c
      function dPaup(j,k)
      implicit none
      integer j, k
      real*8 dPaup, pjk, rjk
      include 'coms.f'
      
      dPaup = -(Pau0/dpPau*pjk(j,k)*
     *                    exp(-(0.5*rjk(j,k)*rjk(j,k)/drPau))*
     *                    exp(-(0.5*pjk(j,k)*pjk(j,k)/dpPau)))
      return
      end
  
c Derivative (r) for Pauli Potential
c
      function dPaur(j,k)
      implicit none
      integer j, k
      real*8 dPaur, pjk, rjk
      include 'coms.f'
      
      dPaur = -(Pau0/drPau*rjk(j,k)*
     *                   exp(-(0.5*rjk(j,k)*rjk(j,k)/drPau))*
     *                   exp(-(0.5*pjk(j,k)*pjk(j,k)/dpPau)))
      return
      end

      function rjk(j,k)
      implicit none 
      integer j, k
      real*8 rjk
      include 'coms.f'

      rjk = sqrt((rx(j)-rx(k))**2+(ry(j)-ry(k))**2+(rz(j)-rz(k))**2)           
      return
      end

      function pjk(j,k)
      implicit none 
      integer j, k
      real*8 pjk
      include 'coms.f'

      pjk = sqrt((px(j)-px(k))**2+(py(j)-py(k))**2+(pz(j)-pz(k))**2)           
      return
      end

      function iPau(j,k)
      implicit none 
      integer j, k
      logical iPau
      include 'coms.f'

      iPau = .false.
      if (iso3(j).eq.iso3(k).and.ityp(j).eq.ityp(k)) iPau = .true.
      return
      end







c $Id: proppot.f,v 1.9 1997/08/25 13:37:55 konopka Exp $
c here are all Jens Konopka routines dealing with propagation and forces

      function timing(itest)
      integer itest, timing 
c      integer tin, tout
c, mclock
c      save tin
      
      timing=0
c AIX
c      if(itest.ne.0) then
c        tout = mclock()
c        timing = tout-tin 
c        tin = tout
c      else
c        tin = mclock()
c      end if

c PC
c      if(itest.ne.0) then
c        call dclock@(tout)
c        timing = tout-tin 
c        tin = tout
c      else
c        call dclock@(tin)
c      end if

      return
      end   



 
