!=============================================================================================
program main
  use mod_bandstructure
  use mod_phonon_dispersion
  implicit none
  include 'mpif.h'
  character string*200
  character suffix*200
  character filename1*200
  character filename2*200
  integer, parameter :: root = 0
  integer, parameter :: CuLedge = 1
  integer, parameter :: OKedge = 2
  integer, parameter :: nthreads = 4
  integer, parameter :: use_band = 2
  integer, parameter :: use_Z = 1
  integer, parameter :: nfwrk1 = 60
  integer, parameter :: nfwrk2 = 61
  integer, parameter :: nfwrk3 = 62
  integer, parameter :: nfwrk4 = 63
  integer, parameter :: Nk = 64
  integer, parameter :: numomega = 100
  integer, parameter :: ibare = 1
  integer, parameter :: ibr = 2
  integer, parameter :: ib1g = 3
  integer, parameter :: ia1g = 4
  integer, parameter :: iapex = 5
  integer, parameter :: iac = 6
  integer myrank, mpistat, numproc
  integer edge
  integer orb
  integer i, j, k
  integer nomega
  integer nkx, nky, npx, npy, nqx, nqy
  double precision, dimension(0:nk,0:numomega,1:6) :: RIXS
  double precision, parameter :: omegascale = 0.1d0
  double precision, parameter :: Temp = 10.0d0
  double precision, parameter :: kb  = 8.617e-5
  double precision, parameter :: beta = 1.0d0/(kb*Temp)
  double precision, parameter :: eps9 = 1e-9
  double precision Qcdw 
  double precision deltacdw
  double precision sigmacdw
  double precision acdw
  double precision ampcdw 
  double precision etacdw, etacdw0 
  double precision gam
  double precision wb1g, wa1g, wbr, wapex, wac
  double precision filling
  double precision zfac, zfacp
  double precision arg
  double precision denom
  double precision fermi_fac
  double precision domega, omega, omegap
  double precision orb_fac, orb_facp
  double precision g0br, g0a1g, g0b1g, g0apex
  double precision gb1g, ga1g, gbr, gapex, gac
  double precision gb1gp, ga1gp, gbrp, gapexp, gacp
  double precision pi
  double precision starttime, endtime
  double precision E2p3d, win, wout
  double precision kx, ky, px, py, qx, qy
  double precision Mq
  double precision Nb
  double precision gcdwbr
  double precision gcdwapex
  double precision gcdwb1g
  double precision gcdwa1g
  double precision gcdwac
  double precision epsbx
  double precision epsby
  double precision epscx
  double precision epscy
  double precision, parameter :: gammae = 0.025d0
  double precision, dimension(-Nk:Nk) :: dk
  double precision, dimension(1:nbands) :: Phik
  double precision, dimension(1:nbands) :: Phip
  double precision, dimension(-Nk:2*Nk,-Nk:2*Nk) :: cosx_LUT, sinx_LUT, cosy_LUT, siny_LUT
  double precision, dimension(-Nk:2*Nk,-Nk:2*Nk) :: Ek_LUT
  double precision, dimension(-Nk:2*Nk,-Nk:2*Nk) :: Zk_LUT
  double precision, dimension(-Nk:2*Nk,-Nk:2*Nk) :: Fermik_LUT
  double precision, dimension(-Nk:2*Nk,-Nk:2*Nk) :: Deltak_LUT
  double precision, dimension(-Nk:Nk,-Nk:Nk) :: epsbx_lut
  double precision, dimension(-Nk:Nk,-Nk:Nk) :: epsby_lut
  double precision, dimension(-Nk:Nk,-Nk:Nk) :: epscx_lut
  double precision, dimension(-Nk:Nk,-Nk:Nk) :: epscy_lut
  double precision, dimension(-Nk:2*Nk,-Nk:2*Nk,1:nbands) :: Phi_LUT
  double complex, dimension(0:Nk,0:numomega) :: chi
  double complex ctmp
  double complex green, greenp, cfac, cfacp
  double complex, dimension(1:6) :: D
  double complex, dimension(0:Nk,0:numomega) :: Db1g, Dbr, Da1g, Dapex, Dac
  double complex, dimension(0:Nk,0:numomega) :: Db1g0, Dbr0, Da1g0, Dapex0, Dac0
  double complex, parameter :: eye = (0.0d0,1.0d0)
  double complex, parameter :: cone = (1.0d0,0.0d0)
 
  call mpi_init(mpistat)
  call mpi_comm_rank(mpi_comm_world,myrank,mpistat)
  call mpi_comm_size(mpi_comm_world,numproc,mpistat)
  print*, 'Process', myrank, ' of ', numproc, ' started.'
 
  edge = OKedge
  if(edge.eq.CuLedge)then
    orb = Cudx2y2
    gam = 0.3d0
    if(myrank.eq.root) print*, 'Computing spectra at the Cu L-edge.'
  elseif(edge.eq.OKedge)then
    orb = Ox
    gam = 0.15d0
    if(myrank.eq.root) print*, 'Computing spectra at the O  K-edge.'
  endif

  !Set the phonon dispersion
  wb1g = wb1g0
  wa1g = wa1g0
  wbr = wbr0
  wapex = wapex0

  !set the prefactor for the e-ph coupling
  g0br = 2.9*0.04433*sqrt(0.066/wbr)*(3.5*1.6/1.85);
  g0b1g = 2.9*3.56*0.04433*sqrt(0.066/wb1g);
  g0a1g = 2.9*3.56*0.04433*sqrt(0.066/wa1g);
  g0apex = 16.0d0*0.04433*sqrt(0.066/wapex);  

  pi = 2.0d0*asin(1.0d0)
  Qcdw = 0.5d0*pi
  acdw = 0.30d0 !0.250d0
  deltacdw = 0.012d0
  etacdw0 = 0.01d0
  gcdwbr = 0.008d0; gcdwapex = 0.005d0; gcdwa1g = 0.005d0; gcdwb1g = 0.005d0; gcdwac = 0.001d0

  call get_suffix(use_band,use_Z,suffix)

  domega = omegascale/dfloat(numomega)
  !Set up the integration weights appropriate for 2D simpson's rule
  do nkx = -Nk,Nk-1,2
    dk(nkx) = 2.0d0
    dk(nkx+1) = 4.0d0
  enddo
  dk(-Nk) = 1.0d0
  dk( Nk) = 1.0d0
  dk = dk/dfloat(3*2*Nk)    !Note: the factor of 2*pi is cancelled by the prefactor for bz integrals
  
  !First we construct look up tables for everything to save time
  filling = 0.0d0
  do nkx = -nk,2*Nk
    kx = pi*dfloat(nkx)/dfloat(Nk)
    do nky = -nk,2*Nk
      ky = pi*dfloat(nky)/dfloat(nk)
      EK_lut(nkx,nky) = energy(kx,ky,use_band)
      Zk_lut(nkx,nky) = quasiparticle_wgt(kx,ky,use_Z)
      call get_phi(kx,ky,band_to_compute,phik)
      Phi_lut(nkx,nky,:) = phik
      Fermik_lut(nkx,nky) = Fermi(Ek_lut(nkx,nky),beta)
      Deltak_LUT(nkx,nky) = Lorentzian(Ek_lut(nkx,nky),gammae)
  
      cosx_Lut(nkx,nky) = cos(kx/2.0d0)
      cosy_Lut(nkx,nky) = cos(ky/2.0d0)
      sinx_Lut(nkx,nky) = sin(kx/2.0d0)
      siny_Lut(nkx,nky) = sin(ky/2.0d0)
  
      if(nkx.le.Nk.and.nky.le.Nk)then
        filling = filling + 2.0d0*Fermi(Ek_lut(nkx,nky),beta)*dk(nkx)*dk(nky)

        if(kx.eq.0.0d0.and.ky.eq.0.0d0)then
          kx = 1e-6
          ky = 1e-6
        endif
        Mq = sqrt(cos(kx/2.0d0)*cos(kx/2.0d0) + cos(ky/2.0d0)*cos(ky/2.0d0));
        epscx_lut(nkx,nky) = cos(ky/2.0d0)/Mq;
        epscy_lut(nkx,nky) = cos(kx/2.0d0)/Mq;
        Mq = sqrt(sin(kx/2.0d0)*sin(kx/2.0d0) + sin(ky/2.0d0)*sin(ky/2.0d0));
        epsbx_lut(nkx,nky) = sin(kx/2.0d0)/Mq;
        epsby_lut(nkx,nky) = sin(ky/2.0d0)/Mq;
      endif
    enddo
  enddo
  if(myrank.eq.root) print*, 'Total electronic filling = ', filling

  ! Working here, insert code for the chi from the paper.
  print*, 'Building chi.'
  chi = dcmplx(0.0d0,0.0d0)
  do nqx = 0,Nk
    nqy = 0
    qx = dfloat(nqx)*pi/dfloat(Nk)
    qy = dfloat(nqy)*pi/dfloat(Nk)
    do nomega =  0,numomega
      omega = dfloat(nomega)*domega 
      chi(nqx,nomega) = 0.0d0
      etacdw = 0.01d0+omega
      do nqy = 0,0
        chi(nqx,nomega) = chi(nqx,nomega) + cone/(deltacdw*deltacdw + acdw*acdw*(qx-Qcdw)**2.0d0 - (omega+eye*etacdw)**2.0d0)
      enddo
    enddo 
  enddo

  do nqx = 0,Nk
    nqy = 0
    qx = dfloat(nqx)*pi/dfloat(Nk)
    qy = dfloat(nqy)*pi/dfloat(Nk)
    wb1g = get_wb1g(qx,qy)
    wa1g = get_wa1g(qx,qy)
    wbr = get_wbr(qx,qy)
    wapex = get_wapex(qx,qy)    
    wac = get_wac(qx,qy) 
   
    do nomega =  0,numomega
      omega = dfloat(nomega)*domega   
 
      Dbr0(nqx,nomega) = cone/(omega - wbr + eye*0.005d0) - cone/(omega + wbr + eye*0.005d0)
      Dapex0(nqx,nomega) = cone/(omega - wapex + eye*0.005d0) - cone/(omega + wapex + eye*0.005d0)
      Db1g0(nqx,nomega) = cone/(omega - wb1g + eye*0.005d0) - cone/(omega + wb1g + eye*0.005d0)
      Da1g0(nqx,nomega) = cone/(omega - wa1g + eye*0.005d0) - cone/(omega + wa1g + eye*0.005d0)
      Dac0(nqx,nomega) = cone/(omega - wac + eye*0.005d0) - cone/(omega + wac + eye*0.005d0)

      Dbr(nqx,nomega) = Dbr0(nqx,nomega)/(cone + gcdwbr*gcdwbr*Dbr0(nqx,nomega)*chi(nqx,nomega))
      Dapex(nqx,nomega) = Dapex0(nqx,nomega)/(cone + gcdwapex*gcdwapex*Dapex0(nqx,nomega)*chi(nqx,nomega))
      Da1g(nqx,nomega) = Da1g0(nqx,nomega)/(cone + gcdwa1g*gcdwa1g*Da1g0(nqx,nomega)*chi(nqx,nomega))
      Db1g(nqx,nomega) = Db1g0(nqx,nomega)/(cone + gcdwb1g*gcdwb1g*Db1g0(nqx,nomega)*chi(nqx,nomega))
      Dac(nqx,nomega) = Dac0(nqx,nomega)/(cone + gcdwac*gcdwac*Dac(nqx,nomega)*chi(nqx,nomega))
    enddo
  enddo
  
  if(myrank.eq.root)then
    print*, 'Outputting chi(q,w), and D(q,w).'
    open(file='chi.dat',unit=nfwrk1,action='write')
    open(file='ImDChi.dat',unit=nfwrk2,action='write')
    do nqx = 0,Nk
      nqy = 0
      qx = dfloat(nqx)*pi/dfloat(Nk)
      do nomega =  0,numomega
        omega = dfloat(nomega)*domega 
        if(nomega.ne.0)then 
          nb = 2.0d0/(1.0d0-exp(-omega*beta))
        else 
          nb = 2.0d0
        endif
        write(unit=nfwrk1,fmt=*) qx/pi, omega, nb*real(chi(nqx,nomega)), nb*imag(chi(nqx,nomega))
        write(unit=nfwrk2,fmt=*) qx/pi, omega,-imag(Dapex(nqx,nomega)), -imag(Dapex0(nqx,nomega)), &
                                              -imag(Da1g(nqx,nomega)), -imag(Da1g0(nqx,nomega)),  &
                                              -imag(Db1g(nqx,nomega)), -imag(Db1g0(nqx,nomega)),  &
                                              -imag(Dbr(nqx,nomega)), -imag(Dbr0(nqx,nomega)), &  
                                              -imag(Dac(nqx,nomega)), -imag(Dac0(nqx,nomega))  
      enddo
      write(unit=nfwrk1,fmt=*) ' '
      write(unit=nfwrk2,fmt=*) ' '
    enddo
    close(unit=nfwrk1) 
    close(unit=nfwrk2)
  endif 


  goto 999
  E2p3d = -932.0d0*0.0d0
  win = -E2p3d
  RIXS = 0.0d0 

  if(myrank.eq.root) print*, 'Computing RIXS spectra.'
  call cpu_time(starttime)
  do nqx = myrank,Nk,numproc
    print*, myrank, nqx, nk
    nqy = 0
    qx = dfloat(nqx)*pi/dfloat(Nk)
    qy = dfloat(nqy)*pi/dfloat(Nk)
          
    wb1g = get_wb1g(qx,qy)
    wa1g = get_wa1g(qx,qy)
    wbr = get_wbr(qx,qy)
    wapex = get_wapex(qx,qy)

  
    do nomega = 0,numomega
      omega = dfloat(nomega)*domega
      wout = win - omega
  
      D(iapex) = Dapex(nqx,nomega)
      D(ia1g) = Da1g(nqx,nomega)
      D(ib1g) = Db1g(nqx,nomega)
      D(ibr) = Dbr(nqx,nomega)
      D(iac) = Dac(nqx,nomega)
 
      do npx = -Nk,Nk
        do npy = -Nk,Nk
          zfac = Zk_lut(npx,npy) 
          zfacp = Zk_lut(npx+nqx,npy+nqy)
          denom = (ek_lut(npx,npy) - wout - E2p3d)**2.0d0 + gam*gam
          arg = omega + ek_lut(npx,npy) - ek_lut(npx+nqx,npy+nqy)
          fermi_fac = fermik_lut(npx,npy) - fermik_lut(npx+nqx,npy+nqy)  
          RIXS(nqx,nomega,ibare) = RIXS(Nqx,nomega,ibare) + zfac*zfacp &  
                                  *phi_lut(npx,npy,orb)*phi_lut(npx,npy,orb) &
                                  *phi_lut(npx+nqx,npy+nqy,orb)*phi_lut(npx+nqx,npy+nqy,orb) &
                                  *fermi_fac*Lorentzian(arg,gammae)*dk(npx)*dk(npy)/denom
        enddo
      enddo
  
      epscx = epscx_lut(nqx,nqy)
      epscy = epscy_lut(nqx,nqy)
      epsbx = epsbx_lut(nqx,nqy)
      epsby = epsby_lut(nqx,nqy)
 
      do npx = -nk,nk
        px = dfloat(npx)*pi/dfloat(nk)
        do npy = -nk,nk
          py = dfloat(npy)*pi/dfloat(nk)

          zfac = Zk_lut(npx,npy)*Zk_lut(npx+nqx,npy+nqy)

          orb_fac = phi_lut(npx,npy,orb)*phi_lut(npx,npy,orb) &  
                   *phi_lut(npx+nqx,npy+nqy,orb)*phi_lut(npx+nqx,npy+nqy,orb)
          cfac = (1.0d0-fermik_lut(npx,npy))/(ek_lut(npx,npy)-wout-E2p3d-eye*gam) & 
               - (1.0d0-fermik_lut(npx+nqx,npy+nqy))/(ek_lut(npx+nqx,npy+nqy)-win-E2p3d-eye*gam)
          !Note that I have -G^* where as Tom writes down G
          !green = cone/(ek_lut(npx,npy)-ek_lut(npx+nqx,npy+nqy)+omega+eye*gammae)
          green =-cone/(ek_lut(npx,npy)-ek_lut(npx+nqx,npy+nqy)+omega-eye*gammae)

          phik = phi_lut(npx,npy,:)
          phip = phi_lut(npx+nqx,npy+nqy,:)
!          gapex = g0apex*0.25d0*(cosx_lut(npx,npy)-cosy_lut(npx,npy))*  &  
!                                (cosx_lut(npx+nqx,npy+nqy)-cosy_lut(npx+nqy,npy+nqy))
!          Gb1g = g0b1g*(epscx*Phik(Ox)*Phip(Ox)-epscy*Phik(Oy)*Phip(Oy));
!          Ga1g = g0a1g*(epscx*Phik(Ox)*Phip(Ox)+epscy*Phik(Oy)*Phip(Oy));
!          Gbr = g0br*epsbx*(cosx_lut(npx+nqx,npy+nqy)*Phip(Cudx2y2)*Phik(Ox) & 
!                           -cosx_lut(npx,npy)*Phik(Cudx2y2)*Phip(Ox)) &
!              - g0br*epsby*(cosy_lut(npx+nqx,npy+nqy)*Phip(Cudx2y2)*Phik(Oy) & 
!                           -cosy_lut(npx,npy)*Phik(Cudx2y2)*Phip(Oy));

          gbr = 0.5d0*sqrt(sinx_lut(nqx,nqy)*sinx_lut(nqx,nqy) + siny_lut(nqx,nqy)*siny_lut(nqx,nqy))
          gb1g = 0.5d0*(sinx_lut(npx,npy)*sinx_lut(npx+nqx,npy+nqy)*cosy_lut(nqx,nqy) & 
               - siny_lut(npx,npy)*siny_lut(npx+nqx,npy+nqy)*cosx_lut(nqx,nqy)) 
          ga1g = 0.50d0*(sinx_lut(npx,npy)*sinx_lut(npx+nqx,npy+nqy)*cosy_lut(nqx,nqy) & 
               + siny_lut(npx,npy)*siny_lut(npx+nqx,npy+nqy)*cosx_lut(nqx,nqy))  
          gapex = 0.25d0*(cosx_lut(npx,npy)-cosy_lut(npx,npy))*  &  
                         (cosx_lut(npx+nqx,npy+nqy)-cosy_lut(npx+nqy,npy+nqy))
          gac = 1.0d0  
 
          do nkx = -nk,nk
            kx = dfloat(nkx)*pi/dfloat(nk)
            do nky = -nk,nk
              ky = dfloat(nky)*pi/dfloat(nk)

              zfacp = Zk_lut(nkx,nky)*Zk_lut(nkx+nqx,nky+nqy)
              orb_facp = phi_lut(nkx,nky,orb)*phi_lut(nkx,nky,orb) &  
                        *phi_lut(nkx+nqx,nky+nqy,orb)*phi_lut(nkx+nqx,nky+nqy,orb)
              cfacp = (1.0d0-fermik_lut(nkx,nky))/(ek_lut(nkx,nky)-wout-E2p3d+eye*gam) & 
                    - (1.0d0-fermik_lut(nkx+nqx,nky+nqy))/(ek_lut(nkx+nqx,nky+nqy)-win-E2p3d+eye*gam)
              greenp = cone/(ek_lut(nkx,nky)-ek_lut(nkx+nqx,nky+nqy)+omega+eye*gammae)
         
              phik = phi_lut(nkx,nky,:)
              phip = phi_lut(nkx+nqx,nky+nqy,:)
              !gapexp = g0apex*0.25d0*(cosx_lut(nkx,nky)-cosy_lut(nkx,nky))*  &  
              !                       (cosx_lut(nkx+nqx,nky+nqy)-cosy_lut(nkx+nqy,nky+nqy))
              !Gb1gp = g0b1g*(epscx*Phik(Ox)*Phip(Ox)-epscy*Phik(Oy)*Phip(Oy));
              !Ga1gp = g0a1g*(epscx*Phik(Ox)*Phip(Ox)+epscy*Phik(Oy)*Phip(Oy));
              !Gbrp = g0br*epsbx*(cosx_lut(nkx+nqx,nky+nqy)*Phip(Cudx2y2)*Phik(Ox) & 
              !                  -cosx_lut(nkx,nky)*Phik(Cudx2y2)*Phip(Ox)) &
              !     - g0br*epsby*(cosy_lut(nkx+nqx,nky+nqy)*Phip(Cudx2y2)*Phik(Oy) & 
              !                  -cosy_lut(nkx,nky)*Phik(Cudx2y2)*Phip(Oy));

              gbrp = 0.50d0*sqrt(sinx_lut(nqx,nqy)*sinx_lut(nqx,nqy) + siny_lut(nqx,nqy)*siny_lut(nqx,nqy))
              gb1gp = 0.50d0*(sinx_lut(nkx,nky)*sinx_lut(nkx+nqx,nky+nqy)*cosy_lut(nqx,nqy) & 
                    - siny_lut(nkx,nky)*siny_lut(nkx+nqx,nky+nqy)*cosx_lut(nqx,nqy))  
              ga1gp = 0.5d0*(sinx_lut(nkx,nky)*sinx_lut(nkx+nqx,nky+nqy)*cosy_lut(nqx,nqy) & 
                    + siny_lut(nkx,nky)*siny_lut(nkx+nqx,nky+nqy)*cosx_lut(nqx,nqy))  
              gapexp = 0.25d0*(cosx_lut(nkx,nky)-cosy_lut(nkx,nky))*  &  
                         (cosx_lut(nkx+nqx,nky+nqy)-cosy_lut(nkx+nqy,nky+nqy))
              gacp = 1.0d0                   
                          
              ctmp = zfac*orb_fac*cfac*gbr*imag(Green*D(ibr)*greenp)*zfacp*cfacp*orb_facp*gbrp/pi 
              RIXS(nqx,nomega,ibr) = RIXS(nqx,nomega,ibr) + real(ctmp)*dk(npx)*dk(npy)*dk(nkx)*dk(nky) 
                                             
              ctmp = zfac*orb_fac*cfac*gb1g*imag(Green*D(ib1g)*greenp)*zfacp*cfacp*orb_facp*gb1gp/pi 
              RIXS(nqx,nomega,ib1g) = RIXS(nqx,nomega,ib1g) + real(ctmp)*dk(npx)*dk(npy)*dk(nkx)*dk(nky) 
                                             
              ctmp = zfac*orb_fac*cfac*ga1g*imag(Green*D(ia1g)*greenp)*zfacp*cfacp*orb_facp*ga1gp/pi 
              RIXS(nqx,nomega,ia1g) = RIXS(nqx,nomega,ia1g) + real(ctmp)*dk(npx)*dk(npy)*dk(nkx)*dk(nky) 
                                             
              ctmp = zfac*orb_fac*cfac*gapex*imag(Green*D(iapex)*greenp)*zfacp*cfacp*orb_facp*gapexp/pi 
              RIXS(nqx,nomega,iapex) = RIXS(nqx,nomega,iapex) + real(ctmp)*dk(npx)*dk(npy)*dk(nkx)*dk(nky) 
              
              ctmp = zfac*orb_fac*cfac*gac*imag(Green*D(iac)*greenp)*zfacp*cfacp*orb_facp*gacp/pi 
              RIXS(nqx,nomega,iac) = RIXS(nqx,nomega,iac) + real(ctmp)*dk(npx)*dk(npy)*dk(nkx)*dk(nky) 
            enddo
          enddo 
  
        enddo
      enddo
    enddo
  enddo
  call MPI_ALLREDUCE(MPI_IN_PLACE,RIXS,size(RIXS),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpistat)
  call cpu_time(endtime)
  if(myrank.eq.root) print*, 'Time (s) = ', endtime-starttime
  
  ! Output the RIXS results

  
  if(myrank.eq.root)then
    if(edge.eq.CuLedge)then
      filename1 = 'rixs_CuLedge' // trim(suffix)
    elseif(edge.eq.OKedge)then
      filename1 = 'rixs_OKedge' // trim(suffix)
    endif

    open(file=trim(filename1),unit=nfwrk1,action='write')
    do nqx = 0,Nk
      qx = dfloat(nqx)*pi/dfloat(Nk)
      do nomega = 0,numomega
        omega = dfloat(nomega)*domega
        write(unit=nfwrk1,fmt=*) qx/pi, omega, (RIXS(nqx,nomega,i),i=1,6)
      enddo
      write(unit=nfwrk1,fmt=*) ' '
    enddo
    close(unit=nfwrk1)
  endif 

999 print*, 'Skipped.'
  call mpi_barrier(mpi_comm_world,mpistat)
  call mpi_finalize(mpistat)

  stop
!=============================================================================================
! End of the main program.
!=============================================================================================
contains
  double precision function fermi(E,beta)
  double precision, intent(in) :: E, Beta
  fermi = 1.0d0/(exp(E*beta)+1.0d0)
  end function fermi

  double precision function Lorentzian(E,gam)
  double precision, intent(in) :: E, gam
  double precision pi
  pi = 2.0d0*asin(1.0d0)
  Lorentzian = (gam/pi)/(E*E + gam*gam)
  end function Lorentzian

  subroutine test_imaginary_part(G,string)
  double complex, intent(in) :: G
  character string*200
  if(imag(G).gt.0.0d0)then
    print*, 'Error: ', string, ' has a positive real part.'
  endif
  end subroutine test_imaginary_part 

end program main
