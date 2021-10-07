module mod_tight_binding
  implicit none
  integer, public, parameter          :: band_to_compute = 3
  integer, public, parameter          :: nbands = 3
  integer, public, parameter          :: Cudx2y2 = 1
  integer, public, parameter          :: Ox = 2
  integer, public, parameter          :: Oy = 3
  double precision, public, parameter :: mu = 1.0325d0
  double precision, public, parameter :: Z = 3.0d0
  double precision, public, parameter :: ed = 0.0d0
  double precision, public, parameter :: ep =-3.6d0
  double precision, public, parameter :: tpd = 1.3d0
  double precision, public, parameter :: tpp =-0.65d0
contains
  subroutine diagonalize_H(kx,ky,E,H)
    integer info
    integer, parameter :: lwork = 3*nbands-1
    double precision, dimension(1:lwork) :: work
    double precision, intent(in) :: kx, ky
    double precision, dimension(1:nbands), intent(out) :: E
    double precision, dimension(1:nbands,1:nbands), intent(out) :: H
    call get_H(kx,ky,H)
    call dsyev('V','U',nbands,H,nbands,E,work,lwork,info)
    E = E/Z
    if(H(1,band_to_compute).lt.0.0d0)then
      H(:,band_to_compute) = -H(:,band_to_compute)
    endif
  end subroutine diagonalize_H

  subroutine get_H(kx,ky,H)
    double precision, intent(in) ::  kx, ky
    double precision, dimension(1:nbands,1:nbands), intent(out) :: H
    H = 0.0d0
    H(Cudx2y2,Cudx2y2) = ed-mu;
    H(Ox,Ox) = ep - mu + 4.0d0*tpp*sin(kx/2.0d0)*sin(kx/2.0d0);
    H(Oy,Oy) = ep - mu + 4.0d0*tpp*sin(ky/2.0d0)*sin(ky/2.0d0);
    H(Cudx2y2,Ox) = 2.0d0*tpd*sin(kx/2.0d0);
    H(Ox,Cudx2y2) = 2.0d0*tpd*sin(kx/2.0d0);
    H(Cudx2y2,Oy) =-2.0d0*tpd*sin(ky/2.0d0);
    H(Oy,Cudx2y2) =-2.0d0*tpd*sin(ky/2.0d0);
    H(Ox,Oy) = 4.0d0*tpp*sin(kx/2.0d0)*sin(ky/2.0d0);
    H(Oy,Ox) = 4.0d0*tpp*sin(kx/2.0d0)*sin(ky/2.0d0);
  end subroutine get_H
end module mod_tight_binding
!=============================================================================================
program main
  use mod_tight_binding
  implicit none
  character string*200
  character filename*200
  integer, parameter :: CuLedge = 1
  integer, parameter :: OKedge = 2
  integer, parameter :: nthreads = 4
  integer, parameter :: nfwrk1 = 60
  integer, parameter :: nfwrk2 = 61
  integer, parameter :: nfwrk3 = 62
  integer, parameter :: Nk = 24
  integer, parameter :: numomega = 100
  integer, parameter :: ibare = 1
  integer, parameter :: ibr = 2
  integer, parameter :: ib1g = 3
  integer, parameter :: ia1g = 4
  integer, parameter :: iapex = 5
  integer edge
  integer orb
  integer i, j, k
  integer nomega
  integer nkx, nky, npx, npy, nqx, nqy
  double precision, dimension(0:nk,0:numomega,1:5) :: RIXS
  double precision, parameter :: omegascale = 0.2d0
  double precision, parameter :: Temp = 10.0d0
  double precision, parameter :: kb  = 8.617e-5
  double precision, parameter :: beta = 1.0d0/(kb*Temp)
  double precision, parameter :: gam = 0.5d0 !core hole lifetime
  double precision wb1g, wa1g, wbr, wapex
  double precision filling
  double precision arg
  double precision denom
  double precision fermi_fac
  double precision domega, omega
  double precision orb_fac, orb_facp
  double precision gb1g, ga1g, gbr, gapex
  double precision gb1gp, ga1gp, gbrp, gapexp
  double precision pi
  double precision starttime, endtime
  double precision E2p3d, win, wout
  double precision kx, ky, px, py, qx, qy
  double precision, parameter :: gammae = 0.025d0
  double precision, dimension(-Nk:Nk) :: dk
  double precision, dimension(1:nbands) :: E
  double precision, dimension(1:nbands,1:nbands) :: H
  double precision, dimension(-Nk:2*Nk,-Nk:2*Nk) :: cosx_LUT, sinx_LUT, cosy_LUT, siny_LUT
  double precision, dimension(-Nk:2*Nk,-Nk:2*Nk) :: Ek_LUT
  double precision, dimension(-Nk:2*Nk,-Nk:2*Nk) :: Fermik_LUT
  double precision, dimension(-Nk:2*Nk,-Nk:2*Nk) :: Deltak_LUT
  double precision, dimension(-Nk:2*Nk,-Nk:2*Nk,1:nbands) :: Phi_LUT
  double complex ctmp
  double complex green, greenp, cfac, cfacp
  double complex Db1g, Dbr, Da1g, Dapex
  double complex, parameter :: eye = (0.0d0,1.0d0)
  double complex, parameter :: cone = (1.0d0,0.0d0)
  
  edge = OKedge
  if(edge.eq.CuLedge)then
    orb = Cudx2y2
    print*, 'Computing spectra at the Cu L-edge.'
  elseif(edge.eq.OKedge)then
    orb = Ox
    print*, 'Computing spectra at the O  K-edge.'
  endif

  pi = 2.0d0*asin(1.0d0)
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
  open(file='bands.dat',unit=nfwrk1,action='write')
  open(file='phi_functions.dat',unit=nfwrk2,action='write')
  do nkx = -nk,2*Nk
    kx = pi*dfloat(nkx)/dfloat(Nk)
    do nky = -nk,2*Nk
      ky = pi*dfloat(nky)/dfloat(nk)
      call diagonalize_H(kx,ky,E,H)
      Ek_lut(nkx,nky) = E(band_to_compute)
      Phi_lut(nkx,nky,:) = H(:,band_to_compute)
      Fermik_lut(nkx,nky) = Fermi(Ek_lut(nkx,nky),beta)
      Deltak_LUT(nkx,nky) = Lorentzian(Ek_lut(nkx,nky),gammae)
  
      cosx_Lut(nkx,nky) = cos(kx/2.0d0)
      cosy_Lut(nkx,nky) = cos(ky/2.0d0)
      sinx_Lut(nkx,nky) = sin(kx/2.0d0)
      siny_Lut(nkx,nky) = sin(ky/2.0d0)
  
      if(nkx.le.Nk.and.nky.le.Nk)then
        filling = filling + 2.0d0*Fermi(Ek_lut(nkx,nky),beta)*dk(nkx)*dk(nky)
        write(unit=nfwrk1,fmt=*) kx/pi,ky/pi,Ek_lut(nkx,nky), Fermik_LUT(nkx,nky), deltak_LUT(nkx,nky)
        write(unit=nfwrk2,fmt=*) kx/pi,ky/pi,(Phi_lut(nkx,nky,i),i=1,3)
      endif
    enddo
    write(unit=nfwrk1,fmt=*) ' '
    write(unit=nfwrk2,fmt=*) ' '
  enddo
  close(unit=nfwrk1)
  close(unit=nfwrk2)
  print*, 'Total electronic filling = ', filling
  
  
  E2p3d = -932.0d0*0.0d0
  win = -E2p3d
  call cpu_time(starttime)
  open(file='phonon_dos.dat',unit=nfwrk3,action='write')
  do nqx = 0,Nk
    nqy = 0
    qx = dfloat(nqx)*pi/dfloat(Nk)
    qy = dfloat(nqy)*pi/dfloat(Nk)
  
    wbr = 0.085*(1.0d0-0.18d0*(sin(qx/2.0d0)*sin(qx/2.0d0)+sin(qy/2.0d0)*sin(qy/2.0d0)))
    wapex = 0.085d0
    wb1g = 0.035d0
    wa1g = 0.040d0 
  
    do nomega = 0,numomega
      omega = dfloat(nomega)*domega
      wout = win - omega
      RIXS(nqx,nomega,:) = 0.0d0
  
      Dbr = cone/(omega - wbr + eye*0.005d0) - cone/(omega + wbr + eye*0.005d0)
      Dapex = cone/(omega - wapex + eye*0.005d0) - cone/(omega + wapex + eye*0.005d0)
      Db1g = cone/(omega - wb1g + eye*0.005d0) - cone/(omega + wb1g + eye*0.005d0)
      Da1g = cone/(omega - wa1g + eye*0.005d0) - cone/(omega + wa1g + eye*0.005d0)
  
      write(unit=nfwrk3,fmt=*) qx/pi, omega, -imag(Dbr)/pi, -imag(Dapex)/pi, -imag(Db1g)/pi, -imag(Da1g)/pi
  
      do npx = -Nk,Nk
        do npy = -Nk,Nk
          denom = (ek_lut(npx,npy) - wout - E2p3d)**2.0d0 + gam*gam
          arg = omega + ek_lut(npx,npy) - ek_lut(npx+nqx,npy+nqy)
          fermi_fac = fermik_lut(npx,npy) - fermik_lut(npx+nqx,npy+nqy)  
          RIXS(nqx,nomega,ibare) = RIXS(Nqx,nomega,ibare) & 
                                 + phi_lut(npx,npy,orb)*phi_lut(npx,npy,orb) &
                                  *phi_lut(npx+nqx,npy+nqy,orb)*phi_lut(npx+nqx,npy+nqy,orb) &
                                  *fermi_fac*Lorentzian(arg,gammae)*dk(npx)*dk(npy)/denom
        enddo
      enddo
  
      RIXS(nqx,nomega,ibr) = RIXS(nqx,nomega,ibare)
      RIXS(nqx,nomega,iapex) = RIXS(nqx,nomega,ibare)
      RIXS(nqx,nomega,ib1g) = RIXS(nqx,nomega,ibare)
      RIXS(nqx,nomega,ia1g) = RIXS(nqx,nomega,ibare)
  
      do npx = -nk,nk
        px = dfloat(npx)*pi/dfloat(nk)
        do npy = -nk,nk
          py = dfloat(npy)*pi/dfloat(nk)
          !orb_fac = 1.0d0
          orb_fac = phi_lut(npx,npy,orb)*phi_lut(npx,npy,orb) &  
                   *phi_lut(npx+nqx,npy+nqy,orb)*phi_lut(npx+nqx,npy+nqy,orb)
          cfac = (1.0d0-fermik_lut(npx,npy))/(ek_lut(npx,npy)-wout-E2p3d-eye*gam) & 
               - (1.0d0-fermik_lut(npx+nqx,npy+nqy))/(ek_lut(npx+nqx,npy+nqy)-win-E2p3d-eye*gam)
          !Note that I have -G^* where as Tom writes down G
          !green = cone/(ek_lut(npx,npy)-ek_lut(npx+nqx,npy+nqy)+omega+eye*gammae)
          green =-cone/(ek_lut(npx,npy)-ek_lut(npx+nqx,npy+nqy)+omega-eye*gammae)
          !gbr=1.0d0
          !gb1g=1.0d0
          !ga1g=1.0d0
          !gapex=1.0d0
          gbr = 0.5d0*sqrt(sinx_lut(nqx,nqy)*sinx_lut(nqx,nqy) + siny_lut(nqx,nqy)*siny_lut(nqx,nqy))
          gb1g = 0.5d0*(sinx_lut(npx,npy)*sinx_lut(npx+nqx,npy+nqy)*cosy_lut(nqx,nqy) & 
               - siny_lut(npx,npy)*siny_lut(npx+nqx,npy+nqy)*cosx_lut(nqx,nqy)) 
          ga1g = 0.50d0*(sinx_lut(npx,npy)*sinx_lut(npx+nqx,npy+nqy)*cosy_lut(nqx,nqy) & 
               + siny_lut(npx,npy)*siny_lut(npx+nqx,npy+nqy)*cosx_lut(nqx,nqy))  
          gapex = 0.25d0*(cosx_lut(npx,npy)-cosy_lut(npx,npy))*  &  
                         (cosx_lut(npx+nqx,npy+nqy)-cosy_lut(npx+nqy,npy+nqy))
  
          do nkx = -nk,nk
            kx = dfloat(nkx)*pi/dfloat(nk)
            do nky = -nk,nk
              ky = dfloat(nky)*pi/dfloat(nk)
          !    orb_facp = 1.0d0       
              orb_facp = phi_lut(nkx,nky,orb)*phi_lut(nkx,nky,orb) &  
                        *phi_lut(nkx+nqx,nky+nqy,orb)*phi_lut(nkx+nqx,nky+nqy,orb)
              cfacp = (1.0d0-fermik_lut(nkx,nky))/(ek_lut(nkx,nky)-wout-E2p3d+eye*gam) & 
                    - (1.0d0-fermik_lut(nkx+nqx,nky+nqy))/(ek_lut(nkx+nqx,nky+nqy)-win-E2p3d+eye*gam)
              greenp = cone/(ek_lut(nkx,nky)-ek_lut(nkx+nqx,nky+nqy)+omega+eye*gammae)
              gbrp = 0.50d0*sqrt(sinx_lut(nqx,nqy)*sinx_lut(nqx,nqy) + siny_lut(nqx,nqy)*siny_lut(nqx,nqy))
              gb1gp = 0.50d0*(sinx_lut(nkx,nky)*sinx_lut(nkx+nqx,nky+nqy)*cosy_lut(nqx,nqy) & 
                    - siny_lut(nkx,nky)*siny_lut(nkx+nqx,nky+nqy)*cosx_lut(nqx,nqy))  
              ga1gp = 0.5d0*(sinx_lut(nkx,nky)*sinx_lut(nkx+nqx,nky+nqy)*cosy_lut(nqx,nqy) & 
                    + siny_lut(nkx,nky)*siny_lut(nkx+nqx,nky+nqy)*cosx_lut(nqx,nqy))  
              gapexp = 0.25d0*(cosx_lut(nkx,nky)-cosy_lut(nkx,nky))*  &  
                         (cosx_lut(nkx+nqx,nky+nqy)-cosy_lut(nkx+nqy,nky+nqy))
          !gbrp=1.0d0
          !gb1gp=1.0d0
          !ga1gp=1.0d0
          !gapexp=1.0d0
                                             
              !ctmp = orb_fac*cfac*gbr*imag(Green*Dbr*greenp)*cfacp*orb_facp*gbrp/pi 
              ctmp = orb_fac*cfac*gbr*imag(Green*Dbr*greenp)*cfacp*orb_facp*gbrp/pi 
              RIXS(nqx,nomega,ibr) = RIXS(nqx,nomega,ibr) + real(ctmp)*dk(npx)*dk(npy)*dk(nkx)*dk(nky) 
                                             
              ctmp = orb_fac*cfac*gb1g*imag(Green*Db1g*greenp)*cfacp*orb_facp*gb1gp/pi 
              RIXS(nqx,nomega,ib1g) = RIXS(nqx,nomega,ib1g) + real(ctmp)*dk(npx)*dk(npy)*dk(nkx)*dk(nky) 
                                             
              ctmp = orb_fac*cfac*ga1g*imag(Green*Da1g*greenp)*cfacp*orb_facp*ga1gp/pi 
              RIXS(nqx,nomega,ia1g) = RIXS(nqx,nomega,ia1g) + real(ctmp)*dk(npx)*dk(npy)*dk(nkx)*dk(nky) 
                                             
              ctmp = orb_fac*cfac*gapex*imag(Green*Dapex*greenp)*cfacp*orb_facp*gapexp/pi 
              RIXS(nqx,nomega,iapex) = RIXS(nqx,nomega,iapex) + real(ctmp)*dk(npx)*dk(npy)*dk(nkx)*dk(nky) 
            enddo
          enddo 
  
        enddo
      enddo
    enddo
    write(unit=nfwrk3,fmt=*) ' '
  enddo
  close(unit=nfwrk3)
  call cpu_time(endtime)
  print*, 'Time (s) = ', endtime-starttime
  
  ! Output the RIXS results
  if(edge.eq.CuLedge)then
    filename = 'rixs_CuLedge.dat'
  elseif(edge.eq.OKedge)then
    filename = 'rixs_OKedge.dat'
  endif
  

  open(file=trim(filename),unit=nfwrk1,action='write')
  do nqx = 0,Nk
    qx = dfloat(nqx)*pi/dfloat(Nk)
    do nomega = 0,numomega
      omega = dfloat(nomega)*domega
      write(unit=nfwrk1,fmt=*) qx/pi, omega, (RIXS(nqx,nomega,i),i=1,5)
    enddo
    write(unit=nfwrk1,fmt=*) ' '
  enddo
  close(unit=nfwrk1)
  
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

  double precision function vert_br(kx,ky,qx,qy) 
  double precision, intent(in) :: kx, ky, qx, qy
  vert_br = sqrt(sin(qx/2.0d0)*sin(qx/2.0d0) + sin(qy/2.0d0)*sin(qy/2.0d0))
  end function vert_br

  double precision function vert_apex(kx,ky,qx,qy)
  double precision, intent(in) :: kx, ky, qx, qy
  double precision px, py
  px = kx + qx
  py = ky + qy
  vert_apex = 0.25d0*(cos(kx/2.0d0)-cos(ky/2.0d0))*(cos(px/2.0d0)-cos(py/2.0d0))
  end function vert_apex
 
  double precision function vert_a1g(kx,ky,qx,qy)
  double precision, intent(in) :: kx, ky, qx, qy
  double precision px, py
  px = kx + qx
  py = ky + qy
  vert_a1g = sin(kx/2.0d0)*sin(px/2.0d0)*cos(qy/2.0d0) & 
           + sin(ky/2.0d0)*sin(py/2.0d0)*cos(qx/2.0d0)
  end function vert_a1g

  double precision function vert_b1g(kx,ky,qx,qy)
  double precision, intent(in) :: kx, ky, qx, qy
  double precision px, py
  px = kx + qx
  py = ky + qy
  vert_b1g = sin(kx/2.0d0)*sin(px/2.0d0)*cos(qy/2.0d0) & 
           - sin(ky/2.0d0)*sin(py/2.0d0)*cos(qx/2.0d0)
  end function vert_b1g

  subroutine test_imaginary_part(G,string)
  double complex, intent(in) :: G
  character string*200
  if(imag(G).gt.0.0d0)then
    print*, 'Error: ', string, ' has a positive real part.'
  endif
  end subroutine test_imaginary_part 

end program main
