!=============================================================================================
module mod_bandstructure
  implicit none
  integer, parameter, public :: nbands = 3
  integer, parameter, public :: band_to_compute = 3
  integer, parameter, public :: Cudx2y2 = 1
  integer, parameter, public :: Ox = 2
  integer, parameter, public :: Oy = 3
  integer, parameter, private :: Zeq1 = 1
  integer, parameter, private :: Znorman = 2
  integer, parameter, private :: BandTB = 1
  integer, parameter, private :: BandARPES = 2
  integer, parameter, private :: BandTBtimeZ = 3
  integer, parameter, private :: Zeff = 1.0d0/3.0d0
  double precision, private :: kf
  !Parameters for the tight-binding model
  double precision, parameter, private :: mu = 2.16d0
  double precision, parameter, private :: ed = 0.0d0
  double precision, parameter, private :: ep =-0.9d0
  double precision, parameter, private :: tpd = 1.6d0
  double precision, parameter, private :: tppi =-1.0d0
  double precision, parameter, private :: tppd = 0d0
  double precision, parameter, private :: tpp = tppi + tppd
contains
  subroutine get_suffix(use_band,use_z,suffix)
    integer, intent(in) :: use_band, use_z
    character, intent(out) :: suffix*200
    character bandstring*50, Zstring*50
    if(use_band.eq.BandTB)then
      Bandstring = '_TBbands'
    elseif(use_band.eq.BandARPES)then
      Bandstring = '_ARPESbands'
    elseif(use_band.eq.BandTBtimeZ)then
      Bandstring= '_TBxZ'
    endif
    if(use_Z.eq.Zeq1)then
      Zstring = '_Zeq1.dat'
    elseif(use_Z.eq.Znorman)then
      Zstring = '_ZNorman.dat'
    endif 
    suffix = trim(bandstring) // trim(Zstring)
    print*, 'Using file suffix = ', trim(suffix)
  end subroutine get_suffix
   
  double precision function energy(kx,ky,use_band)
    integer, intent(in) :: use_band
    double precision, intent(in) :: kx, ky
    if(use_band.eq.BandARPES)then
      energy = ARPES_bands(kx,ky)
    elseif(use_band.eq.BandTB)then
      energy = TB_bands(kx,ky,band_to_compute)
    elseif(use_band.eq.BandTBtimeZ)then
      energy = TB_bands(kx,ky,band_to_compute)*Zeff 
    endif
  end function energy

  double precision function ARPES_Bands(kx,ky)
    double precision, intent(in) :: kx, ky
    double precision, parameter :: t0 = 0.1305d0-0.023d0   
    double precision, parameter :: t1 =-0.5951d0
    double precision, parameter :: t2 = 0.1636d0
    double precision, parameter :: t3 =-0.0519d0
    double precision, parameter :: t4 =-0.1117d0
    double precision, parameter :: t5 = 0.0510d0 
    ARPES_Bands = t0 + 0.5*t1*(cos(kx)+cos(ky)) + t2*cos(kx)*cos(ky) &
                         + 0.5*t3*(cos(2*kx)+cos(2*ky)) &
                         + 0.5*t4*(cos(2*kx)*cos(ky) + cos(kx)*cos(2*ky)) &
                         + t5*cos(2*kx)*cos(2*ky)
  end function ARPES_Bands
  
  subroutine get_phi(kx,ky,band,phi)
    integer, parameter :: lwork = 3*nbands-1
    integer info
    integer, intent(in) :: band
    double precision, intent(in) :: kx, ky
    double precision, dimension(1:nbands), intent(out) :: phi 
    double precision, dimension(1:nbands,1:nbands) :: H
    double precision, dimension(1:nbands) :: E
    double precision, dimension(1:lwork) :: work 
    call get_H(kx,ky,H)
    call dsyev('V','U',nbands,H,nbands,E,work,lwork,info)
    phi(:) = H(:,band)
    if(phi(1).lt.0.0d0)then
      phi = -phi
    endif
  end subroutine get_phi

  double precision function TB_bands(kx,ky,band)
    integer, parameter :: lwork = 3*nbands-1
    integer info
    integer, intent(in) :: band
    double precision, intent(in) :: kx, ky 
    double precision, dimension(1:nbands,1:nbands) :: H
    double precision, dimension(1:nbands) :: E
    double precision, dimension(1:lwork) :: work 
    call get_H(kx,ky,H)
    call dsyev('N','U',nbands,H,nbands,E,work,lwork,info)
    TB_bands = E(band)
  end function TB_bands

  subroutine get_H(kx,ky,H)
    double precision, intent(in) :: kx, ky
    double precision, dimension(1:nbands,1:nbands), intent(out) :: H
    H = 0.0d0
    H(Cudx2y2,Cudx2y2) = ed-mu
    H(Ox,Ox) = ep - mu + 4*tppi*sin(kx/2)*sin(kx/2)
    H(Oy,Oy) = ep - mu + 4*tppi*sin(ky/2)*sin(ky/2)
    H(Cudx2y2,Ox) = 2*tpd*sin(kx/2)  
    H(Ox,Cudx2y2) = 2*tpd*sin(kx/2)
    H(Cudx2y2,Oy) =-2*tpd*sin(ky/2)  
    H(Oy,Cudx2y2) =-2*tpd*sin(ky/2)
    H(Ox,Oy) = 4*tpp*sin(kx/2)*sin(ky/2) 
    H(Oy,Ox) = 4*tpp*sin(kx/2)*sin(ky/2) 
  end subroutine get_H

  double precision function quasiparticle_wgt(kx_in,ky_in,use_Z)
    integer, intent(in) :: use_Z
    double precision, intent(in) :: kx_in, ky_in
    double precision kx, ky
    double precision pi, theta_min, theta_max, theta
    pi = 2.0d0*asin(1.0d0)
    if(use_Z.eq.Zeq1)then
      quasiparticle_wgt = 1.0d0
    elseif(use_z.eq.Znorman)then
      kx = kx_in
      ky = ky_in
      if(kx.gt.pi)then
        kx = kx - 2.0d0*pi
      endif
      if(ky.gt.pi)then
        ky = ky - 2.0d0*pi
      endif
      kx = abs(kx)
      ky = abs(ky)
      theta_min = atan(kf/pi)
      theta_max = atan(pi/kf)
      if(kx.eq.0.0d0.and.ky.eq.0.0d0)then
        kx = 1e-7;
        ky = 1e-7;
      endif
      
      theta = atan(ky/kx)
      if(theta.ge.pi/4.0d0)then
        quasiparticle_wgt = (1.0d0-cos(theta-theta_max))/(1.0d0-cos(pi/4.0d0-theta_max))
      else
        quasiparticle_wgt = (1.0d0-cos(theta-theta_min))/(1.0d0-cos(pi/4.0d0-theta_min))
      endif
    endif
  end function quasiparticle_wgt 

  subroutine get_kf(use_band)
    integer, intent(in) :: use_band
    double precision kL, kR, kM, pi, EL, ER, EM
    pi = 2.0d0*asin(1.0d0)
    kL = 0.0d0
    kR = pi
    kM = 0.5d0*(kR+kL)
    EM = energy(kM,pi,use_band)
    do while(abs(EM).ge.1e-9)
      EL = energy(kL,pi,use_band)
      EM = energy(kM,pi,use_band)
      ER = energy(kR,pi,use_band)
      if(EM.gt.0.0d0)then
        kr = km
      elseif(Em.lt.0.0d0)then
        kl = km
      endif
      km = 0.5d0*(kr+kl)
    end do
    kf = KM 
    print*, 'kf/pi = ', kf/pi
  end subroutine get_kf  

end module mod_bandstructure
