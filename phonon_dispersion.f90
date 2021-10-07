module mod_phonon_dispersion
  implicit none
  double precision, parameter, public :: wb1g0 = 0.045d0
  double precision, parameter, public :: wa1g0 = 0.055d0
  double precision, parameter, public :: wapex0 = 0.055d0
  double precision, parameter, public :: wbr0 = 0.090d0
  double precision, parameter, public :: wbr1 = 0.30d0 
  double precision, parameter, public :: wac0 = 0.018d0
contains
  double precision function get_wb1g(qx,qy)
    double precision, intent(in) :: qx
    double precision, intent(out) :: qy
    get_wb1g = wb1g0
  end function 

  double precision function get_wa1g(qx,qy)
    double precision, intent(in) :: qx
    double precision, intent(out) :: qy
    get_wa1g = wa1g0
  end function 

  double precision function get_wapex(qx,qy)
    double precision, intent(in) :: qx
    double precision, intent(out) :: qy
    get_wapex = wapex0
  end function 

  double precision function get_wbr(qx,qy)
    double precision, intent(in) :: qx
    double precision, intent(out) :: qy
    get_wbr = wbr0*(1.0d0-wbr1*(sin(qx/2.0d0)*sin(qx/2.0d0)+sin(qy/2.0d0)*sin(qy/2.0d0)))
  end function 

  double precision function get_wac(qx,qy)
    double precision, intent(in) :: qx
    double precision, intent(out) :: qy
    get_wac = wac0*sqrt(sin(qx/2.0d0)*sin(qx/2.0d0)+sin(qy/2.0d0)*sin(qy/2.0d0))
  end function 

end module mod_phonon_dispersion
