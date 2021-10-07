program main
  use mod_bandstructure
  use mod_phonon_dispersion
  implicit none
  integer, parameter :: nk = 16
  integer, parameter :: numw = 1200
  integer, parameter :: use_band = 2
  integer, parameter :: use_Z = 2
  integer nkx, nky, npx, npy, nqx, nqy     !counters for momentum
  integer i, j, k                          !Other counters
  integer nw, nwp
  double precision, parameter :: kb = 8.617e-5
  double precision, parameter :: temp = 20.0d0
  double precision, parameter :: beta = 1.0d0/(temp*kb)
  double precision, parameter :: sig = 0.025d0
  double precision, parameter :: eta = 0.0005d0
  double precision Mq
  double precision epsbx
  double precision epsby
  double precision epscx
  double precision epscy
  double precision wb1g, wa1g, wbr, wapex
  double precision gb1g, ga1g, gbr, gapex
  double precision g0b1g, g0a1g, g0br, g0apex
  double precision g2b1g, g2a1g, g2br, g2apex
  double precision Zk, Zp
  double precision pi
  double precision filling
  double precision dw, w, wp, den
  double precision wscale
  double precision pre
  double precision kx,ky,px,py,qx,qy
  double precision starttime
  double precision endtime
  double precision rtmp
  double precision delta 
  double precision Ek, Ep, Fermik, Fermip
  double precision bose_br, bose_b1g, bose_a1g, bose_apex
  double precision lambda_br, lambda_b1g, lambda_a1g, lambda_apex
  double precision, dimension(-Nk:Nk,-Nk:Nk) :: epsbx_lut
  double precision, dimension(-Nk:Nk,-Nk:Nk) :: epsby_lut
  double precision, dimension(-Nk:Nk,-Nk:Nk) :: epscx_lut
  double precision, dimension(-Nk:Nk,-Nk:Nk) :: epscy_lut
  double precision, dimension(-Nk:2*Nk,-Nk:2*Nk) :: dk_lut
  double precision, dimension(-Nk:2*Nk,-Nk:2*Nk) :: Ek_LUT,Zk_LUT
  double precision, dimension(-Nk:2*Nk,-Nk:2*Nk,1:nbands) :: Phi_LUT
  double precision, dimension(-Nk:2*Nk,-Nk:2*Nk) :: fermi_LUT, delta_LUT
  double precision, dimension(1:nbands) :: phik
  double precision, dimension(1:nbands) :: phip
  double precision, dimension(0:Nk,-numw:numw) :: SE2_br
  double precision, dimension(0:Nk,-numw:numw) :: SE2_b1g
  double precision, dimension(0:Nk,-numw:numw) :: SE2_a1g
  double precision, dimension(0:Nk,-numw:numw) :: SE2_apex
  double precision, dimension(0:Nk,-numw:numw) :: SE1_br
  double precision, dimension(0:Nk,-numw:numw) :: SE1_b1g
  double precision, dimension(0:Nk,-numw:numw) :: SE1_a1g
  double precision, dimension(0:Nk,-numw:numw) :: SE1_apex
  double precision, dimension(0:Nk,0:numw) :: Pol1_br
  double precision, dimension(0:Nk,0:numw) :: Pol1_b1g
  double precision, dimension(0:Nk,0:numw) :: Pol1_a1g
  double precision, dimension(0:Nk,0:numw) :: Pol1_apex
  double precision, dimension(0:Nk,0:numw) :: Pol2_br
  double precision, dimension(0:Nk,0:numw) :: Pol2_b1g
  double precision, dimension(0:Nk,0:numw) :: Pol2_a1g
  double precision, dimension(0:Nk,0:numw) :: Pol2_apex
  double precision, dimension(0:Nk,-numw:numw) :: A_br
  double precision, dimension(0:Nk,-numw:numw) :: A_b1g
  double precision, dimension(0:Nk,-numw:numw) :: A_a1g
  double precision, dimension(0:Nk,-numw:numw) :: A_apex
  double complex Pol
  double complex SE
  double complex ctmp
  double complex, parameter :: eye = (0.0d0,1.0d0)
  double complex, parameter :: cone = (1.0d0,0.0d0)
  double complex, dimension(0:Nk,0:numw) :: D_br
  double complex, dimension(0:Nk,0:numw) :: D_b1g
  double complex, dimension(0:Nk,0:numw) :: D_a1g
  double complex, dimension(0:Nk,0:numw) :: D_apex
  character suffix*200, filename1*200, filename2*200, filename3*200
  character filename4*200

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
  pre = 1.0d0/dfloat(4*nk*nk)
  call get_suffix(use_band,use_Z,suffix)
  print*, 'Computing the polarizability phonon dispersion along (H,0).'

  !If we are using Mike Norman's PG model, we need to determine kf
  call get_kf(use_band)
    
  !Build the look-up tables
  filename1 = 'Ek' // trim(suffix)
  filename2 = 'Zk' // trim(suffix)
  filename3 = 'Phi' // trim(suffix)
  open(file=filename1,unit=50,action='write')
  open(file=filename2,unit=51,action='write')
  open(file=filename3,unit=52,action='write')
  filling = 0.0d0
  do nkx = -Nk,2*Nk
    kx = dfloat(nkx)*pi/dfloat(Nk)
    do nky = -Nk,2*Nk
      ky  = dfloat(nky)*pi/dfloat(Nk)
      EK_lut(nkx,nky) = energy(kx,ky,use_band)
      Zk_lut(nkx,nky) = quasiparticle_wgt(kx,ky,use_Z)
      call get_phi(kx,ky,band_to_compute,phik)
      Phi_lut(nkx,nky,:) = phik(:) 
      Fermi_lut(nkx,nky) = fermi(Ek_lut(nkx,nky),beta)
      delta_lut(nkx,nky) = lorentzian(Ek_lut(nkx,nky),sig)
      dk_lut(nkx,nky) = 0.5d0*(cos(kx/2.0d0)-cos(ky/2.0d0))
      if(nkx.le.Nk.and.nky.le.Nk)then
        write(50,*) kx/pi,ky/pi,Ek_lut(nkx,nky)
        write(51,*) kx/pi,ky/pi,Zk_lut(nkx,nky)
        write(52,*) kx/pi,ky/pi,(Phi_lut(nkx,nky,i),i=1,nbands)

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
     
      if(nkx.lt.Nk.and.nky.lt.Nk)then
        filling = filling+2.0d0*pre*fermi(Ek_lut(nkx,nky),beta)
      endif
    enddo
    write(50,*) ' '
    write(51,*) ' '
    write(52,*) ' '
  enddo
  close(unit=50)
  close(unit=51)
  close(unit=52)

  print*, 'Computing lambda.'
  lambda_b1g = 0.0d0 
  lambda_a1g = 0.0d0
  lambda_br = 0.0d0
  lambda_apex = 0.0d0
  call cpu_time(starttime)
  pre = 2.0d0/dfloat(4*nk*nk)
  rtmp = 0.0d0
  do nkx = 0,Nk
    kx = dfloat(nkx)*pi/dfloat(Nk)
    do nky = 0,nkx
      ky = dfloat(nky)*pi/dfloat(Nk)
      rtmp = rtmp + 8.0d0*delta_lut(nkx,nky)

      do nqx = -NK+1,Nk
        qx = dfloat(nqx)*pi/dfloat(nk)
        px = kx+qx
        npx = nkx + nqx
        do nqy = -Nk+1,Nk
          qy = dfloat(nqy)*pi/dfloat(nk)
          npy = nky+nqy
          py = ky + qy

          phik(:) = Phi_lut(nkx,nky,:)
          phip(:) = Phi_lut(npx,npy,:)
          epscx = epscx_lut(nqx,nqy)
          epscy = epscy_lut(nqx,nqy)
          epsbx = epsbx_lut(nqx,nqy)
          epsby = epsby_lut(nqx,nqy)
  
          wb1g = get_wb1g(qx,qy)
          wa1g = get_wa1g(qx,qy)
          wbr = get_wbr(qx,qy)
          wapex = get_wapex(qx,qy)

          gapex = g0apex*dk_lut(nkx,nky)*dk_lut(npy,npy)
          Gb1g = g0b1g*(epscx*Phik(Ox)*Phip(Ox)-epscy*Phik(Oy)*Phip(Oy));
          Ga1g = g0a1g*(epscx*Phik(Ox)*Phip(Ox)+epscy*Phik(Oy)*Phip(Oy));
          Gbr = g0br*epsbx*(cos(px/2)*Phip(Cudx2y2)*Phik(Ox)-cos(kx/2)*Phik(Cudx2y2)*Phip(Ox)) &
              - g0br*epsby*(cos(py/2)*Phip(Cudx2y2)*Phik(Oy)-cos(ky/2)*Phik(Cudx2y2)*Phip(Oy)); 
     
          lambda_br = lambda_br + delta_lut(nkx,nky)*delta_lut(npx,npy)*gbr*gbr/wbr
          lambda_apex = lambda_apex + delta_lut(nkx,nky)*delta_lut(npx,npy)*gapex*gapex/wapex
          lambda_b1g = lambda_b1g + delta_lut(nkx,nky)*delta_lut(npx,npy)*gb1g*gb1g/wb1g
          lambda_a1g = lambda_a1g + delta_lut(nkx,nky)*delta_lut(npx,npy)*ga1g*ga1g/wa1g
        enddo
      enddo
    enddo
  enddo
  call cpu_time(endtime)
  pre = pre/rtmp
  print*, 'Lambda_br = ', 8.0d0*pre*lambda_br
  print*, 'Lambda_apex = ', 8.0d0*pre*lambda_apex
  print*, 'Lambda_b1g = ', 8.0d0*pre*lambda_b1g
  print*, 'Lambda_a1g = ', 8.0d0*pre*lambda_a1g
  print*, 'Time = ', endtime-starttime

  wscale = 1.2d0 
  dw = wscale/dfloat(numw)
  print*, 'Wscale = ', wscale
  print*, 'dw = ', dw
  print*, 'Filling = ', filling 
  print*, 'Computing the polarizability.'

  pre = 1.0d0/dfloat(4*nk*nk)
  Pol2_apex = 0.0d0
  Pol2_br = 0.0d0
  Pol2_b1g = 0.0d0
  Pol2_a1g = 0.0d0

  call cpu_time(starttime)
  do nqx = 0,nk
    qx = pi*dfloat(nqx)/nk
    nqy = 0
    qy = 0.0d0

    if(nqx.eq.0)then
      qx = 1e-6
    endif

    epscx = epscx_lut(nqx,nqy)
    epscy = epscy_lut(nqx,nqy)
    epsbx = epsbx_lut(nqx,nqy)
    epsby = epsby_lut(nqx,nqy)
    
    wb1g = get_wb1g(qx,qy)
    wa1g = get_wa1g(qx,qy)
    wbr = get_wbr(qx,qy)
    wapex = get_wapex(qx,qy)
   
    do nkx = -Nk+1,Nk
      kx = dfloat(nkx)*pi/dfloat(Nk)
      npx = nkx + nqx
      px = kx + qx
      do nky = -nk+1,Nk
        ky = dfloat(nky)*pi/dfloat(Nk)
        npy = nky + nqy
        py =  ky + qy

        Ek = Ek_lut(nkx,nky) 
        Zk = Zk_lut(nkx,nky)
        Fermik = Fermi_lut(nkx,nky)
        phik(:) = Phi_lut(nkx,nky,:)

        Ep = Ek_lut(npx,npy)
        Zp = Zk_lut(npx,npy)
        Fermip = Fermi_lut(npx,npy)
        phip(:) = Phi_lut(npx,npy,:)

        gapex = g0apex*dk_lut(nkx,nky)*dk_lut(npy,npy)
        Gb1g = g0b1g*(epscx*Phik(Ox)*Phip(Ox)-epscy*Phik(Oy)*Phip(Oy));
        Ga1g = g0a1g*(epscx*Phik(Ox)*Phip(Ox)+epscy*Phik(Oy)*Phip(Oy));
        Gbr = g0br*epsbx*(cos(px/2)*Phip(Cudx2y2)*Phik(Ox)-cos(kx/2)*Phik(Cudx2y2)*Phip(Ox)) &
            - g0br*epsby*(cos(py/2)*Phip(Cudx2y2)*Phik(Oy)-cos(ky/2)*Phik(Cudx2y2)*Phip(Oy)); 

        g2apex = gapex*gapex
        g2b1g = gb1g*gb1g
        g2a1g = ga1g*ga1g
        g2br = gbr*gbr

        rtmp = 2.0d0*pre*pi*Zk*Zp*(Fermik-Fermip)
        do nw = 0,numw
          w = dfloat(nw)*dw
          delta = lorentzian(w + Ek - Ep,sig) 
          Pol2_br(nqx,nw)  = Pol2_br(nqx,nw) - rtmp*g2br*delta
          Pol2_b1g(nqx,nw)  = Pol2_b1g(nqx,nw) - rtmp*g2b1g*delta
          Pol2_a1g(nqx,nw)  = Pol2_a1g(nqx,nw) - rtmp*g2a1g*delta
          Pol2_apex(nqx,nw)  = Pol2_apex(nqx,nw) - rtmp*g2apex*delta
        enddo
      enddo
    enddo
  enddo

  !Now compute the real part
  Pol1_apex = 0.0d0
  Pol1_br = 0.0d0
  Pol1_b1g = 0.0d0
  Pol1_a1g = 0.0d0
  do nqx = 0,Nk 
    nqy = 0.0d0
    qx = pi*dfloat(nqx)/Nk
    qy = pi*dfloat(nqy)/Nk

    wb1g = get_wb1g(qx,qy)
    wa1g = get_wa1g(qx,qy)
    wbr = get_wbr(qx,qy)
    wapex = get_wapex(qx,qy)

    do nw = 0,numw-1
      w = dfloat(nw)*dw
      do nwp = 0,numw
        wp = dfloat(nwp)*dw
        den = wp*wp-w*w
        if(den.ne.0.0d0)then
          rtmp = (2.0d0*dw/pi)*(1.0d0/den)
          Pol1_apex(nqx,nw) = Pol1_apex(nqx,nw) + rtmp*(wp*Pol2_apex(nqx,nwp)-w*Pol2_apex(nqx,nw))
          Pol1_a1g(nqx,nw) = Pol1_a1g(nqx,nw) + rtmp*(wp*Pol2_a1g(nqx,nwp)-w*Pol2_a1g(nqx,nw))
          Pol1_b1g(nqx,nw) = Pol1_b1g(nqx,nw) + rtmp*(wp*Pol2_b1g(nqx,nwp)-w*Pol2_b1g(nqx,nw))
          Pol1_br(nqx,nw) = Pol1_br(nqx,nw) + rtmp*(wp*Pol2_br(nqx,nwp)-w*Pol2_br(nqx,nw))
        endif        
      enddo
      rtmp = log(abs(wscale-w)/abs(wscale+w))/pi
      Pol1_apex(nqx,nw) = Pol1_apex(nqx,nw) + rtmp*Pol2_apex(nqx,nw)
      Pol1_a1g(nqx,nw) = Pol1_a1g(nqx,nw) + rtmp*Pol2_a1g(nqx,nw)
      Pol1_b1g(nqx,nw) = Pol1_b1g(nqx,nw) + rtmp*Pol2_b1g(nqx,nw)
      Pol1_br(nqx,nw) = Pol1_br(nqx,nw) + rtmp*Pol2_br(nqx,nw)

      Pol = dcmplx(Pol1_apex(nqx,nw),Pol2_apex(nqx,nw))
      ctmp = cone/(w*w  - wapex*wapex + eye*eta - 2.0d0*wapex*Pol)
      D_apex(nqx,nw) = 2.0d0*wapex*ctmp

      Pol = dcmplx(Pol1_a1g(nqx,nw),Pol2_a1g(nqx,nw))
      ctmp = cone/(w*w  - wa1g*wa1g + eye*eta - 2.0d0*wa1g*Pol)
      D_a1g(nqx,nw) = 2.0d0*wa1g*ctmp

      Pol = dcmplx(Pol1_b1g(nqx,nw),Pol2_b1g(nqx,nw))
      ctmp = cone/(w*w  - wb1g*wb1g + eye*eta - 2.0d0*wb1g*Pol)
      D_b1g(nqx,nw) = 2.0d0*wb1g*ctmp
    
      Pol = dcmplx(Pol1_br(nqx,nw),Pol2_br(nqx,nw))
      ctmp = cone/(w*w  - wbr*wbr + eye*eta - 2.0d0*wbr*Pol)
      D_br(nqx,nw) = 2.0d0*wbr*ctmp
    enddo 
  enddo
  call cpu_time(endtime)
  print*, 'Time = ', (endtime-starttime)

  filename1 = 'Pol1' // trim(suffix)
  filename2 = 'Pol2' // trim(suffix)
  filename3 = 'ReDqw' // trim(suffix)
  filename4 = 'ImDqw' // trim(suffix)
  open(file=filename1,unit=50,action='write')
  open(file=filename2,unit=51,action='write')
  open(file=filename3,unit=52,action='write')
  open(file=filename4,unit=53,action='write')
  do nqx = 0,Nk
    qx = dfloat(nqx)*pi/dfloat(Nk)
    do nw = 0,numw
      w = dw*dfloat(nw)
      write(50,*) qx/pi, w, Pol1_apex(nqx,nw), Pol1_a1g(nqx,nw), Pol1_b1g(nqx,nw), Pol1_br(nqx,nw)
      write(51,*) qx/pi, w, Pol2_apex(nqx,nw), Pol2_a1g(nqx,nw), Pol2_b1g(nqx,nw), Pol2_br(nqx,nw)
      if(w.le.0.2d0)then
        write(52,*) qx/pi, w, real(D_apex(nqx,nw)), real(D_a1g(nqx,nw)),&
                              real(D_b1g(nqx,nw)), real(D_br(nqx,nw))
        write(53,*) qx/pi, w, imag(D_apex(nqx,nw)), imag(D_a1g(nqx,nw)), &
                              imag(D_b1g(nqx,nw)), imag(D_br(nqx,nw))
      endif
    enddo
    write(50,*) ' '
    write(51,*) ' '
    write(52,*) ' '
    write(53,*) ' '
  enddo
  close(unit=50)
  close(unit=51)
  close(unit=52)
  close(unit=53)

  print*, 'Computing the electron spectral function.'
  call cpu_time(starttime)
  SE2_br = 0.0d0
  SE2_a1g = 0.0d0
  SE2_b1g = 0.0d0
  SE2_apex = 0.0d0
  SE1_br = 0.0d0
  SE1_a1g = 0.0d0
  SE1_b1g = 0.0d0
  SE1_apex = 0.0d0
  pre = pi/dfloat(4*nk*nk)
  do nkx = 0,Nk
    kx = dfloat(nkx)*pi/dfloat(nk)
    ky = kx
    nky = nkx

    !sum over q
    do nqx = -Nk+1,Nk
      qx = dfloat(nqx)*pi/dfloat(Nk)
      npx = nkx+nqx
      px = kx + qx
      do nqy = -Nk+1,Nk
        qy =  dfloat(nqy)*pi/dfloat(Nk)
        npy = nky + nqy
        py = ky + qy 
   
        wb1g = get_wb1g(qx,qy)
        wa1g = get_wa1g(qx,qy)
        wbr = get_wbr(qx,qy)
        wapex = get_wapex(qx,qy)

        bose_br = bose(wbr,beta)
        bose_b1g = bose(wb1g,beta)
        bose_a1g = bose(wa1g,beta)
        bose_apex = bose(wapex,beta)

        Ep = Ek_lut(npx,npy)
        Fermip = Fermi_lut(npx,npy)
        Zp = Zk_lut(npx,npy)

        epscx = epscx_lut(nqx,nqy)
        epscy = epscy_lut(nqx,nqy)
        epsbx = epsbx_lut(nqx,nqy)
        epsby = epsby_lut(nqx,nqy)

        gapex = g0apex*dk_lut(nkx,nky)*dk_lut(npy,npy)
        Gb1g = g0b1g*(epscx*Phik(Ox)*Phip(Ox)-epscy*Phik(Oy)*Phip(Oy));
        Ga1g = g0a1g*(epscx*Phik(Ox)*Phip(Ox)+epscy*Phik(Oy)*Phip(Oy));
        Gbr = g0br*epsbx*(cos(px/2)*Phip(Cudx2y2)*Phik(Ox)-cos(kx/2)*Phik(Cudx2y2)*Phip(Ox)) &
            - g0br*epsby*(cos(py/2)*Phip(Cudx2y2)*Phik(Oy)-cos(ky/2)*Phik(Cudx2y2)*Phip(Oy)); 

        g2apex = gapex*gapex
        g2b1g = gb1g*gb1g
        g2a1g = ga1g*ga1g
        g2br = gbr*gbr
        
        do nw = -numw,numw
          w = dfloat(nw)*dw
          SE2_br(nkx,nw) = SE2_br(nkx,nw) - pre*Zp*g2br*(bose_br+fermip)*lorentzian(w+wbr-Ep,sig) & 
                                          - pre*Zp*g2br*(bose_br+1.0d0-fermip)*lorentzian(w-wbr-Ep,sig)
          SE2_apex(nkx,nw) = SE2_apex(nkx,nw) - pre*Zp*g2apex*(bose_apex+fermip)*lorentzian(w+wapex-Ep,sig) & 
                                              - pre*Zp*g2apex*(bose_apex+1.0d0-fermip)*lorentzian(w-wapex-Ep,sig)
          SE2_b1g(nkx,nw) = SE2_b1g(nkx,nw) - pre*Zp*g2b1g*(bose_b1g+fermip)*lorentzian(w+wb1g-Ep,sig) & 
                                            - pre*Zp*g2b1g*(bose_b1g+1.0d0-fermip)*lorentzian(w-wb1g-Ep,sig)
          SE2_a1g(nkx,nw) = SE2_a1g(nkx,nw) - pre*Zp*g2a1g*(bose_a1g+fermip)*lorentzian(w+wa1g-Ep,sig) & 
                                            - pre*Zp*g2a1g*(bose_a1g+1.0d0-fermip)*lorentzian(w-wa1g-Ep,sig)
        enddo
      enddo
    enddo
  enddo

  !Compute the real part via the KK relations 
  do nkx = 0,Nk
    do nw = -numw+1,numw-1
      w = dfloat(nw)*dw
      do nwp = -numw,numw
        wp = dfloat(nwp)*dw
        den = wp-w
        rtmp = (dw/pi)/den
        if(den.ne.0.0d0)then
          SE1_apex(nkx,nw) =  SE1_apex(nkx,nw) + rtmp*(SE2_apex(nkx,nwp)-SE2_apex(nkx,nw))
          SE1_a1g(nkx,nw) =  SE1_a1g(nkx,nw) + rtmp*(SE2_a1g(nkx,nwp)-SE2_a1g(nkx,nw))
          SE1_b1g(nkx,nw) =  SE1_b1g(nkx,nw) + rtmp*(SE2_b1g(nkx,nwp)-SE2_b1g(nkx,nw))
          SE1_br(nkx,nw) =  SE1_br(nkx,nw) + rtmp*(SE2_br(nkx,nwp)-SE2_br(nkx,nw))
        endif
      enddo
      rtmp = dw*log(abs(wscale-w)/abs(wscale+w))/pi
      SE1_apex(nkx,nw) = SE1_apex(nkx,nw) + rtmp*SE2_apex(nkx,nw)
      SE1_a1g(nkx,nw) = SE1_a1g(nkx,nw) + rtmp*SE2_a1g(nkx,nw)
      SE1_b1g(nkx,nw) = SE1_b1g(nkx,nw) + rtmp*SE2_b1g(nkx,nw)
      SE1_br(nkx,nw) = SE1_br(nkx,nw) + rtmp*SE2_br(nkx,nw)
    enddo
  enddo

  do nkx = 0,Nk
    kx = dfloat(nkx)*pi/dfloat(nk)
    ky = kx
    nky = nkx
    Ek = Ek_lut(nkx,nky)
    rtmp = 2.0d0*Zk_lut(nkx,nky)*fermi_lut(nkx,nky)/pi
    do nw = -numw,numw
      w = dfloat(nw)*dw

      SE = dcmplx(SE1_br(nkx,nw),SE2_br(nkx,nw) - sig)
      A_br(nkx,nw) = -rtmp*imag(cone/(w - Ek - SE))
    
      SE = dcmplx(SE1_a1g(nkx,nw),SE2_a1g(nkx,nw) - sig)
      A_a1g(nkx,nw) = -rtmp*imag(cone/(w - Ek - SE))
    
      SE = dcmplx(SE1_b1g(nkx,nw),SE2_b1g(nkx,nw) - sig)
      A_b1g(nkx,nw) = -rtmp*imag(cone/(w - Ek - SE))
    
      SE = dcmplx(SE1_apex(nkx,nw),SE2_apex(nkx,nw) - sig)
      A_apex(nkx,nw) = -rtmp*imag(cone/(w - Ek - SE))
    enddo
  enddo

  call cpu_time(endtime)
  print*, 'Time = ', endtime-starttime
  
  filename1 = 'SE1' // trim(suffix)
  filename2 = 'SE2' // trim(suffix)
  filename3 = 'Akw' // trim(suffix)
  open(file=filename1,unit=50,action='write')
  open(file=filename2,unit=51,action='write')
  open(file=filename3,unit=52,action='write')
  do nkx = 0,Nk
    kx = dfloat(nkx)*pi/dfloat(Nk)
    do nw = -numw,numw
      w = dw*dfloat(nw)
      write(50,*) kx/pi, w, SE1_apex(nkx,nw), SE1_a1g(nkx,nw), SE1_b1g(nkx,nw), SE1_br(nkx,nw)
      write(51,*) kx/pi, w, SE2_apex(nkx,nw), SE2_a1g(nkx,nw), SE2_b1g(nkx,nw), SE2_br(nkx,nw)

      if(w.gt.-0.6d0.and.w.le.0.05d0)then
        write(52,*) kx/pi, w, A_apex(nkx,nw), A_a1g(nkx,nw), A_b1g(nkx,nw), A_br(nkx,nw)
      endif
    enddo
    write(50,*) ' '
    write(51,*) ' '
    write(52,*) ' '
  enddo
  close(unit=50)
  close(unit=51)
  close(unit=52)

  stop
contains 
  double precision function lorentzian(arg,gam)
    double precision, intent(in) :: arg,gam
    double precision pi
    pi = 2.0d0*asin(1.0d0)
    lorentzian = (gam/pi)/(arg*arg + gam*gam)
  end function lorentzian

  double precision function fermi(E,beta)
    double precision, intent(in) :: E, beta
    fermi = 1.0d0/(exp(E*beta)+1.0d0)
  end function fermi

  double precision function bose(E,beta)
    double precision, intent(in) :: E, beta
    if(E.le.0.0d0)then
      bose = 0.0d0
    else
      bose = 1.0d0/(exp(E*beta)-1.0d0)
    endif
  end function bose
end program main
