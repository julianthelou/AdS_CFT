! NPScalars
! NP_calc_psi_gfs.F90 : Actual calculation of the NP scalar grid functions
!                         (Psi4 for now)
!
!=============================================================================

#define REAL_PART 1
#define IMAGINARY_PART 2

!need to determine input arguments
SUBROUTINE CalcPsi4(Q, gradQ, psi0, psi4, phi0, phi1, phi2)
  USE typesDef, ONLY : nVar, nDim
  implicit none
  REAL, INTENT(IN)  :: Q(nVar),  gradQ(nVar,nDim) ! CCZ4 state vector and its first derivative
  REAL, INTENT(OUT) :: psi0(2), psi4(2), phi0(2), phi1(2), phi2(2)  ! Complex numbers as vector, shape (re,im)

  REAL :: gd(3,3), gu(3,3), kd(3,3), detgd, mu, g_cov(3,3), Aex(3,3), dAex(3,3,3)
  REAL :: d1_gd(3,3,3), d1_gu(3,3,3), d2_gd(3,3,3,3), d1_kd(3,3,3),      &
            cd1_kd(3,3,3)
  REAL :: cf1(3,3,3), cf2(3,3,3), ri(3,3), d1_cf2(3,3,3,3)
  REAL :: dx(3), xx(3)
  REAL :: u_vec(3), v_vec(3), w_vec(3), ud_vec(3), dotp1, dotp2
  REAL :: eps_lc_d(3,3,3), eps_lc_u(3,3,3), elec(3,3), mag(3,3), electr, magtr
  REAL :: W, traceK, dtraceK(3), phi, PP(3), dPP(3,3), DD(3,3,3), dDD(3,3,3,3)
 

  INTEGER :: i, j, k, m, n, p, a, b, c, d

  xx(1) = Q(61)
  xx(2) = Q(62)
  xx(3) = Q(63)

  phi = Q(55)             !conformal factor from the state vector
  W = EXP(-2. * phi) 
  PP = Q(56:58)           !first derivative of the conformal factor
  dPP(:, 1) = gradQ(56,:) !second derivatives of the conformal factor
  dPP(:, 2) = gradQ(57,:)
  dPP(:, 3) = gradQ(58,:) 

  !=== Initialize grid functions as zero ===
  !Need to determine if we calculate these or psi0 and psi4, Phi0, 1, 2
  !psi4re = 0
  !psi4im = 0
  psi4 = 0

  !psi0 and NPScalar terms, not currently calculated
  psi0 = 0
  phi0 = 0
  phi1 = 0
  phi2 = 0

  g_cov(1,1) = Q(1)
  g_cov(1,2) = Q(2)
  g_cov(1,3) = Q(3)
  g_cov(2,1) = Q(2)
  g_cov(2,2) = Q(4)
  g_cov(2,3) = Q(5)
  g_cov(3,1) = Q(3)
  g_cov(3,2) = Q(5)
  g_cov(3,3) = Q(6)

  traceK = Q(54)
  dtraceK(:) = gradQ(54,:) 


  Aex(1,1) = Q(7) 
  Aex(1,2) = Q(8) 
  Aex(1,3) = Q(9) 
  Aex(2,1) = Q(8) 
  Aex(2,2) = Q(10) 
  Aex(2,3) = Q(11) 
  Aex(3,1) = Q(9) 
  Aex(3,2) = Q(11) 
  Aex(3,3) = Q(12)

  dAex(:,1,1) = gradQ(7,:)
  dAex(:,1,2) = gradQ(8,:)
  dAex(:,1,3) = gradQ(9,:)
  dAex(:,2,1) = gradQ(8,:)
  dAex(:,2,2) = gradQ(10,:)
  dAex(:,2,3) = gradQ(11,:)
  dAex(:,3,1) = gradQ(9,:)
  dAex(:,3,2) = gradQ(11,:)
  dAex(:,3,3) = gradQ(12,:)

  DD(1,1,1)=Q(36) 
  DD(1,1,2)=Q(37) 
  DD(1,1,3)=Q(38) 
  DD(1,2,1)=Q(37) 
  DD(1,2,2)=Q(39) 
  DD(1,2,3)=Q(40)
  DD(1,3,1)=Q(38) 
  DD(1,3,2)=Q(40) 
  DD(1,3,3)=Q(41)
  ! 
  DD(2,1,1)=Q(42) 
  DD(2,1,2)=Q(43) 
  DD(2,1,3)=Q(44) 
  DD(2,2,1)=Q(43) 
  DD(2,2,2)=Q(45) 
  DD(2,2,3)=Q(46)
  DD(2,3,1)=Q(44) 
  DD(2,3,2)=Q(46) 
  DD(2,3,3)=Q(47)
  !
  DD(3,1,1)=Q(48) 
  DD(3,1,2)=Q(49) 
  DD(3,1,3)=Q(50) 
  DD(3,2,1)=Q(49) 
  DD(3,2,2)=Q(51) 
  DD(3,2,3)=Q(52)
  DD(3,3,1)=Q(50) 
  DD(3,3,2)=Q(52) 
  DD(3,3,3)=Q(53)

  !DD = 1/2 * deriv(conformal gamma_ij)
  DD = 2*DD

  dDD(:,1,1,1)=gradQ(36,:) 
  dDD(:,1,1,2)=gradQ(37,:) 
  dDD(:,1,1,3)=gradQ(38,:) 
  dDD(:,1,2,1)=gradQ(37,:) 
  dDD(:,1,2,2)=gradQ(39,:) 
  dDD(:,1,2,3)=gradQ(40,:)
  dDD(:,1,3,1)=gradQ(38,:) 
  dDD(:,1,3,2)=gradQ(40,:) 
  dDD(:,1,3,3)=gradQ(41,:)
  dDD(:,2,1,1)=gradQ(42,:) 
  dDD(:,2,1,2)=gradQ(43,:) 
  dDD(:,2,1,3)=gradQ(44,:) 
  dDD(:,2,2,1)=gradQ(43,:) 
  dDD(:,2,2,2)=gradQ(45,:) 
  dDD(:,2,2,3)=gradQ(46,:)
  dDD(:,2,3,1)=gradQ(44,:) 
  dDD(:,2,3,2)=gradQ(46,:) 
  dDD(:,2,3,3)=gradQ(47,:) 
  dDD(:,3,1,1)=gradQ(48,:) 
  dDD(:,3,1,2)=gradQ(49,:) 
  dDD(:,3,1,3)=gradQ(50,:) 
  dDD(:,3,2,1)=gradQ(49,:) 
  dDD(:,3,2,2)=gradQ(51,:) 
  dDD(:,3,2,3)=gradQ(52,:)
  dDD(:,3,3,1)=gradQ(50,:) 
  dDD(:,3,3,2)=gradQ(52,:) 
  dDD(:,3,3,3)=gradQ(53,:)

  dDD = 2. * dDD

  !metric values from state vector are conformal quantities, so we'll transform them back to non-conformal quantities
  gd = W * g_cov
  kd = W *(Aex + 1./3. * traceK * g_cov)
    
  !-------------- Invert metric-original code----
  detgd =        gd(1,1) * gd(2,2) * gd(3,3)                                &
           + 2 * gd(1,2) * gd(1,3) * gd(2,3)                                &
           -     gd(1,1) * gd(2,3) ** 2                                     &
           -     gd(2,2) * gd(1,3) ** 2                                     &
           -     gd(3,3) * gd(1,2) ** 2
  gu(1,1) = (gd(2,2) * gd(3,3) - gd(2,3) ** 2     ) / detgd
  gu(2,2) = (gd(1,1) * gd(3,3) - gd(1,3) ** 2     ) / detgd
  gu(3,3) = (gd(1,1) * gd(2,2) - gd(1,2) ** 2     ) / detgd
  gu(1,2) = (gd(1,3) * gd(2,3) - gd(1,2) * gd(3,3)) / detgd
  gu(1,3) = (gd(1,2) * gd(2,3) - gd(1,3) * gd(2,2)) / detgd
  gu(2,3) = (gd(1,3) * gd(1,2) - gd(2,3) * gd(1,1)) / detgd
  gu(2,1) = gu(1,2)
  gu(3,1) = gu(1,3)
  gu(3,2) = gu(2,3)
  !----------------------------------------------

  !--------new derivatives based on state vector and conformal & transformed quantities
  do k = 1, 3
  do j = 1, 3
  do i = 1, 3
    d1_gd(k,i,j) =   W * (DD(k, i, j)                                       &
                   - 2. * g_cov(i,j) * PP(k))

    d1_kd(k,i,j) =   W * (dAex(k, i, j)                                     &
                   + 1./3. * (g_cov(i,j) * dtraceK(k)                       &
                   + traceK * DD(k,i,j)))                                   &
                   - 2. * kd(i,j) * PP(k)
    
    do m = 1, 3
      d2_gd(m,k,i,j) = - 2. * PP(m) * d1_gd(k,i,j)                          &
                       + W * (dDD(m,k,i,j)                                  &
                       - 2. * (DD(m,i,j) * PP(k) + g_cov(i,j) * dPP(m,k)))
    end do
  end do
  end do
  end do
  
  
  !------------------------------------------

  !------------ Christoffel symbols ---------
  cf1 = 0
  do a = 1, 3
    do b = 1, 3
      do c = 1, 3
        cf1(a,b,c) = 0.5d0 * (d1_gd(a,b,c) + d1_gd(a,c,b) - d1_gd(b,c,a))
      end do
    end do
  end do

  cf2 = 0
  do a = 1, 3
    do b = 1, 3
      do c = 1, 3
        do m = 1, 3
          cf2(a,b,c) = cf2(a,b,c) + gu(a,m) * cf1(m,b,c)
        end do
      end do
    end do
  end do
  !------------------------------------------


  !------------ covariant derivs ------------
  cd1_kd = d1_kd
  do a = 1, 3
    do b = 1, 3
      do c = 1, 3
        do m = 1, 3
          cd1_kd(a,b,c) = cd1_kd(a,b,c) - cf2(m,a,c) * kd(m,b)           &
                                        - cf2(m,b,c) * kd(a,m)
        end do
      end do
    end do
  end do
  !------------------------------------------


  !----- d1 of inverse metric and cf2 -------
  d1_gu = 0
  do a = 1, 3
    do b = 1, 3
      do c = 1, 3
        do m = 1, 3
          do n = 1, 3
            d1_gu(a,b,c) = d1_gu(a,b,c) - gu(a,m) * gu(b,n) * d1_gd(m,n,c)
          end do
        end do
      end do
    end do
  end do

  d1_cf2 = 0
  do a = 1, 3
    do b = 1, 3
      do c = 1, 3
        do d = 1, 3
          do m = 1, 3
            d1_cf2(a,b,c,d) = d1_cf2(a,b,c,d) + gu(a,m) * (d2_gd(c,m,b,d)&
                              + d2_gd(m,b,c,d) - d2_gd(b,c,m,d)) / 2     &
                              + d1_gu(a,m,d) * cf1(m,b,c)
          end do
        end do
      end do
    end do
  end do
  !------------------------------------------


  !------------ Ricci tensor ----------------
  ri = 0
  do a = 1, 3
    do b = 1, 3
      do m = 1, 3
        ri(a,b) = ri(a,b) + d1_cf2(m,b,a,m) - d1_cf2(m,m,a,b)
        do n = 1, 3
          ri(a,b) = ri(a,b) + cf2(m,m,n) * cf2(n,b,a)                    &
                            - cf2(m,b,n) * cf2(n,m,a)
        end do
      end do
    end do
  end do
  !------------------------------------------


  !------------ Levi-Civita tensor ----------
  eps_lc_u        = 0
  eps_lc_u(1,2,3) = 1
  eps_lc_u(2,3,1) = 1
  eps_lc_u(3,1,2) = 1
  eps_lc_u(3,2,1) = -1
  eps_lc_u(2,1,3) = -1
  eps_lc_u(1,3,2) = -1
  eps_lc_u = eps_lc_u / sqrt(detgd)

  eps_lc_d = eps_lc_u * detgd
  !------------------------------------------


  !------------ Orthonormal basis -----------
  ! All points on the z-axis are pathological, since the triad vectors
  ! in the phi and theta direction are not well-defined. Take points
  ! just a little off, say at x = +epsilon.

  if( xx(1)**2 + xx(2)**2 < 1.0d-12 ) xx(1) = xx(1) + 1.0d-10
  u_vec(:) = xx(:)
  v_vec(:) = (/ xx(1)*xx(3), xx(2)*xx(3), -xx(1)**2 - xx(2)**2 /)
  w_vec(:) = (/ -xx(2), xx(1), 0.0d0 /)

  ! Orthonormalization
  dotp1 =   gd(1,1) * u_vec(1) * u_vec(1) + gd(1,2) * u_vec(1) * u_vec(2)&
          + gd(1,3) * u_vec(1) * u_vec(3) + gd(2,1) * u_vec(2) * u_vec(1)&
          + gd(2,2) * u_vec(2) * u_vec(2) + gd(2,3) * u_vec(2) * u_vec(3)&
          + gd(3,1) * u_vec(3) * u_vec(1) + gd(3,2) * u_vec(3) * u_vec(2)&
          + gd(3,3) * u_vec(3) * u_vec(3)
  u_vec = u_vec / sqrt(dotp1)

  dotp1 =   gd(1,1) * u_vec(1) * v_vec(1) + gd(1,2) * u_vec(1) * v_vec(2)&
          + gd(1,3) * u_vec(1) * v_vec(3) + gd(2,1) * u_vec(2) * v_vec(1)&
          + gd(2,2) * u_vec(2) * v_vec(2) + gd(2,3) * u_vec(2) * v_vec(3)&
          + gd(3,1) * u_vec(3) * v_vec(1) + gd(3,2) * u_vec(3) * v_vec(2)&
          + gd(3,3) * u_vec(3) * v_vec(3)
  v_vec = v_vec - dotp1 * u_vec

  dotp1 =   gd(1,1) * v_vec(1) * v_vec(1) + gd(1,2) * v_vec(1) * v_vec(2)&
          + gd(1,3) * v_vec(1) * v_vec(3) + gd(2,1) * v_vec(2) * v_vec(1)&
          + gd(2,2) * v_vec(2) * v_vec(2) + gd(2,3) * v_vec(2) * v_vec(3)&
          + gd(3,1) * v_vec(3) * v_vec(1) + gd(3,2) * v_vec(3) * v_vec(2)&
          + gd(3,3) * v_vec(3) * v_vec(3)
  v_vec = v_vec / sqrt(dotp1)

  dotp1 =   gd(1,1) * u_vec(1) * w_vec(1) + gd(1,2) * u_vec(1) * w_vec(2)&
          + gd(1,3) * u_vec(1) * w_vec(3) + gd(2,1) * u_vec(2) * w_vec(1)&
          + gd(2,2) * u_vec(2) * w_vec(2) + gd(2,3) * u_vec(2) * w_vec(3)&
          + gd(3,1) * u_vec(3) * w_vec(1) + gd(3,2) * u_vec(3) * w_vec(2)&
          + gd(3,3) * u_vec(3) * w_vec(3)

  dotp2 =   gd(1,1) * v_vec(1) * w_vec(1) + gd(1,2) * v_vec(1) * w_vec(2)&
          + gd(1,3) * v_vec(1) * w_vec(3) + gd(2,1) * v_vec(2) * w_vec(1)&
          + gd(2,2) * v_vec(2) * w_vec(2) + gd(2,3) * v_vec(2) * w_vec(3)&
          + gd(3,1) * v_vec(3) * w_vec(1) + gd(3,2) * v_vec(3) * w_vec(2)&
          + gd(3,3) * v_vec(3) * w_vec(3)
  w_vec = w_vec - dotp1 * u_vec - dotp2 * v_vec

  dotp1 =   gd(1,1) * w_vec(1) * w_vec(1) + gd(1,2) * w_vec(1) * w_vec(2)&
          + gd(1,3) * w_vec(1) * w_vec(3) + gd(2,1) * w_vec(2) * w_vec(1)&
          + gd(2,2) * w_vec(2) * w_vec(2) + gd(2,3) * w_vec(2) * w_vec(3)&
          + gd(3,1) * w_vec(3) * w_vec(1) + gd(3,2) * w_vec(3) * w_vec(2)&
          + gd(3,3) * w_vec(3) * w_vec(3)
  w_vec = w_vec / sqrt(dotp1)

  ud_vec = matmul( gd, u_vec )
  !------------------------------------------


  !------------ EB part of Weyl -------------
  elec = ri
  mag  = 0
  do a = 1, 3
    do b = 1, 3
      do m = 1, 3
        do n = 1, 3
          elec(a,b) = elec(a,b) + gu(m,n) * (kd(a,b) * kd(m,n)           &
                                             - kd(a,m) * kd(b,n))
          do p = 1, 3
            mag(a,b) = mag(a,b) + gd(b,p) * eps_lc_u(p,m,n) * cd1_kd(n,a,m)   &
                                + gd(a,p) * eps_lc_u(p,m,n) * cd1_kd(n,b,m)
          end do
        end do
      end do
    end do
  end do
  mag = 0.5 * mag

  ! construct tracefree part of the electric and magnetic parts of the Weyl tensor
  electr = 0.0
  magtr  = 0.0
  do m = 1, 3
    do n = 1, 3
      electr = electr + gu(m,n) * elec(m,n)
      magtr  = magtr  + gu(m,n) * mag(m,n)
    end do
  end do

  elec = elec - gd * electr / 3.0
  mag  = mag  - gd * magtr  / 3.0

  !------------------------------------------


  !------------ Psi4 ------------------------
  !psi4re(i,j,k) = 0
  !psi4im(i,j,k) = 0
  do m = 1, 3
    do n = 1, 3
      !psi4re(i,j,k) = psi4re(i,j,k)                                      &
      psi4(REAL_PART) = psi4(REAL_PART)                                  &
                      + 0.5 * (   elec(m,n) * (   v_vec(m) * v_vec(n)    &
                                                - w_vec(m) * w_vec(n) )  &
                                - mag(m,n)  * (   v_vec(m) * w_vec(n)    &
                                                + w_vec(m) * v_vec(n) ) )
      !psi4im(i,j,k) = psi4im(i,j,k)                                      &
      psi4(IMAGINARY_PART) = psi4(IMAGINARY_PART)                        &
                      + 0.5 * (  -elec(m,n) * (   v_vec(m) * w_vec(n)    &
                                                + w_vec(m) * v_vec(n) )  &
                                + mag(m,n)  * (   w_vec(m) * w_vec(n)    &
                                                - v_vec(m) * v_vec(n) ) )
    end do
  end do
  !------------------------------------------

end subroutine CalcPsi4
