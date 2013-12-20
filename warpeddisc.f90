!
! Routines to solve diffusion equation for warped accretion discs
! in the manner described by Pringle (1992) and Ogilvie (1999)
!
! Many aspects of the code from an original code
! by Jim Pringle, with modifications by Giuseppe Lodato
!
! This version written by Daniel Price Dec. 2013
!
module warped_disc
 implicit none
 integer, parameter :: db = kind(0.d0)

contains
!
! set up arrays containing radial grid
!
pure subroutine set_grid(R_in,R_out,r,dr,rsqrt,r32,nr2)
 integer, intent(in) :: nr2
 real(kind=db), intent(in) :: R_in,R_out
 real(kind=db), dimension(nr2), intent(out) :: r,dr,rsqrt,r32
 real(kind=db) :: dri
 integer :: i

 dri   = (R_out - R_in)/real(nr2 - 2)
 do i=2,nr2
    r(i) = R_in + (i-2)*dri
 enddo
 do i=3,nr2-1
    dr(i) = r(i+1) - r(i-1)
 enddo
 do i=2,nr2
    rsqrt(i) = sqrt(r(i))
    r32(i)   = r(i)**1.5
 enddo

end subroutine set_grid
!
! set up initial conditions for surface density (sigma)
! and the unit vector of angular momentum (unitl)
!
pure subroutine initial_conds(r,p_index,r_in,ampl,sigma,unitl,bigL,nr,nr2)
 integer, intent(in) :: nr,nr2
 real(kind=db), dimension(nr2),   intent(in)  :: r
 real(kind=db),                   intent(in)  :: p_index,r_in,ampl
 real(kind=db), dimension(nr),    intent(out) :: sigma
 real(kind=db), dimension(3,nr),  intent(out) :: bigL
 real(kind=db), dimension(3,nr2), intent(out) :: unitl
 real(kind=db) :: ri,r1,r2,rwarp,sini
 integer :: i
 real(kind=db), parameter :: pi = 3.14159265358979d0

 do i=1,nr
    ri = r(2*i)
    sigma(i) = ri**(-p_index)*(1. - sqrt(r_in/ri))
    
    r1 = 3.5
    r2 = 6.5
    rwarp = 0.5*(r1 + r2)
    
    if (ri < r1) then
       sini = 0.
    elseif (ri < r2) then
       sini = ampl*0.5*(1. + sin(pi/(r2 - r1)*(ri - rwarp)))
    else
       sini = ampl
    endif
    unitl(1,2*i) = sini
    unitl(2,2*i) = 0.
    unitl(3,2*i) = sqrt(1. - sini**2)
    
    bigL(:,i) = sigma(i)*sqrt(ri)*unitl(:,2*i)
 enddo

end subroutine initial_conds
!
! copy values between cells in staggered grid
!
pure subroutine update_variables(r,dr,bigL,sigma,psi,unitl,dldr,nr,nr2)
 integer, intent(in) :: nr,nr2
 real(kind=db), dimension(nr2),   intent(in)  :: r,dr
 real(kind=db), dimension(3,nr),  intent(in)  :: bigL
 real(kind=db), dimension(nr),    intent(out) :: sigma,psi
 real(kind=db), dimension(3,nr2), intent(out) :: unitl,dldr
 real(kind=db) :: ri,bigLsqi,unitlmag
 integer :: i,k
!
! update unit vector l and surface density Sigma
! from the angular momentum vector \vec{L} = Sigma*sqrt(GMr)*\vec{l}
!
 do i=1,nr
    ri = r(2*i)
    bigLsqi  = dot_product(bigL(:,i),bigL(:,i))
    sigma(i) = sqrt(bigLsqi/ri)
    if (bigLsqi > tiny(bigLsqi)) then
       unitl(:,2*i) = bigL(:,i)/sqrt(bigLsqi)    
    else
       unitl(:,2*i) = 0.
    endif
 enddo
!
! compute the unit vector at intermediate grid points
! by averaging between full grid points
!
 do i=2,nr
    k = 2*i - 1
    unitl(:,k) = 0.5*(unitl(:,k-1) + unitl(:,k+1))
    unitlmag = sqrt(dot_product(unitl(:,k),unitl(:,k)))
    if (unitlmag > tiny(unitlmag)) then
       unitl(:,k) = unitl(:,k)/unitlmag
    endif
 enddo
!
! compute the radial derivative of tilt vector at every grid and mid-grid point
!
 do i=3,2*nr-1
    dldr(:,i) = (unitl(:,i+1) - unitl(:,i-1))/dr(i)
 enddo
! set to zero for boundaries
 dldr(:,2) = 0.
 dldr(:,nr2) = 0.
 dldr(:,3) = 0.
 dldr(:,nr2-1) = 0.
!
! compute the square of the unit tilt vector at full grid points
!
 do i=1,nr
    ri = r(2*i)
    psi(i) = ri*sqrt(dot_product(dldr(:,2*i),dldr(:,2*i)))
 enddo

end subroutine update_variables
!
! compute the diffusion/precession coefficients from Ogilvie (1999) theory
!
subroutine get_coeffs(alpha,psi,nu1,nu2,nu3,nr)
 integer, intent(in) :: nr
 real(kind=db), intent(in) :: alpha
 real(kind=db), dimension(nr), intent(in)  :: psi
 real(kind=db), dimension(nr), intent(out) :: nu1,nu2,nu3
 real(kind=db) :: alphab,alpha1,alpha2,alpha3,psii
 integer :: i
!
! compute the nu1, nu2 and nu3 diffusion/precession coefficients
! 
 !$omp parallel do default(none) schedule(dynamic) &
 !$omp private(i,alphab,alpha1,alpha2,alpha3,psii) &
 !$omp shared(nr,psi,nu1,nu2,nu3,alpha) 
 do i=1,nr
    alphab = 5./3.*alpha
    psii = max(psi(i),1.d-9)
    if (psii > 0.1) then
       call get_alpha_gordon(alpha,alphab,psii,alpha1,alpha2,alpha3)
    else
       ! Eq. 136 in Ogilvie 1999
       alpha1 = alpha - 2./3.*(1. - 17.*alpha**2 + 21.*alpha**4)/(4.*alpha*(4. + alpha**2))*psii**2
       ! Eq. 145 in Ogilvie 1999 (see Eq. 8 in Lodato & Price (2010))
       alpha2 = 1./(2.*alpha)*(4.*(1. + 7.*alpha**2))/(4. + alpha**2)
       ! Eq. 12 in Ogilvie & Dubus (2001)
       alpha3 = 1.5*(1. - 2.*alpha**2)/(4. + alpha**2)
    endif
    nu1(i) = alpha1/alpha
    nu2(i) = alpha2/alpha
    nu3(i) = alpha3/alpha
 enddo
 !$omp end parallel do

end subroutine get_coeffs

pure subroutine get_vel(r,nu1,nu2,dldr,vel2,nr,nr2)
 integer, intent(in) :: nr,nr2
 real(kind=db), dimension(nr2),   intent(in) :: r
 real(kind=db), dimension(nr),    intent(in) :: nu1,nu2
 real(kind=db), dimension(3,nr2), intent(in) :: dldr
 real(kind=db), dimension(nr),    intent(out) :: vel2
 integer :: i
 real(kind=db) :: vel1i
!
! compute the advective velocity
!
 do i=1,nr
    vel1i   = -nu2(i)*(r(2*i)**2)*dot_product(dldr(:,2*i),dldr(:,2*i))
    vel2(i) = 1.5*nu1(i) + vel1i
 enddo

end subroutine get_vel
!
! return timestep constrained by Courant condition 
! from advection and diffusion terms
!
pure subroutine get_timestep(dt,courant_fac,dr,sigma,vel2,nu1,nu2,nu3,nr,nr2)
 integer, intent(in) :: nr,nr2
 real(kind=db), intent(out) :: dt
 real(kind=db), intent(in) :: courant_fac
 real(kind=db), dimension(nr),  intent(in) :: sigma,vel2,nu1,nu2,nu3
 real(kind=db), dimension(nr2), intent(in) :: dr
 real(kind=db) :: dti,absv,dri
 integer :: i
 
 dt = huge(dt)
 over_grid: do i=2,nr-1
!
!   skip if no mass
!
    if (sigma(i) < tiny(sigma)) cycle over_grid
!
!   timestep from advection term
!
    dri  = dr(2*i)
    absv = abs(vel2(i))
    if (absv > tiny(vel2)) then
       dti = dri/absv
       if (dti < dt) dt = dti
    endif
!
!   timestep from diffusion term
!
    dti = dri**2/(nu1(i) + nu2(i) + nu3(i))
    if (dti < dt) dt = dti

 enddo over_grid
 dt = courant_fac*dt
 
end subroutine get_timestep
!
! advance one timestep: this is the main routine specifying the equations
!
pure subroutine advance(dt,r,rsqrt,r32,dr,sigma,nu1,nu2,nu3,vel2,dldr,unitl,bigL,nr,nr2)
 integer, intent(in) :: nr,nr2
 real(kind=db), intent(in) :: dt
 real(kind=db), dimension(nr2), intent(in) :: r,rsqrt,r32,dr
 real(kind=db), dimension(nr),  intent(in) :: sigma,nu1,nu2,nu3,vel2
 real(kind=db), dimension(3,nr2), intent(in)    :: unitl,dldr
 real(kind=db), dimension(3,nr),  intent(inout) :: bigL
 real(kind=db), dimension(3,nr) :: bigLnew
 integer :: i,k
 real(kind=db), dimension(3) :: term1,term2,term3,term4,lcrossdldrp1,lcrossdldrm1
 real(kind=db) :: ri,dri,c1,c2,c3,c4

 do i=2,nr-1
    k = 2*i
    ri  = r(k)
    dri = dr(k)
!
! diffusion terms
!
    if (i.eq.2) then
       term1(:) = r(k+1)*unitl(:,k+1)*(nu1(i+1)*sigma(i+1)*rsqrt(k+2) - &
                                       nu1(i)*sigma(i)*rsqrt(k))/dr(k+1) &
! these extra bits give a zero viscous torque boundary
                 -r(k-1)*unitl(:,k-1)*(nu1(i)*sigma(i)*rsqrt(k) - &
                                       nu1(i-1)*sigma(i-1)*rsqrt(k-2))/dr(k-1)
    elseif (i.eq.nr-1) then
       term1(:) = -r(k-1)*unitl(:,k-1)*(nu1(i)*sigma(i)*rsqrt(k) - &
                                        nu1(i-1)*sigma(i-1)*rsqrt(k-2))/dr(k-1)
    else
       term1(:) = r(k+1)*unitl(:,k+1)*(nu1(i+1)*sigma(i+1)*rsqrt(k+2) - &
                                       nu1(i)*sigma(i)*rsqrt(k))/dr(k+1) &
                 -r(k-1)*unitl(:,k-1)*(nu1(i)*sigma(i)*rsqrt(k) - &
                                       nu1(i-1)*sigma(i-1)*rsqrt(k-2))/dr(k-1)
    endif
    term1(:) = term1(:)*3./ri/dri
!
! advection term
!
    if (vel2(i) < 0.) then
       term2(:) = -(vel2(i+1)*bigL(:,i+1) - vel2(i)*bigL(:,i))/dr(k+1)    
    else
       term2(:) = -(vel2(i)*bigL(:,i) - vel2(i-1)*bigL(:,i-1))/dr(k-1)
    endif
    term2(:) = term2(:)/ri
!
! diffusion term
!
    c1 = 0.5*(nu2(i+1)*sigma(i+1) + nu2(i)*sigma(i))
    c2 = 0.5*(nu2(i)*sigma(i) + nu2(i-1)*sigma(i-1))
    term3(:) = c1*r32(k+1)*dldr(:,k+1) - c2*r32(k-1)*dldr(:,k-1)
    term3(:) = term3(:)*0.5/ri/dri
!
! precession term
!
    c3 = 0.5*(nu3(i+1)*sigma(i+1) + nu3(i)*sigma(i))
    c4 = 0.5*(nu3(i)*sigma(i) + nu3(i-1)*sigma(i-1))
    call cross_product3D(unitl(:,k+1),dldr(:,k+1),lcrossdldrp1)
    call cross_product3D(unitl(:,k-1),dldr(:,k-1),lcrossdldrm1)
    term4(:) = c3*r32(k+1)*lcrossdldrp1(:) - c4*r32(k-1)*lcrossdldrm1(:)
    term4(:) = term4(:)/ri/dri
!
! update the angular momentum vector
!   
    bigLnew(:,i) = bigL(:,i) + dt*(term1(:) + term2(:) + term3(:) + term4(:))
 enddo

 bigL = bigLnew
!
! boundary conditions: bigL = 0 at i=1 and i=nr
!
 bigL(:,1)  = 0.
 bigL(:,nr) = 0.

end subroutine advance
!
! take the cross product of two vectors
!
pure subroutine cross_product3D(veca,vecb,vecc)
 real(kind=db), dimension(3), intent(in) :: veca,vecb
 real(kind=db), dimension(3), intent(out) :: vecc
 
 vecc(1) = veca(2)*vecb(3) - veca(3)*vecb(2)
 vecc(2) = veca(3)*vecb(1) - veca(1)*vecb(3)
 vecc(3) = veca(1)*vecb(2) - veca(2)*vecb(1)

end subroutine cross_product3D
!
! decide whether or not to write a file this timestep or not
!
subroutine decide_to_write_output(time,fac,dt_out,dt,do_output)
 real(kind=db), intent(in)    :: time,fac,dt_out
 real(kind=db), intent(inout) :: dt
 logical,       intent(out)   :: do_output
 real(kind=db) :: tvisc,tviscnext,t_next_output
 integer :: next_output
 
 do_output = .false.

 tvisc = time/fac
 tviscnext = (time + dt)/fac
 next_output = int(tvisc/dt_out) + 1
 t_next_output = next_output*dt_out

 if (tviscnext > t_next_output) then
    do_output = .true.
    dt = (t_next_output - tvisc)*fac + epsilon(0.d0)
    if (dt < 0.) print*,' error: dt = ',dt,' in decide_to_write_output'
 endif

end subroutine decide_to_write_output
!
! output results to file
!
subroutine write_output(nout,time,r,sig0,sigma,unitl,psi,nr,nr2)
 integer, intent(in) :: nr,nr2
 integer, intent(inout) :: nout
 real(kind=db), intent(in) :: time,sig0
 real(kind=db), dimension(nr2),   intent(in) :: r
 real(kind=db), dimension(nr),    intent(in) :: sigma,psi
 real(kind=db), dimension(3,nr2), intent(in) :: unitl
 integer :: i,lu,ierr
 character(len=8) :: filename
 
 write(filename,"(a,i3.3)") 'dump',nout
 open(newunit=lu,file=filename,status='replace',iostat=ierr)
 if (ierr /= 0) then
    print*,' error writing to '//trim(filename)
    return
 endif
 
 print "(a,f12.5,a)",' t = ',time,' writing to '//trim(filename)
 write(lu,"(6(a10,1x))") 'r','sigma','lx','ly','lz','psi'
 write(lu,*) time
 do i=1,nr
    write(lu,"(6(es13.6,1x))") r(2*i),sig0*sigma(i),unitl(1,2*i),unitl(2,2*i),unitl(3,2*i),psi(i)
 enddo
 close(unit=lu)
 nout = nout + 1
 
end subroutine write_output
!
! main driver routine
!
subroutine evolve(t_end,dt_out)
 real(kind=db), intent(in) :: t_end,dt_out
 integer, parameter :: nr = 201
 integer, parameter :: nr2 = 2*nr
 real(kind=db), dimension(nr2) :: r,dr,rsqrt,r32
 real(kind=db), dimension(nr)  :: sigma,nu1,nu2,nu3,vel2,psi
 real(kind=db), dimension(3,nr)  :: bigL
 real(kind=db), dimension(3,nr2) :: unitl,dldr
 real(kind=db) :: R_in,R_out,warp_ampl,p_index,courant_fac,alpha
 real(kind=db) :: time,tvisc,sig0,dt,factor,hoverr
 logical :: do_output
 integer :: nout
 
 R_in = 0.5d0
 R_out = 10.5d0
 warp_ampl = 0.05d0
 p_index = 1.5d0
 alpha = 0.05d0
 courant_fac = 0.1  ! timestep constraint
 sig0   = 2.679581e-6/0.305 ! normalisation of surface density, for printing only
 hoverr = 0.02 ! H/R at r=1: ONLY used to convert the time, see below
 factor = alpha*hoverr**2
 nout = 0

 call set_grid(R_in,R_out,r,dr,rsqrt,r32,nr2)
 call initial_conds(r,p_index,R_in,warp_ampl,sigma,unitl,bigL,nr,nr2)
 call update_variables(r,dr,bigL,sigma,psi,unitl,dldr,nr,nr2)
 call get_coeffs(alpha,psi,nu1,nu2,nu3,nr)
 call get_vel(r,nu1,nu2,dldr,vel2,nr,nr2)

 time = 0.
 tvisc = 0.
 call write_output(nout,tvisc,r,sig0,sigma,unitl,psi,nr,nr2)

 do while (tvisc < t_end)
    call get_timestep(dt,courant_fac,dr,sigma,vel2,nu1,nu2,nu3,nr,nr2)
    call decide_to_write_output(time,factor,dt_out,dt,do_output)
    call advance(dt,r,rsqrt,r32,dr,sigma,nu1,nu2,nu3,vel2,dldr,unitl,bigL,nr,nr2)

    call update_variables(r,dr,bigL,sigma,psi,unitl,dldr,nr,nr2)
    call get_coeffs(alpha,psi,nu1,nu2,nu3,nr)
    call get_vel(r,nu1,nu2,dldr,vel2,nr,nr2)

    time = time + dt
    tvisc = time/factor

    if (do_output) call write_output(nout,tvisc,r,sig0,sigma,unitl,psi,nr,nr2)
 enddo

end subroutine evolve

end module warped_disc
