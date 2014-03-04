!
! Warp diffusion code (these are the main driving routines)
! Daniel Price
!
module rungrid
 implicit none
contains

subroutine run_grid
 use warped_disc, only:evolve,db
 real(kind=db) :: time,dtout,psi,alpha
 integer, parameter :: nalpha = 50, npsi = 15
 real(kind=db) :: alphamin,dalpha,psimin,dpsi
 character(len=7) :: psidir
 character(len=17) :: alphadir
 character(len=30) :: dir
 logical :: fail,iexist
 integer :: j,i,nr

 alphamin = 0.01d0
 dalpha   = 0.01d0
 psimin = 0.1d0
 dpsi = 0.1d0
 time = 1000.
 dtout = 5.
 nr = 201

 do i=1,npsi
    psi = psimin + (i-1)*dpsi
    write(psidir,"(a,f3.1)") 'psi-',psi

    do j=1,nalpha
       alpha = alphamin + (j-1)*dalpha
       write(alphadir,"(a,f5.3)") 'warp_alphaSS',alpha
       
       dir = psidir//'/'//alphadir//'/'
       !
       !--only redo the calculation if not already done
       !
       inquire(file=trim(dir)//'angm001',exist=iexist)
       if (.not.iexist) then
          call system('mkdir -p '//trim(dir))
          print*,psi,alpha,trim(dir)
          call evolve(time,dtout,psi,alpha,trim(dir)//'angm',nr,fail)
          if (fail) print*,' FAILED'
       else
          print*,' skipping '//trim(dir)//': angm200 exists'
       endif
    enddo
 enddo
 
 end subroutine run_grid

end module rungrid


program disc
 use warped_disc, only:evolve,db
 use rungrid,     only:run_grid
 implicit none
 real(kind=db) :: time,dtout,psi,alpha
 logical :: fail,grid
 integer :: nr
 
 time = 1000.d0
 dtout = 5.d0
 psi = 0.8d0
 alpha = 0.29d0
 grid = .false.
 nr = 401
 
 if (grid) then
    call run_grid()
 else
    print*,' running psi = ',psi,' alpha = ',alpha
    call evolve(time,dtout,psi,alpha,'dump',nr,fail)
    if (fail) print*,' FAILED'
 endif

end program disc
