!
! Program to solve diffusion equation for accretion discs
! in the manner described by Pringle (1992) and Ogilvie (1999)
!
! Many aspects of the code from an original code by Jim Pringle
!
! This version by Daniel Price Dec. 2013
!
program disc
 use warped_disc, only:evolve,db
 real(kind=db) :: time,dtout
 
 time = 1000.
 dtout = 5.

 call evolve(time,dtout)

end program disc
