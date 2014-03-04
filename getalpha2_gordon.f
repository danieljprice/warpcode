      module gordon
      
      public :: get_alpha_gordon
      
      private
      
      contains

      subroutine get_alpha_gordon(alpha,alphab,psi,alpha1,alpha2,alpha3,
     &                            fail)

c Integrates the dimensionless ODEs for a nonlinear warp in a polytropic
c disc.
c Computes the dimensionless coefficients Q_1, Q_2 and Q_3.
c Compiles a table with |psi| between 0.00 and 5.00.

c Independent variables: phi (azimuthal angle)

c Dependent variables: y(1): f2-1.0d0 (representing h,
c                                      part proportional to z^2)
c                      y(2): f3 (representing v_r)
c                      y(3): f4 (representing v_theta)
c                      y(4): f5 (representing v_phi)
c                      y(5): (d/d lambda1)y(1)
c                      y(6): (d/d lambda1)y(2)
c                      y(7): (d/d lambda1)y(3)
c                      y(8): (d/d lambda1)y(4)
c                      y(9): (d/d lambda2)y(1)
c                     y(10): (d/d lambda2)y(2)
c                     y(11): (d/d lambda2)y(3)
c                     y(12): (d/d lambda2)y(4)
c                     y(13): (d/d lambda3)y(1)
c                     y(14): (d/d lambda3)y(2)
c                     y(15): (d/d lambda3)y(3)
c                     y(16): (d/d lambda3)y(4)
c                     y(17): (d/d lambda4)y(1)
c                     y(18): (d/d lambda4)y(2)
c                     y(19): (d/d lambda4)y(3)
c                     y(20): (d/d lambda4)y(4)
c                     y(21): f6 (representing azimuthal dependence of I)
c                     y(22): integral of y(21)
c                     y(23): integral for Q_1
c                     y(24): integral for Q_2
c                     y(25): integral for Q_3
c                     y(26): integral for check on Q_1 (real)
c                     y(27): integral for check on Q_1 (imaginary)
c                     y(28): integral for check on Q_2

c Parameters: mu(1): |psi|           (amplitude of warp)
c             mu(2): kappa^2/omega^2 (dimensionless epicyclic freq^2)
c             mu(3): gamma           (polytropic exponent)
c             mu(4): alpha           (shear viscosity coefficient)
c             mu(5): alpha_b         (bulk viscosity coefficient)

c Eigenvalues: lambda(1): f1(0)
c              lambda(2): f2(0)
c              lambda(3): f3(0)
c              lambda(4): f4(0)

c Boundary conditions: periodic on (0,2 pi).

c f1 and f6 are decoupled and have arbitrary normalizations
c f1(0)=f6(0)=1.  However, f6 must be renormalized to ensure
c <f6>=1, which is necessary for the evaluation of Q_1, Q_2 and Q_3.
c Hence the need for y(22).

      implicit none
      real*8, intent(in)   :: psi,alpha,alphab
      real*8, intent(out)  :: alpha1,alpha2,alpha3
      logical, intent(out) :: fail

      integer nmu
      parameter (nmu=5)
      real*8 mu(nmu)
      character*15 filename

      integer i,imax,nmvary,ndmu,indx(4)
      real*8 x1,x2,y1(28),y2(28),lambda(4),acc,mm(4),m(4,4),del(4),
     +  delta,dmu,mu0,d
      complex*16 temp
      logical converged

      parameter (imax=50,acc=1.0d-8)

      fail = .false.
      mu(1) = psi
      !print *,'Enter psi'
      !read *,mu(1)

      !write(filename,"('psi',f5.3,'.dat')") mu(1)
     
      !open (unit=12,file=filename,form='formatted')

c      mu(1)=0.005d0
      mu(2)=1.0d0
      mu(3)=1.0d0
      mu(4)=0.01d0
      mu(5)=alphab
      lambda(1)=0.0001d0
      lambda(2)=0.05d0
      lambda(3)=-0.0008d0
      lambda(4)=0.004d0
      nmvary=4
      ndmu=99
      dmu=0.01d0

      temp=dcmplx(-2.0d0*mu(4),1.0d0+(7.0d0-mu(2))*mu(4)*
     +  mu(4))/dcmplx(2.0d0*(1.0d0-mu(2))-2.0d0*mu(4)*mu(4),4.0d0*mu(4))
c      write (12,'(4d20.12)') 0.0d0,-0.5d0*(4.0d0-mu(2))*mu(4),
c     +  dreal(temp),dimag(temp)

      
      mu(nmvary) = alpha
      !mu0=mu(nmvary)

      !do j=0,ndmu

       ! mu(nmvary)=mu0+real(j)*dmu

        i=0
        converged=.false.
        do while (i.lt.imax .and. .not.converged)
        !do i=1,imax
          i = i + 1

          x1=0.0d0
          x2=8.0d0*datan(1.0d0)
          y1(1)=lambda(1)
          y1(2)=lambda(2)
          y1(3)=lambda(3)
          y1(4)=lambda(4)
          y1(5)=1.0d0
          y1(6)=0.0d0
          y1(7)=0.0d0
          y1(8)=0.0d0
          y1(9)=0.0d0
          y1(10)=1.0d0
          y1(11)=0.0d0
          y1(12)=0.0d0
          y1(13)=0.0d0
          y1(14)=0.0d0
          y1(15)=1.0d0
          y1(16)=0.0d0
          y1(17)=0.0d0
          y1(18)=0.0d0
          y1(19)=0.0d0
          y1(20)=1.0d0
          y1(21)=1.0d0
          y1(22)=0.0d0
          y1(23)=0.0d0
          y1(24)=0.0d0
          y1(25)=0.0d0
          y1(26)=0.0d0
          y1(27)=0.0d0
          y1(28)=0.0d0
          call odeint(x1,x2,y1,y2,lambda,nmu,mu,.false.,fail)
          if (fail) then
            write (*,'(a13)') 'odeint failed'
            return
          endif
          mm(1)=y2(1)-lambda(1)
          mm(2)=y2(2)-lambda(2)
          mm(3)=y2(3)-lambda(3)
          mm(4)=y2(4)-lambda(4)
          m(1,1)=y2(5)-1.0d0
          m(2,1)=y2(6)
          m(3,1)=y2(7)
          m(4,1)=y2(8)
          m(1,2)=y2(9)
          m(2,2)=y2(10)-1.0d0
          m(3,2)=y2(11)
          m(4,2)=y2(12)
          m(1,3)=y2(13)
          m(2,3)=y2(14)
          m(3,3)=y2(15)-1.0d0
          m(4,3)=y2(16)
          m(1,4)=y2(17)
          m(2,4)=y2(18)
          m(3,4)=y2(19)
          m(4,4)=y2(20)-1.0d0
          call ludcmp(m,4,4,indx,d)
          del(1)=-mm(1)
          del(2)=-mm(2)
          del(3)=-mm(3)
          del(4)=-mm(4)
          call lubksb(m,4,4,indx,del)
          delta=dmax1(dabs(del(1)),dabs(del(2)),dabs(del(3)),
     +      dabs(del(4)))
          lambda(1)=lambda(1)+del(1)
          lambda(2)=lambda(2)+del(2)
          lambda(3)=lambda(3)+del(3)
          lambda(4)=lambda(4)+del(4)
          if (delta.lt.acc) then
!            write (12,'(4d20.12)') mu(nmvary),y2(23)/y2(22),
!     +        y2(24)/y2(22),y2(25)/y2(22)
!            goto 1
             alpha1 = -2./3.*y2(23)/y2(22)
             alpha2 = 2.*y2(24)/y2(22)
             alpha3 = y2(25)/y2(22)
             converged = .true.
          endif
        enddo
        if (.not.converged) then
           write (*,*) 'Failed to converge: alpha=',alpha,' psi=',psi
           fail = .true.
           return
        endif

! 1    enddo

!      close (unit=12)

      end


      SUBROUTINE odeint(x1,x2,y1,y2,lambda,nmu,mu,dispef,fail)
      implicit none
      LOGICAL dispef,fail
      INTEGER nmu
      REAL*8 x1,x2,y1(28),y2(28),lambda(4),mu(nmu)
C
      LOGICAL last,good,revers
      INTEGER stepmx,i,j
      REAL*8 eps,hmax,x,h,hnext,y(28)
      PARAMETER (stepmx=10000,eps=1.0D-8)
      IF (x2.LT.x1) THEN
        revers=.TRUE.
      ELSE
        revers=.FALSE.
      ENDIF
      x=x1
      h=(x2-x1)*0.001D0
      hmax=h*10.0D0
      fail=.FALSE.
      last=.FALSE.
      DO i=1,28
        y(i)=y1(i)
      ENDDO
      IF (dispef) THEN
        WRITE (*,'(5D20.12)') x,y(1),y(2),y(3),y(4)
      ENDIF
      DO i=1,stepmx
        IF ((((x+h).GT.x2).AND.(.NOT.revers)).OR.
     +    (((x+h).LT.x2).AND.revers)) THEN
          h=x2-x
          last=.TRUE.
        ENDIF
        CALL qstep(x,y,lambda,nmu,mu,h,hnext,eps,hmax,good,
     +    fail)
        IF (fail) RETURN
        IF (last.AND.good) THEN
          DO j=1,28
            y2(j)=y(j)
          ENDDO
          IF (dispef) THEN
            WRITE (*,'(5D20.12)') x,y(1),y(2),y(3),y(4)
          ENDIF
          RETURN
        ENDIF
        last=.FALSE.
        h=hnext
        IF (dispef.AND.good) THEN
          WRITE (*,'(5D20.12)') x,y(1),y(2),y(3),y(4)
        ENDIF
      ENDDO
      WRITE (*,'(A40)') 'odeint: maximum number of steps exceeded'
      fail=.TRUE.
      RETURN
      END

      SUBROUTINE qstep(x,y,lambda,nmu,mu,htry,hnext,eps,hmax,good,fail)
      implicit none
      LOGICAL good,fail
      INTEGER nmu,ierr
      REAL*8 x,y(28),lambda(4),mu(nmu),htry,hnext,eps,hmax
C
      INTEGER i
      REAL*8 safety,errmax,errmin,h,error,xnext,ynext(28),yerror(28)
      PARAMETER (safety=0.9D0,errmax=1.845281D3,errmin=1.889568D-4)
      h=htry
      good=.TRUE.
      fail=.false.
 1    CALL ckstep(x,y,lambda,nmu,mu,h,ynext,yerror,ierr)
      if (ierr /= 0) then
         fail = .true.
         print*,' got ckstep failure',ierr
         return
      endif
      error=(DMAX1(DABS(yerror(1)),DABS(yerror(2)),DABS(yerror(3)),
     +  DABS(yerror(4)))/
     +  DMAX1(DABS(y(1)),DABS(y(2)),DABS(y(3)),DABS(y(4))))/eps
      IF (error.GT.1.0D0) THEN
        IF (error.GT.errmax) THEN
          h=h*0.2D0
        ELSE
          h=safety*h*(error**(-0.2D0))
        ENDIF
        good=.FALSE.
        xnext=x+h
        IF (xnext.EQ.x) THEN
          WRITE(*,'(A25)') 'qstep: stepsize underflow'
          fail=.TRUE.
          RETURN
        ELSE
          GOTO 1
        ENDIF
      ELSE
        IF (error.LT.errmin) THEN
          hnext=5.0D0*h
        ELSE
          hnext=safety*h*(error**(-0.2D0))
        ENDIF
        IF (DABS(hnext).GT.DABS(hmax)) hnext=hmax
        x=x+h
        DO i=1,28
          y(i)=ynext(i)
        ENDDO
        RETURN
      ENDIF
      END

      SUBROUTINE ckstep(x,y,lambda,nmu,mu,h,ynext,yerror,ierr)
      implicit none
      INTEGER nmu
      REAL*8 x,y(28),lambda(4),mu(nmu),h,ynext(28),yerror(28)
      integer, intent(out) :: ierr
C
      INTEGER i
      REAL*8 a2,a3,a4,a5,a6,b21,b31,b32,b41,b42,b43,b51,b52,b53,b54,b61,
     +  b62,b63,b64,b65,c1,c3,c4,c6,d1,d3,d4,d5,d6,
     +  ytemp(28),f1(28),f2(28),f3(28),f4(28),f5(28),f6(28)
      PARAMETER (a2=0.2D0,a3=0.3D0,a4=0.6D0,a5=1.0D0,a6=0.875D0,
     +  b21=0.2D0,b31=3.0D0/40.0D0,b32=9.0D0/40.0D0,b41=0.3D0,
     +  b42=-0.9D0,b43=1.2D0,b51=-11.0D0/54.0D0,b52=2.5D0,
     +  b53=-70.0D0/27.0D0,b54=35.0D0/27.0D0,b61=1631.0D0/55296.0D0,
     +  b62=175.0D0/512.0D0,b63=575.0D0/13824.0D0,
     +  b64=44275.0D0/110592.0D0,b65=253.0D0/4096.0D0,c1=37.0D0/378.0D0,
     +  c3=250.0D0/621.0D0,c4=125.0D0/594.0D0,c6=512.0D0/1771.0D0,
     +  d1=c1-2825.0D0/27648.0D0,d3=c3-18575.0D0/48384.0D0,
     +  d4=c4-13525.0D0/55296.0D0,d5=-277.0D0/14336.0D0,d6=c6-0.25D0)
      CALL derivs(x,y,lambda,nmu,mu,f1,ierr)
      if (ierr /= 0) return
      DO i=1,28
        ytemp(i)=y(i)+h*b21*f1(i)
      ENDDO
      CALL derivs(x+a2*h,ytemp,lambda,nmu,mu,f2,ierr)
      if (ierr /= 0) return
      DO i=1,28
        ytemp(i)=y(i)+h*(b31*f1(i)+b32*f2(i))
      ENDDO
      CALL derivs(x+a3*h,ytemp,lambda,nmu,mu,f3,ierr)
      if (ierr /= 0) return
      DO i=1,28
        ytemp(i)=y(i)+h*(b41*f1(i)+b42*f2(i)+b43*f3(i))
      ENDDO
      CALL derivs(x+a4*h,ytemp,lambda,nmu,mu,f4,ierr)
      if (ierr /= 0) return
      DO i=1,28
        ytemp(i)=y(i)+h*(b51*f1(i)+b52*f2(i)+b53*f3(i)+b54*f4(i))
      ENDDO
      CALL derivs(x+a5*h,ytemp,lambda,nmu,mu,f5,ierr)
      if (ierr /= 0) return
      DO i=1,28
        ytemp(i)=y(i)+h*(b61*f1(i)+b62*f2(i)+b63*f3(i)+b64*f4(i)+
     +    b65*f5(i))
      ENDDO
      CALL derivs(x+a6*h,ytemp,lambda,nmu,mu,f6,ierr)
      if (ierr /= 0) return
      DO i=1,28
        ynext(i)=y(i)+h*(c1*f1(i)+c3*f3(i)+c4*f4(i)+c6*f6(i))
        yerror(i)=h*(d1*f1(i)+d3*f3(i)+d4*f4(i)+d5*f5(i)+d6*f6(i))
      ENDDO
      RETURN
      END

      SUBROUTINE derivs(x,y,lambda,nmu,mu,f,ierr)
      implicit none
      INTEGER nmu
      REAL*8 x,y(28),lambda(4),mu(nmu),f(28)
      integer, intent(out) :: ierr
C
      ierr = 0
      if ((1.0d0+y(1)).lt.0.0d0) then
        write (*,'(a23)') 'Warning: disc explodes!'
        ierr = 1
        return
      endif
      f(1)=(mu(3)+1.0d0)*y(3)*(1.0d0+y(1))
      f(2)=y(3)*y(2)+2.0d0*y(4)+(1.0d0+(mu(5)+(mu(4)/3.0d0))*y(3))*
     +  (1.0d0+y(1))*mu(1)*dcos(x)-
     +  mu(4)*(1.0d0+y(1))*y(2)*(1.0d0+mu(1)*mu(1)*dcos(x)*dcos(x))-
     +  mu(4)*(1.0d0+y(1))*mu(1)*dsin(x)
      f(3)=-f(2)*mu(1)*dcos(x)+2.0d0*y(2)*mu(1)*dsin(x)+
     +  y(3)*(y(3)+y(2)*mu(1)*dcos(x))-y(1)-
     +  ((mu(5)+(mu(4)/3.0d0))*y(3))*(1.0d0+y(1))-
     +  mu(4)*(1.0d0+y(1))*(y(3)+y(2)*mu(1)*dcos(x))*
     +  (1.0d0+mu(1)*mu(1)*dcos(x)*dcos(x))+
     +  mu(4)*(1.0d0+y(1))*mu(1)*mu(1)*dcos(x)*dsin(x)
      f(4)=y(3)*y(4)-0.5d0*mu(2)*y(2)-
     +  mu(4)*(1.0d0+y(1))*y(4)*(1.0d0+mu(1)*mu(1)*dcos(x)*dcos(x))+
     +  0.5d0*(4.0d0-mu(2))*mu(4)*(1.0d0+y(1))*mu(1)*dcos(x)
      f(5)=(mu(3)+1.0d0)*(y(7)*(1.0d0+y(1))+y(3)*y(5))
      f(6)=y(7)*y(2)+y(3)*y(6)+2.0d0*y(8)+
     +  (1.0d0+(mu(5)+(mu(4)/3.0d0))*y(3))*y(5)*mu(1)*dcos(x)+
     +  ((mu(5)+(mu(4)/3.0d0))*y(7))*(1.0d0+y(1))*mu(1)*dcos(x)-
     +  mu(4)*(y(5)*y(2)+(1.0d0+y(1))*y(6))*(1.0d0+mu(1)*mu(1)*dcos(x)*
     +  dcos(x))-mu(4)*y(5)*mu(1)*dsin(x)
      f(7)=-f(6)*mu(1)*dcos(x)+2.0d0*y(6)*mu(1)*dsin(x)+
     +  y(7)*(y(3)+y(2)*mu(1)*dcos(x))+y(3)*(y(7)+y(6)*mu(1)*dcos(x))-
     +  (1.0d0+(mu(5)+(mu(4)/3.0d0))*y(3))*y(5)-
     +  ((mu(5)+(mu(4)/3.0d0))*y(7))*(1.0d0+y(1))-
     +  mu(4)*y(5)*(y(3)+y(2)*mu(1)*dcos(x))*
     +  (1.0d0+mu(1)*mu(1)*dcos(x)*dcos(x))-
     +  mu(4)*(1.0d0+y(1))*(y(7)+y(6)*mu(1)*dcos(x))*
     +  (1.0d0+mu(1)*mu(1)*dcos(x)*dcos(x))+
     +  mu(4)*y(5)*mu(1)*mu(1)*dcos(x)*dsin(x)
      f(8)=y(7)*y(4)+y(3)*y(8)-0.5d0*mu(2)*y(6)-
     +  mu(4)*(y(5)*y(4)+(1.0d0+y(1))*y(8))*(1.0d0+mu(1)*mu(1)*dcos(x)*
     +  dcos(x))+0.5d0*(4.0d0-mu(2))*mu(4)*y(5)*mu(1)*dcos(x)
      f(9)=(mu(3)+1.0d0)*(y(11)*(1.0d0+y(1))+y(3)*y(9))
      f(10)=y(11)*y(2)+y(3)*y(10)+2.0d0*y(12)+
     +  (1.0d0+(mu(5)+(mu(4)/3.0d0))*y(3))*y(9)*mu(1)*dcos(x)+
     +  ((mu(5)+(mu(4)/3.0d0))*y(11))*(1.0d0+y(1))*mu(1)*dcos(x)-
     +  mu(4)*(y(9)*y(2)+(1.0d0+y(1))*y(10))*
     +  (1.0d0+mu(1)*mu(1)*dcos(x)*dcos(x))-
     +  mu(4)*y(9)*mu(1)*dsin(x)
      f(11)=-f(10)*mu(1)*dcos(x)+2.0d0*y(10)*mu(1)*dsin(x)+
     +  y(11)*(y(3)+y(2)*mu(1)*dcos(x))+y(3)*
     +  (y(11)+y(10)*mu(1)*dcos(x))-
     +  (1.0d0+(mu(5)+(mu(4)/3.0d0))*y(3))*y(9)-
     +  ((mu(5)+(mu(4)/3.0d0))*y(11))*(1.0d0+y(1))-
     +  mu(4)*y(9)*(y(3)+y(2)*mu(1)*dcos(x))*
     +  (1.0d0+mu(1)*mu(1)*dcos(x)*dcos(x))-
     +  mu(4)*(1.0d0+y(1))*(y(11)+y(10)*mu(1)*dcos(x))*
     +  (1.0d0+mu(1)*mu(1)*dcos(x)*dcos(x))+
     +  mu(4)*y(9)*mu(1)*mu(1)*dcos(x)*dsin(x)
      f(12)=y(11)*y(4)+y(3)*y(12)-0.5d0*mu(2)*y(10)-
     +  mu(4)*(y(9)*y(4)+(1.0d0+y(1))*y(12))*
     +  (1.0d0+mu(1)*mu(1)*dcos(x)*dcos(x))+
     +  0.5d0*(4.0d0-mu(2))*mu(4)*y(9)*mu(1)*dcos(x)
      f(13)=(mu(3)+1.0d0)*(y(15)*(1.0d0+y(1))+y(3)*y(13))
      f(14)=y(15)*y(2)+y(3)*y(14)+2.0d0*y(16)+
     +  (1.0d0+(mu(5)+(mu(4)/3.0d0))*y(3))*y(13)*mu(1)*dcos(x)+
     +  ((mu(5)+(mu(4)/3.0d0))*y(15))*(1.0d0+y(1))*mu(1)*dcos(x)-
     +  mu(4)*(y(13)*y(2)+(1.0d0+y(1))*y(14))*
     +  (1.0d0+mu(1)*mu(1)*dcos(x)*dcos(x))-
     +  mu(4)*y(13)*mu(1)*dsin(x)
      f(15)=-f(14)*mu(1)*dcos(x)+2.0d0*y(14)*mu(1)*dsin(x)+
     +  y(15)*(y(3)+y(2)*mu(1)*dcos(x))+y(3)*
     +  (y(15)+y(14)*mu(1)*dcos(x))-
     +  (1.0d0+(mu(5)+(mu(4)/3.0d0))*y(3))*y(13)-
     +  ((mu(5)+(mu(4)/3.0d0))*y(15))*(1.0d0+y(1))-
     +  mu(4)*y(13)*(y(3)+y(2)*mu(1)*dcos(x))*
     +  (1.0d0+mu(1)*mu(1)*dcos(x)*dcos(x))-
     +  mu(4)*(1.0d0+y(1))*(y(15)+y(14)*mu(1)*dcos(x))*
     +  (1.0d0+mu(1)*mu(1)*dcos(x)*dcos(x))+
     +  mu(4)*y(13)*mu(1)*mu(1)*dcos(x)*dsin(x)
      f(16)=y(15)*y(4)+y(3)*y(16)-0.5d0*mu(2)*y(14)-
     +  mu(4)*(y(13)*y(4)+(1.0d0+y(1))*y(16))*
     +  (1.0d0+mu(1)*mu(1)*dcos(x)*dcos(x))+
     +  0.5d0*(4.0d0-mu(2))*mu(4)*y(13)*mu(1)*dcos(x)
      f(17)=(mu(3)+1.0d0)*(y(19)*(1.0d0+y(1))+y(3)*y(17))
      f(18)=y(19)*y(2)+y(3)*y(18)+2.0d0*y(20)+
     +  (1.0d0+(mu(5)+(mu(4)/3.0d0))*y(3))*y(17)*mu(1)*dcos(x)+
     +  ((mu(5)+(mu(4)/3.0d0))*y(19))*(1.0d0+y(1))*mu(1)*dcos(x)-
     +  mu(4)*(y(17)*y(2)+(1.0d0+y(1))*y(18))*
     +  (1.0d0+mu(1)*mu(1)*dcos(x)*dcos(x))-
     +  mu(4)*y(17)*mu(1)*dsin(x)
      f(19)=-f(18)*mu(1)*dcos(x)+2.0d0*y(18)*mu(1)*dsin(x)+
     +  y(19)*(y(3)+y(2)*mu(1)*dcos(x))+y(3)*
     +  (y(19)+y(18)*mu(1)*dcos(x))-
     +  (1.0d0+(mu(5)+(mu(4)/3.0d0))*y(3))*y(17)-
     +  ((mu(5)+(mu(4)/3.0d0))*y(19))*(1.0d0+y(1))-
     +  mu(4)*y(17)*(y(3)+y(2)*mu(1)*dcos(x))*
     +  (1.0d0+mu(1)*mu(1)*dcos(x)*dcos(x))-
     +  mu(4)*(1.0d0+y(1))*(y(19)+y(18)*mu(1)*dcos(x))*
     +  (1.0d0+mu(1)*mu(1)*dcos(x)*dcos(x))+
     +  mu(4)*y(17)*mu(1)*mu(1)*dcos(x)*dsin(x)
      f(20)=y(19)*y(4)+y(3)*y(20)-0.5d0*mu(2)*y(18)-
     +  mu(4)*(y(17)*y(4)+(1.0d0+y(1))*y(20))*
     +  (1.0d0+mu(1)*mu(1)*dcos(x)*dcos(x))+
     +  0.5d0*(4.0d0-mu(2))*mu(4)*y(17)*mu(1)*dcos(x)
      f(21)=-2.0d0*y(3)*y(21)
      f(22)=y(21)
      f(23)=y(21)*(-0.5d0*(4.0d0-mu(2))*mu(4)*(1.0d0+y(1))-y(2)*y(4)+
     +  mu(4)*(1.0d0+y(1))*y(4)*mu(1)*dcos(x))
      f(24)=(1.0d0/mu(1))*y(21)*(dcos(x)*y(2)-dsin(x)*
     +  (-y(2)*(y(3)+y(2)*mu(1)*dcos(x))+mu(4)*(1.0d0+y(1))*mu(1)*
     +  dcos(x)*(y(3)+y(2)*mu(1)*dcos(x))-mu(4)*(1.0d0+y(1))*y(2)-
     +  mu(4)*(1.0d0+y(1))*mu(1)*dsin(x)))
      f(25)=(1.0d0/mu(1))*y(21)*(dsin(x)*y(2)+dcos(x)*
     +  (-y(2)*(y(3)+y(2)*mu(1)*dcos(x))+mu(4)*(1.0d0+y(1))*mu(1)*
     +  dcos(x)*(y(3)+y(2)*mu(1)*dcos(x))-mu(4)*(1.0d0+y(1))*y(2)-
     +  mu(4)*(1.0d0+y(1))*mu(1)*dsin(x)))
      f(26)=(1.0d0/mu(1))*y(21)*(dcos(x)*(-0.5d0*mu(2)*y(2)-y(4)*
     +  (y(3)+y(2)*mu(1)*dcos(x))-mu(4)*(1.0d0+y(1))*y(4))-dsin(x)*
     +  (y(4)+y(2)*
     +  y(4)*mu(1)*dsin(x)-mu(4)*(1.0d0+y(1))*mu(1)*dsin(x)*(-0.5d0*
     +  (4.0d0-mu(2))+y(4)*mu(1)*dcos(x))))
      f(27)=(1.0d0/mu(1))*y(21)*(dsin(x)*(-0.5d0*mu(2)*y(2)-y(4)*
     +  (y(3)+y(2)*mu(1)*dcos(x))-mu(4)*(1.0d0+y(1))*y(4))+dcos(x)*
     +  (y(4)+y(2)*
     +  y(4)*mu(1)*dsin(x)-mu(4)*(1.0d0+y(1))*mu(1)*dsin(x)*(-0.5d0*
     +  (4.0d0-mu(2))+y(4)*mu(1)*dcos(x))))
      f(28)=(1.0d0/(mu(1)*mu(1)))*y(21)*((y(3)+y(2)*mu(1)*dcos(x))*
     +  (1.0d0+y(2)*mu(1)*dsin(x))+mu(4)*(1.0d0+y(1))*y(2)*mu(1)*
     +  dsin(x)-mu(4)*
     +  (1.0d0+y(1))*mu(1)*mu(1)*dcos(x)*dsin(x)*(y(3)+y(2)*mu(1)*
     +  dcos(x))+
     +  mu(4)*(1.0d0+y(1))*mu(1)*mu(1)*dsin(x)*dsin(x))
      RETURN
      END


      subroutine ludcmp(a,n,np,indx,d)

      implicit none

      integer n,np,indx(n),nmax
      real*8 d,a(np,np),tiny

      parameter (nmax=500,tiny=1.0d-20)

      integer i,imax,j,k
      real*8 aamax,dum,sum,vv(nmax)
      d=1.0d0
      do i=1,n
        aamax=0.0d0
        do j=1,n
          if (dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
        enddo
        if (aamax.eq.0.0d0) then
          write (*,'(a23)') 'ludcmp: singular matrix'
          stop
        endif
        vv(i)=1.0d0/aamax
      enddo
      do j=1,n
        do i=1,j-1
          sum=a(i,j)
          do k=1,i-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
        enddo
        aamax=0.0d0
        do i=j,n
          sum=a(i,j)
          do k=1,j-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
          dum=vv(i)*dabs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
        enddo
        if (j.ne.imax) then
          do k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          enddo
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if (a(j,j).eq.0.0d0) a(j,j)=tiny
        if (j.ne.n) then
          dum=1.0d0/a(j,j)
          do i=j+1,n
            a(i,j)=a(i,j)*dum
          enddo
        endif
      enddo

      return
      end


      subroutine lubksb(a,n,np,indx,b)

      implicit none

      integer n,np,indx(n)
      real*8 a(np,np),b(n)

      integer i,ii,j,ll
      real*8 sum

      ii=0
      do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0) then
          do j=ii,i-1
            sum=sum-a(i,j)*b(j)
          enddo
        else if (sum.ne.0.0d0) then
          ii=1
        endif
        b(i)=sum
      enddo
      do i=n,1,-1
        sum=b(i)
        do j=i+1,n
          sum=sum-a(i,j)*b(j)
        enddo
        b(i)=sum/a(i,i)
      enddo

      return
      end

      end module gordon
