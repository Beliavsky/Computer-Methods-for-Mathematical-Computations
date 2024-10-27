PROGRAM ForsytheSamples
! ---------------------------------------------------------------------------
! PURPOSE - Test the procedures in ForsytheMalcolmMoler. Combined programs
!   from the original reference, "Computer Methods for Mathematical
!   Computations", by George E. Forsythe, Michael A. Malcolm, and 
!   Cleve B. Moler. Prentice-Hall, 1977.
! Writes a file called samples.dbg with the output from the routines.

! AUTHORS - George E. Forsythe, Michael A. Malcolm, & 
!    Cleve B. Moler, Stanford University (1977)
!    Ralph L. Carmichael, Public Domain Aeronautical Software

! REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
!   1977   0.1   FMM   Publication of book
! 10Apr92  1.0   RLC   Original coding for F90
! 01Jun97  1.1   RLC   Compatible with Elf90
! 04Apr02  1.2   RLC   Combined all sample programs into one.
! 05Apr02  1.3   RLC   Put the test functions in a module
! 18Apr02  1.4   RLC   Adjusted error settings in Sample6
! 27Apr02  1.5   RLC   Changed calling sequence to FMMspline,Rkf45
! 28Jun02  1.6   RLC   Adjusted calling sequence of Decomp
! 22Aug02  1.7   RLC   Using the FMM module with both single and double precision
! 11Oct02  1.8   RLC   Some final renaming of variables and files  

USE ForsytheMalcolmMoler
IMPLICIT NONE

  INTEGER, PARAMETER :: DBG=3
  CHARACTER (LEN=*), PARAMETER :: GREETING = 'Samples - test the FMM routines'
  CHARACTER (LEN=*), PARAMETER :: AUTHOR= &
    'Ralph L. Carmichael, Public Domain Aeronautical Software'
  CHARACTER (LEN=*), PARAMETER :: VERSION='1.8 (11Oct02)'

  INTEGER,PARAMETER:: SP = SELECTED_REAL_KIND(6,20)
  INTEGER,PARAMETER:: DP = SELECTED_REAL_KIND(10,50) 
!----------------------------------------------------------------------------  
  WRITE(*,*) GREETING
  WRITE(*,*) AUTHOR
  WRITE(*,*) "VERSION "//VERSION
  OPEN(UNIT=DBG,FILE='samples.dbg',STATUS='REPLACE',ACTION='WRITE')
  WRITE(*,*) 'Uses module ForsytheMalcolmMoler, version '//FMM_VERSION
  WRITE(DBG,*) 'Uses module ForsytheMalcolmMoler, version '//FMM_VERSION

  CALL Sample3single()
  CALL Sample3Double()
  CALL Sample4single()
  CALL Sample4double()
  CALL Sample5single()
  CALL Sample5double()
  CALL Sample6single()
  CALL Sample6Double()
  CALL Sample7single()
  CALL Sample7double()
  CALL Sample8single()
  CALL Sample8double()
  CALL Sample9single()
  CALL Sample9double()
  
  STOP

CONTAINS

!+
Subroutine Sample3double()
! ---------------------------------------------------------------------------
! PURPOSE - Test the Decomp and Solve subroutines

   INTEGER, PARAMETER :: NDIM = 10
   REAL(DP)    :: cond,condp1
   INTEGER:: errCode
   INTEGER :: i,j
   INTEGER :: n
   REAL(DP), DIMENSION(NDIM,NDIM) :: a
   REAL(DP), DIMENSION(NDIM)      :: b
   INTEGER, DIMENSION(NDIM)   :: ipvt
!----------------------------------------------------------------------------
  n=3
  a(1,1)=10
  a(2,1)=-3
  a(3,1)=5
  a(1,2)=-7
  a(2,2)=2
  a(3,2)=-1
  a(1,3)=0
  a(2,3)=6
  a(3,3)=5

  DO i=1,n
    WRITE(DBG, '(1X,10F5.0)' ) (a(i,j),j=1,n)
  END DO

  CALL Decomp(a(1:3,1:3), ipvt, errCode, cond)
  IF (errCode /= 0) THEN
    WRITE(DBG,*) "Error in Decomp. The matrix is singular."
    RETURN
  ELSE    
    WRITE(DBG,*) 'condition number=', cond
    condp1=cond+1.0
    IF (condp1 <= cond) THEN
      WRITE(DBG,*) 'Matrix is singular to working precision'
      RETURN
    END IF
  END IF
  
  b(1)=7
  b(2)=4
  b(3)=6
  DO i=1,n
    WRITE(DBG, '(1X,8F10.5)' ) b(i)
  END DO
  CALL Solve(a(1:n,1:n), b(1:n), ipvt(1:n) )
  DO I=1,N
    WRITE(DBG, '(1X,8F10.5)' ) b(i)
  END DO
  WRITE(*,*)   'Successful completion of Sample3'
  WRITE(DBG,*) 'Successful completion of Sample3'
  RETURN
END Subroutine Sample3double   ! ==================================================

!+
Subroutine Sample3single()
! ---------------------------------------------------------------------------
! PURPOSE - Test the Decomp and Solve subroutines

   INTEGER, PARAMETER :: NDIM = 10
   REAL(SP)    :: cond,condp1
   INTEGER:: errCode
   INTEGER :: i,j
   INTEGER :: n
   REAL(SP), DIMENSION(NDIM,NDIM) :: a
   REAL(SP), DIMENSION(NDIM)      :: b
   INTEGER, DIMENSION(NDIM)   :: ipvt
!----------------------------------------------------------------------------
  n=3
  a(1,1)=10
  a(2,1)=-3
  a(3,1)=5
  a(1,2)=-7
  a(2,2)=2
  a(3,2)=-1
  a(1,3)=0
  a(2,3)=6
  a(3,3)=5

  DO i=1,n
    WRITE(DBG, '(1X,10F5.0)' ) (a(i,j),j=1,n)
  END DO

  CALL Decomp(a(1:3,1:3), ipvt, errCode, cond)
  IF (errCode /= 0) THEN
    WRITE(DBG,*) "Error in Decomp. The matrix is singular."
    RETURN
  ELSE    
    WRITE(DBG,*) 'condition number=', cond
    condp1=cond+1.0
    IF (condp1 <= cond) THEN
      WRITE(DBG,*) 'Matrix is singular to working precision'
      RETURN
    END IF
  END IF
  
  b(1)=7
  b(2)=4
  b(3)=6
  DO i=1,n
    WRITE(DBG, '(1X,8F10.5)' ) b(i)
  END DO
  CALL Solve(a(1:n,1:n), b(1:n), ipvt(1:n) )
  DO I=1,N
    WRITE(DBG, '(1X,8F10.5)' ) b(i)
  END DO
  WRITE(*,*)   'Successful completion of Sample3'
  WRITE(DBG,*) 'Successful completion of Sample3'
  RETURN
END Subroutine Sample3single   ! ==================================================

!+
Subroutine Sample4double()
! ---------------------------------------------------------------------------
! PURPOSE - Test the FMMspline and Seval subroutines
!
  REAL(DP) :: s,u
  INTEGER :: i,n
  REAL(DP), DIMENSION(10) :: x,y
  REAL(DP), DIMENSION(10) :: b,c,d
!----------------------------------------------------------------------------
  n=10
  DO i=1,n
    x(i)=REAL(i)
    y(i)=x(i)**3
  END DO

  WRITE(DBG,*) 'Preparing to call Spline with n=', n
  CALL FMMspline(x(1:n),y(1:n),b,c,d)

  u=2.5
  s=Seval(u,x,y,b,c,d)
  WRITE(DBG, *) u,s
  WRITE(*,*)   'Successful completion of Sample4'
  WRITE(DBG,*) 'Successful completion of Sample4'
  RETURN
END Subroutine Sample4double ! ----------------------------------------------------

!+
Subroutine Sample4single()
! ---------------------------------------------------------------------------
! PURPOSE - Test the FMMspline and Seval subroutines
!
  REAL(SP) :: s,u
  INTEGER :: i,n
  REAL(SP), DIMENSION(10) :: x,y
  REAL(SP), DIMENSION(10) :: b,c,d
!----------------------------------------------------------------------------
  n=10
  DO i=1,n
    x(i)=REAL(i)
    y(i)=x(i)**3
  END DO

  WRITE(DBG,*) 'Preparing to call Spline with n=', n
  CALL FMMspline(x(1:n),y(1:n),b,c,d)

  u=2.5
  s=Seval(u,x,y,b,c,d)
  WRITE(DBG, *) u,s
  WRITE(*,*)   'Successful completion of Sample4'
  WRITE(DBG,*) 'Successful completion of Sample4'
  RETURN
END Subroutine Sample4single ! ----------------------------------------------------

!+
Subroutine Sample5double()
! ---------------------------------------------------------------------------
! PURPOSE - Test Quanc8
USE TestFunctions,ONLY: F5double

  REAL(DP) :: a = 0.0
  REAL(DP) :: b = 2.0
  REAL(DP),PARAMETER :: abserr = 0.0
  REAL(DP) :: relerr = 1E-10
  REAL(DP) :: result,errest,flag
  INTEGER :: nofun
!----------------------------------------------------------------------------
  CALL Quanc8(F5double,a,b,abserr,relerr,result,errest,nofun,flag)
  WRITE(DBG, '(A,F15.9)' ) ' result=', result
  WRITE(DBG, '(A,E10.2)' ) ' error estimate=', errest
  WRITE(DBG,*) nofun, " function evaluations"
  IF (FLAG /= 0.0) WRITE(DBG,*) 'WARNING. Result may be unreliable. Flag=',flag
  WRITE(*,*)   'Successful completion of Sample5'
  WRITE(DBG,*) 'Successful completion of Sample5'

  RETURN
END Subroutine Sample5double   ! ==================================================

!+
Subroutine Sample5single()
! ---------------------------------------------------------------------------
! PURPOSE - Test Quanc8
USE TestFunctions,ONLY: F5single

  REAL(SP) :: a = 0.0
  REAL(SP) :: b = 2.0
  REAL(SP),PARAMETER :: abserr = 0.0
  REAL(SP) :: relerr = 1E-10
  REAL(SP) :: result,errest,flag
  INTEGER :: nofun
!----------------------------------------------------------------------------
  CALL Quanc8(F5single,a,b,abserr,relerr,result,errest,nofun,flag)
  WRITE(DBG, '(A,F15.9)' ) ' result=', result
  WRITE(DBG, '(A,E10.2)' ) ' error estimate=', errest
  WRITE(DBG,*) nofun, " function evaluations"
  IF (FLAG /= 0.0) WRITE(DBG,*) 'WARNING. Result may be unreliable. Flag=',flag
  WRITE(*,*)   'Successful completion of Sample5'
  WRITE(DBG,*) 'Successful completion of Sample5'

  RETURN
END Subroutine Sample5single   ! ==================================================


!+
Subroutine Sample6double()
! ---------------------------------------------------------------------------
! PURPOSE - Test RKF45
USE TestFunctions,ONLY: alfasqD,OrbitDouble
   INTEGER, PARAMETER :: NEQN = 4
   INTEGER :: iflag = 1
   REAL(DP) :: t
   REAL(DP) :: tout
   REAL(DP) :: tprint = 0.5
   REAL(DP) :: tfinal = 12.0
   REAL(DP) :: alfa = 3.141592653589/4.0
   REAL(DP) :: ecc
   REAL(DP) :: relerr = 1.0E-5
   REAL(DP) :: abserr = 1.0E-5
   REAL(DP), DIMENSION(NEQN) :: y
   REAL(DP), DIMENSION(6*NEQN+3) :: work
   INTEGER, DIMENSION(5) :: iwork
!----------------------------------------------------------------------------
  WRITE(DBG,*) 'Sample6 Results'
  ecc=0.25
  alfasqD=alfa*alfa
  y(1)=1.0-ecc
  y(2)=0.
  y(3)=0.
  y(4)=alfa*SQRT((1.0+ecc)/(1.0-ecc))
  t=0.0
  tout=t

  DO
    CALL RKF45(OrbitDouble,y,t,tout,relerr,abserr,iflag,work,iwork)
    WRITE(DBG, '(F5.1,2F15.9)' ) t,y(1),y(2)
    SELECT CASE (iflag)
      CASE (1)
        WRITE(DBG,*) 'IMPROPER CALL'
        STOP
      CASE(3)
        WRITE(DBG,*) 'Tolerances reset due to iflag=3 ', relerr, abserr
      CASE(4)
        WRITE(DBG,*) 'Warning.  Many steps...'
      CASE(5)
        abserr=1.0E-9
        WRITE(DBG,*) 'Tolerances reset due to iflag=5', relerr, abserr
      CASE(6)
        relerr=10.0*relerr
        WRITE(DBG,*) 'Tolerances reset due to iflag=6 ', relerr, abserr
      CASE(7)
        WRITE(DBG,*) 'WARNING. Much output'
        iflag=2
      CASE DEFAULT
        tout=t+tprint
        IF (t >= tfinal) EXIT
    END SELECT
  END DO
  WRITE(*,*)   'Successful completion of Sample6'
  WRITE(DBG,*) 'Successful completion of Sample6'

  RETURN
END Subroutine Sample6double   ! --------------------------------------------------

!+
Subroutine Sample6single()
! ---------------------------------------------------------------------------
! PURPOSE - Test RKF45
USE TestFunctions,ONLY: alfasqS,OrbitSingle
   INTEGER, PARAMETER :: NEQN = 4
   INTEGER :: iflag = 1
   REAL(SP) :: t
   REAL(SP) :: tout
   REAL(SP) :: tprint = 0.5
   REAL(SP) :: tfinal = 12.0
   REAL(SP) :: alfa = 3.141592653589/4.0
   REAL(SP) :: ecc
   REAL(SP) :: relerr = 1.0E-5
   REAL(SP) :: abserr = 1.0E-5
   REAL(SP), DIMENSION(NEQN) :: y
   REAL(SP), DIMENSION(6*NEQN+3) :: work
   INTEGER, DIMENSION(5) :: iwork
!----------------------------------------------------------------------------
  WRITE(DBG,*) 'Sample6 Results'
  ecc=0.25
  alfasqS=alfa*alfa
  y(1)=1.0-ecc
  y(2)=0.
  y(3)=0.
  y(4)=alfa*SQRT((1.0+ecc)/(1.0-ecc))
  t=0.0
  tout=t

  DO
    CALL RKF45(OrbitSingle,y,t,tout,relerr,abserr,iflag,work,iwork)
    WRITE(DBG, '(F5.1,2F15.9)' ) t,y(1),y(2)
    SELECT CASE (iflag)
      CASE (1)
        WRITE(DBG,*) 'IMPROPER CALL'
        STOP
      CASE(3)
        WRITE(DBG,*) 'Tolerances reset due to iflag=3 ', relerr, abserr
      CASE(4)
        WRITE(DBG,*) 'Warning.  Many steps...'
      CASE(5)
        abserr=1.0E-9
        WRITE(DBG,*) 'Tolerances reset due to iflag=5', relerr, abserr
      CASE(6)
        relerr=10.0*relerr
        WRITE(DBG,*) 'Tolerances reset due to iflag=6 ', relerr, abserr
      CASE(7)
        WRITE(DBG,*) 'WARNING. Much output'
        iflag=2
      CASE DEFAULT
        tout=t+tprint
        IF (t >= tfinal) EXIT
    END SELECT
  END DO
  WRITE(*,*)   'Successful completion of Sample6'
  WRITE(DBG,*) 'Successful completion of Sample6'

  RETURN
END Subroutine Sample6single   ! --------------------------------------------------

!+
Subroutine Sample7double()
! ---------------------------------------------------------------------------
! PURPOSE - Test Zeroin
USE TestFunctions,ONLY: F7double
  REAL(DP),PARAMETER:: a = 2.0
  REAL(DP),PARAMETER:: b = 3.0
  REAL(DP),PARAMETER:: TOL = 1E-10
  REAL(DP):: z
!----------------------------------------------------------------------------
  z=Zeroin(A,B,F7double,TOL)
  WRITE(DBG,*) 'Sample7 result=', z
  WRITE(*,*)   'Successful completion of Sample7'
  WRITE(DBG,*) 'Successful completion of Sample7'
  RETURN
END Subroutine Sample7double   ! --------------------------------------------------

!+
Subroutine Sample7single()
! ---------------------------------------------------------------------------
! PURPOSE - Test Zeroin
USE TestFunctions,ONLY: F7single
  REAL(SP),PARAMETER:: a = 2.0
  REAL(SP),PARAMETER:: b = 3.0
  REAL(SP),PARAMETER:: TOL = 1E-10
  REAL(SP):: z
!----------------------------------------------------------------------------
  z=Zeroin(A,B,F7single,TOL)
  WRITE(DBG,*) 'Sample7 result=', z
  WRITE(*,*)   'Successful completion of Sample7'
  WRITE(DBG,*) 'Successful completion of Sample7'
  RETURN
END Subroutine Sample7single   ! --------------------------------------------------

!+
Subroutine Sample8double()
! ---------------------------------------------------------------------------
! PURPOSE - Test Fmin
USE TestFunctions,ONLY:F8double

  REAL(DP) :: a = 0.0
  REAL(DP) :: b = 1.0
  REAL(DP),PARAMETER :: TOL = 1.0E-6
  REAL(DP) :: z
!----------------------------------------------------------------------------
  z=Fmin(a,b,F8double,TOL)
  WRITE(DBG,*) 'Sample8 result=', z
  WRITE(*,*)   'Successful completion of Sample8'
  WRITE(DBG,*) 'Successful completion of Sample8'

  RETURN
END Subroutine Sample8double   ! --------------------------------------------------

!+
Subroutine Sample8single()
! ---------------------------------------------------------------------------
! PURPOSE - Test Fmin
USE TestFunctions,ONLY:F8single

  REAL(SP) :: a = 0.0
  REAL(SP) :: b = 1.0
  REAL(SP),PARAMETER :: TOL = 1.0E-6
  REAL(SP) :: z
!----------------------------------------------------------------------------
  z=Fmin(a,b,F8single,TOL)
  WRITE(DBG,*) 'Sample8 result=', z
  WRITE(*,*)   'Successful completion of Sample8'
  WRITE(DBG,*) 'Successful completion of Sample8'

  RETURN
END Subroutine Sample8single   ! --------------------------------------------------

!+
Subroutine Sample9double()
! ---------------------------------------------------------------------------
! PURPOSE - Test SVD

  INTEGER:: i,ierr,j
  INTEGER:: m=5
  INTEGER:: n=3
  REAL(DP),DIMENSION(5,3):: a,u
  REAL(DP),DIMENSION(3,3):: v
  REAL(DP),DIMENSION(3):: sigma
!----------------------------------------------------------------------------
  DO i=1,m
    DO j=1,n
      a(i,j)=i+(j-1)*m
    END DO
  END DO

  CALL SVD(a, sigma, .TRUE., u, .TRUE., v, ierr)
  IF (IERR /= 0)  WRITE(*,*) 'TROUBLE!  ierr=', ierr
  WRITE(DBG,*) 'sigma'
  DO J=1,N
    WRITE(DBG,*) sigma(J)
  END DO
  WRITE(DBG,*) 'u'
  DO I=1,M
    WRITE(DBG,*) (u(i,j),j=1,n)
  END DO
  WRITE(DBG,*) 'v'
  DO I=1,N
    WRITE(DBG,*) (v(i,j),j=1,n)
  END DO
  WRITE(*,*)   'Successful completion of Sample9'
  WRITE(DBG,*) 'Successful completion of Sample9'

  RETURN
END Subroutine Sample9double   ! --------------------------------------------------

!+
Subroutine Sample9single()
! ---------------------------------------------------------------------------
! PURPOSE - Test SVD

  INTEGER:: i,ierr,j
  INTEGER:: m=5
  INTEGER:: n=3
  REAL(SP),DIMENSION(5,3):: a,u
  REAL(SP),DIMENSION(3,3):: v
  REAL(SP),DIMENSION(3):: sigma
!----------------------------------------------------------------------------
  DO i=1,m
    DO j=1,n
      a(i,j)=i+(j-1)*m
    END DO
  END DO

  CALL SVD(a, sigma, .TRUE., u, .TRUE., v, ierr)
  IF (IERR /= 0)  WRITE(*,*) 'TROUBLE!  ierr=', ierr
  WRITE(DBG,*) 'sigma'
  DO J=1,N
    WRITE(DBG,*) sigma(J)
  END DO
  WRITE(DBG,*) 'u'
  DO I=1,M
    WRITE(DBG,*) (u(i,j),j=1,n)
  END DO
  WRITE(DBG,*) 'v'
  DO I=1,N
    WRITE(DBG,*) (v(i,j),j=1,n)
  END DO
  WRITE(*,*)   'Successful completion of Sample9'
  WRITE(DBG,*) 'Successful completion of Sample9'

  RETURN
END Subroutine Sample9single   ! --------------------------------------------------


END Program ForsytheSamples   ! =============================================
