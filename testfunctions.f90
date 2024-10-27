MODULE TestFunctions
! ---------------------------------------------------------------------------
! PURPOSE - Collect some fumctions used by the sample programs
!----------------------------------------------------------------------------
IMPLICIT NONE

  INTEGER,PRIVATE,PARAMETER:: SP = SELECTED_REAL_KIND(6,20)
  INTEGER,PRIVATE,PARAMETER:: DP = SELECTED_REAL_KIND(10,50)


  REAL(SP),PUBLIC:: alfasqS   ! used by Orbit (single precision)
  REAL(DP),PUBLIC:: alfasqD   ! used by Orbit (double precision)

  INTERFACE F5
    MODULE PROCEDURE F5single, F5double
  END INTERFACE

  INTERFACE Orbit
    MODULE PROCEDURE OrbitSingle, OrbitDouble
  END INTERFACE

  INTERFACE F7
    MODULE PROCEDURE F7single, F7double
  END INTERFACE

  INTERFACE F8
    MODULE PROCEDURE F8single, F8double
  END INTERFACE

CONTAINS

!+
FUNCTION F5double(x) RESULT (F5Result)                            ! Used by Sample5
! ---------------------------------------------------------------------------
  REAL(DP), INTENT(IN) :: x
  REAL(DP) :: F5Result
!----------------------------------------------------------------------------
  IF (x == 0.0) THEN
    F5Result=1.0
  ELSE
    F5Result=SIN(x)/x
  END IF
  RETURN
END Function F5double   ! ---------------------------------------------------------

!+
FUNCTION F5single(x) RESULT(F5result)
! ---------------------------------------------------------------------------
  REAL(SP), INTENT(IN) :: x
  REAL(SP) :: F5Result
!----------------------------------------------------------------------------
  IF (x == 0.0) THEN
    F5Result=1.0
  ELSE
    F5Result=SIN(x)/x
  END IF
  RETURN
END Function F5single   ! ---------------------------------------------------------

!+
SUBROUTINE OrbitDouble (t,y,yp)                                   ! Used by Sample6
! ---------------------------------------------------------------------------
  REAL(DP), INTENT(IN) :: t
  REAL(DP), DIMENSION(:), INTENT(IN) :: y
  REAL(DP), DIMENSION(:), INTENT(OUT) :: yp
  REAL(DP) :: r
!----------------------------------------------------------------------------
  r=y(1)*y(1)+y(2)*y(2)
  r=r*SQRT(r)/alfasqD     ! alfasqD is a module variable, used here and Sample6
  yp(1)=y(3)
  yp(2)=y(4)
  yp(3)=-y(1)/r
  yp(4)=-y(2)/r
  RETURN 
END Subroutine OrbitDouble   ! ----------------------------------------------------

!+
SUBROUTINE OrbitSingle (t,y,yp)                                   ! Used by Sample6
! ---------------------------------------------------------------------------
  REAL(SP), INTENT(IN) :: t
  REAL(SP), DIMENSION(:), INTENT(IN) :: y
  REAL(SP), DIMENSION(:), INTENT(OUT) :: yp
  REAL(SP) :: r
!----------------------------------------------------------------------------
  r=y(1)*y(1)+y(2)*y(2)
  r=r*SQRT(r)/alfasqS     ! alfasqS is a module variable, used here and Sample6
  yp(1)=y(3)
  yp(2)=y(4)
  yp(3)=-y(1)/r
  yp(4)=-y(2)/r
  RETURN 
END Subroutine OrbitSingle   ! ----------------------------------------------------

!+
FUNCTION F7Double(x) RESULT (F7Return)                            ! Used by Sample7
! ---------------------------------------------------------------------------
  REAL(DP), INTENT(IN) :: x
  REAL(DP) :: F7Return
!----------------------------------------------------------------------------
  F7Return=x*(x*x-2.0)-5.0
  RETURN
END Function F7double   ! ---------------------------------------------------------

!+
FUNCTION F7single(x) RESULT (F7Return)                            ! Used by Sample7
! ---------------------------------------------------------------------------
  REAL(SP), INTENT(IN) :: x
  REAL(SP) :: F7Return
!----------------------------------------------------------------------------
  F7Return=x*(x*x-2.0)-5.0
  RETURN
END Function F7single   ! ---------------------------------------------------------

!+
FUNCTION F8double(x)   RESULT(F8Return)                           ! Used by Sample8
! ---------------------------------------------------------------------------
   REAL(DP), INTENT(IN) :: x
   REAL(DP):: F8Return
   REAL(DP),PARAMETER:: TWO=2.0, FIVE=5.0
!----------------------------------------------------------------------------
  F8Return=x*(x*x-TWO)-FIVE
  RETURN
END Function F8double   ! ---------------------------------------------------------

!+
FUNCTION F8single(x)   RESULT(F8Return)                           ! Used by Sample8
! ---------------------------------------------------------------------------
   REAL(SP), INTENT(IN) :: x
   REAL(SP):: F8Return
   REAL(SP),PARAMETER:: TWO=2.0, FIVE=5.0
!----------------------------------------------------------------------------
  F8Return=x*(x*x-TWO)-FIVE
  RETURN
END Function F8single   ! ---------------------------------------------------------

END Module TestFunctions   ! ================================================
