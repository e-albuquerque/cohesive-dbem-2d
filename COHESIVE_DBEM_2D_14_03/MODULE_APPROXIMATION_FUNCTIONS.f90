MODULE MODULE_APPROXIMATION_FUNCTIONS
! -------------------------------------------------------------------------------------!
!   Description.                                                                       !
!                                                                                      !
! -------------------------------------------------------------------------------------!
!   This module contains the following subroutines:                                    !
!                                                                                      !                                              
!   [1]. SHAPE_FUNCTIONS (XSI,ELEM,PHI)                                                !
!   [2]. SHAPE_FUNCTIONS_DERIVATIVES (XSI,ELEM,dPHI)                                   !
!   [3]. BASIS_FUNCTIONS (XSI,ELEM,PHI)                                                !
!   [4]. BASIS_FUNCTIONS_DERIVATIVES (XSI,ELEM,dPHI)                                   !
! -------------------------------------------------------------------------------------!

! -------------------------------------------------------------------------------------!
!   USED DATA:                                                                         !
! -------------------------------------------------------------------------------------!
USE MODULE_DATA

CONTAINS

! ========================================================================================!
![1]                                                                                      !
SUBROUTINE SHAPE_FUNCTIONS (XSI,ELEM,PHI)                                                 !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!        
!
	IMPLICIT NONE 
!
    INTEGER,INTENT(IN)::ELEM
    REAL(KIND=wp),INTENT(IN)::XSI
    REAL(KIND=wp),INTENT(OUT)::PHI(ORDER(ELEM)+1)
!    
	INTEGER::I,J
	REAL(KIND=wp)::AUX1,XSI_I,XSI_J
!
	DO I=1,ORDER(ELEM)+1
		AUX1=1.D0
		XSI_I=XSI_NODES(ELEM,I)
		DO J=1,ORDER(ELEM)+1
			IF(I.NE.J) THEN
                XSI_J=XSI_NODES(ELEM,J)
				AUX1=AUX1*((XSI-XSI_J)/(XSI_I-XSI_J))
			ENDIF
		ENDDO
		PHI(I)=AUX1
	ENDDO
!       
END SUBROUTINE SHAPE_FUNCTIONS

! ========================================================================================!
![2]                                                                                      !
SUBROUTINE SHAPE_FUNCTIONS_DERIVATIVES (XSI,ELEM,dPHI)                                    !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!        
!
	IMPLICIT NONE 
!
    INTEGER,INTENT(IN)::ELEM
    REAL(KIND=wp),INTENT(IN)::XSI
    REAL(KIND=wp),INTENT(OUT)::dPHI(ORDER(ELEM)+1)
    
!    
	INTEGER::I,J,K,APONTA    		 
	REAL(KIND=wp)::AUX1,AUX2,AUX3,AUX4,XSI_I,XSI_J
!
	dPHI=0.D0
	DO I=1,ORDER(ELEM)+1
		AUX3=0.D0
		APONTA=0
        XSI_I=XSI_NODES(ELEM,I)
		DO K=1,ORDER(ELEM)
			AUX1=1.D0
			APONTA=APONTA+1
			IF(APONTA.EQ.I) THEN
				APONTA=APONTA+1
			ENDIF
!
			DO J=1,ORDER(ELEM)+1
				IF(I.NE.J) THEN
                    XSI_J=XSI_NODES(ELEM,J)
					IF(J.EQ.APONTA) THEN
						AUX4=1.D0
					ELSE
						AUX4=XSI-XSI_J
					ENDIF
					AUX2=((AUX4)/(XSI_I-XSI_J))
					AUX1=AUX1*AUX2
				ENDIF
			ENDDO
			AUX3=AUX3+AUX1
		ENDDO
		dPHI(I)=AUX3
    ENDDO
!
END SUBROUTINE SHAPE_FUNCTIONS_DERIVATIVES  

! ========================================================================================!
![3]                                                                                      !
SUBROUTINE BASIS_FUNCTIONS (XSI,ELEM,PHI)                                                 !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!        
!
	IMPLICIT NONE 
!
    INTEGER,INTENT(IN)::ELEM
    REAL(KIND=wp),INTENT(IN)::XSI
    REAL(KIND=wp),INTENT(OUT)::PHI(ORDER(ELEM)+1)
!    
	INTEGER::I,J		 
	REAL(KIND=wp)::AUX1,XSI_I,XSI_J
!
	DO I=1,ORDER(ELEM)+1
		AUX1=1.D0
		XSI_I=XSI_COLLOCATIONPOINTS(ELEM,I)
		DO J=1,ORDER(ELEM)+1
			IF(I.NE.J) THEN
                XSI_J=XSI_COLLOCATIONPOINTS(ELEM,J)
				AUX1=AUX1*((XSI-XSI_J)/(XSI_I-XSI_J))
			ENDIF
		ENDDO
		PHI(I)=AUX1
	ENDDO
!   
END SUBROUTINE BASIS_FUNCTIONS

! ========================================================================================!
![4]                                                                                      !
SUBROUTINE BASIS_FUNCTIONS_DERIVATIVES (XSI,ELEM,dPHI)                                    !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!        
!
	IMPLICIT NONE 
!
    INTEGER,INTENT(IN)::ELEM
    REAL(KIND=wp),INTENT(IN)::XSI
    REAL(KIND=wp),INTENT(OUT)::dPHI(ORDER(ELEM)+1)
    
!    
	INTEGER::I,J,K,APONTA    		 
	REAL(KIND=wp)::AUX1,AUX2,AUX3,AUX4,XSI_I,XSI_J
!
	dPHI=0.D0
	DO I=1,ORDER(ELEM)+1
		AUX3=0.D0
		APONTA=0
        XSI_I=XSI_COLLOCATIONPOINTS(ELEM,I)
		DO K=1,ORDER(ELEM)
			AUX1=1.D0
			APONTA=APONTA+1
			IF(APONTA.EQ.I) THEN
				APONTA=APONTA+1
			ENDIF
!
			DO J=1,ORDER(ELEM)+1
				IF(I.NE.J) THEN
                    XSI_J=XSI_COLLOCATIONPOINTS(ELEM,J)
					IF(J.EQ.APONTA) THEN
						AUX4=1.D0
					ELSE
						AUX4=XSI-XSI_J
					ENDIF
					AUX2=((AUX4)/(XSI_I-XSI_J))
					AUX1=AUX1*AUX2
				ENDIF
			ENDDO
			AUX3=AUX3+AUX1
		ENDDO
		dPHI(I)=AUX3
    ENDDO
!
END SUBROUTINE BASIS_FUNCTIONS_DERIVATIVES
!
END MODULE MODULE_APPROXIMATION_FUNCTIONS