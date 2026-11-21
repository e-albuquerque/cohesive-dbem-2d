MODULE MODULE_GEOMETRIC_MAPPINGS
! -------------------------------------------------------------------------------------!
!   Description.                                                                       !
!                                                                                      !
! -------------------------------------------------------------------------------------!
!   This module contains the following subroutines:                                    !
!                                                                                      !                                              
!   [1]. COORDINATES (XSI,ELEM,VALUES,X,Y)                                             !
!   [2]. JACOBIAN_AND_NORMAL_VECTOR (XSI,ELEM,VALUES,JAC,n)                            !
! -------------------------------------------------------------------------------------!

! -------------------------------------------------------------------------------------!
!   USED DATA:                                                                         !
! -------------------------------------------------------------------------------------!
USE MODULE_DATA

CONTAINS

! ========================================================================================!
![1]                                                                                      !
SUBROUTINE COORDINATES(XSI,ELEM,VALUES,X,Y)                                              !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!        
!
IMPLICIT NONE 
!
INTEGER,INTENT(IN)::ELEM
REAL(KIND=wp),INTENT(IN)::XSI,VALUES((ORDER(ELEM)+1),2)
REAL(KIND=wp),INTENT(OUT)::X,Y
!    
INTEGER::I,J  		 
REAL(KIND=wp)::AUX1,XSI_I,XSI_J
!
X=0.D0
Y=0.D0
DO I=1,ORDER(ELEM)+1
	AUX1=1.D0
	XSI_I=XSI_NODES(ELEM,I)
	DO J=1,ORDER(ELEM)+1
		IF(I.NE.J) THEN
            XSI_J=XSI_NODES(ELEM,J)
			AUX1=AUX1*((XSI-XSI_J)/(XSI_I-XSI_J))
		ENDIF
    ENDDO  
	X=X+AUX1*VALUES(I,1)
	Y=Y+AUX1*VALUES(I,2)
ENDDO
!         
END SUBROUTINE COORDINATES

! ========================================================================================!
![2]                                                                                      !
SUBROUTINE JACOBIAN_AND_NORMAL_VECTOR(XSI,ELEM,VALUES,JAC,n)                             !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!        
!
IMPLICIT NONE 
!
INTEGER,INTENT(IN)::ELEM
REAL(KIND=wp),INTENT(IN)::XSI,VALUES((ORDER(ELEM)+1),2)
REAL(KIND=wp),INTENT(OUT)::JAC,n(2)
!    
INTEGER::I,J,K,APONTA  		 
REAL(KIND=wp)::AUX1,AUX2,AUX3,AUX4,XSI_I,XSI_J,T(2)
!
T(1)=0.D0
T(2)=0.D0
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
	T(1)=T(1)+AUX3*VALUES(I,1)
	T(2)=T(2)+AUX3*VALUES(I,2)      
ENDDO
JAC=DSQRT(T(1)*T(1)+T(2)*T(2))
n(1)=T(2)/JAC
n(2)=-T(1)/JAC  
!
END SUBROUTINE JACOBIAN_AND_NORMAL_VECTOR    
!
END MODULE MODULE_GEOMETRIC_MAPPINGS