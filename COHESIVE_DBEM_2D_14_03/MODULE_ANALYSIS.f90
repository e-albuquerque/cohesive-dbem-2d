MODULE MODULE_ANALYSIS
! -------------------------------------------------------------------------------------!
!   Description.                                                                       !
!                                                                                      !
! -------------------------------------------------------------------------------------!
!   This module contains the following subroutines:                                    !
!                                                                                      !
!   [1]. H_G_BEM_MATRICES                                                              !
!   [2]. LINEAR_ANALYSIS                                                               !
!   [3]. COHESIVE_ANALYSIS                                                             !
!   [4]. ANALYSIS                                                                      !
! -------------------------------------------------------------------------------------!

! -------------------------------------------------------------------------------------!
!   USED MODULE:                                                                       !
! -------------------------------------------------------------------------------------!
USE MODULE_DATA
USE MODULE_NUMERICAL_INTEGRATION
USE MODULE_SOLVER

CONTAINS

! ========================================================================================!
![1]                                                                                      !
SUBROUTINE H_G_MATRICES                                                                   !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!        
!
IMPLICIT NONE 
!
INTEGER::I,J,K,JJ,KK,ELEM0
REAL(KIND=wp),ALLOCATABLE::ge(:,:),he(:,:)
REAL(KIND=wp)::DIST,DIST_MIN,free_term(2,2),TOL
!
ALLOCATE(H(2*N_NODES,2*N_NODES),G(2*N_NODES,2*N_NODES))
H=0.D0
G=0.D0
TOL=1.D-10
!
DO I=1,N_NODES
!        
    JJ=0
	DO K=1,N_ELEM
		DO KK=1,(ORDER(K)+1)
			IF(CONECTIVITY(K,KK).EQ.I) THEN
				JJ=K
			ENDIF
        ENDDO
    ENDDO    
!        
    DO J=1,N_ELEM
        DO K=1,ORDER(J)+1
            IF(CONECTIVITY(J,K).EQ.I)THEN
                ELEM0=J
            ENDIF
        ENDDO
    ENDDO
!        
    DO J=1,N_ELEM
        IF(SUBREGION(JJ).EQ.SUBREGION(J)) THEN
            DIST_MIN=1.D0
            DO K=1,ORDER(J)+1
                DIST=DSQRT((COORDCOLOCATIONPOINTS(CONECTIVITY(J,K),1)-COORDCOLOCATIONPOINTS(I,1))**2.D0+&
                    (COORDCOLOCATIONPOINTS(CONECTIVITY(J,K),2)-COORDCOLOCATIONPOINTS(I,2))**2.D0)    
                IF(DIST.LE.DIST_MIN)THEN
                    DIST_MIN=DIST
                ENDIF    
            ENDDO
            ALLOCATE(he(2,2*(ORDER(J)+1)),ge(2,2*(ORDER(J)+1)))
            ge=0.d0
            he=0.d0
            IF(DIST_MIN.GT.TOL)THEN
                IF(COLLOCATION_IE(I).EQ."D")THEN
                    CALL NONSINGULARINTEGRATIONS_DBIE(I,J,ge,he)
                ELSE 
                    CALL NONSINGULARINTEGRATIONS_TBIE(I,J,ELEM0,ge,he)                
                ENDIF
            ELSE
                IF(ABS(CP_DISTANCING).LE.1.D-14)THEN
                    IF(COLLOCATION_IE(I).EQ."D")THEN
                        CALL SINGULARINTEGRATIONS_DBIE(I,J,ge,he) 
                    ELSE 
                        CALL SINGULARINTEGRATIONS_TBIE(I,J,ge,he)                     
                    ENDIF
                ELSE
                    IF(COLLOCATION_IE(I).EQ."D")THEN
                        CALL NONSINGULARINTEGRATIONS_DBIE(I,J,ge,he)
                    ELSE
                        CALL NONSINGULARINTEGRATIONS_TBIE(I,J,ELEM0,ge,he)                                    
                    ENDIF
                ENDIF
            END IF 
            DO  K=1,ORDER(J)+1  
                G(2*I-1,2*CONECTIVITY(J,K)-1)=G(2*I-1,2*CONECTIVITY(J,K)-1)+ge(1,2*K-1)
                G(2*I-1,2*CONECTIVITY(J,K))=G(2*I-1,2*CONECTIVITY(J,K))+ge(1,2*K)
                G(2*I,2*CONECTIVITY(J,K)-1)=G(2*I,2*CONECTIVITY(J,K)-1)+ge(2,2*K-1)
                G(2*I,2*CONECTIVITY(J,K))=G(2*I,2*CONECTIVITY(J,K))+ge(2,2*K)
!                
                H(2*I-1,2*CONECTIVITY(J,K)-1)=H(2*I-1,2*CONECTIVITY(J,K)-1)+he(1,2*K-1)
                H(2*I-1,2*CONECTIVITY(J,K))=H(2*I-1,2*CONECTIVITY(J,K))+he(1,2*K)
                H(2*I,2*CONECTIVITY(J,K)-1)=H(2*I,2*CONECTIVITY(J,K)-1)+he(2,2*K-1)
                H(2*I,2*CONECTIVITY(J,K))=H(2*I,2*CONECTIVITY(J,K))+he(2,2*K)
            ENDDO
            DEALLOCATE(ge,he)
        ENDIF
    ENDDO  
ENDDO
!    
DO I=1,N_NODES
    IF(COLLOCATION_IE(I).EQ."D")THEN
        free_term=0.D0
        DO J=1,N_NODES
            free_term(1,1)=free_term(1,1)-H(2*I-1,2*J-1)
            free_term(1,2)=free_term(1,2)-H(2*I-1,2*J)
            free_term(2,1)=free_term(2,1)-H(2*I,2*J-1)
            free_term(2,2)=free_term(2,2)-H(2*I,2*J)
        ENDDO
        H(2*I-1,2*I-1)=H(2*I-1,2*I-1)+free_term(1,1)
        H(2*I,2*I)=H(2*I,2*I)+free_term(2,2)
    ELSE IF(COLLOCATION_IE(I).EQ."T")THEN
        free_term=0.D0
        DO J=1,N_NODES
            free_term(1,1)=free_term(1,1)-H(2*I-1,2*J-1)
            free_term(1,2)=free_term(1,2)-H(2*I-1,2*J)
            free_term(2,1)=free_term(2,1)-H(2*I,2*J-1)
            free_term(2,2)=free_term(2,2)-H(2*I,2*J)
        ENDDO
        H(2*I-1,2*I-1)=H(2*I-1,2*I-1)+free_term(1,1)
        H(2*I,2*I)=H(2*I,2*I)+free_term(2,2)
        IF(ABS(CP_DISTANCING).LE.1.D-14)THEN
            free_term=0.D0
            free_term(1,1)=0.5D0
            free_term(2,2)=0.5D0
            G(2*I-1,2*I-1)=G(2*I-1,2*I-1)-free_term(1,1)
            G(2*I,2*I)=G(2*I,2*I)-free_term(2,2)
        ENDIF
    ENDIF
ENDDO
!        
END SUBROUTINE H_G_MATRICES     

! ========================================================================================!
![2]                                                                                      !
SUBROUTINE LINEAR_ANALYSIS                                                                !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!        
!
IMPLICIT NONE 
!
CALL H_G_MATRICES
CALL SOLVE_LINEAR_SYSTEM
! 
END SUBROUTINE LINEAR_ANALYSIS 

! ========================================================================================!
![3]                                                                                      !
SUBROUTINE COHESIVE_ANALYSIS                                                              !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!        
!
IMPLICIT NONE
REAL(KIND=wp)::TI_HG,TF_HG,DT_HG
!
OPEN(40,FILE='7_OUTPUT_ANALYSIS_TIME.txt',STATUS='UNKNOWN')
WRITE(40,*) "BUILDING H AND G MATRICES"

CALL CPU_TIME(TI_HG)
CALL H_G_MATRICES
CALL CPU_TIME(TF_HG)
DT_HG=TF_HG-TI_HG

WRITE(40,500)DT_HG
CALL SOLVE_NONLINEAR_SYSTEM
!
500 FORMAT(3x,E14.7)
!
END SUBROUTINE COHESIVE_ANALYSIS 

! ========================================================================================!
![4]                                                                                      !
SUBROUTINE ANALYSIS                                                                       !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!        
!        
IMPLICIT NONE 
!
IF(ANALYSIS_TYPE.EQ."LE")THEN
    CALL LINEAR_ANALYSIS
ELSE IF(ANALYSIS_TYPE.EQ."CP")THEN
    CALL COHESIVE_ANALYSIS
END IF
!      
END SUBROUTINE ANALYSIS
!
END MODULE MODULE_ANALYSIS