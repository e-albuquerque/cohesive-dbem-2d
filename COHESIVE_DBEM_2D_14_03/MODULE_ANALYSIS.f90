MODULE MODULE_ANALYSIS
! -------------------------------------------------------------------------------------!
!   Description.                                                                       !
!                                                                                      !
! -------------------------------------------------------------------------------------!
!   This module contains the following subroutines:                                    !
!                                                                                      !
!   [1]. H_G_MATRICES                                                                  !
!   [2]. LINEAR_ANALYSIS                                                               !
!   [3]. COHESIVE_ANALYSIS                                                             !
!   [4]. ANALYSIS                                                                      !
! -------------------------------------------------------------------------------------!

USE MODULE_DATA
USE MODULE_NUMERICAL_INTEGRATION
USE MODULE_SOLVER
USE omp_lib

CONTAINS

! ========================================================================================!
![1]                                                                                      !
SUBROUTINE H_G_MATRICES                                                                   !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
IMPLICIT NONE

INTEGER::I,J,K,JJ,KK,ELEM0,N_COLS
INTEGER::NODE_ELEM(N_NODES)
REAL(KIND=wp)::ge(2,2*(MAX_ORDER+1)),he(2,2*(MAX_ORDER+1))
REAL(KIND=wp)::DIST,DIST_MIN,free_term(2,2),TOL

ALLOCATE(H(2*N_NODES,2*N_NODES),G(2*N_NODES,2*N_NODES))
H = 0.D0
G = 0.D0
TOL = 1.D-10

! Precompute element index associated with each collocation node (O(N_ELEM) lookup table)
NODE_ELEM = 0
DO K=1,N_ELEM
    DO KK=1,(ORDER(K)+1)
        NODE_ELEM(CONECTIVITY(K,KK)) = K
    ENDDO
ENDDO

! Parallel assembly of influence matrices H and G across collocation nodes
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(I, J, K, JJ, ELEM0, N_COLS, DIST, DIST_MIN, ge, he) &
!$OMP SCHEDULE(DYNAMIC)
DO I=1,N_NODES
    JJ = NODE_ELEM(I)
    ELEM0 = JJ

    DO J=1,N_ELEM
        IF(SUBREGION(JJ).EQ.SUBREGION(J)) THEN
            DIST_MIN = 1.D0
            DO K=1,ORDER(J)+1
                DIST = DSQRT((COORDCOLOCATIONPOINTS(CONECTIVITY(J,K),1)-COORDCOLOCATIONPOINTS(I,1))**2.D0+&
                    (COORDCOLOCATIONPOINTS(CONECTIVITY(J,K),2)-COORDCOLOCATIONPOINTS(I,2))**2.D0)
                IF(DIST.LE.DIST_MIN)THEN
                    DIST_MIN = DIST
                ENDIF
            ENDDO

            N_COLS = 2*(ORDER(J)+1)
            ge(:,1:N_COLS) = 0.D0
            he(:,1:N_COLS) = 0.D0

            IF(DIST_MIN.GT.TOL)THEN
                IF(COLLOCATION_IE(I).EQ."D")THEN
                    CALL NONSINGULARINTEGRATIONS_DBIE(I,J,ge(:,1:N_COLS),he(:,1:N_COLS))
                ELSE
                    CALL NONSINGULARINTEGRATIONS_TBIE(I,J,ELEM0,ge(:,1:N_COLS),he(:,1:N_COLS))
                ENDIF
            ELSE
                IF(ABS(CP_DISTANCING).LE.1.D-14)THEN
                    IF(COLLOCATION_IE(I).EQ."D")THEN
                        CALL SINGULARINTEGRATIONS_DBIE(I,J,ge(:,1:N_COLS),he(:,1:N_COLS))
                    ELSE
                        CALL SINGULARINTEGRATIONS_TBIE(I,J,ge(:,1:N_COLS),he(:,1:N_COLS))
                    ENDIF
                ELSE
                    IF(COLLOCATION_IE(I).EQ."D")THEN
                        CALL NONSINGULARINTEGRATIONS_DBIE(I,J,ge(:,1:N_COLS),he(:,1:N_COLS))
                    ELSE
                        CALL NONSINGULARINTEGRATIONS_TBIE(I,J,ELEM0,ge(:,1:N_COLS),he(:,1:N_COLS))
                    ENDIF
                ENDIF
            END IF

            DO K=1,ORDER(J)+1
                G(2*I-1,2*CONECTIVITY(J,K)-1) = G(2*I-1,2*CONECTIVITY(J,K)-1) + ge(1,2*K-1)
                G(2*I-1,2*CONECTIVITY(J,K))   = G(2*I-1,2*CONECTIVITY(J,K))   + ge(1,2*K)
                G(2*I,2*CONECTIVITY(J,K)-1)   = G(2*I,2*CONECTIVITY(J,K)-1)   + ge(2,2*K-1)
                G(2*I,2*CONECTIVITY(J,K))     = G(2*I,2*CONECTIVITY(J,K))     + ge(2,2*K)

                H(2*I-1,2*CONECTIVITY(J,K)-1) = H(2*I-1,2*CONECTIVITY(J,K)-1) + he(1,2*K-1)
                H(2*I-1,2*CONECTIVITY(J,K))   = H(2*I-1,2*CONECTIVITY(J,K))   + he(1,2*K)
                H(2*I,2*CONECTIVITY(J,K)-1)   = H(2*I,2*CONECTIVITY(J,K)-1)   + he(2,2*K-1)
                H(2*I,2*CONECTIVITY(J,K))     = H(2*I,2*CONECTIVITY(J,K))     + he(2,2*K)
            ENDDO
        ENDIF
    ENDDO
ENDDO
!$OMP END PARALLEL DO

! Rigid-body motion diagonal free terms
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J, free_term)
DO I=1,N_NODES
    IF(COLLOCATION_IE(I).EQ."D")THEN
        free_term = 0.D0
        DO J=1,N_NODES
            free_term(1,1) = free_term(1,1) - H(2*I-1,2*J-1)
            free_term(1,2) = free_term(1,2) - H(2*I-1,2*J)
            free_term(2,1) = free_term(2,1) - H(2*I,2*J-1)
            free_term(2,2) = free_term(2,2) - H(2*I,2*J)
        ENDDO
        H(2*I-1,2*I-1) = H(2*I-1,2*I-1) + free_term(1,1)
        H(2*I,2*I)     = H(2*I,2*I)     + free_term(2,2)
    ELSE IF(COLLOCATION_IE(I).EQ."T")THEN
        free_term = 0.D0
        DO J=1,N_NODES
            free_term(1,1) = free_term(1,1) - H(2*I-1,2*J-1)
            free_term(1,2) = free_term(1,2) - H(2*I-1,2*J)
            free_term(2,1) = free_term(2,1) - H(2*I,2*J-1)
            free_term(2,2) = free_term(2,2) - H(2*I,2*J)
        ENDDO
        H(2*I-1,2*I-1) = H(2*I-1,2*I-1) + free_term(1,1)
        H(2*I,2*I)     = H(2*I,2*I)     + free_term(2,2)
        IF(ABS(CP_DISTANCING).LE.1.D-14)THEN
            G(2*I-1,2*I-1) = G(2*I-1,2*I-1) - 0.5D0
            G(2*I,2*I)     = G(2*I,2*I)     - 0.5D0
        ENDIF
    ENDIF
ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE H_G_MATRICES

! ========================================================================================!
![2]                                                                                      !
SUBROUTINE LINEAR_ANALYSIS                                                                !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
IMPLICIT NONE

CALL H_G_MATRICES
CALL SOLVE_LINEAR_SYSTEM

END SUBROUTINE LINEAR_ANALYSIS

! ========================================================================================!
![3]                                                                                      !
SUBROUTINE COHESIVE_ANALYSIS                                                              !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
IMPLICIT NONE
REAL(KIND=wp)::TI_HG,TF_HG,DT_HG

OPEN(40,FILE='7_OUTPUT_ANALYSIS_TIME.txt',STATUS='UNKNOWN')
WRITE(40,*) "BUILDING H AND G MATRICES"

TI_HG = omp_get_wtime()
CALL H_G_MATRICES
TF_HG = omp_get_wtime()
DT_HG = TF_HG - TI_HG

WRITE(40,500)DT_HG
CALL SOLVE_NONLINEAR_SYSTEM

500 FORMAT(3x,E14.7)

END SUBROUTINE COHESIVE_ANALYSIS

! ========================================================================================!
![4]                                                                                      !
SUBROUTINE ANALYSIS                                                                       !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
IMPLICIT NONE

IF(ANALYSIS_TYPE.EQ."LE")THEN
    CALL LINEAR_ANALYSIS
ELSE IF(ANALYSIS_TYPE.EQ."CP")THEN
    CALL COHESIVE_ANALYSIS
END IF

END SUBROUTINE ANALYSIS
!
END MODULE MODULE_ANALYSIS
