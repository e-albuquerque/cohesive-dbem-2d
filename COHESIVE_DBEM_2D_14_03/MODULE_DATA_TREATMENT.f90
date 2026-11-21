MODULE MODULE_DATA_TREATMENT
! -------------------------------------------------------------------------------------!
!   Description.                                                                       !
!                                                                                      !
! -------------------------------------------------------------------------------------!
!   This module contains the following subroutines:                                    !
!                                                                                      !
!   [1]. COLLOCATION_POINT_DISTANCING                                                  !
!   [2]. COLLOCATION_POINTS_COORDINATES                                                !
!   [3]. COHESIVE_SURFACES_DEFINITIONS                                                 !
!   [4]. PRESCRIBED_VALUES_FOR_U_AND_T                                                 !
!   [5]. TREATMENT_OF_DATA                                                             !
! -------------------------------------------------------------------------------------!

USE MODULE_DATA
USE MODULE_GEOMETRIC_MAPPINGS

CONTAINS

! ========================================================================================!
![1]
SUBROUTINE COLLOCATION_POINT_DISTANCING
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    IMPLICIT NONE

    INTEGER::J,K
    REAL(KIND=wp)::L,n(2),JAC
    REAL(KIND=8)::VALUES(MAX_ORDER+1,2)

    DO J=1,N_ELEM
        L=DSQRT((COORDNODES(CONECTIVITY(J,ORDER(J)+1),1)-COORDNODES(CONECTIVITY(J,1),1))**2.D0+&
            (COORDNODES(CONECTIVITY(J,ORDER(J)+1),2)-COORDNODES(CONECTIVITY(J,1),2))**2.D0)
        DO K=1,ORDER(J)+1
            VALUES(K,1)=COORDNODES(CONECTIVITY(J,K),1)
            VALUES(K,2)=COORDNODES(CONECTIVITY(J,K),2)
        END DO
        DO K=1,ORDER(J)+1
            CALL JACOBIAN_AND_NORMAL_VECTOR (XSI_COLLOCATIONPOINTS(J,K),J,VALUES(1:ORDER(J)+1,:),JAC,n)
            COORDCOLOCATIONPOINTS(CONECTIVITY(J,K),1)=COORDCOLOCATIONPOINTS(CONECTIVITY(J,K),1)+&
                CP_DISTANCING*L*n(1)
            COORDCOLOCATIONPOINTS(CONECTIVITY(J,K),2)=COORDCOLOCATIONPOINTS(CONECTIVITY(J,K),2)+&
                CP_DISTANCING*L*n(2)
        END DO
    END DO

END SUBROUTINE COLLOCATION_POINT_DISTANCING

! ========================================================================================!
![2]
SUBROUTINE COLLOCATION_POINTS_COORDINATES
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    IMPLICIT NONE

    INTEGER::I,J,K,CONT,ELEM,NODE_LOCAL
    REAL(KIND=wp)::X,Y
    REAL(KIND=8)::VALUES(MAX_ORDER+1,2)

!   COMPUTING THE COLLOCATION POINTS PHYSICAL AND PARAMETRIC COORDINATES
    COORDCOLOCATIONPOINTS=COORDNODES
    XSI_COLLOCATIONPOINTS=XSI_NODES

    DO I=1,N_NODES
        CONT=0
        ELEM=1
        NODE_LOCAL=1
        DO J=1,N_ELEM
            IF(CONECTIVITY(J,1).EQ.I) THEN
                CONT=CONT+1
                ELEM=J
                NODE_LOCAL=1
            ELSE IF (CONECTIVITY(J,ORDER(J)+1).EQ.I) THEN
                CONT=CONT+1
                ELEM=J
                NODE_LOCAL=ORDER(J)+1
            END IF
        ENDDO
        IF (CONT.EQ.1)THEN
            DO K=1,ORDER(ELEM)+1
                VALUES(K,1)=COORDNODES(CONECTIVITY(ELEM,K),1)
                VALUES(K,2)=COORDNODES(CONECTIVITY(ELEM,K),2)
            ENDDO
            IF(ELEM_TYPE(ELEM).EQ."D")THEN
                IF(NODE_LOCAL.EQ.1) THEN
                    XSI_COLLOCATIONPOINTS(ELEM,NODE_LOCAL)=-1.D0+(2.D0/DBLE(ORDER(ELEM)))/3.D0
                    CALL COORDINATES (XSI_COLLOCATIONPOINTS(ELEM,NODE_LOCAL),ELEM,VALUES(1:ORDER(ELEM)+1,:),X,Y)
                    COORDCOLOCATIONPOINTS(CONECTIVITY(ELEM,NODE_LOCAL),1)=X
                    COORDCOLOCATIONPOINTS(CONECTIVITY(ELEM,NODE_LOCAL),2)=Y
                ELSE IF (NODE_LOCAL.EQ.ORDER(ELEM)+1) THEN
                    XSI_COLLOCATIONPOINTS(ELEM,NODE_LOCAL)=1.D0-(2.D0/DBLE(ORDER(ELEM)))/3.D0
                    CALL COORDINATES (XSI_COLLOCATIONPOINTS(ELEM,NODE_LOCAL),ELEM,VALUES(1:ORDER(ELEM)+1,:),X,Y)
                    COORDCOLOCATIONPOINTS(CONECTIVITY(ELEM,NODE_LOCAL),1)=X
                    COORDCOLOCATIONPOINTS(CONECTIVITY(ELEM,NODE_LOCAL),2)=Y
                END IF
            ELSE IF (ELEM_TYPE(ELEM).EQ."LD")THEN
                IF(NODE_LOCAL.EQ.1) THEN
                    XSI_COLLOCATIONPOINTS(ELEM,NODE_LOCAL)=-1.D0+(2.D0/DBLE(ORDER(ELEM)))/3.D0
                    CALL COORDINATES (XSI_COLLOCATIONPOINTS(ELEM,NODE_LOCAL),ELEM,VALUES(1:ORDER(ELEM)+1,:),X,Y)
                    COORDCOLOCATIONPOINTS(CONECTIVITY(ELEM,NODE_LOCAL),1)=X
                    COORDCOLOCATIONPOINTS(CONECTIVITY(ELEM,NODE_LOCAL),2)=Y
                ENDIF
            ELSE IF (ELEM_TYPE(ELEM).EQ."RD")THEN
                IF (NODE_LOCAL.EQ.ORDER(ELEM)+1) THEN
                    XSI_COLLOCATIONPOINTS(ELEM,NODE_LOCAL)=1.D0-(2.D0/DBLE(ORDER(ELEM)))/3.D0
                    CALL COORDINATES (XSI_COLLOCATIONPOINTS(ELEM,NODE_LOCAL),ELEM,VALUES(1:ORDER(ELEM)+1,:),X,Y)
                    COORDCOLOCATIONPOINTS(CONECTIVITY(ELEM,NODE_LOCAL),1)=X
                    COORDCOLOCATIONPOINTS(CONECTIVITY(ELEM,NODE_LOCAL),2)=Y
                END IF
            ENDIF
        ENDIF
    END DO

    CALL COLLOCATION_POINT_DISTANCING

    OPEN(13,FILE='6_OUTPUT_COLLOCATION_POINTS_COORDINATES.txt',STATUS='UNKNOWN')
    DO I=1,N_NODES
        WRITE(13,113) COORDCOLOCATIONPOINTS(I,1), COORDCOLOCATIONPOINTS(I,2)
    ENDDO
    113 FORMAT(F14.9,16x,F14.9)
    CLOSE(13)

END SUBROUTINE COLLOCATION_POINTS_COORDINATES

! ========================================================================================!
![3]
SUBROUTINE COHESIVE_SURFACES_DEFINITIONS
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    IMPLICIT NONE

    INTEGER::I,J,K
    REAL(KIND=wp)::D

    ALLOCATE (COHESIVE_INTERFACES_Ip_Im(N_NODES,3),COHESIVE_INTERFACES_Im_Ip(N_NODES,3))
    COHESIVE_INTERFACES_Ip_Im=0
    COHESIVE_INTERFACES_Im_Ip=0

!   DEFINING THE COHESIVE INTERFACES PAIR OF NODES
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J, K, D)
    DO I=1,N_NODES
        IF(COLLOCATION_TYPE(I).EQ."Ip")THEN
            DO J=1,N_ELEM
                DO K=1,ORDER(J)+1
                    IF(COLLOCATION_TYPE(CONECTIVITY(J,K)).EQ."Im")THEN
                        D=DSQRT((COORDCOLOCATIONPOINTS(I,1)-COORDCOLOCATIONPOINTS(CONECTIVITY(J,K),1))**2.D0+&
                            (COORDCOLOCATIONPOINTS(I,2)-COORDCOLOCATIONPOINTS(CONECTIVITY(J,K),2))**2.D0)
                        IF(D.LT.1.D-10)THEN
                            COHESIVE_INTERFACES_Ip_Im(I,2)=J
                            COHESIVE_INTERFACES_Ip_Im(I,3)=K
                        ENDIF
                    ENDIF
                ENDDO
            ENDDO
            COHESIVE_INTERFACES_Ip_Im(I,1)=CONECTIVITY(COHESIVE_INTERFACES_Ip_Im(I,2),COHESIVE_INTERFACES_Ip_Im(I,3))
        ELSE IF(COLLOCATION_TYPE(I).EQ."Im")THEN
            DO J=1,N_ELEM
                DO K=1,ORDER(J)+1
                    IF(COLLOCATION_TYPE(CONECTIVITY(J,K)).EQ."Ip")THEN
                        D=DSQRT((COORDCOLOCATIONPOINTS(I,1)-COORDCOLOCATIONPOINTS(CONECTIVITY(J,K),1))**2.D0+&
                            (COORDCOLOCATIONPOINTS(I,2)-COORDCOLOCATIONPOINTS(CONECTIVITY(J,K),2))**2.D0)
                        IF(D.LT.1.D-10)THEN
                            COHESIVE_INTERFACES_Im_Ip(I,2)=J
                            COHESIVE_INTERFACES_Im_Ip(I,3)=K
                        ENDIF
                    ENDIF
                ENDDO
            ENDDO
            COHESIVE_INTERFACES_Im_Ip(I,1)=CONECTIVITY(COHESIVE_INTERFACES_Im_Ip(I,2),COHESIVE_INTERFACES_Im_Ip(I,3))
        ENDIF
    ENDDO
    !$OMP END PARALLEL DO

END SUBROUTINE COHESIVE_SURFACES_DEFINITIONS

! ========================================================================================!
![4]
SUBROUTINE PRESCRIBED_VALUES_FOR_U_AND_T
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    IMPLICIT NONE

    INTEGER::J,K
    REAL(KIND=wp)::L,S

    U=0.D0
    T=0.D0

    DO J=1,N_ELEM
        L=DSQRT((COORDNODES(CONECTIVITY(J,ORDER(J)+1),1)-COORDNODES(CONECTIVITY(J,1),1))**2.D0+&
            (COORDNODES(CONECTIVITY(J,ORDER(J)+1),2)-COORDNODES(CONECTIVITY(J,1),2))**2.D0)
        DO K=1,ORDER(J)+1
            S=0.5D0*L*(XSI_COLLOCATIONPOINTS(J,K)+1.D0)
            IF(BC(2*CONECTIVITY(J,K)-1).EQ.1) THEN
                U(2*CONECTIVITY(J,K)-1)=BC_COEFFICIENTS(2*J-1,1)*S**2.D0+BC_COEFFICIENTS(2*J-1,2)*S+BC_COEFFICIENTS(2*J-1,3)
            ELSE
                T(2*CONECTIVITY(J,K)-1)=BC_COEFFICIENTS(2*J-1,1)*S**2.D0+BC_COEFFICIENTS(2*J-1,2)*S+BC_COEFFICIENTS(2*J-1,3)
            END IF
            IF(BC(2*CONECTIVITY(J,K)).EQ.1) THEN
                U(2*CONECTIVITY(J,K))=BC_COEFFICIENTS(2*J,1)*S**2.D0+BC_COEFFICIENTS(2*J,2)*S+BC_COEFFICIENTS(2*J,3)
            ELSE
                T(2*CONECTIVITY(J,K))=BC_COEFFICIENTS(2*J,1)*S**2.D0+BC_COEFFICIENTS(2*J,2)*S+BC_COEFFICIENTS(2*J,3)
            END IF
        END DO
    END DO

END SUBROUTINE PRESCRIBED_VALUES_FOR_U_AND_T

! ========================================================================================!
![5]
SUBROUTINE TREATMENT_OF_DATA
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    IMPLICIT NONE

    CALL COLLOCATION_POINTS_COORDINATES

    CALL COHESIVE_SURFACES_DEFINITIONS

    CALL PRESCRIBED_VALUES_FOR_U_AND_T

END SUBROUTINE TREATMENT_OF_DATA
!
END MODULE MODULE_DATA_TREATMENT
