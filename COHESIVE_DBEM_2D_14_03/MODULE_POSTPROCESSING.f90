MODULE MODULE_POSTPROCESSING
! -------------------------------------------------------------------------------------!
!   Description.                                                                       !
!                                                                                      !
! -------------------------------------------------------------------------------------!
!   This module contains the following subroutines:                                    !
!                                                                                      !
!   [1]. POSTPROCESSING_INTERNAL_DISPLACEMENTS                                         !
!   [2]. POSTPROCESSING_INTERNAL_STRESSES_AND_STRAINS                                  !
!   [3]. POSTPROCESSING_BOUNDARY_DISPLACEMENTS_AND_TRACTIONS                           !
!   [4]. POSTPROCESSING_BOUNDARY_STRESSES_AND_STRAINS                                  !
!   [5]. PRINCIPAL_STRESSES                                                            !
!   [6]. PARAVIEW_OUTPUT                                                               !
!   [7]. BOUNDARY_OUTPUT                                                               !
!   [8]. POSTPROCESSING                                                                !
! -------------------------------------------------------------------------------------!

USE MODULE_DATA
USE omp_lib
USE MODULE_GEOMETRIC_MAPPINGS
USE MODULE_APPROXIMATION_FUNCTIONS

CONTAINS

! ========================================================================================!
![1]
SUBROUTINE POSTPROCESSING_INTERNAL_DISPLACEMENTS
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    IMPLICIT NONE

    INTEGER::II,JJ,KK,I,J,N_LOC,N_COLS
    REAL(KIND=8)::VALUES(MAX_ORDER+1,2),PHI(MAX_ORDER+1)
    REAL(KIND=8)::he(2,2*(MAX_ORDER+1)),ge(2,2*(MAX_ORDER+1)),MAT_PHI(2,2*(MAX_ORDER+1))
    REAL(KIND=wp)::pi,r(2),llrll,drdx(2),drdn,X,Y,JAC,JW,nv(2),Ustar(2,2),Tstar(2,2),IDENTITY(2,2)
    REAL(KIND=wp)::SHEAR_MODULUS,LONGITUDINAL_MODULUS,Nu

    IDENTITY = 0.D0
    IDENTITY(1,1) = 1.D0
    IDENTITY(2,2) = 1.D0
    pi = DACOS(-1.D0)
    ALLOCATE(HH(2*NI,2*N_NODES),GG(2*NI,2*N_NODES),UI(2*NI))
    HH = 0.D0
    GG = 0.D0

    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP PRIVATE(II, JJ, KK, I, J, N_LOC, N_COLS, VALUES, PHI, he, ge, MAT_PHI, &
    !$OMP         X, Y, JAC, JW, nv, r, llrll, drdx, drdn, LONGITUDINAL_MODULUS, Nu, SHEAR_MODULUS, Ustar, Tstar) &
    !$OMP SCHEDULE(DYNAMIC)
    DO II=1,NI
        DO JJ=1,N_ELEM
            IF(SUBREGION(JJ).EQ.SUBREGION_POINTS(NB+II)) THEN
                N_LOC  = ORDER(JJ) + 1
                N_COLS = 2*N_LOC
                DO KK=1,N_LOC
                    VALUES(KK,1) = COORDNODES(CONECTIVITY(JJ,KK),1)
                    VALUES(KK,2) = COORDNODES(CONECTIVITY(JJ,KK),2)
                ENDDO
                LONGITUDINAL_MODULUS = Young_Modulus(SUBREGION(JJ))
                Nu = Poisson(SUBREGION(JJ))
                SHEAR_MODULUS = LONGITUDINAL_MODULUS/(2.D0*(1.D0+Nu))

                he(:,1:N_COLS) = 0.D0
                ge(:,1:N_COLS) = 0.D0
                MAT_PHI(:,1:N_COLS) = 0.D0

                DO KK=1,N_GAUSS
                    CALL COORDINATES (COORD_GAUSS(KK),JJ,VALUES(1:N_LOC,:),X,Y)
                    CALL JACOBIAN_AND_NORMAL_VECTOR (COORD_GAUSS(KK),JJ,VALUES(1:N_LOC,:),JAC,nv)
                    CALL BASIS_FUNCTIONS (COORD_GAUSS(KK),JJ,PHI(1:N_LOC))
                    DO I=1,N_LOC
                        MAT_PHI(1,2*I-1) = PHI(I)
                        MAT_PHI(2,2*I)   = PHI(I)
                    ENDDO
                    r(1) = X - COORDINTERNALPOINTS(II,1)
                    r(2) = Y - COORDINTERNALPOINTS(II,2)
                    llrll = DSQRT(r(1)*r(1) + r(2)*r(2))
                    drdx = r / llrll
                    drdn = drdx(1)*nv(1) + drdx(2)*nv(2)

                    DO I=1,2
                        DO J=1,2
                            Ustar(I,J) = 1.D0/(8.D0*pi*SHEAR_MODULUS*(1.D0-Nu))*((3.D0-4.D0*Nu)*DLOG(1.D0/llrll)*&
                                IDENTITY(I,J)+drdx(I)*drdx(J))
                            Tstar(I,J) = -1.D0/(4.D0*pi*(1.D0-Nu)*llrll)*(drdn*((1.D0-2.D0*Nu)*IDENTITY(I,J)+&
                                2.D0*drdx(I)*drdx(J))+(1.D0-2.D0*Nu)*(drdx(J)*nv(I)-drdx(I)*nv(J)))
                        ENDDO
                    ENDDO
                    JW = JAC * WEIGHT_GAUSS(KK)
                    ge(:,1:N_COLS) = ge(:,1:N_COLS) + MATMUL(Ustar,MAT_PHI(:,1:N_COLS))*JW
                    he(:,1:N_COLS) = he(:,1:N_COLS) + MATMUL(Tstar,MAT_PHI(:,1:N_COLS))*JW
                ENDDO
                DO KK=1,N_LOC
                    HH(2*II-1,2*CONECTIVITY(JJ,KK)-1) = HH(2*II-1,2*CONECTIVITY(JJ,KK)-1) + he(1,2*KK-1)
                    HH(2*II-1,2*CONECTIVITY(JJ,KK))   = HH(2*II-1,2*CONECTIVITY(JJ,KK))   + he(1,2*KK)
                    HH(2*II,2*CONECTIVITY(JJ,KK)-1)   = HH(2*II,2*CONECTIVITY(JJ,KK)-1)   + he(2,2*KK-1)
                    HH(2*II,2*CONECTIVITY(JJ,KK))     = HH(2*II,2*CONECTIVITY(JJ,KK))     + he(2,2*KK)

                    GG(2*II-1,2*CONECTIVITY(JJ,KK)-1) = GG(2*II-1,2*CONECTIVITY(JJ,KK)-1) + ge(1,2*KK-1)
                    GG(2*II-1,2*CONECTIVITY(JJ,KK))   = GG(2*II-1,2*CONECTIVITY(JJ,KK))   + ge(1,2*KK)
                    GG(2*II,2*CONECTIVITY(JJ,KK)-1)   = GG(2*II,2*CONECTIVITY(JJ,KK)-1)   + ge(2,2*KK-1)
                    GG(2*II,2*CONECTIVITY(JJ,KK))     = GG(2*II,2*CONECTIVITY(JJ,KK))     + ge(2,2*KK)
                ENDDO
            ENDIF
        ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    UI = MATMUL(GG,T) - MATMUL(HH,U)
    DEALLOCATE (HH,GG)

END SUBROUTINE POSTPROCESSING_INTERNAL_DISPLACEMENTS

! ========================================================================================!
![2]
SUBROUTINE POSTPROCESSING_INTERNAL_STRESSES_AND_STRAINS
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    IMPLICIT NONE

    INTEGER::I,J,K,II,JJ,KK,N_LOC,N_COLS
    REAL(KIND=8)::VALUES(MAX_ORDER+1,2),PHI(MAX_ORDER+1)
    REAL(KIND=8)::he(4,2*(MAX_ORDER+1)),ge(4,2*(MAX_ORDER+1)),MAT_PHI(2,2*(MAX_ORDER+1))
    REAL(KIND=wp)::pi,r(2),llrll,drdx(2),drdn,X,Y,JAC,JW,nv(2),Sstar(2,2,2),Dstar(2,2,2),S(4,2),D(4,2),IDENTITY(3,3)
    REAL(KIND=wp)::SHEAR_MODULUS,LONGITUDINAL_MODULUS,Nu

    IDENTITY = 0.D0
    IDENTITY(1,1) = 1.D0
    IDENTITY(2,2) = 1.D0
    IDENTITY(3,3) = 1.D0
    pi = DACOS(-1.D0)
    ALLOCATE(HHH(4*NI,2*N_NODES),GGG(4*NI,2*N_NODES),SIGMAI(4*NI),SIGMA3D(N_POINTS,3,3),EPSILON3D(N_POINTS,3,3))
    HHH = 0.D0
    GGG = 0.D0

    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP PRIVATE(II, JJ, KK, I, J, K, N_LOC, N_COLS, VALUES, PHI, he, ge, MAT_PHI, &
    !$OMP         X, Y, JAC, JW, nv, r, llrll, drdx, drdn, LONGITUDINAL_MODULUS, Nu, SHEAR_MODULUS, Dstar, Sstar, D, S) &
    !$OMP SCHEDULE(DYNAMIC)
    DO II=1,NI
        DO JJ=1,N_ELEM
            IF(SUBREGION(JJ).EQ.SUBREGION_POINTS(NB+II)) THEN
                N_LOC  = ORDER(JJ) + 1
                N_COLS = 2*N_LOC
                DO KK=1,N_LOC
                    VALUES(KK,1) = COORDNODES(CONECTIVITY(JJ,KK),1)
                    VALUES(KK,2) = COORDNODES(CONECTIVITY(JJ,KK),2)
                ENDDO
                LONGITUDINAL_MODULUS = Young_Modulus(SUBREGION(JJ))
                Nu = Poisson(SUBREGION(JJ))
                SHEAR_MODULUS = LONGITUDINAL_MODULUS/(2.D0*(1.D0+Nu))

                he(:,1:N_COLS) = 0.D0
                ge(:,1:N_COLS) = 0.D0
                MAT_PHI(:,1:N_COLS) = 0.D0

                DO KK=1,N_GAUSS
                    CALL COORDINATES (COORD_GAUSS(KK),JJ,VALUES(1:N_LOC,:),X,Y)
                    CALL JACOBIAN_AND_NORMAL_VECTOR (COORD_GAUSS(KK),JJ,VALUES(1:N_LOC,:),JAC,nv)
                    CALL BASIS_FUNCTIONS (COORD_GAUSS(KK),JJ,PHI(1:N_LOC))
                    DO I=1,N_LOC
                        MAT_PHI(1,2*I-1) = PHI(I)
                        MAT_PHI(2,2*I)   = PHI(I)
                    ENDDO
                    r(1) = X - COORDINTERNALPOINTS(II,1)
                    r(2) = Y - COORDINTERNALPOINTS(II,2)
                    llrll = DSQRT(r(1)*r(1) + r(2)*r(2))
                    drdx = r / llrll
                    drdn = drdx(1)*nv(1) + drdx(2)*nv(2)

                    DO I=1,2
                        DO J=1,2
                            DO K=1,2
                                Dstar(I,J,K) = 1.D0/(4.D0*pi*(1.D0-Nu)*llrll)*((1.D0-2.D0*Nu)*(IDENTITY(I,J)*drdx(K)+&
                                    IDENTITY(I,K)*drdx(J)-IDENTITY(J,K)*drdx(I))+2.D0*drdx(I)*drdx(J)*drdx(K))

                                Sstar(I,J,K) = 2.D0*drdn*((1.D0-2.D0*Nu)*IDENTITY(I,J)*drdx(K)+&
                                    Nu*(IDENTITY(I,K)*drdx(J)+IDENTITY(J,K)*drdx(I))-4.D0*drdx(I)*drdx(J)*drdx(K))
                                Sstar(I,J,K) = Sstar(I,J,K)+2.D0*Nu*(nv(I)*drdx(J)*drdx(K)+nv(J)*drdx(I)*drdx(K))-&
                                    (1.D0-4.D0*Nu)*nv(K)*IDENTITY(I,J)+(1.D0-2.D0*Nu)*(nv(I)*IDENTITY(J,K)+&
                                    nv(J)*IDENTITY(I,K)+2.D0*nv(K)*drdx(I)*drdx(J))
                                Sstar(I,J,K) = 2.D0*SHEAR_MODULUS/(4.D0*pi*(1.D0-Nu)*llrll*llrll)*Sstar(I,J,K)
                            ENDDO
                        ENDDO
                    ENDDO
                    D(1,1)=Dstar(1,1,1); D(1,2)=Dstar(2,1,1)
                    D(2,1)=Dstar(1,2,1); D(2,2)=Dstar(2,2,1)
                    D(3,1)=Dstar(1,1,2); D(3,2)=Dstar(2,1,2)
                    D(4,1)=Dstar(1,2,2); D(4,2)=Dstar(2,2,2)

                    S(1,1)=Sstar(1,1,1); S(1,2)=Sstar(2,1,1)
                    S(2,1)=Sstar(1,2,1); S(2,2)=Sstar(2,2,1)
                    S(3,1)=Sstar(1,1,2); S(3,2)=Sstar(2,1,2)
                    S(4,1)=Sstar(1,2,2); S(4,2)=Sstar(2,2,2)

                    JW = JAC * WEIGHT_GAUSS(KK)
                    ge(:,1:N_COLS) = ge(:,1:N_COLS) + MATMUL(D,MAT_PHI(:,1:N_COLS))*JW
                    he(:,1:N_COLS) = he(:,1:N_COLS) + MATMUL(S,MAT_PHI(:,1:N_COLS))*JW
                ENDDO
                DO KK=1,N_LOC
                    GGG(4*II-3,2*CONECTIVITY(JJ,KK)-1) = GGG(4*II-3,2*CONECTIVITY(JJ,KK)-1) + ge(1,2*KK-1)
                    GGG(4*II-3,2*CONECTIVITY(JJ,KK))   = GGG(4*II-3,2*CONECTIVITY(JJ,KK))   + ge(1,2*KK)
                    GGG(4*II-2,2*CONECTIVITY(JJ,KK)-1) = GGG(4*II-2,2*CONECTIVITY(JJ,KK)-1) + ge(2,2*KK-1)
                    GGG(4*II-2,2*CONECTIVITY(JJ,KK))   = GGG(4*II-2,2*CONECTIVITY(JJ,KK))   + ge(2,2*KK)
                    GGG(4*II-1,2*CONECTIVITY(JJ,KK)-1) = GGG(4*II-1,2*CONECTIVITY(JJ,KK)-1) + ge(3,2*KK-1)
                    GGG(4*II-1,2*CONECTIVITY(JJ,KK))   = GGG(4*II-1,2*CONECTIVITY(JJ,KK))   + ge(3,2*KK)
                    GGG(4*II,2*CONECTIVITY(JJ,KK)-1)   = GGG(4*II,2*CONECTIVITY(JJ,KK)-1)   + ge(4,2*KK-1)
                    GGG(4*II,2*CONECTIVITY(JJ,KK))     = GGG(4*II,2*CONECTIVITY(JJ,KK))     + ge(4,2*KK)

                    HHH(4*II-3,2*CONECTIVITY(JJ,KK)-1) = HHH(4*II-3,2*CONECTIVITY(JJ,KK)-1) + he(1,2*KK-1)
                    HHH(4*II-3,2*CONECTIVITY(JJ,KK))   = HHH(4*II-3,2*CONECTIVITY(JJ,KK))   + he(1,2*KK)
                    HHH(4*II-2,2*CONECTIVITY(JJ,KK)-1) = HHH(4*II-2,2*CONECTIVITY(JJ,KK)-1) + he(2,2*KK-1)
                    HHH(4*II-2,2*CONECTIVITY(JJ,KK))   = HHH(4*II-2,2*CONECTIVITY(JJ,KK))   + he(2,2*KK)
                    HHH(4*II-1,2*CONECTIVITY(JJ,KK)-1) = HHH(4*II-1,2*CONECTIVITY(JJ,KK)-1) + he(3,2*KK-1)
                    HHH(4*II-1,2*CONECTIVITY(JJ,KK))   = HHH(4*II-1,2*CONECTIVITY(JJ,KK))   + he(3,2*KK)
                    HHH(4*II,2*CONECTIVITY(JJ,KK)-1)   = HHH(4*II,2*CONECTIVITY(JJ,KK)-1)   + he(4,2*KK-1)
                    HHH(4*II,2*CONECTIVITY(JJ,KK))     = HHH(4*II,2*CONECTIVITY(JJ,KK))     + he(4,2*KK)
                ENDDO
            ENDIF
        ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    SIGMAI = MATMUL(GGG,T) - MATMUL(HHH,U)
    DEALLOCATE (HHH,GGG)
    SIGMA3D = 0.D0
    DO I=1,NI
        II = NB + I
        LONGITUDINAL_MODULUS = Young_Modulus(SUBREGION_POINTS(II))
        Nu = Poisson(SUBREGION_POINTS(II))
        SHEAR_MODULUS = LONGITUDINAL_MODULUS/(2.D0*(1.D0+Nu))
        SIGMA3D(II,1,1) = SIGMAI(4*I-3)
        SIGMA3D(II,1,2) = SIGMAI(4*I-2)
        SIGMA3D(II,2,1) = SIGMAI(4*I-1)
        SIGMA3D(II,2,2) = SIGMAI(4*I)
        IF(PLANE_STATE.EQ."PD")THEN
            SIGMA3D(II,3,3) = Nu*(SIGMA3D(II,1,1)+SIGMA3D(II,2,2))
        ENDIF

    !   GENERALIZED HOOKE'S LAW FOR THE COMPUTATION OF THE STRAIN TENSOR
        DO J=1,3
            DO K=1,3
                EPSILON3D(II,J,K) = (1.D0/(2.D0*SHEAR_MODULUS))*SIGMA3D(II,J,K)-(Nu/LONGITUDINAL_MODULUS)*(SIGMA3D(II,1,1)+&
                    SIGMA3D(II,2,2)+SIGMA3D(II,3,3))*IDENTITY(J,K)
            ENDDO
        ENDDO
    ENDDO

END SUBROUTINE POSTPROCESSING_INTERNAL_STRESSES_AND_STRAINS

! ========================================================================================!
![3]
SUBROUTINE POSTPROCESSING_BOUNDARY_DISPLACEMENTS_AND_TRACTIONS
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    IMPLICIT NONE

    INTEGER::J,II,JJ
    REAL(KIND=8)::UVALUES(MAX_ORDER+1,2),TVALUES(MAX_ORDER+1,2),PHI(MAX_ORDER+1)
    REAL(KIND=wp)::U_COLLOCPOINTS(2*N_NODES),T_COLLOCPOINTS(2*N_NODES)

    U_COLLOCPOINTS = U
    T_COLLOCPOINTS = T
    U = 0.D0
    T = 0.D0

    DO J=1,N_ELEM
        DO II=1,ORDER(J)+1
            UVALUES(II,1) = U_COLLOCPOINTS(2*CONECTIVITY(J,II)-1)
            UVALUES(II,2) = U_COLLOCPOINTS(2*CONECTIVITY(J,II))
            TVALUES(II,1) = T_COLLOCPOINTS(2*CONECTIVITY(J,II)-1)
            TVALUES(II,2) = T_COLLOCPOINTS(2*CONECTIVITY(J,II))
        ENDDO
        DO II=1,ORDER(J)+1
            CALL BASIS_FUNCTIONS(XSI_NODES(J,II),J,PHI(1:ORDER(J)+1))
            U(2*CONECTIVITY(J,II)-1) = 0.D0
            U(2*CONECTIVITY(J,II))   = 0.D0
            T(2*CONECTIVITY(J,II)-1) = 0.D0
            T(2*CONECTIVITY(J,II))   = 0.D0
            DO JJ=1,ORDER(J)+1
                U(2*CONECTIVITY(J,II)-1) = U(2*CONECTIVITY(J,II)-1) + PHI(JJ)*UVALUES(JJ,1)
                U(2*CONECTIVITY(J,II))   = U(2*CONECTIVITY(J,II))   + PHI(JJ)*UVALUES(JJ,2)
                T(2*CONECTIVITY(J,II)-1) = T(2*CONECTIVITY(J,II)-1) + PHI(JJ)*TVALUES(JJ,1)
                T(2*CONECTIVITY(J,II))   = T(2*CONECTIVITY(J,II))   + PHI(JJ)*TVALUES(JJ,2)
            ENDDO
        ENDDO
    ENDDO

END SUBROUTINE POSTPROCESSING_BOUNDARY_DISPLACEMENTS_AND_TRACTIONS

! ========================================================================================!
![4]
SUBROUTINE POSTPROCESSING_BOUNDARY_STRESSES_AND_STRAINS
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    IMPLICIT NONE

    INTEGER::I,J,II,JJ,POINT1,POINT2
    REAL(KIND=8),ALLOCATABLE::VALUES(:,:),U_BAR(:,:),X_BAR(:,:),dPHI(:),EPSILONXX_BAR(:)
    REAL(KIND=8),ALLOCATABLE::T_BAR(:,:),SIGMAXX_BAR(:),SIGMAXY_BAR(:),SIGMAYY_BAR(:)
    REAL(KIND=8),ALLOCATABLE::SIGMA_BAR(:,:,:),SIGMA(:,:,:)
    REAL(KIND=wp)::nv(2),JAC,IDENTITY(3,3),dUBARdXSI,dXBARdXSI,DIST
    REAL(KIND=wp)::R(2,2),Rt(2,2),U_NODAL(2),T_NODAL(2),X(2)
    REAL(KIND=wp)::SHEAR_MODULUS,LONGITUDINAL_MODULUS,Nu

    IDENTITY = 0.D0
    IDENTITY(1,1) = 1.D0
    IDENTITY(2,2) = 1.D0
    IDENTITY(3,3) = 1.D0
    ALLOCATE(SIGMAB(4*NB))
    SIGMAB = 0.D0
    POINT1 = 1
    POINT2 = 1

    DO J=1,N_ELEM
        LONGITUDINAL_MODULUS = Young_Modulus(SUBREGION(J))
        Nu = Poisson(SUBREGION(J))
        SHEAR_MODULUS = LONGITUDINAL_MODULUS/(2.D0*(1.D0+Nu))
        ALLOCATE(U_BAR(ORDER(J)+1,2),T_BAR(ORDER(J)+1,2),X_BAR(ORDER(J)+1,2),EPSILONXX_BAR(ORDER(J)+1),&
            SIGMAXX_BAR(ORDER(J)+1),SIGMAXY_BAR(ORDER(J)+1),SIGMAYY_BAR(ORDER(J)+1),SIGMA_BAR(ORDER(J)+&
            1,2,2),SIGMA(ORDER(J)+1,2,2))
        ALLOCATE(VALUES((ORDER(J)+1),2),dPHI(ORDER(J)+1))
        DO II=1,ORDER(J)+1
            VALUES(II,1) = COORDNODES(CONECTIVITY(J,II),1)
            VALUES(II,2) = COORDNODES(CONECTIVITY(J,II),2)
        ENDDO
        DO II=1,ORDER(J)+1
            CALL JACOBIAN_AND_NORMAL_VECTOR (XSI_NODES(J,II),J,VALUES,JAC,nv)
            R(1,1) = nv(2)
            R(1,2) = nv(1)
            R(2,1) = -nv(1)
            R(2,2) = nv(2)
            Rt = TRANSPOSE(R)

            X(1) = COORDNODES(CONECTIVITY(J,II),1)
            X(2) = COORDNODES(CONECTIVITY(J,II),2)
            X_BAR(II,:) = MATMUL(Rt,X)
            U_NODAL(1) = U(2*CONECTIVITY(J,II)-1)
            U_NODAL(2) = U(2*CONECTIVITY(J,II))
            U_BAR(II,:) = MATMUL(Rt,U_NODAL)
            T_NODAL(1) = T(2*CONECTIVITY(J,II)-1)
            T_NODAL(2) = T(2*CONECTIVITY(J,II))
            T_BAR(II,:) = MATMUL(Rt,T_NODAL)
            SIGMAXY_BAR(II) = T_BAR(II,1)
            SIGMAYY_BAR(II) = T_BAR(II,2)
        ENDDO
        DO II=1,ORDER(J)+1
            EPSILONXX_BAR = 0.D0
            CALL SHAPE_FUNCTIONS_DERIVATIVES (XSI_NODES(J,II),J,dPHI)
            dXBARdXSI = 0.D0
            dUBARdXSI = 0.D0
            DO JJ=1,ORDER(J)+1
                dXBARdXSI = dXBARdXSI + dPHI(JJ)*X_BAR(JJ,1)
                dUBARdXSI = dUBARdXSI + dPHI(JJ)*U_BAR(JJ,1)
            ENDDO
            EPSILONXX_BAR(II) = dUBARdXSI / dXBARdXSI
            IF(PLANE_STATE.EQ."PT")THEN
                SIGMAXX_BAR(II) = LONGITUDINAL_MODULUS*EPSILONXX_BAR(II) + Nu*SIGMAYY_BAR(II)
            ELSE IF(PLANE_STATE.EQ."PD")THEN
                SIGMAXX_BAR(II) = (LONGITUDINAL_MODULUS*EPSILONXX_BAR(II) + Nu*SIGMAYY_BAR(II))/(1.D0-Nu**2.D0)
            ENDIF
            SIGMA_BAR(II,1,1) = SIGMAXX_BAR(II)
            SIGMA_BAR(II,1,2) = SIGMAXY_BAR(II)
            SIGMA_BAR(II,2,1) = SIGMAXY_BAR(II)
            SIGMA_BAR(II,2,2) = SIGMAYY_BAR(II)
            SIGMA(II,:,:) = MATMUL(R,SIGMA_BAR(II,:,:))
            SIGMA(II,:,:) = MATMUL(SIGMA(II,:,:),Rt)
            DIST = 0.D0
            IF(COLLOCATION_TYPE(CONECTIVITY(J,II)).EQ."B")THEN
                DO I=1,NB
                    DIST = DSQRT((COORDBOUNDARYPOINTS(I,1)-COORDNODES(CONECTIVITY(J,II),1))**2.D0+&
                        (COORDBOUNDARYPOINTS(I,2)-COORDNODES(CONECTIVITY(J,II),2))**2.D0)
                    IF(DIST.LT.1.D-10)THEN
                        POINT1 = I
                    ENDIF
                ENDDO
            ELSE
                DO I=1,NB
                    DIST = DSQRT((COORDBOUNDARYPOINTS(I,1)-COORDNODES(CONECTIVITY(J,II),1))**2.D0+&
                        (COORDBOUNDARYPOINTS(I,2)-COORDNODES(CONECTIVITY(J,II),2))**2.D0)
                    IF(DIST.LT.1.D-10)THEN
                        POINT1 = I
                    ENDIF
                ENDDO
                DO I=1,NB
                    DIST = DSQRT((COORDBOUNDARYPOINTS(I,1)-COORDNODES(CONECTIVITY(J,II),1))**2.D0+&
                        (COORDBOUNDARYPOINTS(I,2)-COORDNODES(CONECTIVITY(J,II),2))**2.D0)
                    IF(I.NE.POINT1.AND.DIST.LT.1.D-10)THEN
                        POINT2 = I
                    ENDIF
                ENDDO
            ENDIF
            IF(COLLOCATION_TYPE(CONECTIVITY(J,II)).EQ."B")THEN
                SIGMAB(4*POINT1-3) = SIGMA(II,1,1)
                SIGMAB(4*POINT1-2) = SIGMA(II,1,2)
                SIGMAB(4*POINT1-1) = SIGMA(II,2,1)
                SIGMAB(4*POINT1)   = SIGMA(II,2,2)
            ELSE
                SIGMAB(4*POINT1-3) = SIGMA(II,1,1)
                SIGMAB(4*POINT1-2) = SIGMA(II,1,2)
                SIGMAB(4*POINT1-1) = SIGMA(II,2,1)
                SIGMAB(4*POINT1)   = SIGMA(II,2,2)

                SIGMAB(4*POINT2-3) = SIGMA(II,1,1)
                SIGMAB(4*POINT2-2) = SIGMA(II,1,2)
                SIGMAB(4*POINT2-1) = SIGMA(II,2,1)
                SIGMAB(4*POINT2)   = SIGMA(II,2,2)
            ENDIF
        ENDDO
        DEALLOCATE(U_BAR,T_BAR,X_BAR,EPSILONXX_BAR,SIGMAXX_BAR,SIGMAXY_BAR,SIGMAYY_BAR,SIGMA_BAR,SIGMA)
        DEALLOCATE(VALUES,dPHI)
    ENDDO

    DO II=1,NB
        LONGITUDINAL_MODULUS = Young_Modulus(SUBREGION_POINTS(II))
        Nu = Poisson(SUBREGION_POINTS(II))
        SHEAR_MODULUS = LONGITUDINAL_MODULUS/(2.D0*(1.D0+Nu))

        SIGMA3D(II,1,1) = SIGMAB(4*II-3)
        SIGMA3D(II,1,2) = SIGMAB(4*II-2)
        SIGMA3D(II,2,1) = SIGMAB(4*II-1)
        SIGMA3D(II,2,2) = SIGMAB(4*II)
        IF(PLANE_STATE.EQ."PD")THEN
            SIGMA3D(II,3,3) = Nu*(SIGMA3D(II,1,1)+SIGMA3D(II,2,2))
        ENDIF

    !   GENERALIZED HOOKE'S LAW FOR THE COMPUTATION OF THE STRAIN TENSOR
        DO I=1,3
            DO J=1,3
                EPSILON3D(II,I,J) = (1.D0/(2.D0*SHEAR_MODULUS))*SIGMA3D(II,I,J)-(Nu/LONGITUDINAL_MODULUS)*&
                    (SIGMA3D(II,1,1)+SIGMA3D(II,2,2)+SIGMA3D(II,3,3))*IDENTITY(I,J)
            ENDDO
        ENDDO
    ENDDO

END SUBROUTINE POSTPROCESSING_BOUNDARY_STRESSES_AND_STRAINS

! ========================================================================================!
![5]
SUBROUTINE PRINCIPAL_STRESSES
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    IMPLICIT NONE
    INTEGER::I
    REAL(KIND=wp)::S1plusS2_2,S1minusS2_2

    ALLOCATE(SIGMAP(N_POINTS,2))
    SIGMAP = 0.D0

    DO I=1,N_POINTS
        S1plusS2_2  = (SIGMA3D(I,1,1)+SIGMA3D(I,2,2))/2.D0
        S1minusS2_2 = (SIGMA3D(I,1,1)-SIGMA3D(I,2,2))/2.D0

        SIGMAP(I,1) = S1plusS2_2 + DSQRT(S1minusS2_2**2.D0+SIGMA3D(I,1,2)**2.D0)
        SIGMAP(I,2) = S1plusS2_2 - DSQRT(S1minusS2_2**2.D0+SIGMA3D(I,1,2)**2.D0)
    ENDDO

END SUBROUTINE PRINCIPAL_STRESSES

! ========================================================================================!
![6]
SUBROUTINE PARAVIEW_OUTPUT
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    IMPLICIT NONE
    INTEGER::I,J,II

    OPEN(2,FILE='OUTPUT.vtk',STATUS='UNKNOWN')
    WRITE(2,1)'vtk DataFile Version 2.0'
    WRITE(2,'(a)')'Animation data'
    WRITE(2,'(a)')'ASCII'
    WRITE(2,'(a)')
    WRITE(2,'(a)')'DATASET UNSTRUCTURED_GRID'

    WRITE(2,'(a,1x,i0,1x,a)')'POINTS ',N_POINTS,' float'
    DO I=1,NB
        WRITE(2,'(3(f0.10,1x))') COORDBOUNDARYPOINTS(I,1),COORDBOUNDARYPOINTS(I,2),0.00000000000000D0
    END DO
    DO I=1,NI
        WRITE(2,'(3(f0.10,1x))') COORDINTERNALPOINTS(I,1),COORDINTERNALPOINTS(I,2),0.00000000000000D0
    END DO
    WRITE(2,'(a)')

    WRITE(2,'(a,1x,i0,1x,i0)')'CELLS ',N_POSPROCESSINGELEM,5*N_POSPROCESSINGELEM
    DO I=1,N_POSPROCESSINGELEM
        WRITE(2,'(100(I0,1X))') 3,(POSTPROCESSING_CONECTIVITY(I,J)-1,J=1,3)
    ENDDO
    WRITE(2,'(a)')

    WRITE(2,'(a,1x,i0)')'CELL_TYPES ',N_POSPROCESSINGELEM
    DO I=1,N_POSPROCESSINGELEM
        WRITE(2,'(i0)') 69
    ENDDO
    WRITE(2,'(a)')

    WRITE(2,'(a,1x,i0)')'POINT_DATA ',N_POINTS
    WRITE(2,'(a,1x,i0)')'SCALARS DISPLACEMENT FLOAT ',3
    WRITE(2,'(a)')'LOOKUP_TABLE default'
    DO I=1,NB
        WRITE(2,'(3(f0.10,1x))') U(2*I-1),U(2*I),0.D0
    END DO
    DO I=1,NI
        II = 2*I - 1
        WRITE(2,'(3(f0.10,1x))') UI(II),UI(II+1),0.D0
    END DO
    WRITE(2,'(a)')

    WRITE(2,'(a,1x,i0)')'SCALARS STRESS FLOAT ',6
    WRITE(2,'(a)')'LOOKUP_TABLE default'
    DO I=1,N_POINTS
        WRITE(2,'(6(f0.10,1x))') SIGMA3D(I,1,1),SIGMA3D(I,2,2),SIGMA3D(I,3,3),SIGMA3D(I,1,2),0.D0,0.D0
    END DO
    WRITE(2,'(a)')

    WRITE(2,'(a,1x,i0)')'SCALARS STRAINS FLOAT ',6
    WRITE(2,'(a)')'LOOKUP_TABLE default'
    DO I=1,N_POINTS
        WRITE(2,'(6(f0.10,1x))') EPSILON3D(I,1,1),EPSILON3D(I,2,2),EPSILON3D(I,3,3),EPSILON3D(I,1,2),0.D0,0.D0
    END DO
    WRITE(2,'(a)')

    WRITE(2,'(a,1x,i0)')'SCALARS PRINCIPAL STRESS FLOAT ',2
    WRITE(2,'(a)')'LOOKUP_TABLE default'
    DO I=1,N_POINTS
        WRITE(2,'(2(f0.10,1x))') SIGMAP(I,1),SIGMAP(I,2)
    END DO
    WRITE(2,'(a)')

    CLOSE(2)

1   FORMAT('#',1x,a)

END SUBROUTINE PARAVIEW_OUTPUT

! ========================================================================================!
![7]
SUBROUTINE BOUNDARY_OUTPUT
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    IMPLICIT NONE

    INTEGER::I,J,II,JJ,N_LOC
    REAL(KIND=8)::XVALUES(MAX_ORDER+1,2),UVALUES(MAX_ORDER+1,2),TVALUES(MAX_ORDER+1,2),PHI(MAX_ORDER+1)
    REAL(KIND=wp)::U_COLLOCPOINTS(2*N_NODES),T_COLLOCPOINTS(2*N_NODES),X_COLLOCPOINTS(N_NODES,2),X(N_NODES,2)

    X_COLLOCPOINTS = COORDCOLOCATIONPOINTS
    U_COLLOCPOINTS = U
    T_COLLOCPOINTS = T
    U = 0.D0
    T = 0.D0

    DO J=1,N_ELEM
        N_LOC = ORDER(J) + 1
        DO II=1,N_LOC
            XVALUES(II,1) = X_COLLOCPOINTS(CONECTIVITY(J,II),1)
            XVALUES(II,2) = X_COLLOCPOINTS(CONECTIVITY(J,II),2)
            UVALUES(II,1) = U_COLLOCPOINTS(2*CONECTIVITY(J,II)-1)
            UVALUES(II,2) = U_COLLOCPOINTS(2*CONECTIVITY(J,II))
            TVALUES(II,1) = T_COLLOCPOINTS(2*CONECTIVITY(J,II)-1)
            TVALUES(II,2) = T_COLLOCPOINTS(2*CONECTIVITY(J,II))
        ENDDO
        DO II=1,N_LOC
            CALL BASIS_FUNCTIONS(XSI_NODES(J,II),J,PHI(1:N_LOC))
            X(CONECTIVITY(J,II),1)   = 0.D0
            X(CONECTIVITY(J,II),2)   = 0.D0
            U(2*CONECTIVITY(J,II)-1) = 0.D0
            U(2*CONECTIVITY(J,II))   = 0.D0
            T(2*CONECTIVITY(J,II)-1) = 0.D0
            T(2*CONECTIVITY(J,II))   = 0.D0
            DO JJ=1,N_LOC
                X(CONECTIVITY(J,II),1)   = X(CONECTIVITY(J,II),1)   + PHI(JJ)*XVALUES(JJ,1)
                X(CONECTIVITY(J,II),2)   = X(CONECTIVITY(J,II),2)   + PHI(JJ)*XVALUES(JJ,2)
                U(2*CONECTIVITY(J,II)-1) = U(2*CONECTIVITY(J,II)-1) + PHI(JJ)*UVALUES(JJ,1)
                U(2*CONECTIVITY(J,II))   = U(2*CONECTIVITY(J,II))   + PHI(JJ)*UVALUES(JJ,2)
                T(2*CONECTIVITY(J,II)-1) = T(2*CONECTIVITY(J,II)-1) + PHI(JJ)*TVALUES(JJ,1)
                T(2*CONECTIVITY(J,II))   = T(2*CONECTIVITY(J,II))   + PHI(JJ)*TVALUES(JJ,2)
            ENDDO
        ENDDO
    ENDDO

    OPEN(10,FILE='11_OUTPUT_BOUNDARY_STRUCTURAL_RESPONSE.txt',STATUS='UNKNOWN')
    WRITE(10,*) "     X                  Y                  Ux                  Uy                  Tx           Ty"
    DO I=1,N_NODES
        WRITE(10,100) X(I,1),X(I,2),U(2*I-1),U(2*I),T(2*I-1),T(2*I)
    ENDDO
    100 FORMAT(F14.9,16x,F14.9,16x,F14.9,16x,F14.9,16x,F14.9,16x,F14.9)
    CLOSE(10)

END SUBROUTINE BOUNDARY_OUTPUT

! ========================================================================================!
![8]
SUBROUTINE POSTPROCESSING
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    IMPLICIT NONE

    IF(N_POSPROCESSINGELEM.GT.0)THEN
        CALL POSTPROCESSING_INTERNAL_DISPLACEMENTS
        CALL POSTPROCESSING_INTERNAL_STRESSES_AND_STRAINS
        CALL POSTPROCESSING_BOUNDARY_DISPLACEMENTS_AND_TRACTIONS
        CALL POSTPROCESSING_BOUNDARY_STRESSES_AND_STRAINS
        CALL PRINCIPAL_STRESSES
        CALL PARAVIEW_OUTPUT
    ENDIF

    CALL BOUNDARY_OUTPUT

END SUBROUTINE POSTPROCESSING
!
END MODULE MODULE_POSTPROCESSING
