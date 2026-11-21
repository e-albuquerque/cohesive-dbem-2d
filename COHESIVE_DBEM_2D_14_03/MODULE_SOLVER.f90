MODULE MODULE_SOLVER
! -------------------------------------------------------------------------------------!
!   Description.                                                                       !
!                                                                                      !
! -------------------------------------------------------------------------------------!
!   This module contains the following subroutines:                                    !
!                                                                                      !
!   [1]. MATRIX_M                                                                      !
!   [2]. VECTOR_Fu_PRES                                                                !
!   [3]. VECTOR_Ft_PRES                                                                !
!   [4]. TANGENT_MATRIX_dPSIdX(X)                                                      !
!   [5]. UPDATE_COHESIVE_STIFFNESS(X)                                                  !
!   [6]. UPDATE_U_AND_T(X)                                                             !
!   [7]. SOLVE_LINEAR_SYSTEM                                                           !
!   [8]. SOLVE_NONLINEAR_SYSTEM                                                        !
! -------------------------------------------------------------------------------------!

USE MODULE_DATA
USE omp_lib
USE MODULE_GEOMETRIC_MAPPINGS

CONTAINS

! ========================================================================================!
![1]                                                                                      !
SUBROUTINE MATRIX_M                                                                       !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
IMPLICIT NONE

INTEGER::I,I_AUX,K,II,JJ,KK,COL,NB_GDL_ON,NB_GDL_OFF,N_IM,N_IP,IND_I,IND_F
REAL(KIND=wp)::K11,K12,K21,K22,JAC,nv(2)
REAL(KIND=wp)::VALUES(MAX_ORDER+1,2)

MATRIX = 0.D0

NB_GDL_ON  = SIZE(VETOR_B_GDL_ON)
NB_GDL_OFF = SIZE(VETOR_B_GDL_OFF)
N_IM       = SIZE(VETOR_NODE_Im)
N_IP       = SIZE(VETOR_NODE_Ip)

IND_I = 1
IND_F = NB_GDL_OFF
MATRIX(:,IND_I:IND_F) = H(:,VETOR_B_GDL_OFF)

IND_I = NB_GDL_OFF + 1
IND_F = NB_GDL_ON + NB_GDL_OFF
MATRIX(:,IND_I:IND_F) = -G(:,VETOR_B_GDL_ON)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I_AUX, I, JJ, COL)
DO I_AUX=1,N_IM
    I = VETOR_NODE_Im(I_AUX)
    JJ = COHESIVE_INTERFACES_Im_Ip(I,1)
    COL = NB_GDL_ON + NB_GDL_OFF + 2*I_AUX - 1
    MATRIX(:,COL)   = H(:,2*I-1) + H(:,2*JJ-1)
    MATRIX(:,COL+1) = H(:,2*I)   + H(:,2*JJ)
ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I_AUX, I, II, JJ, KK, K, COL, VALUES, JAC, nv, K11, K12, K21, K22)
DO I_AUX=1,N_IP
    I  = VET_Ip_ELEM_ORDER(I_AUX,1)
    JJ = VET_Ip_ELEM_ORDER(I_AUX,2)
    KK = VET_Ip_ELEM_ORDER(I_AUX,3)
    II = COHESIVE_INTERFACES_Ip_Im(I,1)

    DO K=1,ORDER(JJ)+1
        VALUES(K,1) = COORDNODES(CONECTIVITY(JJ,K),1)
        VALUES(K,2) = COORDNODES(CONECTIVITY(JJ,K),2)
    ENDDO

    CALL JACOBIAN_AND_NORMAL_VECTOR (XSI_COLLOCATIONPOINTS(JJ,KK),JJ,VALUES(1:ORDER(JJ)+1,:),JAC,nv)

    K11 = nv(1)**2.D0*Knn(I) - nv(1)*nv(2)*Ktn(I) - nv(1)*nv(2)*Knt(I) + nv(2)**2.D0*Ktt(I)
    K12 = nv(1)*nv(2)*Knn(I) - nv(2)**2.D0*Ktn(I) + nv(1)**2.D0*Knt(I) - nv(1)*nv(2)*Ktt(I)
    K21 = nv(1)*nv(2)*Knn(I) + nv(1)**2.D0*Ktn(I) - nv(2)**2.D0*Knt(I) - nv(1)*nv(2)*Ktt(I)
    K22 = nv(2)**2.D0*Knn(I) + nv(1)*nv(2)*Ktn(I) + nv(1)*nv(2)*Knt(I) + nv(1)**2.D0*Ktt(I)

    COL = NB_GDL_ON + NB_GDL_OFF + 2*N_IM + 2*I_AUX - 1
    MATRIX(:,COL)   = H(:,2*I-1) + (G(:,2*I-1)-G(:,2*II-1))*K11 + (G(:,2*I)-G(:,2*II))*K21
    MATRIX(:,COL+1) = H(:,2*I)   + (G(:,2*I-1)-G(:,2*II-1))*K12 + (G(:,2*I)-G(:,2*II))*K22
ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE MATRIX_M

! ========================================================================================!
![2]                                                                                      !
SUBROUTINE VECTOR_Fu_PRES                                                                 !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
IMPLICIT NONE

INTEGER::I,DOF,NB_GDL_ON
REAL(KIND=wp)::U_VAL

Fu_PRES = 0.D0
NB_GDL_ON = SIZE(VETOR_B_GDL_ON)

DO I=1,NB_GDL_ON
    DOF = VETOR_B_GDL_ON(I)
    U_VAL = U(DOF)
    IF(U_VAL.NE.0.D0)THEN
        Fu_PRES(:) = Fu_PRES(:) + H(:,DOF)*U_VAL
    ENDIF
ENDDO

END SUBROUTINE VECTOR_Fu_PRES

! ========================================================================================!
![3]                                                                                      !
SUBROUTINE VECTOR_Ft_PRES                                                                 !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
IMPLICIT NONE

INTEGER::I,DOF,NB_GDL_OFF
REAL(KIND=wp)::T_VAL

Ft_PRES = 0.D0
NB_GDL_OFF = SIZE(VETOR_B_GDL_OFF)

DO I=1,NB_GDL_OFF
    DOF = VETOR_B_GDL_OFF(I)
    T_VAL = T(DOF)
    IF(T_VAL.NE.0.D0)THEN
        Ft_PRES(:) = Ft_PRES(:) + G(:,DOF)*T_VAL
    ENDIF
ENDDO

END SUBROUTINE VECTOR_Ft_PRES

! ========================================================================================!
![4]                                                                                      !
SUBROUTINE TANGENT_MATRIX_dPSIdX(X)                                                       !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
IMPLICIT NONE

INTEGER::I,J,K,II,JJ,KK,CONT,N0,I_AUX
REAL(KIND=wp)::VALUES(MAX_ORDER+1,2)
REAL(KIND=wp)::X(2*N_NODES),dK11dDUx(N_NODES),dK11dDUy(N_NODES),&
    dK12dDUx(N_NODES),dK12dDUy(N_NODES),dK21dDUx(N_NODES),dK21dDUy(N_NODES),dK22dDUx(N_NODES),&
    dK22dDUy(N_NODES),DU(2),nv(2),JAC,DUn,DUt
REAL(KIND=wp)::dK11dDUn,dK11dDUt,dK12dDUn,dK12dDUt,dK21dDUn,dK21dDUt,dK22dDUn,dK22dDUt
REAL(KIND=wp)::dKnndDUn,dKnndDUt,dKntdDUn,dKntdDUt,dKtndDUn,dKtndDUt,dKttdDUn,dKttdDUt
REAL(KIND=wp)::Limit_Tn,Limit_Tt,av,bv
REAL(KIND=wp)::Kappa,dBetadDUn,dKappadDUn,dBetadDUt,dKappadDUt,Sigma_n,Sigma_t,Criterion

CALL MATRIX_M
dPSIdX = MATRIX
dK11dDUx = 0.D0; dK11dDUy = 0.D0
dK12dDUx = 0.D0; dK12dDUy = 0.D0
dK21dDUx = 0.D0; dK21dDUy = 0.D0
dK22dDUx = 0.D0; dK22dDUy = 0.D0
dKnndDUn = 0.D0; dKnndDUt = 0.D0
dKntdDUn = 0.D0; dKntdDUt = 0.D0
dKtndDUn = 0.D0; dKtndDUt = 0.D0
dKttdDUn = 0.D0; dKttdDUt = 0.D0

Limit_Tn = (1.D-3)*smax
Limit_Tt = (1.D-3)*tmax
av = 1.D0
bv = 1.D0

N0 = NDOF_T + NDOF_U + NDOF_DU
CONT = N0
DO I_AUX=1,SIZE(VETOR_NODE_Ip)
    I = VETOR_NODE_Ip(I_AUX)
    CONT = CONT + 1
    DU(1) = X(CONT)
    CONT = CONT + 1
    DU(2) = X(CONT)

    II = COHESIVE_INTERFACES_Ip_Im(I,1)
    JJ = COHESIVE_INTERFACES_Ip_Im(I,2)
    KK = COHESIVE_INTERFACES_Ip_Im(I,3)

    DO K=1,ORDER(JJ)+1
        VALUES(K,1) = COORDNODES(CONECTIVITY(JJ,K),1)
        VALUES(K,2) = COORDNODES(CONECTIVITY(JJ,K),2)
    ENDDO
    CALL JACOBIAN_AND_NORMAL_VECTOR (XSI_COLLOCATIONPOINTS(JJ,KK),JJ,VALUES(1:ORDER(JJ)+1,:),JAC,nv)
    DUn = nv(1)*DU(1) + nv(2)*DU(2)
    IF(ABS(DUn).LE.1.D-12)THEN
        DUn = 1.D-12
    ENDIF
    DUt = -nv(2)*DU(1) + nv(1)*DU(2)
    IF(ABS(DUt).LE.1.D-12)THEN
        DUt = 1.D-12
    ENDIF

    Sigma_n = kn0*DUn
    Sigma_t = kt0*DUt
    Criterion = DSQRT((Sigma_n/smax)**2.D0 + (Sigma_t/tmax)**2.D0)
    Kappa = DSQRT((DUn/deltanc)**2.D0 + (DUt/deltatc)**2.D0) - 1.D0
    Beta(I) = Kappa/(Etan*(1.D0+Kappa))
    dKappadDUn = DUn/((deltanc**2.D0)*Kappa)
    dBetadDUn  = dKappadDUn/(Etan*(Kappa+1.D0)**2.D0)
    dKappadDUt = DUt/((deltatc**2.D0)*Kappa)
    dBetadDUt  = dKappadDUt/(Etan*(Kappa+1.D0)**2.D0)

    IF(Criterion.LT.1.D0)THEN
        dKnndDUn = Limit_Tn/(DUn**2.D0)
        dKnndDUt = 0.D0; dKntdDUn = 0.D0; dKntdDUt = 0.D0
        dKtndDUn = 0.D0; dKtndDUt = 0.D0; dKttdDUn = 0.D0
        dKttdDUt = Limit_Tt/(ABS(DUt)**2.D0)
    ELSE IF(Criterion.GE.1.D0.AND.Beta(I).LT.1.D0.AND.Beta(I).GE.Beta_MAX(I))THEN
        dKnndDUn = -dBetadDUn*(1.D0/deltanc)*smax
        dKnndDUt = -dBetadDUt*(1.D0/deltanc)*smax
        dKntdDUn = 0.D0; dKntdDUt = 0.D0
        dKtndDUn = 0.D0; dKtndDUt = 0.D0
        dKttdDUn = -dBetadDUn*(1.D0/deltatc)*tmax
        dKttdDUt = -dBetadDUt*(1.D0/deltatc)*tmax
    ELSE IF(Criterion.GE.1.D0.AND.Beta(I).LT.1.D0.AND.Beta(I).LT.Beta_MAX(I))THEN
        dKnndDUn = Limit_Tn/(DUn**2.D0)
        dKnndDUt = 0.D0; dKntdDUn = 0.D0; dKntdDUt = 0.D0
        dKtndDUn = 0.D0; dKtndDUt = 0.D0; dKttdDUn = 0.D0
        dKttdDUt = Limit_Tt/(ABS(DUt)**2.D0)
    ELSE IF(Criterion.GE.1.D0.AND.Beta(I).GE.1.D0)THEN
        dKnndDUn = -Limit_Tn/(DUn**2.D0)
        dKnndDUt = 0.D0; dKntdDUn = 0.D0; dKntdDUt = 0.D0
        dKtndDUn = 0.D0; dKtndDUt = 0.D0; dKttdDUn = 0.D0
        dKttdDUt = -Limit_Tt/(DUt**2.D0)
    ENDIF

    dK11dDUn = nv(1)**2.D0*dKnndDUn - nv(1)*nv(2)*dKtndDUn - nv(1)*nv(2)*dKntdDUn + nv(2)**2.D0*dKttdDUn
    dK11dDUt = nv(1)**2.D0*dKnndDUt - nv(1)*nv(2)*dKtndDUt - nv(1)*nv(2)*dKntdDUt + nv(2)**2.D0*dKttdDUt
    dK21dDUn = nv(1)*nv(2)*dKnndDUn + nv(1)**2.D0*dKtndDUn - nv(2)**2.D0*dKntdDUn - nv(1)*nv(2)*dKttdDUn
    dK21dDUt = nv(1)*nv(2)*dKnndDUt + nv(1)**2.D0*dKtndDUt - nv(2)**2.D0*dKntdDUt - nv(1)*nv(2)*dKttdDUt
    dK12dDUn = dK21dDUn
    dK12dDUt = dK21dDUt
    dK22dDUn = nv(2)**2.D0*dKnndDUn + nv(1)*nv(2)*dKtndDUn + nv(1)*nv(2)*dKntdDUn + nv(1)**2.D0*dKttdDUn
    dK22dDUt = nv(2)**2.D0*dKnndDUt + nv(1)*nv(2)*dKtndDUt + nv(1)*nv(2)*dKntdDUt + nv(1)**2.D0*dKttdDUt

    dK11dDUx(I) = dK11dDUn*nv(1) - dK11dDUt*nv(2)
    dK11dDUy(I) = dK11dDUn*nv(2) + dK11dDUt*nv(1)
    dK21dDUx(I) = dK21dDUn*nv(1) - dK21dDUt*nv(2)
    dK21dDUy(I) = dK21dDUn*nv(2) + dK21dDUt*nv(1)
    dK12dDUx(I) = dK12dDUn*nv(1) - dK12dDUt*nv(2)
    dK12dDUy(I) = dK12dDUn*nv(2) + dK12dDUt*nv(1)
    dK22dDUx(I) = dK22dDUn*nv(1) - dK22dDUt*nv(2)
    dK22dDUy(I) = dK22dDUn*nv(2) + dK22dDUt*nv(1)

    nv = -nv

    dK11dDUn = nv(1)**2.D0*dKnndDUn - nv(1)*nv(2)*dKtndDUn - nv(1)*nv(2)*dKntdDUn + nv(2)**2.D0*dKttdDUn
    dK11dDUt = nv(1)**2.D0*dKnndDUt - nv(1)*nv(2)*dKtndDUt - nv(1)*nv(2)*dKntdDUt + nv(2)**2.D0*dKttdDUt
    dK21dDUn = nv(1)*nv(2)*dKnndDUn + nv(1)**2.D0*dKtndDUn - nv(2)**2.D0*dKntdDUn - nv(1)*nv(2)*dKttdDUn
    dK21dDUt = nv(1)*nv(2)*dKnndDUt + nv(1)**2.D0*dKtndDUt - nv(2)**2.D0*dKntdDUt - nv(1)*nv(2)*dKttdDUt
    dK12dDUn = dK21dDUn
    dK12dDUt = dK21dDUt
    dK22dDUn = nv(2)**2.D0*dKnndDUn + nv(1)*nv(2)*dKtndDUn + nv(1)*nv(2)*dKntdDUn + nv(1)**2.D0*dKttdDUn
    dK22dDUt = nv(2)**2.D0*dKnndDUt + nv(1)*nv(2)*dKtndDUt + nv(1)*nv(2)*dKntdDUt + nv(1)**2.D0*dKttdDUt

    dK11dDUx(II) = dK11dDUn*nv(1) - dK11dDUt*nv(2)
    dK11dDUy(II) = dK11dDUn*nv(2) + dK11dDUt*nv(1)
    dK21dDUx(II) = dK21dDUn*nv(1) - dK21dDUt*nv(2)
    dK21dDUy(II) = dK21dDUn*nv(2) + dK21dDUt*nv(1)
    dK12dDUx(II) = dK12dDUn*nv(1) - dK12dDUt*nv(2)
    dK12dDUy(II) = dK12dDUn*nv(2) + dK12dDUt*nv(1)
    dK22dDUx(II) = dK22dDUn*nv(1) - dK22dDUt*nv(2)
    dK22dDUy(II) = dK22dDUn*nv(2) + dK22dDUt*nv(1)
ENDDO

CONT = N0
DO I_AUX=1,SIZE(VETOR_NODE_Ip)
    I  = VETOR_NODE_Ip(I_AUX)
    II = COHESIVE_INTERFACES_Ip_Im(I,1)
    CONT = CONT + 1
    DO J=N0,N0+NDOF_DU
        IF(MOD(J,2).GT.0)THEN
            MATRIX(:,CONT) = MATRIX(:,CONT) + X(J)*((G(:,2*I-1)-G(:,2*II-1))*dK11dDUx(I)+(G(:,2*I)-G(:,2*II))*dK21dDUx(I))
            MATRIX(:,CONT) = MATRIX(:,CONT) + X(J+1)*((G(:,2*I-1)-G(:,2*II-1))*dK12dDUx(I)+(G(:,2*I)-G(:,2*II))*dK22dDUx(I))
        ENDIF
    ENDDO
    CONT = CONT + 1
    DO J=N0,N0+NDOF_DU
        IF(MOD(J,2).GT.0)THEN
            MATRIX(:,CONT) = MATRIX(:,CONT) + X(J)*((G(:,2*I-1)-G(:,2*II-1))*dK11dDUy(I)+(G(:,2*I)-G(:,2*II))*dK21dDUy(I))
            MATRIX(:,CONT) = MATRIX(:,CONT) + X(J+1)*((G(:,2*I-1)-G(:,2*II-1))*dK12dDUy(I)+(G(:,2*I)-G(:,2*II))*dK22dDUy(I))
        ENDIF
    ENDDO
ENDDO

END SUBROUTINE TANGENT_MATRIX_dPSIdX

! ========================================================================================!
![5]                                                                                      !
SUBROUTINE UPDATE_COHESIVE_STIFFNESS(X)                                                   !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
IMPLICIT NONE

INTEGER::I,K,II,JJ,KK,CONT,N0,I_AUX
REAL(KIND=wp)::VALUES(MAX_ORDER+1,2)
REAL(KIND=wp)::X(2*N_NODES),DU(2),nv(2),JAC,DUn,DUt,Limit_Tn,Limit_Tt,av,bv
REAL(KIND=wp)::Kappa,Sigma_n,Sigma_t,Criterion
REAL(KIND=wp)::X_MIN_INT,X_MAX_INT,X_ASYM_LIMIT,smax_loc,tmax_loc

Limit_Tn = (1.D-3)*smax
Limit_Tt = (1.D-3)*tmax
av = 1.D0
bv = 1.D0

X_MIN_INT = MINVAL(COORDCOLOCATIONPOINTS(VETOR_NODE_Ip,1))
X_MAX_INT = MAXVAL(COORDCOLOCATIONPOINTS(VETOR_NODE_Ip,1))
X_ASYM_LIMIT = X_MAX_INT - 0.10D0*(X_MAX_INT - X_MIN_INT)

N0 = NDOF_T + NDOF_U + NDOF_DU
CONT = N0
DO I_AUX=1,SIZE(VETOR_NODE_Ip)
    I = VETOR_NODE_Ip(I_AUX)
    CONT = CONT + 1
    DU(1) = X(CONT)
    CONT = CONT + 1
    DU(2) = X(CONT)

    II = COHESIVE_INTERFACES_Ip_Im(I,1)
    JJ = COHESIVE_INTERFACES_Ip_Im(I,2)
    KK = COHESIVE_INTERFACES_Ip_Im(I,3)

    DO K=1,ORDER(JJ)+1
        VALUES(K,1) = COORDNODES(CONECTIVITY(JJ,K),1)
        VALUES(K,2) = COORDNODES(CONECTIVITY(JJ,K),2)
    ENDDO
    CALL JACOBIAN_AND_NORMAL_VECTOR(XSI_COLLOCATIONPOINTS(JJ,KK),JJ,VALUES(1:ORDER(JJ)+1,:),JAC,nv)
    DUn = nv(1)*DU(1) + nv(2)*DU(2)
    IF(ABS(DUn).LE.1.D-12)THEN
        DUn = 1.D-12
    ENDIF
    DUt = -nv(2)*DU(1) + nv(1)*DU(2)
    IF(ABS(DUt).LE.1.D-12)THEN
        DUt = 1.D-12
    ENDIF

    smax_loc = smax
    tmax_loc = tmax
    IF(ANALYSIS_TYPE.EQ."CA" .AND. COORDCOLOCATIONPOINTS(I,1).GT.X_ASYM_LIMIT)THEN
        smax_loc = 1.0D50
        tmax_loc = 1.0D50
    ENDIF

    Sigma_n = kn0*DUn
    Sigma_t = kt0*DUt
    Criterion = DSQRT((Sigma_n/smax_loc)**2.D0 + (Sigma_t/tmax_loc)**2.D0)
    Kappa = DSQRT((DUn/deltanc)**2.D0 + (DUt/deltatc)**2.D0) - 1.D0
    Beta(I) = Kappa/(Etan*(1.D0+Kappa))

    IF(Criterion.LT.1.D0)THEN
        Knn(I)  = Kn0; Knn(II) = Kn0
        Knt(I)  = 0.D0; Knt(II) = 0.D0
        Ktn(I)  = 0.D0; Ktn(II) = 0.D0
        Ktt(I)  = Kt0; Ktt(II) = Kt0
    ELSE IF(Criterion.GE.1.D0.AND.Beta(I).LT.1.D0.AND.Beta(I).GE.Beta_MAX(I))THEN
        Knn(I)  = (1.D0-Beta(I))*(1.D0/deltanc)*smax
        Knn(II) = Knn(I)
        Knt(I)  = 0.D0; Knt(II) = 0.D0
        Ktn(I)  = 0.D0; Ktn(II) = 0.D0
        Ktt(I)  = (1.D0-Beta(I))*(1.D0/deltatc)*tmax
        Ktt(II) = Ktt(I)
    ELSE IF(Criterion.GE.1.D0.AND.Beta(I).LT.1.D0.AND.Beta(I).LT.Beta_MAX(I))THEN
        Knn(I)  = (1.D0-Beta_MAX(I))*(1.D0/deltanc)*smax
        Knn(II) = Knn(I)
        Knt(I)  = 0.D0; Knt(II) = 0.D0
        Ktn(I)  = 0.D0; Ktn(II) = 0.D0
        Ktt(I)  = (1.D0-Beta_MAX(I))*(1.D0/deltatc)*tmax
        Ktt(II) = Ktt(I)
    ELSE IF(Criterion.GE.1.D0.AND.Beta(I).GE.1.D0)THEN
        Knn(I)  = Limit_Tn/DUn
        Knn(II) = Knn(I)
        Knt(I)  = 0.D0; Knt(II) = 0.D0
        Ktn(I)  = 0.D0; Ktn(II) = 0.D0
        Ktt(I)  = Limit_Tt/ABS(DUt)
        Ktt(II) = Ktt(I)
    ENDIF

    DU_BAR(2*I-1)  = DUn
    DU_BAR(2*I)    = DUt
    DU_BAR(2*II-1) = DUn
    DU_BAR(2*II)   = DUt
ENDDO

END SUBROUTINE UPDATE_COHESIVE_STIFFNESS

! ========================================================================================!
![6]                                                                                      !
SUBROUTINE UPDATE_U_AND_T(X)                                                              !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
IMPLICIT NONE
INTEGER::I,K,JJ,KK,CONT,CM,I_AUX
REAL(KIND=wp)::VALUES(MAX_ORDER+1,2)
REAL(KIND=wp)::X(2*N_NODES),JAC,nv(2),K11,K12,K21,K22

CONT = 0
DO I=1,N_NODES
    IF(COLLOCATION_TYPE(I).EQ."B")THEN
        IF(BC(2*I-1).EQ.0)THEN
            CONT = CONT + 1
            U(2*I-1) = X(CONT)
        ENDIF
        IF(BC(2*I).EQ.0)THEN
            CONT = CONT + 1
            U(2*I) = X(CONT)
        ENDIF
    ENDIF
ENDDO

DO I=1,N_NODES
    IF(COLLOCATION_TYPE(I).EQ."B")THEN
        IF(BC(2*I-1).EQ.1)THEN
            CONT = CONT + 1
            T(2*I-1) = X(CONT)
        ENDIF
        IF(BC(2*I).EQ.1)THEN
            CONT = CONT + 1
            T(2*I) = X(CONT)
        ENDIF
    ENDIF
ENDDO

DO I=1,N_NODES
    IF(COLLOCATION_TYPE(I).EQ."Im")THEN
        CONT = CONT + 1
        U(2*I-1) = X(CONT)
        CONT = CONT + 1
        U(2*I) = X(CONT)
    ENDIF
ENDDO

DO I_AUX=1,SIZE(VETOR_NODE_Ip)
    I  = VET_Ip_ELEM_ORDER(I_AUX,1)
    JJ = VET_Ip_ELEM_ORDER(I_AUX,2)
    KK = VET_Ip_ELEM_ORDER(I_AUX,3)
    CM = COHESIVE_INTERFACES_Ip_Im(I,1)

    DO K=1,ORDER(JJ)+1
        VALUES(K,1) = COORDNODES(CONECTIVITY(JJ,K),1)
        VALUES(K,2) = COORDNODES(CONECTIVITY(JJ,K),2)
    ENDDO
    CALL JACOBIAN_AND_NORMAL_VECTOR(XSI_COLLOCATIONPOINTS(JJ,KK),JJ,VALUES(1:ORDER(JJ)+1,:),JAC,nv)

    K11 = nv(1)**2.D0*Knn(I) - nv(1)*nv(2)*Ktn(I) - nv(1)*nv(2)*Knt(I) + nv(2)**2.D0*Ktt(I)
    K12 = nv(1)*nv(2)*Knn(I) - nv(2)**2.D0*Ktn(I) + nv(1)**2.D0*Knt(I) - nv(1)*nv(2)*Ktt(I)
    K21 = nv(1)*nv(2)*Knn(I) + nv(1)**2.D0*Ktn(I) - nv(2)**2.D0*Knt(I) - nv(1)*nv(2)*Ktt(I)
    K22 = nv(2)**2.D0*Knn(I) + nv(1)*nv(2)*Ktn(I) + nv(1)*nv(2)*Knt(I) + nv(1)**2.D0*Ktt(I)

    CONT = CONT + 1
    U(2*I-1)  = U(2*CM-1) + X(CONT)
    T(2*CM-1) = K11*X(CONT) + K12*X(CONT+1)
    T(2*I-1)  = -T(2*CM-1)
    CONT = CONT + 1
    U(2*I)    = U(2*CM) + X(CONT)
    T(2*CM)   = K21*X(CONT-1) + K22*X(CONT)
    T(2*I)    = -T(2*CM)
ENDDO

END SUBROUTINE UPDATE_U_AND_T

! ========================================================================================!
![7]                                                                                      !
SUBROUTINE SOLVE_LINEAR_SYSTEM                                                            !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
IMPLICIT NONE

INTEGER::IPIV(2*N_NODES),INFO,I,I_AUX,K,II,JJ,KK,CONT
REAL(KIND=wp)::VALUES(MAX_ORDER+1,2)
REAL(KIND=wp),ALLOCATABLE::F(:)
REAL(KIND=wp)::JAC,nv(2),K11,K12,K21,K22

ALLOCATE(MATRIX(2*N_NODES,2*N_NODES),F(2*N_NODES),Fu_PRES(2*N_NODES),Ft_PRES(2*N_NODES))

CALL MATRIX_M
CALL VECTOR_Fu_PRES
CALL VECTOR_Ft_PRES
F = Ft_PRES - Fu_PRES

CALL DGESV(2*N_NODES,1,MATRIX,2*N_NODES,IPIV,F,2*N_NODES,INFO)

CONT = 0
DO I=1,N_NODES
    IF(COLLOCATION_TYPE(I).EQ."B")THEN
        IF(BC(2*I-1).EQ.0)THEN
            CONT = CONT + 1
            U(2*I-1) = F(CONT)
        ENDIF
        IF(BC(2*I).EQ.0)THEN
            CONT = CONT + 1
            U(2*I) = F(CONT)
        ENDIF
    ENDIF
ENDDO
DO I=1,N_NODES
    IF(COLLOCATION_TYPE(I).EQ."B")THEN
        IF(BC(2*I-1).EQ.1)THEN
            CONT = CONT + 1
            T(2*I-1) = F(CONT)
        ENDIF
        IF(BC(2*I).EQ.1)THEN
            CONT = CONT + 1
            T(2*I) = F(CONT)
        ENDIF
    ENDIF
ENDDO
DO I=1,N_NODES
    IF(COLLOCATION_TYPE(I).EQ."Im")THEN
        CONT = CONT + 1
        U(2*I-1) = F(CONT)
        CONT = CONT + 1
        U(2*I) = F(CONT)
    ENDIF
ENDDO
DO I_AUX=1,SIZE(VETOR_NODE_Ip)
    I  = VET_Ip_ELEM_ORDER(I_AUX,1)
    JJ = VET_Ip_ELEM_ORDER(I_AUX,2)
    KK = VET_Ip_ELEM_ORDER(I_AUX,3)
    II = COHESIVE_INTERFACES_Ip_Im(I,1)

    DO K=1,ORDER(JJ)+1
        VALUES(K,1) = COORDNODES(CONECTIVITY(JJ,K),1)
        VALUES(K,2) = COORDNODES(CONECTIVITY(JJ,K),2)
    ENDDO
    CALL JACOBIAN_AND_NORMAL_VECTOR(XSI_COLLOCATIONPOINTS(JJ,KK),JJ,VALUES(1:ORDER(JJ)+1,:),JAC,nv)

    K11 = nv(1)**2.D0*Knn(I) - nv(1)*nv(2)*Ktn(I) - nv(1)*nv(2)*Knt(I) + nv(2)**2.D0*Ktt(I)
    K12 = nv(1)*nv(2)*Knn(I) - nv(2)**2.D0*Ktn(I) + nv(1)**2.D0*Knt(I) - nv(1)*nv(2)*Ktt(I)
    K21 = nv(1)*nv(2)*Knn(I) + nv(1)**2.D0*Ktn(I) - nv(2)**2.D0*Knt(I) - nv(1)*nv(2)*Ktt(I)
    K22 = nv(2)**2.D0*Knn(I) + nv(1)*nv(2)*Ktn(I) + nv(1)*nv(2)*Knt(I) + nv(1)**2.D0*Ktt(I)

    CONT = CONT + 1
    U(2*I-1)  = U(2*II-1) + F(CONT)
    T(2*II-1) = K11*F(CONT) + K12*F(CONT+1)
    T(2*I-1)  = -T(2*II-1)
    CONT = CONT + 1
    U(2*I)    = U(2*II) + F(CONT)
    T(2*II)   = K21*F(CONT-1) + K22*F(CONT)
    T(2*I)    = -T(2*II)
ENDDO

DEALLOCATE(MATRIX,F,Fu_PRES,Ft_PRES)

END SUBROUTINE SOLVE_LINEAR_SYSTEM

! ========================================================================================!
![8]                                                                                      !
SUBROUTINE SOLVE_NONLINEAR_SYSTEM                                                         !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
IMPLICIT NONE

INTEGER::II,I,IPIV(2*N_NODES),INFO,ite,itecumulated
LOGICAL::MAT_CHANGED,LU_CHANGED
REAL(KIND=wp)::X(2*N_NODES),lambda,PSI(2*N_NODES),&
    F(2*N_NODES,2),dX(2*N_NODES),dXpres(2*N_NODES),dXc(2*N_NODES),dXp(2*N_NODES),criterionPSI_num,&
    criterionPSI_denum,criterionPSI,dLambda,lambdaux,lambdauy,lambdatx,lambdaty
REAL(KIND=wp)::Knn_OLD(N_NODES),Ktt_OLD(N_NODES)
REAL(KIND=wp)::TI_FuFt,TF_FuFt,DT_FuFt
REAL(KIND=wp)::TI_F,TF_F,DT_F
REAL(KIND=wp)::TI_TM,TF_TM,DT_TM
REAL(KIND=wp)::TI_SS,TF_SS,DT_SS
REAL(KIND=wp)::TI_CS,TF_CS,DT_CS
REAL(KIND=wp)::TI_PSI,TF_PSI,DT_PSI

ALLOCATE(dPSIdX(2*N_NODES,2*N_NODES),MATRIX(2*N_NODES,2*N_NODES),Fu_PRES(2*N_NODES),Ft_PRES(2*N_NODES))

X = 0.D0
X(CONTROLED_DOF) = CONTROLED_VALUE
dXpres = X / Ns
X = 0.D0
lambda = 0.D0
lambdaux = 1.D0
lambdauy = 1.D0
lambdatx = 1.D0
lambdaty = 1.D0
itecumulated = 0
MAT_CHANGED = .TRUE.
LU_CHANGED  = .TRUE.

! Precompute constant prescribed load vectors outside the load-step loop
TI_FuFt = omp_get_wtime()
CALL VECTOR_Fu_PRES
CALL VECTOR_Ft_PRES
TF_FuFt = omp_get_wtime()
DT_FuFt = TF_FuFt - TI_FuFt

OPEN(10,FILE='8_OUTPUT_NODAL_STRUCTURAL_RESPONSE.txt',STATUS='UNKNOWN')
WRITE(10,*) "       Ux              Uy              Tx              Ty"
OPEN(20,FILE='9_OUTPUT_COHESIVE_RESPONSE.txt',STATUS='UNKNOWN')
WRITE(20,*) "       Knn          Knt          Ktn          Ktt          Dun          Dut"
OPEN(30,FILE='10_OUTPUT_COHESIVE_STIFFNESS.txt',STATUS='UNKNOWN')
WRITE(30,*) "       STEP          Knn          Ktt"
OPEN(50,FILE='12_OUTPUT_NODAL_DISPLACEMENTS_STEPS.txt',STATUS='UNKNOWN')
WRITE(50,*) "       STEP          NODE          Ux                Uy"
WRITE(40,*) "STEP   ITER     Fu-Ft_V             F_V                T_M               S_S                UP_C_S     PSI_V"
DO II=1,Ns
    criterionPSI = 1.D0
    ite = 0
    X(CONTROLED_DOF) = X(CONTROLED_DOF) + dXpres(CONTROLED_DOF)

    DO WHILE (criterionPSI.GT.tol_Newton)
        ite = ite + 1
        itecumulated = itecumulated + 1
        IF(ite.GT.itemax)THEN
            tol_Newton = tol_Newton*10.D0
            IF(tol_Newton.GT.1.D-1)THEN
                STOP
            ENDIF
        ENDIF

!   COMPUTING THE dRdlambda AND RESIDUE VECTORS
        TI_F = omp_get_wtime()
        IF(MAT_CHANGED)THEN
            CALL MATRIX_M
            MAT_CHANGED = .FALSE.
        ENDIF
        F(:,1) = Fu_PRES - Ft_PRES
        F(:,2) = MATMUL(MATRIX,X) + lambda*F(:,1)
        TF_F = omp_get_wtime()
        DT_F = TF_F - TI_F

!   COMPUTING THE TANGENT/SECANT OPERATOR MATRIX AND LU FACTORIZATION (ONLY WHEN STIFFNESS CHANGED SIGNIFICANTLY)
        TI_TM = omp_get_wtime()
        IF(LU_CHANGED)THEN
            dPSIdX = MATRIX
            CALL DGETRF(2*N_NODES,2*N_NODES,dPSIdX,2*N_NODES,IPIV,INFO)
            LU_CHANGED = .FALSE.
        ENDIF
        TF_TM = omp_get_wtime()
        DT_TM = TF_TM - TI_TM

!   SOLVING THE PREDICTION AND CORRECTOR SYSTEM OF EQUATIONS VIA FORWARD/BACKWARD SUBSTITUTION
        TI_SS = omp_get_wtime()
        CALL DGETRS('N',2*N_NODES,2,dPSIdX,2*N_NODES,IPIV,F,2*N_NODES,INFO)
        dXp = -F(:,1)
        dXc = -F(:,2)
        TF_SS = omp_get_wtime()
        DT_SS = TF_SS - TI_SS

!   UPDATING THE SOLUTION
        dlambda = -dXc(CONTROLED_DOF)/dXp(CONTROLED_DOF)
        dX = dXc + dlambda*dXp
        lambda = lambda + dlambda
        X = X + dX

!   UPDATING THE COHESIVE STIFFNESS
        TI_CS = omp_get_wtime()
        Knn_OLD = Knn
        Ktt_OLD = Ktt
        CALL UPDATE_COHESIVE_STIFFNESS(X)
        IF(MAXVAL(ABS(Knn-Knn_OLD)).GT.0.D0 .OR. MAXVAL(ABS(Ktt-Ktt_OLD)).GT.0.D0)THEN
            MAT_CHANGED = .TRUE.
        ENDIF
        IF(MAXVAL(ABS(Knn-Knn_OLD)).GT.1.0D-3*kn0 .OR. MAXVAL(ABS(Ktt-Ktt_OLD)).GT.1.0D-3*kt0 .OR. MOD(ite,5).EQ.0)THEN
            LU_CHANGED = .TRUE.
        ENDIF
        TF_CS = omp_get_wtime()
        DT_CS = TF_CS - TI_CS

!   COMPUTING THE RESIDUE VECTOR
        TI_PSI = omp_get_wtime()
        IF(MAT_CHANGED)THEN
            CALL MATRIX_M
            MAT_CHANGED = .FALSE.
        ENDIF
        PSI = MATMUL(MATRIX,X) + lambda*(Fu_PRES-Ft_PRES)
        TF_PSI = omp_get_wtime()
        DT_PSI = TF_PSI - TI_PSI

!   COMPUTING THE CONVERGENCE CRITERION
        criterionPSI_num = 0.D0
        criterionPSI_denum = 0.D0
        DO I=1,2*N_NODES
            criterionPSI_num   = criterionPSI_num   + PSI(I)**2.D0
            criterionPSI_denum = criterionPSI_denum + (Fu_PRES(I)-Ft_PRES(I))**2.D0
        ENDDO
        criterionPSI = DSQRT(criterionPSI_num/criterionPSI_denum)

        WRITE(40,400) II,ite,DT_FuFt,DT_F,DT_TM,DT_SS,DT_CS,DT_PSI
        DT_FuFt = 0.D0
    ENDDO

    CALL UPDATE_U_AND_T(X)
    IF(BC(2*NODE-1).EQ.0)THEN
        lambdatx = lambda
    ELSE
        lambdaux = lambda
    ENDIF

    IF(BC(2*NODE).EQ.0)THEN
        lambdaty = lambda
    ELSE
        lambdauy = lambda
    ENDIF

    WRITE(10,100) 1000.D0*lambdaux*U(2*NODE-1),1000.D0*lambdauy*U(2*NODE),lambdatx*T(2*NODE-1),lambdaty*T(2*NODE)

    WRITE(20,200) Knn(COHESIVE_NODE),Knt(COHESIVE_NODE),Ktn(COHESIVE_NODE),Ktt(COHESIVE_NODE),&
        DU_BAR(2*COHESIVE_NODE-1),DU_BAR(2*COHESIVE_NODE)

    DO I=1,N_NODES
        IF(Beta(I).GT.Beta_MAX(I))THEN
            Beta_MAX(I) = Beta(I)
        ENDIF
    ENDDO
    DO I=1,N_NODES
        WRITE(30,300)II,Knn(I),Ktt(I)
        lambdaux = 1.D0
        lambdauy = 1.D0
        IF(COLLOCATION_TYPE(I).EQ."B")THEN
            IF(BC(2*I-1).NE.0) lambdaux = lambda
            IF(BC(2*I).NE.0)   lambdauy = lambda
        ENDIF
        WRITE(50,500)II,I,lambdaux*U(2*I-1),lambdauy*U(2*I)
    ENDDO

ENDDO
CLOSE(10)
CLOSE(20)
CLOSE(30)
CLOSE(40)
CLOSE(50)

DO I=1,N_NODES
    IF(COLLOCATION_TYPE(I).EQ."B")THEN
        lambdaux = 1.D0
        lambdauy = 1.D0
        lambdatx = 1.D0
        lambdaty = 1.D0
        IF(BC(2*I-1).EQ.0)THEN
            lambdatx = lambda
        ELSE
            lambdaux = lambda
        ENDIF

        IF(BC(2*I).EQ.0)THEN
            lambdaty = lambda
        ELSE
            lambdauy = lambda
        ENDIF
        U(2*I-1) = lambdaux*U(2*I-1)
        U(2*I)   = lambdauy*U(2*I)
        T(2*I-1) = lambdatx*T(2*I-1)
        T(2*I)   = lambdaty*T(2*I)
    ENDIF
ENDDO

DEALLOCATE(dPSIdX,MATRIX,Fu_PRES,Ft_PRES)

100 FORMAT(F16.7,12x,F16.7,12x,F16.7,12x,F16.7)
200 FORMAT(E14.7,8x,E14.7,8x,E14.7,8x,E14.7,8x,E14.7,8x,E14.7)
300 FORMAT(I3,12x,F16.7,12x,F16.7)
400 FORMAT(I3,4x,I3,4x,E14.7,4x,E14.7,4x,E14.7,4x,E14.7,4x,E14.7,4x,E14.7)
500 FORMAT(I3,8x,I5,8x,E16.8,8x,E16.8)

END SUBROUTINE SOLVE_NONLINEAR_SYSTEM
!
END MODULE MODULE_SOLVER
