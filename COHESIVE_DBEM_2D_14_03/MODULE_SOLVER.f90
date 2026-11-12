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

! -------------------------------------------------------------------------------------!
!   USED DATA:                                                                         !
! -------------------------------------------------------------------------------------!
! USE BOUNDARY_CONDITIONS_DATA
! USE ANALYSIS_PARAMETERS_DATA
! USE MATERIAL_PROPERTIES_DATA

USE MODULE_DATA
!USE OMP_LIB

! -------------------------------------------------------------------------------------!
!   USED MODULE:                                                                       !
! -------------------------------------------------------------------------------------!
USE MODULE_GEOMETRIC_MAPPINGS

CONTAINS

! ========================================================================================!
![1]                                                                                      !
SUBROUTINE MATRIX_M                                                                       !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!        
!
IMPLICIT NONE 
!
INTEGER::I,I_AUX,K,II,JJ,KK,CONT,NB_GDL_ON,NB_GDL_OFF,IND_I,IND_F
REAL(KIND=wp)::K11,K12,K21,K22,JAC,nv(2)
REAL(KIND=wp),ALLOCATABLE,DIMENSION(:,:)::VALUES

MATRIX=0.D0
CONT=0

NB_GDL_ON = SIZE(VETOR_B_GDL_ON)
NB_GDL_OFF = SIZE(VETOR_B_GDL_OFF)

IND_I = 1
IND_F = NB_GDL_OFF
MATRIX(:,IND_I:IND_F) = H(:,VETOR_B_GDL_OFF)

IND_I = NB_GDL_OFF + 1
IND_F = NB_GDL_ON + NB_GDL_OFF
MATRIX(:,IND_I:IND_F) = -G(:,VETOR_B_GDL_ON)

CONT = CONT + NB_GDL_ON + NB_GDL_OFF

DO I_AUX=1,SIZE(VETOR_NODE_Im)
    I = VETOR_NODE_Im(I_AUX)
    JJ=COHESIVE_INTERFACES_Im_Ip(I,1)
    CONT=CONT+1
    MATRIX(:,CONT)=H(:,2*I-1)+H(:,2*JJ-1)
    CONT=CONT+1
    MATRIX(:,CONT)=H(:,2*I)+H(:,2*JJ)
ENDDO

DO I_AUX=1,SIZE(VETOR_NODE_Ip)
    
    I = VET_Ip_ELEM_ORDER(I_AUX,1)
    JJ = VET_Ip_ELEM_ORDER(I_AUX,2)
    KK = VET_Ip_ELEM_ORDER(I_AUX,3)

    II=COHESIVE_INTERFACES_Ip_Im(I,1)

    ALLOCATE(VALUES(ORDER(JJ)+1,2))

    DO K=1,ORDER(JJ)+1
        VALUES(K,1)=COORDNODES(CONECTIVITY(JJ,K),1) 
        VALUES(K,2)=COORDNODES(CONECTIVITY(JJ,K),2) 
    ENDDO

    CALL JACOBIAN_AND_NORMAL_VECTOR (XSI_COLLOCATIONPOINTS(JJ,KK),JJ,VALUES,JAC,nv)
    
    DEALLOCATE(VALUES)  
    K11=nv(1)**2.D0*Knn(I)-nv(1)*nv(2)*Ktn(I)-nv(1)*nv(2)*Knt(I)+nv(2)**2.D0*Ktt(I)
    K12=nv(1)*nv(2)*Knn(I)-nv(2)**2.D0*Ktn(I)+nv(1)**2.D0*Knt(I)-nv(1)*nv(2)*Ktt(I)
    K21=nv(1)*nv(2)*Knn(I)+nv(1)**2.D0*Ktn(I)-nv(2)**2.D0*Knt(I)-nv(1)*nv(2)*Ktt(I)
    K22=nv(2)**2.D0*Knn(I)+nv(1)*nv(2)*Ktn(I)+nv(1)*nv(2)*Knt(I)+nv(1)**2.D0*Ktt(I)
!            
    CONT=CONT+1
    MATRIX(:,CONT)=H(:,2*I-1)+(G(:,2*I-1)-G(:,2*II-1))*K11+(G(:,2*I)-G(:,2*II))*K21
    CONT=CONT+1
    MATRIX(:,CONT)=H(:,2*I)+(G(:,2*I-1)-G(:,2*II-1))*K12+(G(:,2*I)-G(:,2*II))*K22

ENDDO   

END SUBROUTINE MATRIX_M    

! ========================================================================================!
![2]                                                                                      !
SUBROUTINE VECTOR_Fu_PRES                                                                 !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!        
!
IMPLICIT NONE 
!
INTEGER::IND_I,IND_F,NB_GDL_ON
REAL(KIND=wp),ALLOCATABLE::Hbu(:,:),U_PRES(:)  
!
ALLOCATE(Hbu(2*N_NODES,NDOF_U),U_PRES(NDOF_U))
!
Hbu=0.D0
U_PRES=0.D0

NB_GDL_ON = SIZE(VETOR_B_GDL_ON)

IND_I = 1
IND_F = NB_GDL_ON

U_PRES(IND_I:IND_F) = U(VETOR_B_GDL_ON)
Hbu(:,IND_I:IND_F) = H(:,VETOR_B_GDL_ON)

Fu_PRES=MATMUL(Hbu,U_PRES)
!
DEALLOCATE(Hbu,U_PRES)
!        
END SUBROUTINE VECTOR_Fu_PRES   

! ========================================================================================!
![3]                                                                                      !
SUBROUTINE VECTOR_Ft_PRES                                                                 !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!        
!
IMPLICIT NONE 
!
INTEGER::IND_I,IND_F,NB_GDL_OFF
REAL(KIND=wp),ALLOCATABLE::Gbt(:,:),T_PRES(:)  
!
ALLOCATE(Gbt(2*N_NODES,NDOF_T),T_PRES(NDOF_T))
!
Gbt=0.D0
T_PRES=0.D0

NB_GDL_OFF = SIZE(VETOR_B_GDL_OFF)

IND_I = 1
IND_F = NB_GDL_OFF

T_PRES(IND_I:IND_F) = T(VETOR_B_GDL_OFF)
Gbt(:,IND_I:IND_F) = G(:,VETOR_B_GDL_OFF)

Ft_PRES=MATMUL(Gbt,T_PRES)
!
DEALLOCATE(Gbt,T_PRES)
!                
END SUBROUTINE VECTOR_Ft_PRES 

! ========================================================================================!
![4]                                                                                      !
SUBROUTINE TANGENT_MATRIX_dPSIdX(X)                                                       !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!        
!
IMPLICIT NONE 
!
INTEGER::I,J,K,II,JJ,KK,CONT,N0,I_AUX
REAL(KIND=wp),ALLOCATABLE::VALUES(:,:)
REAL(KIND=wp)::X(2*N_NODES),dK11dDUx(N_NODES),dK11dDUy(N_NODES),&
dK12dDUx(N_NODES),dK12dDUy(N_NODES),dK21dDUx(N_NODES),dK21dDUy(N_NODES),dK22dDUx(N_NODES),&
dK22dDUy(N_NODES),DU(2),nv(2),JAC,DUn,DUt
REAL(KIND=wp)::dK11dDUn,dK11dDUt,dK12dDUn,dK12dDUt,dK21dDUn,dK21dDUt,dK22dDUn,dK22dDUt
REAL(KIND=wp)::dKnndDUn,dKnndDUt,dKntdDUn,dKntdDUt,dKtndDUn,dKtndDUt,dKttdDUn,dKttdDUt
REAL(KIND=wp)::Limit_Tn,Limit_Tt,av,bv
REAL(KIND=wp)::Kappa,dBetadDUn,dKappadDUn,dBetadDUt,dKappadDUt,Sigma_n,Sigma_t,Criterion
!
CALL MATRIX_M
dPSIdX=MATRIX
dK11dDUx=0.D0
dK11dDUy=0.D0
dK12dDUx=0.D0
dK12dDUy=0.D0
dK21dDUx=0.D0
dK21dDUy=0.D0
dK22dDUx=0.D0
dK22dDUy=0.D0
Limit_Tn=(1.D-3)*smax
Limit_Tt=(1.D-3)*tmax
av=1.D0
bv=1.D0
!
N0=NDOF_T+NDOF_U+NDOF_DU
CONT=N0
DO I_AUX=1,SIZE(VETOR_NODE_Ip)
    I = VETOR_NODE_Ip(I_AUX)
    CONT=CONT+1
    DU(1)=X(CONT)
    CONT=CONT+1
    DU(2)=X(CONT)
!
    II=COHESIVE_INTERFACES_Ip_Im(I,1)
    JJ=COHESIVE_INTERFACES_Ip_Im(I,2)
    KK=COHESIVE_INTERFACES_Ip_Im(I,3)
!           
    ALLOCATE(VALUES(ORDER(JJ)+1,2))
    DO K=1,ORDER(JJ)+1
        VALUES(K,1)=COORDNODES(CONECTIVITY(JJ,K),1) 
        VALUES(K,2)=COORDNODES(CONECTIVITY(JJ,K),2) 
    ENDDO
    CALL JACOBIAN_AND_NORMAL_VECTOR (XSI_COLLOCATIONPOINTS(JJ,KK),JJ,VALUES,JAC,nv)
    DUn=nv(1)*DU(1)+nv(2)*DU(2)
    !DUn_int = IDINT(DUn)
    IF(ABS(DUn).LE.1.D-12)THEN
        DUn=1.D-12  
    ENDIF
    DUt=-nv(2)*DU(1)+nv(1)*DU(2)
    !DUt_int = IDINT(DUt)
    IF(ABS(DUt).LE.1.D-12)THEN
        DUt=1.D-12    
    ENDIF
!
    Sigma_n=kn0*DUn
    Sigma_t=kt0*DUt
    Criterion=((Sigma_n/smax)**2.D0+(Sigma_t/tmax)**2.D0)**0.5D0
    Kappa=((DUn/deltanc)**2.D0+(DUt/deltatc)**2.D0)**0.5D0-1.D0
    Beta(I)=Kappa/(Etan*(1+Kappa))
    dKappadDUn=DUn/((deltanc**2.D0)*Kappa)
    dBetadDUn=dKappadDUn/(Etan*(Kappa+1.D0)**2.D0)
    dKappadDUt=DUt/((deltatc**2.D0)*Kappa)
    dBetadDUt=dKappadDUt/(Etan*(Kappa+1.D0)**2.D0)
!
    IF(Criterion.LT.1.D0)THEN
        dKnndDUn=Limit_Tn/(DUn**2.D0)
        dKnndDUt=0.D0
        dKntdDUn=0.D0
        dKntdDUt=0.D0
!                
        dKtndDUn=0.D0
        dKtndDUt=0.D0
        dKttdDUn=0.D0
        dKttdDUt=Limit_Tt/(abs(DUt)**2.D0)
    ELSE IF(Criterion.GE.1.D0.AND.Beta(I).LT.1.D0.AND.Beta(I).GE.Beta_MAX(I))THEN             
        dKnndDUn=-dBetadDUn*(1.D0/deltanc)*smax
        dKnndDUt=-dBetadDUt*(1.D0/deltanc)*smax
        dKntdDUn=0.D0
        dKntdDUt=0.D0
!  
        dKtndDUn=0.D0
        dKtndDUt=0.D0
        dKttdDUn=-dBetadDUn*(1.D0/deltatc)*tmax
        dKttdDUt=-dBetadDUt*(1.D0/deltatc)*tmax  
    ELSE IF(Criterion.GE.1.D0.AND.Beta(I).LT.1.D0.AND.Beta(I).LT.Beta_MAX(I))THEN
        dKnndDUn=Limit_Tn/(DUn**2.D0)
        dKnndDUt=0.D0
        dKntdDUn=0.D0
        dKntdDUt=0.D0
!                
        dKtndDUn=0.D0
        dKtndDUt=0.D0
        dKttdDUn=0.D0
        dKttdDUt=Limit_Tt/(abs(DUt)**2.D0)
    ELSE IF(Criterion.GE.1.D0.AND.Beta(I).GE.1.D0)THEN
        dKnndDUn=-Limit_Tn/(DUn**2.D0)
        dKnndDUt=0.D0
        dKntdDUn=0.D0
        dKntdDUt=0.D0 
!
        dKtndDUn=0.D0
        dKtndDUt=0.D0
        dKttdDUn=0.D0
        dKttdDUt=-Limit_Tt/(DUt**2.D0)
    ENDIF               
!
    dK11dDUn=nv(1)**2.D0*dKnndDUn-nv(1)*nv(2)*dKtndDUn-nv(1)*nv(2)*dKntdDUn+nv(2)**2.D0*dKttdDUn
    dK11dDUt=nv(1)**2.D0*dKnndDUt-nv(1)*nv(2)*dKtndDUt-nv(1)*nv(2)*dKntdDUt+nv(2)**2.D0*dKttdDUt
    dK21dDUn=nv(1)*nv(2)*dKnndDUn+nv(1)**2.D0*dKtndDUn-nv(2)**2.D0*dKntdDUn-nv(1)*nv(2)*dKttdDUn
    dK21dDUt=nv(1)*nv(2)*dKnndDUt+nv(1)**2.D0*dKtndDUt-nv(2)**2.D0*dKntdDUt-nv(1)*nv(2)*dKttdDUt
    dK12dDUn=dK21dDUn
    dK12dDUt=dK21dDUt
    dK22dDUn=nv(2)**2.D0*dKnndDUn+nv(1)*nv(2)*dKtndDUn+nv(1)*nv(2)*dKntdDUn+nv(1)**2.D0*dKttdDUn
    dK22dDUt=nv(2)**2.D0*dKnndDUt+nv(1)*nv(2)*dKtndDUt+nv(1)*nv(2)*dKntdDUt+nv(1)**2.D0*dKttdDUt
!
    dK11dDUx(I)=dK11dDUn*nv(1)-dK11dDUt*nv(2)
    dK11dDUy(I)=dK11dDUn*nv(2)+dK11dDUt*nv(1)
    dK21dDUx(I)=dK21dDUn*nv(1)-dK21dDUt*nv(2)
    dK21dDUy(I)=dK21dDUn*nv(2)+dK21dDUt*nv(1) 
    dK12dDUx(I)=dK12dDUn*nv(1)-dK12dDUt*nv(2)
    dK12dDUy(I)=dK12dDUn*nv(2)+dK12dDUt*nv(1)            
    dK22dDUx(I)=dK22dDUn*nv(1)-dK22dDUt*nv(2)
    dK22dDUy(I)=dK22dDUn*nv(2)+dK22dDUt*nv(1)
!            
    nv=-nv  
!            
    dK11dDUn=nv(1)**2.D0*dKnndDUn-nv(1)*nv(2)*dKtndDUn-nv(1)*nv(2)*dKntdDUn+nv(2)**2.D0*dKttdDUn
    dK11dDUt=nv(1)**2.D0*dKnndDUt-nv(1)*nv(2)*dKtndDUt-nv(1)*nv(2)*dKntdDUt+nv(2)**2.D0*dKttdDUt
    dK21dDUn=nv(1)*nv(2)*dKnndDUn+nv(1)**2.D0*dKtndDUn-nv(2)**2.D0*dKntdDUn-nv(1)*nv(2)*dKttdDUn
    dK21dDUt=nv(1)*nv(2)*dKnndDUt+nv(1)**2.D0*dKtndDUt-nv(2)**2.D0*dKntdDUt-nv(1)*nv(2)*dKttdDUt
    dK12dDUn=dK21dDUn
    dK12dDUt=dK21dDUt
    dK22dDUn=nv(2)**2.D0*dKnndDUn+nv(1)*nv(2)*dKtndDUn+nv(1)*nv(2)*dKntdDUn+nv(1)**2.D0*dKttdDUn
    dK22dDUt=nv(2)**2.D0*dKnndDUt+nv(1)*nv(2)*dKtndDUt+nv(1)*nv(2)*dKntdDUt+nv(1)**2.D0*dKttdDUt
!  
    dK11dDUx(II)=dK11dDUn*nv(1)-dK11dDUt*nv(2)
    dK11dDUy(II)=dK11dDUn*nv(2)+dK11dDUt*nv(1)
    dK21dDUx(II)=dK21dDUn*nv(1)-dK21dDUt*nv(2)
    dK21dDUy(II)=dK21dDUn*nv(2)+dK21dDUt*nv(1) 
    dK12dDUx(II)=dK12dDUn*nv(1)-dK12dDUt*nv(2)
    dK12dDUy(II)=dK12dDUn*nv(2)+dK12dDUt*nv(1)            
    dK22dDUx(II)=dK22dDUn*nv(1)-dK22dDUt*nv(2)
    dK22dDUy(II)=dK22dDUn*nv(2)+dK22dDUt*nv(1)
!            
    DEALLOCATE(VALUES)
ENDDO
!    
CONT=N0
DO I_AUX=1,SIZE(VETOR_NODE_Ip)
    I = VETOR_NODE_Ip(I_AUX)
    II=COHESIVE_INTERFACES_Ip_Im(I,1)
    CONT=CONT+1
    DO J=N0,N0+NDOF_DU
        IF(MOD(J,2).GT.0)THEN
            MATRIX(:,CONT)=MATRIX(:,CONT)+X(J)*((G(:,2*I-1)-G(:,2*II-1))*dK11dDUx(I)+(G(:,2*I)-G(:,2*II))*dK21dDUx(I))
            MATRIX(:,CONT)=MATRIX(:,CONT)+X(J+1)*((G(:,2*I-1)-G(:,2*II-1))*dK12dDUx(I)+(G(:,2*I)-G(:,2*II))*dK22dDUx(I))
        ENDIF
    ENDDO
    CONT=CONT+1
    DO J=N0,N0+NDOF_DU
        IF(MOD(J,2).GT.0)THEN
            MATRIX(:,CONT)=MATRIX(:,CONT)+X(J)*((G(:,2*I-1)-G(:,2*II-1))*dK11dDUy(I)+(G(:,2*I)-G(:,2*II))*dK21dDUy(I))
            MATRIX(:,CONT)=MATRIX(:,CONT)+X(J+1)*((G(:,2*I-1)-G(:,2*II-1))*dK12dDUy(I)+(G(:,2*I)-G(:,2*II))*dK22dDUy(I))
        ENDIF
    ENDDO
ENDDO
!     
END SUBROUTINE TANGENT_MATRIX_dPSIdX

! ========================================================================================!
![5]                                                                                      !
SUBROUTINE UPDATE_COHESIVE_STIFFNESS(X)                                                   !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!        
!
IMPLICIT NONE 
!
INTEGER::I,K,II,JJ,KK,CONT,N0,I_AUX
REAL(KIND=wp),ALLOCATABLE::VALUES(:,:)
REAL(KIND=wp)::X(2*N_NODES),DU(2),nv(2),JAC,DUn,DUt,Limit_Tn,Limit_Tt,av,bv
REAL(KIND=wp)::Kappa,Sigma_n,Sigma_t,Criterion
!
Limit_Tn=(1.D-3)*smax
Limit_Tt=(1.D-3)*tmax
av=1.D0
bv=1.D0
!
N0=NDOF_T+NDOF_U+NDOF_DU
CONT=N0
DO I_AUX=1,SIZE(VETOR_NODE_Ip)
    I = VETOR_NODE_Ip(I_AUX)
    CONT=CONT+1
    DU(1)=X(CONT)
    CONT=CONT+1
    DU(2)=X(CONT)
!
    II=COHESIVE_INTERFACES_Ip_Im(I,1)
    JJ=COHESIVE_INTERFACES_Ip_Im(I,2)
    KK=COHESIVE_INTERFACES_Ip_Im(I,3)
!  
    ALLOCATE(VALUES(ORDER(JJ)+1,2))
    DO K=1,ORDER(JJ)+1
        VALUES(K,1)=COORDNODES(CONECTIVITY(JJ,K),1) 
        VALUES(K,2)=COORDNODES(CONECTIVITY(JJ,K),2) 
    ENDDO
    CALL JACOBIAN_AND_NORMAL_VECTOR(XSI_COLLOCATIONPOINTS(JJ,KK),JJ,VALUES,JAC,nv)
    DUn=nv(1)*DU(1)+nv(2)*DU(2)
    IF(ABS(DUn).LE.1.D-12)THEN
        DUn=1.D-12    
    ENDIF
    DUt=-nv(2)*DU(1)+nv(1)*DU(2)
    IF(ABS(DUt).LE.1.D-12)THEN
        DUt=1.D-12    
    ENDIF
!
    Sigma_n=kn0*DUn
    Sigma_t=kt0*DUt
    Criterion=((Sigma_n/smax)**2.D0+(Sigma_t/tmax)**2.D0)**0.5D0
    Kappa=((DUn/deltanc)**2.D0+(DUt/deltatc)**2.D0)**0.5D0-1.D0
    Beta(I)=Kappa/(Etan*(1.D0+Kappa))

    IF(Criterion.LT.1.D0)THEN
        Knn(I)=Kn0
        Knn(II)=Kn0
        Knt(I)=0.D0
        Knt(II)=0.D0
!                
        Ktn(I)=0.D0
        Ktn(II)=0.D0
        Ktt(I)=Kt0
        Ktt(II)=Kt0 
    ELSE IF(Criterion.GE.1.D0.AND.Beta(I).LT.1.D0.AND.Beta(I).GE.Beta_MAX(I))THEN
        Knn(I)=(1-Beta(I))*(1.D0/deltanc)*smax
        Knn(II)=(1-Beta(I))*(1.D0/deltanc)*smax
        Knt(I)=0.D0
        Knt(II)=0.D0
!            
        Ktn(I)=0.D0
        Ktn(II)=0.D0
        Ktt(I)=(1-Beta(I))*(1.D0/deltatc)*tmax
        Ktt(II)=(1-Beta(I))*(1.D0/deltatc)*tmax
    ELSE IF(Criterion.GE.1.D0.AND.Beta(I).LT.1.D0.AND.Beta(I).LT.Beta_MAX(I))THEN
        Knn(I)=(1-Beta_MAX(I))*(1.D0/deltanc)*smax
        Knn(II)=(1-Beta_MAX(I))*(1.D0/deltanc)*smax              
        Knt(I)=0.D0
        Knt(II)=0.D0     
!           
        Ktn(I)=0.D0
        Ktn(II)=0.D0
        Ktt(I)=(1-Beta_MAX(I))*(1.D0/deltatc)*tmax
        Ktt(II)=(1-Beta_MAX(I))*(1.D0/deltatc)*tmax 
    ELSE IF(Criterion.GE.1.D0.AND.Beta(I).GE.1.D0)THEN
        Knn(I)=Limit_Tn/DUn
        Knn(II)=Limit_Tn/DUn
        Knt(I)=0.D0
        Knt(II)=0.D0 
!                
        Ktn(I)=0.D0
        Ktn(II)=0.D0
        Ktt(I)=Limit_Tt/abs(DUt)
        Ktt(II)=Limit_Tt/abs(DUt) 
    ENDIF       
!
    DU_BAR(2*I-1)=DUn
    DU_BAR(2*I)=DUt
    DU_BAR(2*II-1)=DUn
    DU_BAR(2*II)=DUt
!
    DEALLOCATE(VALUES)
ENDDO
!   
END SUBROUTINE UPDATE_COHESIVE_STIFFNESS

! ========================================================================================!
![6]                                                                                      !
SUBROUTINE UPDATE_U_AND_T(X)                                                   !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!        
!
INTEGER::I,J,K,JJ,KK,CONT,CM,I_AUX
REAL(KIND=wp),ALLOCATABLE::VALUES(:,:)
REAL(KIND=wp)::X(2*N_NODES),JAC,nv(2),K11,K12,K21,K22
!
CONT=0
DO I=1,N_NODES
    IF(COLLOCATION_TYPE(I).EQ."B")THEN
        IF(BC(2*I-1).EQ.0)THEN
            CONT=CONT+1
            U(2*I-1)=X(CONT)
        ENDIF
        IF(BC(2*I).EQ.0)THEN
            CONT=CONT+1
            U(2*I)=X(CONT)
        ENDIF
    ENDIF
ENDDO
!    
DO I=1,N_NODES
    IF(COLLOCATION_TYPE(I).EQ."B")THEN
        IF(BC(2*I-1).EQ.1)THEN
            CONT=CONT+1
            T(2*I-1)=X(CONT)
        ENDIF
        IF(BC(2*I).EQ.1)THEN
            CONT=CONT+1
            T(2*I)=X(CONT)
        ENDIF
    ENDIF
ENDDO
!    
DO I=1,N_NODES
    IF(COLLOCATION_TYPE(I).EQ."Im")THEN
        CONT=CONT+1
        U(2*I-1)=X(CONT)
        CONT=CONT+1
        U(2*I)=X(CONT)
    ENDIF   
ENDDO         
!
DO I_AUX=1,SIZE(VETOR_NODE_Ip)
    I = VETOR_NODE_Ip(I_AUX)
        CM=COHESIVE_INTERFACES_Ip_Im(I,1)
!            
        DO J=1,N_ELEM
            DO K=1,ORDER(J)+1
                IF(I.EQ.CONECTIVITY(J,K))THEN
                    JJ=J
                    KK=K
                ENDIF
            ENDDO
        ENDDO
        ALLOCATE(VALUES(ORDER(JJ)+1,2))
        DO K=1,ORDER(JJ)+1
            VALUES(K,1)=COORDNODES(CONECTIVITY(JJ,K),1) 
            VALUES(K,2)=COORDNODES(CONECTIVITY(JJ,K),2) 
        ENDDO
        CALL JACOBIAN_AND_NORMAL_VECTOR(XSI_COLLOCATIONPOINTS(JJ,KK),JJ,VALUES,JAC,nv)
        DEALLOCATE(VALUES) 
!
        K11=nv(1)**2.D0*Knn(I)-nv(1)*nv(2)*Ktn(I)-nv(1)*nv(2)*Knt(I)+nv(2)**2.D0*Ktt(I)
        K12=nv(1)*nv(2)*Knn(I)-nv(2)**2.D0*Ktn(I)+nv(1)**2.D0*Knt(I)-nv(1)*nv(2)*Ktt(I)
        K21=nv(1)*nv(2)*Knn(I)+nv(1)**2.D0*Ktn(I)-nv(2)**2.D0*Knt(I)-nv(1)*nv(2)*Ktt(I)
        K22=nv(2)**2.D0*Knn(I)+nv(1)*nv(2)*Ktn(I)+nv(1)*nv(2)*Knt(I)+nv(1)**2.D0*Ktt(I)            
!            
        CONT=CONT+1
        U(2*I-1)=U(2*CM-1)+X(CONT)
        T(2*CM-1)=K11*X(CONT)+K12*X(CONT+1)   
        T(2*I-1)=-T(2*CM-1)
        CONT=CONT+1
        U(2*I)=U(2*CM)+X(CONT)
        T(2*CM)=K21*X(CONT-1)+K22*X(CONT)
        T(2*I)=-T(2*CM)                
ENDDO  
!      
END SUBROUTINE UPDATE_U_AND_T

! ========================================================================================!
![7]                                                                                      !
SUBROUTINE SOLVE_LINEAR_SYSTEM                                                            !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!        
!
IMPLICIT NONE 
!
INTEGER::IPIV(2*N_NODES),INFO,I,J,K,II,JJ,KK,CONT
REAL(KIND=wp),ALLOCATABLE::VALUES(:,:),F(:)
REAL(KIND=wp)::JAC,nv(2),K11,K12,K21,K22
!
ALLOCATE(MATRIX(2*N_NODES,2*N_NODES),F(2*N_NODES),Fu_PRES(2*N_NODES),Ft_PRES(2*N_NODES))
!
CALL MATRIX_M
CALL VECTOR_Fu_PRES
CALL VECTOR_Ft_PRES
F=Ft_PRES-Fu_PRES                     
!
CALL DGESV(2*N_NODES,1,MATRIX,2*N_NODES,IPIV,F,2*N_NODES,INFO)
!
CONT=0
DO I=1,N_NODES
    IF(COLLOCATION_TYPE(I).EQ."B")THEN
        IF(BC(2*I-1).EQ.0)THEN
            CONT=CONT+1
            U(2*I-1)=F(CONT)
        ENDIF
        IF(BC(2*I).EQ.0)THEN
            CONT=CONT+1
            U(2*I)=F(CONT)
        ENDIF          
    ENDIF
ENDDO
DO I=1,N_NODES
    IF(COLLOCATION_TYPE(I).EQ."B")THEN
        IF(BC(2*I-1).EQ.1)THEN
            CONT=CONT+1
            T(2*I-1)=F(CONT)
        ENDIF
        IF(BC(2*I).EQ.1)THEN
            CONT=CONT+1
            T(2*I)=F(CONT)
        ENDIF          
    ENDIF
ENDDO
DO I=1,N_NODES
    IF(COLLOCATION_TYPE(I).EQ."Im")THEN
        CONT=CONT+1
        U(2*I-1)=F(CONT)
        CONT=CONT+1
        U(2*I)=F(CONT)          
    ENDIF
ENDDO
DO I=1,N_NODES
    IF(COLLOCATION_TYPE(I).EQ."Ip")THEN 
        II=COHESIVE_INTERFACES_Ip_Im(I,1)            
        DO J=1,N_ELEM
            DO K=1,ORDER(J)+1
                IF(I.EQ.CONECTIVITY(J,K))THEN
                    JJ=J
                    KK=K
                ENDIF
            ENDDO
        ENDDO
        ALLOCATE(VALUES(ORDER(JJ)+1,2))
        DO K=1,ORDER(JJ)+1
            VALUES(K,1)=COORDNODES(CONECTIVITY(JJ,K),1) 
            VALUES(K,2)=COORDNODES(CONECTIVITY(JJ,K),2) 
        ENDDO
        CALL JACOBIAN_AND_NORMAL_VECTOR(XSI_COLLOCATIONPOINTS(JJ,KK),JJ,VALUES,JAC,nv)
        DEALLOCATE(VALUES) 
!
        K11=nv(1)**2.D0*Knn(I)-nv(1)*nv(2)*Ktn(I)-nv(1)*nv(2)*Knt(I)+nv(2)**2.D0*Ktt(I)
        K12=nv(1)*nv(2)*Knn(I)-nv(2)**2.D0*Ktn(I)+nv(1)**2.D0*Knt(I)-nv(1)*nv(2)*Ktt(I)
        K21=nv(1)*nv(2)*Knn(I)+nv(1)**2.D0*Ktn(I)-nv(2)**2.D0*Knt(I)-nv(1)*nv(2)*Ktt(I)
        K22=nv(2)**2.D0*Knn(I)+nv(1)*nv(2)*Ktn(I)+nv(1)*nv(2)*Knt(I)+nv(1)**2.D0*Ktt(I)
!            
        CONT=CONT+1
        U(2*I-1)=U(2*II-1)+F(CONT)
        T(2*II-1)=K11*F(CONT)+K12*F(CONT+1)   
        T(2*I-1)=-T(2*II-1)
        CONT=CONT+1
        U(2*I)=U(2*II)+F(CONT)
        T(2*II)=K21*F(CONT-1)+K22*F(CONT)
        T(2*I)=-T(2*II)          
    ENDIF
ENDDO
!      
END SUBROUTINE SOLVE_LINEAR_SYSTEM

! ========================================================================================!
![8]                                                                                      !
SUBROUTINE SOLVE_NONLINEAR_SYSTEM                                                         !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!        
!
IMPLICIT NONE 
!
INTEGER::II,I,IPIV(2*N_NODES),INFO,ite,itecumulated!,thread_id
REAL(KIND=wp)::X(2*N_NODES),lambda,PSI(2*N_NODES),&
F(2*N_NODES,2),dX(2*N_NODES),dXpres(2*N_NODES),dXc(2*N_NODES),dXp(2*N_NODES),criterionPSI_num,&
criterionPSI_denum,criterionPSI,dLambda,lambdaux,lambdauy,lambdatx,lambdaty
REAL(KIND=wp)::TI_FuFt,TF_FuFt,DT_FuFt
REAL(KIND=wp)::TI_F,TF_F,DT_F
REAL(KIND=wp)::TI_TM,TF_TM,DT_TM
REAL(KIND=wp)::TI_SS,TF_SS,DT_SS
REAL(KIND=wp)::TI_CS,TF_CS,DT_CS
REAL(KIND=wp)::TI_PSI,TF_PSI,DT_PSI
!
ALLOCATE(dPSIdX(2*N_NODES,2*N_NODES),MATRIX(2*N_NODES,2*N_NODES),Fu_PRES(2*N_NODES),Ft_PRES(2*N_NODES))
!    
X=0.D0
X(CONTROLED_DOF)=CONTROLED_VALUE
dXpres=X/Ns
X=0.D0
lambda=0.D0
lambdaux=1.D0
lambdauy=1.D0
lambdatx=1.D0
lambdaty=1.D0
itecumulated=0
!
OPEN(10,FILE='8_OUTPUT_NODAL_STRUCTURAL_RESPONSE.txt',STATUS='UNKNOWN')
WRITE(10,*) "       Ux              Uy              Tx              Ty"
OPEN(20,FILE='9_OUTPUT_COHESIVE_RESPONSE.txt',STATUS='UNKNOWN')
WRITE(20,*) "       Knn          Knt          Ktn          Ktt          Dun          Dut"
OPEN(30,FILE='10_OUTPUT_COHESIVE_STIFFNESS.txt',STATUS='UNKNOWN')
WRITE(30,*) "       STEP          Knn          Ktt"
WRITE(40,*) "STEP   ITER     Fu-Ft_V             F_V                T_M               S_S                UP_C_S     PSI_V"
DO II=1,Ns      
    criterionPSI=1.D0        
    ite = 0
    X(CONTROLED_DOF)=X(CONTROLED_DOF)+dXpres(CONTROLED_DOF)
!        
    DO WHILE (criterionPSI.GT.tol_Newton)
        ite=ite+1
        itecumulated=itecumulated+1
        IF(ite.GT.itemax)THEN
            tol_Newton=tol_Newton*10.D0
            IF(tol_Newton.gt.1.D-1)THEN
                STOP    
            ENDIF
        ENDIF
! 
!   COMPUTING THE dRdlambda VECTOR
! 
        !thread_id = OMP_GET_THREAD_NUM()
        CALL CPU_TIME(TI_FuFt)
        !WRITE(*,*)"Hello from process: ", thread_id
        CALL VECTOR_Fu_PRES
        !WRITE(*,*)"Hello from process: ", thread_id
        CALL VECTOR_Ft_PRES
        !WRITE(*,*)"Hello from process: ", thread_id
        F(:,1)=Fu_PRES-Ft_PRES
        CALL CPU_TIME(TF_FuFt)
        DT_FuFt=TF_FuFt-TI_FuFt
        !dPSIdlambda=Fu_PRES-Ft_PRES
   
!   COMPUTING THE RESIDUE VECTOR
!
        CALL CPU_TIME(TI_F)
        CALL MATRIX_M
        F(:,2)=MATMUL(MATRIX,X)+lambda*(Fu_PRES-Ft_PRES)
        CALL CPU_TIME(TF_F)
        DT_F=TF_F-TI_F
! 
!   COMPUTING THE TANGENT OPERATOR MATRIX
! 
        CALL CPU_TIME(TI_TM)
        CALL TANGENT_MATRIX_dPSIdX(X) 
        CALL CPU_TIME(TF_TM)
        DT_TM=TF_TM-TI_TM
! 
!   SOLVING THE PREDICTION AND CORRECTOR SYSTEM OF EQUATIONS
!
        CALL CPU_TIME(TI_SS)
        CALL DGESV(2*N_NODES,2,dPSIdX,2*N_NODES,IPIV,F,2*N_NODES,INFO)  !SOLVING THE SYSTEM OF EQUATIONS FOR PREVISIONAL DISPLACEMENTS
        dXp=-F(:,1)
!  
        dXc=-F(:,2)
        CALL CPU_TIME(TF_SS)
        DT_SS=TF_SS-TI_SS
! 
!   UPDATING THE SOLUTION
!
        dlambda=-dXc(CONTROLED_DOF)/dXp(CONTROLED_DOF)
        dX=dXc+dlambda*dXp
        lambda=lambda+dlambda
        X=X+dX
! 
!   UPDATING THE COHESIVE STIFFNESS
!
        CALL CPU_TIME(TI_CS)
        CALL UPDATE_COHESIVE_STIFFNESS(X)
        CALL CPU_TIME(TF_CS)
        DT_CS=TF_CS-TI_CS
! 
!   COMPUTING THE RESIDUE VECTOR
!
        CALL CPU_TIME(TI_PSI)
        CALL MATRIX_M
        PSI=MATMUL(MATRIX,X)+lambda*(Fu_PRES-Ft_PRES)  
        CALL CPU_TIME(TF_PSI)
        DT_PSI=TF_PSI-TI_PSI
! 
!   COMPUTING THE CONVERGENCE CRITERION
!                 
        criterionPSI_num=0.D0
        criterionPSI_denum=0.D0
        DO I=1,2*N_NODES
            criterionPSI_num=criterionPSI_num+PSI(I)**2.D0       
            criterionPSI_denum=criterionPSI_denum+(Fu_PRES(I)-Ft_PRES(I))**2.D0 
        ENDDO
        criterionPSI=criterionPSI_num/criterionPSI_denum
        criterionPSI=DSQRT(criterionPSI) 
!
        WRITE(40,400) II,ite,DT_FuFt,DT_F,DT_TM,DT_SS,DT_CS,DT_PSI
    ENDDO
!        
    CALL UPDATE_U_AND_T(X)
    IF(BC(2*NODE-1).EQ.0)THEN
        lambdatx=lambda
    ELSE
        lambdaux=lambda
    ENDIF
!        
    IF(BC(2*NODE).EQ.0)THEN
        lambdaty=lambda
    ELSE
        lambdauy=lambda
    ENDIF
!
         
    WRITE(10,100) 1000*lambdaux*U(2*NODE-1),1000*lambdauy*U(2*NODE),lambdatx*T(2*NODE-1),lambdaty*T(2*NODE)
    !
    WRITE(20,200) Knn(COHESIVE_NODE),Knt(COHESIVE_NODE),Ktn(COHESIVE_NODE),Ktt(COHESIVE_NODE),&
    DU_BAR(2*COHESIVE_NODE-1),DU_BAR(2*COHESIVE_NODE)
    !        
    DO I=1,N_NODES
        IF(Beta(I).GT.Beta_MAX(I))THEN
            Beta_MAX(I)=Beta(I)
        ENDIF    
    ENDDO
    DO I=1,N_NODES
        WRITE(30,300)II,Knn(I),Ktt(I)
    ENDDO
!        
ENDDO
CLOSE(10)
CLOSE(20)
CLOSE(30)
CLOSE(40)
!   
DO I=1,N_NODES
    IF(COLLOCATION_TYPE(I).EQ."B")THEN
        lambdaux=1.D0
        lambdauy=1.D0
        lambdatx=1.D0
        lambdaty=1.D0
        IF(BC(2*I-1).EQ.0)THEN
            lambdatx=lambda
        ELSE
            lambdaux=lambda
        ENDIF
    !        
        IF(BC(2*I).EQ.0)THEN
            lambdaty=lambda
        ELSE
            lambdauy=lambda
        ENDIF
        U(2*I-1)=lambdaux*U(2*I-1)
        U(2*I)=lambdauy*U(2*I)
        T(2*I-1)=lambdatx*T(2*I-1)
        T(2*I)=lambdaty*T(2*I)
    ENDIF
ENDDO
!
100	FORMAT(F16.7,12x,F16.7,12x,F16.7,12x,F16.7)
200 FORMAT(E14.7,8x,E14.7,8x,E14.7,8x,E14.7,8x,E14.7,8x,E14.7)
300 FORMAT(I3,12x,F16.7,12x,F16.7)
400	FORMAT(I3,4x,I3,4x,E14.7,4x,E14.7,4x,E14.7,4x,E14.7,4x,E14.7,4x,E14.7)
!   
END SUBROUTINE SOLVE_NONLINEAR_SYSTEM
!
END MODULE MODULE_SOLVER