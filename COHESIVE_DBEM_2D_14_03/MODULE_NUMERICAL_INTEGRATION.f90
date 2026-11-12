MODULE MODULE_NUMERICAL_INTEGRATION
! -------------------------------------------------------------------------------------!
!   Description. Module that contains subroutines...                                   !
!                                                                                      !
! -------------------------------------------------------------------------------------!
!   This module contains the following subroutines:                                    !
!                                                                                      !
!   [1]. GAUSS_LEGENDRE(GaussPos,GaussWgt,xin,xfin,numGauss)                           !
!   [2]. NONSINGULARINTEGRATION_DBIE(NODE,ELEM,ge,he)                                  !
!   [3]. NONSINGULARINTEGRATION_TBIE(NODE,ELEM,ge,he)                                  !
!   [4]. SINGULARINTEGRATION_DBIE(NODE,ELEM,ge,he)                                     !
!   [5]. SINGULARINTEGRATION_TBIE(NODE,ELEM,ge,he)                                     !
! -------------------------------------------------------------------------------------!

! -------------------------------------------------------------------------------------!
!   USED DATA:                                                                         !
! -------------------------------------------------------------------------------------!
USE MODULE_DATA
! USE BOUNDARY_MESH_DATA
! USE MATERIAL_PROPERTIES_DATA
! USE NUMERICAL_INTEGRATION_DATA

! -------------------------------------------------------------------------------------!
!   USED MODULE:                                                                       !
! -------------------------------------------------------------------------------------!
USE MODULE_GEOMETRIC_MAPPINGS
USE MODULE_APPROXIMATION_FUNCTIONS

CONTAINS

! =====================================================================================!
![1]                                                                                   !
SUBROUTINE GAUSS_LEGENDRE(GaussPos,GaussWgt,xin,xfin,numGauss)                         !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!
	!IMPLICIT NONE 
!
    INTEGER,INTENT(IN)::numGauss
    REAL(KIND=8),INTENT(IN)::xin,xfin
    REAL(KIND=8)::pi
    REAL(KIND=8),ALLOCATABLE,INTENT(OUT)::GaussPos(:),GaussWgt(:)
!
    INTEGER::i,j,m,cont=0
    REAL(KIND=8)::eps,xm,xl,p1,p2,p3,pp,z,z1
!
    pi=DACOS(-1.D0)
!
    ALLOCATE(GaussPos(numGauss),GaussWgt(numGauss))
!
    eps = 3e-14
    m = NINT((numGauss + 1)/2.0)
    xm = 0.5*(xfin + xin)
    xl = 0.5*(xfin - xin)
!
    DO i = 1,m
        z = COS(pi*(i - 0.25)/(numGauss + 0.5))
        cont = 0
        DO WHILE (1 == 1)
            cont = cont + 1
            p1 = 1.0
            p2 = 0.0
            DO j = 1,numGauss
                p3 = p2
                p2 = p1
                p1 = ( (2.0*j - 1.0)*z*p2 - (j - 1.0)*p3 )/j
            END DO
            pp = numGauss*(z*p1 - p2)/(z*z - 1.0)
            z1 = z
            z = z1 - p1/pp
!
            IF (ABS(z - z1)<eps) THEN
                EXIT
            END IF
!
            IF (cont==20) EXIT
        END DO
        GaussPos(i) = xm - xl*z
        GaussPos(numGauss + 1 - i) = xm + xl*z
        GaussWgt(i) = 2*xl/((1 - z*z)*pp*pp)
        GaussWgt(numGauss + 1 - i) = GaussWgt(i)
    END DO
!
END SUBROUTINE GAUSS_LEGENDRE

! =====================================================================================!
![2]                                                                                   !
SUBROUTINE NONSINGULARINTEGRATIONS_DBIE(NODE,ELEM,ge,he)                                !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    
    IMPLICIT NONE
    INTEGER::II,I,J,NODE,ELEM
    REAL(KIND=8),ALLOCATABLE,INTENT(INOUT)::ge(:,:),he(:,:)
    REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:):: VALUES
    REAL(KIND=8),ALLOCATABLE,DIMENSION(:)::PHI
    REAL(KIND=8):: X,Y,JAC,nv(2),pi,r(2),drdx(2),drdn,llrll,Ustar(2,2),Tstar(2,2),&
    IDENTITY(2,2)
    REAL(KIND=8),ALLOCATABLE::MAT_PHI(:,:)
    REAL(KIND=8)::SHEAR_MODULUS,LONGITUDINAL_MODULUS,Nu
!   
    DEALLOCATE(ge,he)
    ALLOCATE(ge(2,2*(ORDER(ELEM)+1)),he(2,2*(ORDER(ELEM)+1)))
    ALLOCATE(MAT_PHI(2,2*(ORDER(ELEM)+1)))

    IDENTITY=0.D0
    IDENTITY(1,1)=1.D0
    IDENTITY(2,2)=1.D0
    pi=DACOS(-1.D0)
    ALLOCATE(VALUES((ORDER(ELEM)+1),2),PHI(ORDER(ELEM)+1))    
    DO II=1,ORDER(ELEM)+1
        VALUES(II,1)=COORDNODES(CONECTIVITY(ELEM,II),1)  !THE MATRIX VALUES RECIVE THE X OF VALUE OF THE "II" LOCAL NODE OF THE ELEMENT "ELEM"
        VALUES(II,2)=COORDNODES(CONECTIVITY(ELEM,II),2)  !THE MATRIX VALUES RECIVE THE Y OF VALUE OF THE "II" LOCAL NODE OF THE ELEMENT "ELEM"
    ENDDO 
!
    ge=0.D0
    he=0.D0
    DO II=1,N_GAUSS
        CALL COORDINATES (COORD_GAUSS(II),ELEM,VALUES,X,Y) !CALL COORDINATES AND INFORM XSI OF A NODE OF THE ELEMENT, THE ELEMENT, THE VALUE OF CSI AND IT RETURN THE X AND Y COORDNATES
        CALL JACOBIAN_AND_NORMAL_VECTOR (COORD_GAUSS(II),ELEM,VALUES,JAC,nv) !CALL JACOBIAN AND INFORM XSI OF A NODE OF THE ELEMENT, THE ELEMENT, THE VALUES MATRIX AND IT RETURN THE JACOBIAN AND A NORMAL VECTOR N
        CALL BASIS_FUNCTIONS (COORD_GAUSS(II),ELEM,PHI)
        MAT_PHI=0.D0
        DO I=1,ORDER(ELEM)+1
            MAT_PHI(1,2*I-1)=PHI(I)
            MAT_PHI(2,2*I)=PHI(I)
        ENDDO
        r(1)=(X-COORDCOLOCATIONPOINTS(NODE,1))
        r(2)=(Y-COORDCOLOCATIONPOINTS(NODE,2))
        llrll=DSQRT(r(1)**2.D0+r(2)**2.D0)
        drdx=r/llrll
        drdn=(drdx(1)*nv(1)+drdx(2)*nv(2))
        LONGITUDINAL_MODULUS=Young_Modulus(SUBREGION(ELEM))
        Nu=Poisson(SUBREGION(ELEM))
        SHEAR_MODULUS=LONGITUDINAL_MODULUS/(2*(1+Nu))
        DO I=1,2
            DO J=1,2
                Ustar(I,J)=1.D0/(8.D0*pi*SHEAR_MODULUS*(1-Nu))*((4*Nu-3.D0)*DLOG(llrll)*IDENTITY(I,J)+drdx(I)*drdx(J))
                Tstar(I,J)=-1.D0/(4.D0*pi*(1.D0-Nu)*llrll)*(drdn*((1.D0-2.D0*Nu)*IDENTITY(I,J)+2.D0*drdx(I)*drdx(J))+&
                    (1.D0-2.D0*Nu)*(drdx(J)*nv(I)-drdx(I)*nv(J)))
            ENDDO
        ENDDO
        ge=ge+matmul(Ustar,MAT_PHI)*JAC*WEIGHT_GAUSS(II)
        he=he+matmul(Tstar,MAT_PHI)*JAC*WEIGHT_GAUSS(II)
    ENDDO
    DEALLOCATE(VALUES,PHI)
!    
END SUBROUTINE NONSINGULARINTEGRATIONS_DBIE

! =====================================================================================!
![3]                                                                                   !
SUBROUTINE NONSINGULARINTEGRATIONS_TBIE(NODE,ELEM,ELEM0,ge,he)                                !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!    
    IMPLICIT NONE
    INTEGER::II,I,J,K,NODE,ELEM,ELEM0,LOCAL_NODE0
    REAL(KIND=8),INTENT(INOUT)::ge(2,2*(ORDER(ELEM)+1)),he(2,2*(ORDER(ELEM)+1))
    REAL(KIND=8),ALLOCATABLE,DIMENSION(:)::PHI
    REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)::VALUES
    REAL(KIND=8)::X,Y,JAC,JAC0,nv(2),n0(2),pi,r(2),drdx(2),drdn,llrll,Dstar(2,2,2),&
    Sstar(2,2,2),D(4,2),S(4,2),IDENTITY(2,2),MAT_PHI(2,2*(ORDER(ELEM)+1))
    REAL(KIND=8)::ge_aux(4,2*(ORDER(ELEM)+1)),he_aux(4,2*(ORDER(ELEM)+1))
    REAL(KIND=8)::SHEAR_MODULUS,LONGITUDINAL_MODULUS,Nu
! 
    IDENTITY=0.D0
    IDENTITY(1,1)=1.D0
    IDENTITY(2,2)=1.D0
    pi=DACOS(-1.D0)
!    
    DO K=1,ORDER(ELEM0)+1
        IF(CONECTIVITY(ELEM0,K).EQ.NODE)THEN
            LOCAL_NODE0=K
        ENDIF
    ENDDO
!    
    ALLOCATE(VALUES((ORDER(ELEM)+1),2),PHI(ORDER(ELEM)+1))    
    DO II=1,ORDER(ELEM)+1
        VALUES(II,1)=COORDNODES(CONECTIVITY(ELEM,II),1)  !THE MATRIX VALUES RECIVE THE X OF VALUE OF THE "II" LOCAL NODE OF THE ELEMENT "ELEM"
        VALUES(II,2)=COORDNODES(CONECTIVITY(ELEM,II),2)  !THE MATRIX VALUES RECIVE THE Y OF VALUE OF THE "II" LOCAL NODE OF THE ELEMENT "ELEM"
    ENDDO 
!
    ge_aux=0.D0
    he_aux=0.D0 
    DO II=1,N_GAUSS
        CALL COORDINATES (COORD_GAUSS(II),ELEM,VALUES,X,Y) !CALL COORDINATES AND INFORM XSI OF A NODE OF THE ELEMENT, THE ELEMENT, THE VALUE OF CSI AND IT RETURN THE X AND Y COORDNATES
        CALL JACOBIAN_AND_NORMAL_VECTOR (COORD_GAUSS(II),ELEM,VALUES,JAC,nv) !CALL JACOBIAN AND INFORM XSI OF A NODE OF THE ELEMENT, THE ELEMENT, THE VALUES MATRIX AND IT RETURN THE JACOBIAN AND A NORMAL VECTOR N
        CALL BASIS_FUNCTIONS (COORD_GAUSS(II),ELEM,PHI)
        MAT_PHI=0.D0
        DO I=1,ORDER(ELEM)+1
            MAT_PHI(1,2*I-1)=PHI(I)
            MAT_PHI(2,2*I)=PHI(I)
        ENDDO
        r(1)=(X-COORDCOLOCATIONPOINTS(NODE,1))
        r(2)=(Y-COORDCOLOCATIONPOINTS(NODE,2))
        llrll=DSQRT(r(1)**2.D0+r(2)**2.D0)
        drdx=r/llrll
        drdn=(drdx(1)*nv(1)+drdx(2)*nv(2))
        LONGITUDINAL_MODULUS=Young_Modulus(SUBREGION(ELEM))
        Nu=Poisson(SUBREGION(ELEM))
        SHEAR_MODULUS=LONGITUDINAL_MODULUS/(2*(1+Nu))
        DO I=1,2
            DO J=1,2
                DO K=1,2
                    Dstar(I,J,K)=1.D0/(4.D0*pi*(1-Nu)*llrll)*((1.D0-2*Nu)*(IDENTITY(I,J)*drdx(K)+IDENTITY(I,K)*drdx(J)-&
                        IDENTITY(J,K)*drdx(I))+2*drdx(I)*drdx(J)*drdx(K))                                                        
                    Sstar(I,J,K)=2.D0*drdn*((1.D0-2.D0*Nu)*IDENTITY(I,J)*drdx(K)+Nu*(IDENTITY(I,K)*drdx(J)+&
                        IDENTITY(J,K)*drdx(I))-4.D0*drdx(I)*drdx(J)*drdx(K))
                    Sstar(I,J,K)=Sstar(I,J,K)+2.D0*Nu*(nv(I)*drdx(J)*drdx(K)+nv(J)*drdx(I)*drdx(K))-&
                    (1.D0-4.D0*Nu)*nv(K)*IDENTITY(I,J)+(1.D0-2.D0*Nu)*(nv(I)*IDENTITY(J,K)+nv(J)*IDENTITY(I,K)+&
                        2.D0*nv(K)*drdx(I)*drdx(J))  
                    Sstar(I,J,K)=2.D0*SHEAR_MODULUS/(4.D0*pi*(1.D0-Nu)*llrll**2.D0)*Sstar(I,J,K)  
                ENDDO
            ENDDO
        ENDDO
        D(1,1)=Dstar(1,1,1)
        D(1,2)=Dstar(2,1,1)
        D(2,1)=Dstar(1,2,1)
        D(2,2)=Dstar(2,2,1)
        D(3,1)=Dstar(1,1,2)
        D(3,2)=Dstar(2,1,2)
        D(4,1)=Dstar(1,2,2)
        D(4,2)=Dstar(2,2,2)
!  
        S(1,1)=Sstar(1,1,1)
        S(1,2)=Sstar(2,1,1)
        S(2,1)=Sstar(1,2,1)
        S(2,2)=Sstar(2,2,1)
        S(3,1)=Sstar(1,1,2)
        S(3,2)=Sstar(2,1,2)
        S(4,1)=Sstar(1,2,2)
        S(4,2)=Sstar(2,2,2)
!                
        ge_aux=ge_aux+matmul(D,MAT_PHI)*JAC*WEIGHT_GAUSS(II)
        he_aux=he_aux+matmul(S,MAT_PHI)*JAC*WEIGHT_GAUSS(II)        
    ENDDO
!    
    DEALLOCATE(VALUES,PHI)    
    ALLOCATE(VALUES((ORDER(ELEM0)+1),2),PHI(ORDER(ELEM0)+1))    
    DO II=1,ORDER(ELEM0)+1
        VALUES(II,1)=COORDNODES(CONECTIVITY(ELEM0,II),1)  !THE MATRIX VALUES RECIVE THE X OF VALUE OF THE "II" LOCAL NODE OF THE ELEMENT "ELEM0"
        VALUES(II,2)=COORDNODES(CONECTIVITY(ELEM0,II),2)  !THE MATRIX VALUES RECIVE THE Y OF VALUE OF THE "II" LOCAL NODE OF THE ELEMENT "ELEM0"
    ENDDO 
    CALL JACOBIAN_AND_NORMAL_VECTOR (XSI_COLLOCATIONPOINTS(ELEM0,LOCAL_NODE0),ELEM0,VALUES,JAC0,n0) !
!    
    DO I=1,2*(ORDER(ELEM)+1)
        DO J=1,2
            ge(J,I)=ge_aux(2*J-1,I)*n0(1)+ge_aux(2*J,I)*n0(2)
            he(J,I)=he_aux(2*J-1,I)*n0(1)+he_aux(2*J,I)*n0(2)
        ENDDO
    ENDDO
    DEALLOCATE(VALUES,PHI)
!      
END SUBROUTINE NONSINGULARINTEGRATIONS_TBIE

! =====================================================================================!
![4]                                                                                   !
SUBROUTINE SINGULARINTEGRATIONS_DBIE(NODE,ELEM,ge,he)                                !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!    
    IMPLICIT NONE
    INTEGER::II,I,J,K,NODE,ELEM,LOCAL_NODE
    REAL(KIND=8),INTENT(INOUT)::ge(2,2*(ORDER(ELEM)+1)),he(2,2*(ORDER(ELEM)+1))
    REAL(KIND=8),ALLOCATABLE,DIMENSION(:)::PHI,PHI0,dPHI0
    REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)::VALUES
    REAL(KIND=8)::X,Y,JAC,JAC0,&
    nv(2),n0(2),pi,r(2),drdx(2),drdxlx0(2),drdn,drdn0,llrll,llrllaux,Ustar(2,2),Tstar(2,2),Ustar_aux(2,2),&
    Tstar_aux(2,2),Tstar_bar(2,2),IDENTITY(2,2),MAT_PHI(2,2*(ORDER(ELEM)+1)),MAT_PHI0(2,2*(ORDER(ELEM)+1)),&
    C1,CPV_lnJ0XSI,CPV_1lXSI,XSI0,X0,Y0,DIST
    REAL(KIND=8)::ge_weaksingular(2,2*(ORDER(ELEM)+1)),he_strongsingular(2,2*(ORDER(ELEM)+1)),TOL
    REAL(KIND=8)::SHEAR_MODULUS,LONGITUDINAL_MODULUS,Nu
!   
    IDENTITY=0.D0
    IDENTITY(1,1)=1.D0
    IDENTITY(2,2)=1.D0
    pi=DACOS(-1.D0)
    TOL=1.D-10
! 
    DIST=0.D0
    DO K=1,ORDER(ELEM)+1
        DIST=DSQRT((COORDCOLOCATIONPOINTS(CONECTIVITY(ELEM,K),1)-COORDCOLOCATIONPOINTS(NODE,1))**2.D0+&
            (COORDCOLOCATIONPOINTS(CONECTIVITY(ELEM,K),2)-COORDCOLOCATIONPOINTS(NODE,2))**2.D0)
        IF(DIST.LT.TOL)THEN
            LOCAL_NODE=K
        ENDIF
    ENDDO
!    
    ALLOCATE(VALUES((ORDER(ELEM)+1),2),PHI(ORDER(ELEM)+1),PHI0(ORDER(ELEM)+1),dPHI0(ORDER(ELEM)+1))    
    DO II=1,ORDER(ELEM)+1
        VALUES(II,1)=COORDNODES(CONECTIVITY(ELEM,II),1)  !THE MATRIX VALUES RECIVE THE X OF VALUE OF THE "II" LOCAL NODE OF THE ELEMENT "ELEM"
        VALUES(II,2)=COORDNODES(CONECTIVITY(ELEM,II),2)  !THE MATRIX VALUES RECIVE THE Y OF VALUE OF THE "II" LOCAL NODE OF THE ELEMENT "ELEM"
    ENDDO 
!
    X0=COORDCOLOCATIONPOINTS(NODE,1)
    Y0=COORDCOLOCATIONPOINTS(NODE,2)
    XSI0=XSI_COLLOCATIONPOINTS(ELEM,LOCAL_NODE)
    CALL JACOBIAN_AND_NORMAL_VECTOR (XSI0,ELEM,VALUES,JAC0,n0)
    CALL BASIS_FUNCTIONS (XSI0,ELEM,PHI0)
    MAT_PHI0=0.D0
    DO II=1,ORDER(ELEM)+1
        MAT_PHI0(1,2*II-1)=PHI0(II)
        MAT_PHI0(2,2*II)=PHI0(II)
    ENDDO  
!    
    CALL SHAPE_FUNCTIONS_DERIVATIVES(XSI0,ELEM,dPHI0)
    drdxlx0=0.D0
    DO II=1,ORDER(ELEM)+1
        drdxlx0(1)=drdxlx0(1)+dPHI0(II)*VALUES(II,1)
        drdxlx0(2)=drdxlx0(2)+dPHI0(II)*VALUES(II,2)
    ENDDO
    drdxlx0=drdxlx0/JAC0
!
! NUMERICAL INTEGRATION OF THE REGULAR PARTS
!        
    ge=0.D0
    he=0.D0
    DO II=1,N_GAUSS
        CALL COORDINATES (COORD_GAUSS(II),ELEM,VALUES,X,Y) !CALL COORDINATES AND INFORM XSI OF A NODE OF THE ELEMENT, THE ELEMENT, THE VALUES MATRIX AND IT RETURN THE X AND Y COORDNATES
        CALL JACOBIAN_AND_NORMAL_VECTOR (COORD_GAUSS(II),ELEM,VALUES,JAC,nv) !CALL JACOBIAN AND INFORM XSI OF A NODE OF THE ELEMENT, THE ELEMENT, THE VALUE OF CSI AND IT RETURN THE JACOBIAN AND A NORMAL VECTOR N
        CALL BASIS_FUNCTIONS (COORD_GAUSS(II),ELEM,PHI)
        MAT_PHI=0.D0
        DO I=1,ORDER(ELEM)+1
            MAT_PHI(1,2*I-1)=PHI(I)
            MAT_PHI(2,2*I)=PHI(I)
        ENDDO
!    
        r(1)=(X-X0)
        r(2)=(Y-Y0)
        llrll=DSQRT(r(1)**2.D0+r(2)**2.D0)
        drdx=r/llrll
        drdn=(drdx(1)*nv(1)+drdx(2)*nv(2))  
        LONGITUDINAL_MODULUS=Young_Modulus(SUBREGION(ELEM))
        Nu=Poisson(SUBREGION(ELEM))
        SHEAR_MODULUS=LONGITUDINAL_MODULUS/(2*(1+Nu))
        C1=1.D0/(8.D0*pi*SHEAR_MODULUS*(1-Nu))*(4*Nu-3.D0)
!        
        DO I=1,2
            DO J=1,2
                Ustar(I,J)=1.D0/(8.D0*pi*SHEAR_MODULUS*(1-Nu))*((4*Nu-3.D0)*DLOG(llrll)*IDENTITY(I,J)+drdx(I)*drdx(J))                
                Tstar(I,J)=-1.D0/(4.D0*pi*(1.D0-Nu)*llrll)*(drdn*((1.D0-2.D0*Nu)*IDENTITY(I,J)+&
                    2.D0*drdx(I)*drdx(J))+(1.D0-2.D0*Nu)*(drdx(J)*nv(I)-drdx(I)*nv(J)))
            ENDDO
        ENDDO
!         
        llrllaux=JAC0*(COORD_GAUSS(II)-XSI0)      
        drdn0=0.D0
!        
        DO I=1,2
            DO J=1,2
                Ustar_aux(I,J)=1.D0/(8.D0*pi*SHEAR_MODULUS*(1-Nu))*(4*Nu-3.D0)*DLOG(ABS(llrllaux))*IDENTITY(I,J)               
                Tstar_aux(I,J)=-1.D0/(4.D0*pi*(1.D0-Nu)*llrllaux)*(drdn0*((1.D0-2.D0*Nu)*IDENTITY(I,J)+&
                    2.D0*drdxlx0(I)*drdxlx0(J))+(1.D0-2.D0*Nu)*(drdxlx0(J)*n0(I)-drdxlx0(I)*n0(J)))
            ENDDO
        ENDDO
!        
        ge=ge+(matmul(Ustar,MAT_PHI)-matmul(Ustar_aux,MAT_PHI0))*JAC*WEIGHT_GAUSS(II)
        he=he+(matmul(Tstar,MAT_PHI)-matmul(Tstar_aux,MAT_PHI0))*JAC*WEIGHT_GAUSS(II)
    ENDDO
!
! CAUCHY PRINCIPAL VALUE (CPV) INTEGRATION FOR THE STRONG SINGULAR INTEGRALS 
!   
    IF(ABS(XSI0-1.D0).GT.1.D-14)THEN !IF(XSI0.NE.1.D0)THEN
        CPV_lnJ0XSI=(1.D0+XSI0)*DLOG(JAC0*(1.D0+XSI0))+(1.D0-XSI0)*DLOG(JAC0*(1.D0-XSI0))-2.D0
    ELSE
        CPV_lnJ0XSI=2.D0*DLOG(2.D0*JAC0)-2.D0
    ENDIF
!    
    IF(ABS(XSI0+1.D0).LE.1.D-14)THEN !IF(XSI0.EQ.-1.D0)THEN
        CPV_1lXSI=DLOG(2.D0)
    ELSE IF(ABS(XSI0-1.D0).LE.1.D-14)THEN !ELSE IF(XSI0.EQ.1.D0)THEN
        CPV_1lXSI=-DLOG(2.D0)
    ELSE
        CPV_1lXSI=DLOG(1.D0-XSI0)-DLOG(1.D0+XSI0)        
    ENDIF
!
    DO I=1,2
        DO J=1,2
            Tstar_bar(I,J)=-1.D0/(4.D0*pi*(1.D0-Nu))*(drdn0*((1.D0-2.D0*Nu)*IDENTITY(I,J)+&
                2.D0*drdxlx0(I)*drdxlx0(J))+(1.D0-2.D0*Nu)*(drdxlx0(J)*n0(I)-drdxlx0(I)*n0(J)))
        ENDDO
    ENDDO  
!
    ge_weaksingular=C1*MAT_PHI0*JAC0*CPV_lnJ0XSI
    he_strongsingular=matmul(Tstar_bar,MAT_PHI0)*CPV_1lXSI
 !   
    ge=ge+ge_weaksingular
    he=he+he_strongsingular
!      
    DEALLOCATE(VALUES,PHI)
!  
END SUBROUTINE SINGULARINTEGRATIONS_DBIE

! =====================================================================================!
![5]                                                                                   !
SUBROUTINE SINGULARINTEGRATIONS_TBIE(NODE,ELEM,ge,he)                                !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!    
    IMPLICIT NONE
    INTEGER::II,I,J,K,NODE,ELEM,LOCAL_NODE
    REAL(KIND=8),INTENT(INOUT)::ge(2,2*(ORDER(ELEM)+1)),he(2,2*(ORDER(ELEM)+1))
    REAL(KIND=8),ALLOCATABLE,DIMENSION(:)::PHI,PHI0,dPHI0
    REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)::VALUES
    REAL(KIND=8)::X,Y,JAC,JAC0,&
    nv(2),n0(2),n0_(2),pi,r(2),drdx(2),drdxlx0(2),drdn,drdn0,llrll,llrllaux,Dstar(2,2,2),Sstar(2,2,2),&
    Dstar_aux(2,2,2),Dstar_bar(2,2,2),Sstar_aux(2,2,2),Sstar_bar(2,2,2),D(4,2),S(4,2),D_aux(4,2),S_aux(4,2),&
    D_bar(4,2),S_bar(4,2),IDENTITY(2,2),MAT_PHI(2,2*(ORDER(ELEM)+1)),MAT_PHI0(2,2*(ORDER(ELEM)+1)),&
    MAT_dPHI0(2,2*(ORDER(ELEM)+1)),MAT_PHI_aux(2,2*(ORDER(ELEM)+1)),CPV,HFP,XSI0,X0,Y0,DIST,NORMALOUTWARD_SIGN
    REAL(KIND=8)::ge_aux(4,2*(ORDER(ELEM)+1)),he_aux(4,2*(ORDER(ELEM)+1))
    REAL(KIND=8)::ge_strongsingular(2,2*(ORDER(ELEM)+1)),he_hypersingular(2,2*(ORDER(ELEM)+1)),TOL
    REAL(KIND=8)::SHEAR_MODULUS,LONGITUDINAL_MODULUS,Nu
!   
    IDENTITY=0.D0
    IDENTITY(1,1)=1.D0
    IDENTITY(2,2)=1.D0
    pi=DACOS(-1.D0)
    TOL=1.D-10
!
    LOCAL_NODE=0
    DO K=1,ORDER(ELEM)+1
        IF(CONECTIVITY(ELEM,K).EQ.NODE)THEN
            LOCAL_NODE=K
            NORMALOUTWARD_SIGN=1.D0
        ENDIF
    ENDDO
    IF(LOCAL_NODE.EQ.0)THEN
        DO K=1,ORDER(ELEM)+1
            DIST=DSQRT((COORDCOLOCATIONPOINTS(CONECTIVITY(ELEM,K),1)-COORDCOLOCATIONPOINTS(NODE,1))**2.D0+&
                (COORDCOLOCATIONPOINTS(CONECTIVITY(ELEM,K),2)-COORDCOLOCATIONPOINTS(NODE,2))**2.D0)
            IF(DIST.LT.TOL)THEN
                LOCAL_NODE=K
                NORMALOUTWARD_SIGN=-1.D0
            ENDIF
        ENDDO
    ENDIF
    !ALLOCATE(VALUES((ORDER(ELEM)+1),2),PHI(ORDER(ELEM)+1),PHI0(ORDER(ELEM)+1),dPHI0(ORDER(ELEM)+1))    
    ALLOCATE(VALUES((ORDER(ELEM)+1),2),PHI(ORDER(ELEM)+1),PHI0(ORDER(ELEM)+1),dPHI0(ORDER(ELEM)+1))
    DO II=1,ORDER(ELEM)+1
        VALUES(II,1)=COORDNODES(CONECTIVITY(ELEM,II),1)  !THE MATRIX VALUES RECIVE THE X OF VALUE OF THE "II" LOCAL NODE OF THE ELEMENT "ELEM"
        VALUES(II,2)=COORDNODES(CONECTIVITY(ELEM,II),2)  !THE MATRIX VALUES RECIVE THE Y OF VALUE OF THE "II" LOCAL NODE OF THE ELEMENT "ELEM"
    ENDDO 
    XSI0=XSI_COLLOCATIONPOINTS(ELEM,LOCAL_NODE)
    CALL JACOBIAN_AND_NORMAL_VECTOR (XSI0,ELEM,VALUES,JAC0,n0)
    n0_=n0*NORMALOUTWARD_SIGN
    CALL SHAPE_FUNCTIONS_DERIVATIVES(XSI0,ELEM,dPHI0)
    drdxlx0=0.D0
    DO II=1,ORDER(ELEM)+1
        drdxlx0(1)=drdxlx0(1)+dPHI0(II)*VALUES(II,1)
        drdxlx0(2)=drdxlx0(2)+dPHI0(II)*VALUES(II,2)
    ENDDO
    drdxlx0=drdxlx0/JAC0
    CALL BASIS_FUNCTIONS (XSI0,ELEM,PHI0)
    MAT_PHI0=0.D0
    DO II=1,ORDER(ELEM)+1
        MAT_PHI0(1,2*II-1)=PHI0(II)
        MAT_PHI0(2,2*II)=PHI0(II)
    ENDDO  
    CALL BASIS_FUNCTIONS_DERIVATIVES(XSI0,ELEM,dPHI0)
    MAT_dPHI0=0.D0
    DO II=1,ORDER(ELEM)+1
        MAT_dPHI0(1,2*II-1)=dPHI0(II)
        MAT_dPHI0(2,2*II)=dPHI0(II)
    ENDDO 
!
    X0=COORDCOLOCATIONPOINTS(NODE,1)
    Y0=COORDCOLOCATIONPOINTS(NODE,2)
!
! NUMERICAL INTEGRATION OF THE REGULAR PARTS
! 
    ge_aux=0.D0
    he_aux=0.D0
    DO II=1,N_GAUSS
        CALL COORDINATES (COORD_GAUSS(II),ELEM,VALUES,X,Y) !CALL COORDINATES AND INFORM XSI OF A NODE OF THE ELEMENT, THE ELEMENT, THE VALUES MATRIX AND IT RETURN THE X AND Y COORDNATES
        CALL JACOBIAN_AND_NORMAL_VECTOR (COORD_GAUSS(II),ELEM,VALUES,JAC,nv) !CALL JACOBIAN AND INFORM XSI OF A NODE OF THE ELEMENT, THE ELEMENT, THE VALUE OF CSI AND IT RETURN THE JACOBIAN AND A NORMAL VECTOR N
        CALL BASIS_FUNCTIONS (COORD_GAUSS(II),ELEM,PHI)
        MAT_PHI=0.D0
        DO I=1,ORDER(ELEM)+1
            MAT_PHI(1,2*I-1)=PHI(I)
            MAT_PHI(2,2*I)=PHI(I)
        ENDDO
!    
        r(1)=(X-X0)
        r(2)=(Y-Y0)
        llrll=DSQRT(r(1)**2.D0+r(2)**2.D0)
        drdx=r/llrll
        drdn=(drdx(1)*nv(1)+drdx(2)*nv(2))    
        LONGITUDINAL_MODULUS=Young_Modulus(SUBREGION(ELEM))
        Nu=Poisson(SUBREGION(ELEM))
        SHEAR_MODULUS=LONGITUDINAL_MODULUS/(2*(1+Nu))
!        
        DO I=1,2
            DO J=1,2
                DO K=1,2
                    Dstar(I,J,K)=1.D0/(4.D0*pi*(1-Nu)*llrll)*((1.D0-2*Nu)*(IDENTITY(I,J)*drdx(K)+IDENTITY(I,K)*drdx(J)-&
                        IDENTITY(J,K)*drdx(I))+2*drdx(I)*drdx(J)*drdx(K))                                                      
                    Sstar(I,J,K)=2.D0*drdn*((1.D0-2.D0*Nu)*IDENTITY(I,J)*drdx(K)+Nu*(IDENTITY(I,K)*drdx(J)+&
                        IDENTITY(J,K)*drdx(I))-4.D0*drdx(I)*drdx(J)*drdx(K))
                    Sstar(I,J,K)=Sstar(I,J,K)+2.D0*Nu*(nv(I)*drdx(J)*drdx(K)+nv(J)*drdx(I)*drdx(K))-&
                    (1.D0-4.D0*Nu)*nv(K)*IDENTITY(I,J)+(1.D0-2.D0*Nu)*(nv(I)*IDENTITY(J,K)+nv(J)*IDENTITY(I,K)+&
                        2.D0*nv(K)*drdx(I)*drdx(J))  
                    Sstar(I,J,K)=2.D0*SHEAR_MODULUS/(4.D0*pi*(1.D0-Nu)*llrll**2.D0)*Sstar(I,J,K)  
                ENDDO
            ENDDO
        ENDDO
        D(1,1)=Dstar(1,1,1)
        D(1,2)=Dstar(2,1,1)
        D(2,1)=Dstar(1,2,1)
        D(2,2)=Dstar(2,2,1)
        D(3,1)=Dstar(1,1,2)
        D(3,2)=Dstar(2,1,2)
        D(4,1)=Dstar(1,2,2)
        D(4,2)=Dstar(2,2,2)
!  
        S(1,1)=Sstar(1,1,1)
        S(1,2)=Sstar(2,1,1)
        S(2,1)=Sstar(1,2,1)
        S(2,2)=Sstar(2,2,1)
        S(3,1)=Sstar(1,1,2)
        S(3,2)=Sstar(2,1,2)
        S(4,1)=Sstar(1,2,2)
        S(4,2)=Sstar(2,2,2)
!
        llrllaux=JAC0*(COORD_GAUSS(II)-XSI0) 
        drdn0=0.D0 
        MAT_PHI_aux=0.D0
        DO I=1,ORDER(ELEM)+1
            MAT_PHI_aux(1,2*I-1)=PHI0(I)+dPHI0(I)*(COORD_GAUSS(II)-XSI0)
            MAT_PHI_aux(2,2*I)=PHI0(I)+dPHI0(I)*(COORD_GAUSS(II)-XSI0)
        ENDDO      
!        
        DO I=1,2
            DO J=1,2
                DO K=1,2
                    Dstar_aux(I,J,K)=1.D0/(4.D0*pi*(1-Nu)*llrllaux)*((1.D0-2*Nu)*(IDENTITY(I,J)*drdxlx0(K)+&
                        IDENTITY(I,K)*drdxlx0(J)-IDENTITY(J,K)*drdxlx0(I))+2*drdxlx0(I)*drdxlx0(J)*drdxlx0(K))                                                        
                    Sstar_aux(I,J,K)=2.D0*drdn0*((1.D0-2.D0*Nu)*IDENTITY(I,J)*drdxlx0(K)+Nu*(IDENTITY(I,K)*drdxlx0(J)+&
                        IDENTITY(J,K)*drdxlx0(I))-4.D0*drdxlx0(I)*drdxlx0(J)*drdxlx0(K))
                    Sstar_aux(I,J,K)=Sstar_aux(I,J,K)+2.D0*Nu*(n0(I)*drdxlx0(J)*drdxlx0(K)+n0(J)*drdxlx0(I)*drdxlx0(K))-&
                    (1.D0-4.D0*Nu)*n0(K)*IDENTITY(I,J)+(1.D0-2.D0*Nu)*(n0(I)*IDENTITY(J,K)+n0(J)*IDENTITY(I,K)+&
                        2.D0*n0(K)*drdxlx0(I)*drdxlx0(J))  
                    Sstar_aux(I,J,K)=2.D0*SHEAR_MODULUS/(4.D0*pi*(1.D0-Nu)*llrllaux**2.D0)*Sstar_aux(I,J,K)  
                ENDDO
            ENDDO
        ENDDO
        D_aux(1,1)=Dstar_aux(1,1,1)
        D_aux(1,2)=Dstar_aux(2,1,1)
        D_aux(2,1)=Dstar_aux(1,2,1)
        D_aux(2,2)=Dstar_aux(2,2,1)
        D_aux(3,1)=Dstar_aux(1,1,2)
        D_aux(3,2)=Dstar_aux(2,1,2)
        D_aux(4,1)=Dstar_aux(1,2,2)
        D_aux(4,2)=Dstar_aux(2,2,2)
!  
        S_aux(1,1)=Sstar_aux(1,1,1)
        S_aux(1,2)=Sstar_aux(2,1,1)
        S_aux(2,1)=Sstar_aux(1,2,1)
        S_aux(2,2)=Sstar_aux(2,2,1)
        S_aux(3,1)=Sstar_aux(1,1,2)
        S_aux(3,2)=Sstar_aux(2,1,2)
        S_aux(4,1)=Sstar_aux(1,2,2)
        S_aux(4,2)=Sstar_aux(2,2,2)
!  
        ge_aux=ge_aux+(matmul(D,MAT_PHI)*JAC-matmul(D_aux,MAT_PHI0)*JAC0)*WEIGHT_GAUSS(II)
        he_aux=he_aux+(matmul(S,MAT_PHI)*JAC-matmul(S_aux,MAT_PHI_aux)*JAC0)*WEIGHT_GAUSS(II)          
    ENDDO    
    DO I=1,2*(ORDER(ELEM)+1)
        DO J=1,2
            ge(J,I)=ge_aux(2*J-1,I)*n0_(1)+ge_aux(2*J,I)*n0_(2)
            he(J,I)=he_aux(2*J-1,I)*n0_(1)+he_aux(2*J,I)*n0_(2)
        ENDDO
    ENDDO
!
! NONSINGULAR PART OF Dstar_aux and Sstar_aux
! 
    DO I=1,2
        DO J=1,2
            DO K=1,2
                Dstar_bar(I,J,K)=1.D0/(4.D0*pi*(1-Nu)*JAC0)*((1.D0-2*Nu)*(IDENTITY(I,J)*drdxlx0(K)+&
                    IDENTITY(I,K)*drdxlx0(J)-IDENTITY(J,K)*drdxlx0(I))+2*drdxlx0(I)*drdxlx0(J)*drdxlx0(K))                                                        
                Sstar_bar(I,J,K)=2.D0*drdn0*((1.D0-2.D0*Nu)*IDENTITY(I,J)*drdxlx0(K)+Nu*(IDENTITY(I,K)*drdxlx0(J)+&
                    IDENTITY(J,K)*drdxlx0(I))-4.D0*drdxlx0(I)*drdxlx0(J)*drdxlx0(K))
                Sstar_bar(I,J,K)=Sstar_bar(I,J,K)+2.D0*Nu*(n0(I)*drdxlx0(J)*drdxlx0(K)+n0(J)*drdxlx0(I)*drdxlx0(K))-&
                (1.D0-4.D0*Nu)*n0(K)*IDENTITY(I,J)+(1.D0-2.D0*Nu)*(n0(I)*IDENTITY(J,K)+n0(J)*IDENTITY(I,K)+&
                    2.D0*n0(K)*drdxlx0(I)*drdxlx0(J))  
                Sstar_bar(I,J,K)=2.D0*SHEAR_MODULUS/(4.D0*pi*(1.D0-Nu)*JAC0**2.D0)*Sstar_bar(I,J,K)  
            ENDDO
         ENDDO
    ENDDO
    D_bar(1,1)=Dstar_bar(1,1,1)
    D_bar(1,2)=Dstar_bar(2,1,1)
    D_bar(2,1)=Dstar_bar(1,2,1)
    D_bar(2,2)=Dstar_bar(2,2,1)
    D_bar(3,1)=Dstar_bar(1,1,2)
    D_bar(3,2)=Dstar_bar(2,1,2)
    D_bar(4,1)=Dstar_bar(1,2,2)
    D_bar(4,2)=Dstar_bar(2,2,2)
!  
    S_bar(1,1)=Sstar_bar(1,1,1)
    S_bar(1,2)=Sstar_bar(2,1,1)
    S_bar(2,1)=Sstar_bar(1,2,1)
    S_bar(2,2)=Sstar_bar(2,2,1)
    S_bar(3,1)=Sstar_bar(1,1,2)
    S_bar(3,2)=Sstar_bar(2,1,2)
    S_bar(4,1)=Sstar_bar(1,2,2)
    S_bar(4,2)=Sstar_bar(2,2,2)
!
! CAUCHY PRINCIPAL VALUE (CPV) INTEGRATION FOR THE STRONG SINGULAR INTEGRALS 
!
    IF((XSI0+1.D0).LE.1.D-14)THEN
        CPV=DLOG(2.D0)
    ELSE IF((XSI0-1.D0).LE.1.D-14)THEN
        CPV=-DLOG(2.D0)
    ELSE
        CPV=DLOG(1.D0-XSI0)-DLOG(1.D0+XSI0)
    ENDIF
!
! HADAMARD FINITE PART (HFP) INTEGRATION FOR THE HYPERSINGULAR SINGULAR INTEGRALS 
!           
    IF((XSI0+1.D0).LE.1.D-14)THEN
        HFP=-0.5D0
    ELSE IF((XSI0-1.D0).LE.1.D-14)THEN
        HFP=-0.5D0
    ELSE 
        HFP=(-1.D0)/(1.D0+XSI0)+(-1.D0)/(1.D0-XSI0)
    ENDIF
!
    ge_aux=matmul(D_bar,MAT_PHI0)*JAC0*CPV
    he_aux=matmul(S_bar,MAT_PHI0)*JAC0*HFP+matmul(S_bar,MAT_dPHI0)*JAC0*CPV
!
    DO I=1,2*(ORDER(ELEM)+1)
        DO J=1,2
            ge_strongsingular(J,I)=ge_aux(2*J-1,I)*n0_(1)+ge_aux(2*J,I)*n0_(2)
            he_hypersingular(J,I)=he_aux(2*J-1,I)*n0_(1)+he_aux(2*J,I)*n0_(2)
        ENDDO
    ENDDO
!
    ge=ge+ge_strongsingular
    he=he+he_hypersingular
!    
    DEALLOCATE(VALUES,PHI)
!     
END SUBROUTINE SINGULARINTEGRATIONS_TBIE
!
END MODULE MODULE_NUMERICAL_INTEGRATION