MODULE MODULE_INPUT_DATA
! -------------------------------------------------------------------------------------!
!   Description.                                                                       !
!                                                                                      !
! -------------------------------------------------------------------------------------!
!   This module contains the following subroutines:                                    !
!                                                                                      !
!   [1]. INPUT_QUADRATURE_DATA                                                         !
!   [2]. INPUT_ANALYSIS_PARAMETERS                                                     !
!   [3]. INPUT_BOUNDARY_MESH                                                           !
!   [4]. INPUT_BOUNDARY_CONDITIONS                                                     !
!   [5]. INPUT_MATERIALS                                                               !
!   [6]. INPUT_POSTPROCESSING_MESH                                                     !                                                
!   [7]. INPUT_DATA                                                                    !
! -------------------------------------------------------------------------------------!

! -------------------------------------------------------------------------------------!
!   USED DATA:                                                                         !
! -------------------------------------------------------------------------------------!
USE MODULE_DATA
! USE BOUNDARY_MESH_DATA
! USE POSTPROCESSING_MESH_DATA
! USE BOUNDARY_CONDITIONS_DATA
! USE ANALYSIS_PARAMETERS_DATA
! USE MATERIAL_PROPERTIES_DATA

! -------------------------------------------------------------------------------------!
!   USED MODULE:                                                                       !
! -------------------------------------------------------------------------------------!
!USE NUMERICAL_INTEGRATION_DATA
USE MODULE_NUMERICAL_INTEGRATION

CONTAINS

! ========================================================================================!
![1]
SUBROUTINE INPUT_QUADRATURE_DATA()    
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!   
!
	IMPLICIT NONE 
!
    REAL(KIND=wp)::xi,xf
    REAL(KIND=wp),ALLOCATABLE::GaussPos(:),GaussWgt(:)
!	
	ALLOCATE(COORD_GAUSS(N_GAUSS),WEIGHT_GAUSS(N_GAUSS))
    xi=-1.D0
    xf=1.D0
!	
    CALL GAUSS_LEGENDRE(GaussPos,GaussWgt,xi,xf,N_GAUSS)    

    COORD_GAUSS=GaussPos   
    WEIGHT_GAUSS=GaussWgt 
    DEALLOCATE(GaussPos,GaussWgt)
!
END SUBROUTINE INPUT_QUADRATURE_DATA     
    
! ========================================================================================!
![2]
SUBROUTINE INPUT_ANALYSIS_PARAMETERS    
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!   
!
	IMPLICIT NONE 
! 
    INTEGER::I,II
!
	OPEN(1,file='4_ANALYSIS_PARAMETERS.txt',status='unknown')
!     	 
	READ(1,*)!**************** ANALYSIS PARAMETERS DEFINITIONS *****************
	READ(1,*)!WARNING: RESPECT THIS FILE FORMAT.THE SPACE BETWEEN DATA
	READ(1,*)!AND THE NUMBER FORMAT
	READ(1,*)!
    READ(1,*)!NUMBER OF GAUSS POINTS
	READ(1,*)N_GAUSS
	READ(1,*)!
    READ(1,*)!COLLOCATION POINT DISTANCING
	READ(1,*)CP_DISTANCING
	READ(1,*)!
    READ(1,*)!PLANE STRAIN OR PLANE STRESS CONDITION
	READ(1,*)PLANE_STATE
	READ(1,*)!
    READ(1,*)!ANALYSIS TYPE: LINEAR ELASTIC (LE) OR COHESIVE PROPAGATION (CP)
	READ(1,*)ANALYSIS_TYPE
	READ(1,*)!
    READ(1,*)!NUMBER OF STEPS IN COHESIVE PROPAGATION ANALYSIS
	READ(1,*)Ns
	READ(1,*)!
    READ(1,*)!TOLERANCE FOR CONVERGENCE OF THE NEWTON-RAPHSON
	READ(1,*)tol_Newton
	READ(1,*)!
    READ(1,*)!MAXIMUM NUMBER OF ITERACTIONS
	READ(1,*)itemax
	READ(1,*)!
    READ(1,*)!NODE FOR THE EVALUATION OF THE STRUCTURAL RESPONSE
	READ(1,*)NODE
	READ(1,*)!
    READ(1,*)!NODE FOR THE EVALUATION OF THE COHESIVE RESPONSE
	READ(1,*)COHESIVE_NODE
	READ(1,*)!
!
!	ALLOCATING SOME VARIABLES CONSIDERING THE INFORMED VALUES
!
	ALLOCATE(SURF_TYPE(N_ELEM),SURF_IE(N_ELEM),COLLOCATION_TYPE(N_NODES),COLLOCATION_IE(N_NODES))
!
	READ(1,*)!ELEMENT AND SURFACE TYPE 
	READ(1,*)!ELEMENT      SURFACE TYPE  (B = BOUNDARY Ip = INTERFACE plus Im = INTERFACE minus)    SURFACE IE  (D = Displacement BIE    T = Traction BIE) 
	DO I=1,N_ELEM
		READ(1,*)II,SURF_TYPE(II),SURF_IE(II)
    ENDDO
!
    CLOSE(1) 
!
! COMPUTING GAUSS COORDINATES AND WEIGHTS
!
    CALL INPUT_QUADRATURE_DATA
!
END SUBROUTINE INPUT_ANALYSIS_PARAMETERS    
    
! ========================================================================================!
![3]
SUBROUTINE INPUT_BOUNDARY_MESH    
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!   
!
	IMPLICIT NONE 
! 
    INTEGER::I,J,K,II,MAXORDER,N_DUPLICATED_NODES,CONT,CONT_Im,CONT_Ip
    REAL(KIND=8)::SPACING,Z
    CHARACTER(LEN=1),ALLOCATABLE,DIMENSION(:)::DUPLICATED_NODES
!
	OPEN(1,file='1_BOUNDARY_MESH.txt',status='unknown')
!      	     	 
	READ(1,*)!**************** BOUNDARY MESH DEFINITIONS *****************
	READ(1,*)!WARNING: RESPECT THIS FILE FORMAT.THE SPACE BETWEEN DATA
	READ(1,*)!AND THE NUMBER FORMAT
	READ(1,*)!
	READ(1,*)!NUMBER OF SUBREGIONS
	READ(1,*)N_SUBREGIONS
	READ(1,*)!
	READ(1,*)!NUMBER OF BOUNDARY ELEMENTS
	READ(1,*)N_ELEM
	READ(1,*)!
    READ(1,*)!NUMBER OF NODES 
	READ(1,*)N_NODES
	READ(1,*)!
    READ(1,*)!NUMBER OF DUPLICATE NODES 
	READ(1,*)N_DUPLICATED_NODES
	READ(1,*)!
!
!	ALLOCATING DUPLICATED_NODES
!
    ALLOCATE(DUPLICATED_NODES(N_NODES+N_DUPLICATED_NODES))
    DUPLICATED_NODES='S'
!
	READ(1,*)!DUPLICATED NODES
	READ(1,*)!NODE
!
	DO I=1,N_DUPLICATED_NODES
		READ(1,*)II
        DUPLICATED_NODES(II)='D'
    ENDDO
 !   
    DO I=N_NODES+1,N_NODES+N_DUPLICATED_NODES
        DUPLICATED_NODES(I)='D'
	ENDDO
	READ(1,*)!    
!    
    N_NODES=N_NODES+N_DUPLICATED_NODES
!
!	ALLOCATING OTHER VARIABLES
!
	ALLOCATE(COORDNODES(N_NODES,2),COORDCOLOCATIONPOINTS(N_NODES,2),ORDER(N_ELEM),SUBREGION(N_ELEM),ELEM_TYPE(N_ELEM))
!     
	COORDNODES=0.D0
	ORDER=0
    ELEM_TYPE='C'
!
	READ(1,*)!NODAL COORDINATES
	READ(1,*)!NODE		X		Y		Z
!
	DO I=1,N_NODES-N_DUPLICATED_NODES
		READ(1,*)II,COORDNODES(II,1),COORDNODES(II,2),Z
    ENDDO
	READ(1,*)!
    CONT=0
    DO I=1,N_NODES-N_DUPLICATED_NODES
        IF(DUPLICATED_NODES(I).EQ.'D')THEN
            CONT=CONT+1
            II=N_NODES-N_DUPLICATED_NODES+CONT
			COORDNODES(II,:)=COORDNODES(I,:)
        ENDIF
    ENDDO
!    
	READ(1,*)!INFORM THE POLYNOMIAL ORDER FOR ALL ELEMENTS
	READ(1,*)!ELEMENT	       ORDER	        SUBREGION  
!
	DO I=1,N_ELEM
		READ(1,*)II,ORDER(II),SUBREGION(II)
	ENDDO
	READ(1,*)!
!
	MAXORDER=0
	DO I=1,N_ELEM
		IF(ORDER(I).GE.MAXORDER) THEN
			MAXORDER=ORDER(I)
		ENDIF
    ENDDO
    MAX_ORDER=MAXORDER
!
	ALLOCATE(CONECTIVITY(N_ELEM,(MAXORDER+1)),XSI_NODES(N_ELEM,(MAXORDER+1)),XSI_COLLOCATIONPOINTS(N_ELEM,(MAXORDER+1)))
	CONECTIVITY=0
    XSI_NODES=0.D0
!    
	READ(1,*)!INFORM THE CONECTIVITY FOR ALL ELEMENTS
	READ(1,*)!ELEMENT		NODES
!	
	DO I=1,N_ELEM
		READ(1,*)(CONECTIVITY(I,J),J=1,(ORDER(I)+1))
    ENDDO
!
    CALL INPUT_ANALYSIS_PARAMETERS
!
    CONT=0
    DO I=1,N_NODES-N_DUPLICATED_NODES
        IF(DUPLICATED_NODES(I).EQ.'D')THEN
            CONT=CONT+1
			DO J=1,N_ELEM
				IF(CONECTIVITY(J,1).EQ.I)THEN
					CONECTIVITY(J,1)=N_NODES-N_DUPLICATED_NODES+CONT
				ENDIF
			ENDDO
		ENDIF
    ENDDO
!
!DEFINING THE COLLOCATION_TYPE AND COLLOCATION_IE
!
    DO I=1,N_ELEM
		IF(SURF_TYPE(I).EQ."B")THEN
            DO J=1,(ORDER(I)+1)
                COLLOCATION_TYPE(CONECTIVITY(I,J))="B"    
            ENDDO
        ELSE IF(SURF_TYPE(I).EQ."Ip")THEN
            DO J=1,(ORDER(I)+1)
                COLLOCATION_TYPE(CONECTIVITY(I,J))="Ip"                   
            ENDDO                
        ELSE IF(SURF_TYPE(I).EQ."Im")THEN
            DO J=1,(ORDER(I)+1)
                COLLOCATION_TYPE(CONECTIVITY(I,J))="Im"                  
            ENDDO    
        ENDIF
    ENDDO  
!
	DO I=1,N_NODES
        J=0
		DO K=1,N_ELEM
			DO II=1,(ORDER(K)+1)
				IF(CONECTIVITY(K,II).EQ.I) THEN
					J=K
				ENDIF
            ENDDO
        ENDDO
        COLLOCATION_IE(I)=SURF_IE(J)            
	ENDDO
!    
    DO I=1,N_ELEM
		IF(DUPLICATED_NODES(CONECTIVITY(I,1)).EQ.'D'.AND.DUPLICATED_NODES(CONECTIVITY(I,(ORDER(I)+1))).EQ.'D')THEN
            ELEM_TYPE(I)='D'    
        ELSE IF(DUPLICATED_NODES(CONECTIVITY(I,1)).EQ.'D')THEN 
            ELEM_TYPE(I)='LD'                
        ELSE IF(DUPLICATED_NODES(CONECTIVITY(I,(ORDER(I)+1))).EQ.'D')THEN
            ELEM_TYPE(I)='RD'               
        ENDIF
        IF(COLLOCATION_TYPE(CONECTIVITY(I,1)).NE."B")THEN
            ELEM_TYPE(I)='D'           
        ENDIF
    ENDDO
!    
	CLOSE(1)    
!
!DEFINING THE XSI COORDINATES FOR EACH NODE
!
    DO I=1,N_ELEM
        SPACING=2.D0/ORDER(I)
        DO J=1,(ORDER(I)+1)
            XSI_NODES(I,J)=-1.D0+SPACING*(J-1)
        ENDDO
    ENDDO 
!
!DEFINING THE NUMBER OF COHESIVE PAIR OF NODES
!   
    ALLOCATE(VETOR_NODE_B_TOT(N_NODES),VETOR_NODE_Im_TOT(N_NODES),VETOR_NODE_Ip_TOT(N_NODES))
    N_COHESIVE_NODES=0
    CONT = 0
    CONT_Im = 0
    CONT_Ip = 0
    DO I=1,N_NODES
        IF(COLLOCATION_TYPE(I).EQ."Im")THEN
            CONT_Im = CONT_Im + 1
            VETOR_NODE_Im_TOT(CONT_Im) = I
            N_COHESIVE_NODES=N_COHESIVE_NODES+1
        ELSEIF(COLLOCATION_TYPE(I).EQ."Ip")THEN
            CONT_Ip = CONT_Ip + 1
            VETOR_NODE_Ip_TOT(CONT_Ip) = I
            N_COHESIVE_NODES=N_COHESIVE_NODES+1
        ELSE
            CONT = CONT + 1
            VETOR_NODE_B_TOT(CONT) = I
        ENDIF
    ENDDO     
    ALLOCATE(VETOR_NODE_B(CONT),VETOR_NODE_Im(CONT_Im),VETOR_NODE_Ip(CONT_Ip))
    VETOR_NODE_B(:)=VETOR_NODE_B_TOT(1:CONT)
    VETOR_NODE_Im(:)=VETOR_NODE_Im_TOT(1:CONT_Im)
    VETOR_NODE_Ip(:)=VETOR_NODE_Ip_TOT(1:CONT_Ip)
    DEALLOCATE(VETOR_NODE_B_TOT,VETOR_NODE_Im_TOT,VETOR_NODE_Ip_TOT,DUPLICATED_NODES)
!
END SUBROUTINE INPUT_BOUNDARY_MESH    

! ========================================================================================!
![4]
SUBROUTINE INPUT_BOUNDARY_CONDITIONS    
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!  
!
	IMPLICIT NONE 
! 
    INTEGER::I,IB,J,II,CONT,CONT_ON,CONT_OFF,I_AUX,K
!
	OPEN(1,file='3_BOUNDARY_CONDITIONS.txt',status='unknown')
!      	     	 
	READ(1,*)!**************** BOUNDARY CONDITIONS DEFINITIONS *****************
	READ(1,*)!WARNING: RESPECT THIS FILE FORMAT.THE SPACE BETWEEN DATA
	READ(1,*)!AND THE NUMBER FORMAT
	READ(1,*)!
    READ(1,*)!NUMBER OF ELEMENTS WITH Ux PRESCRIBED
	READ(1,*)NU1
	READ(1,*)!
    READ(1,*)!NUMBER OF ELEMENTS WITH Uy PRESCRIBED
	READ(1,*)NU2
	READ(1,*)!
    READ(1,*)!NUMBER OF ELEMENTS WITH Tx PRESCRIBED
	READ(1,*)NT1
	READ(1,*)!
    READ(1,*)!NUMBER OF ELEMENTS WITH Ty PRESCRIBED
	READ(1,*)NT2
	READ(1,*)!
!
!	ALLOCATING SOME VARIABLES
!
    ALLOCATE(BC_COEFFICIENTS(2*N_ELEM,3),U(2*N_NODES),DU_BAR(2*N_NODES),T(2*N_NODES),BC(2*N_NODES),Beta(N_NODES))
    DU_BAR=0.D0
    Beta=0.D0
!
    BC=0
    BC_COEFFICIENTS=0
    READ(1,*)!BOUNDARY CONDITIONS
	READ(1,*)!II	Aux	  Bux	 Cux   
    DO I=1,NU1
		READ(1,*)II,BC_COEFFICIENTS(2*II-1,1),BC_COEFFICIENTS(2*II-1,2),BC_COEFFICIENTS(2*II-1,3)
        DO J=1,ORDER(II)+1
			BC(2*CONECTIVITY(II,J)-1)=1
		ENDDO    
    ENDDO
	READ(1,*)
    READ(1,*)!II	Auy	  Buy	 Cuy   
    DO I=1,NU2
		READ(1,*)II,BC_COEFFICIENTS(2*II,1),BC_COEFFICIENTS(2*II,2),BC_COEFFICIENTS(2*II,3)
        DO J=1,ORDER(II)+1
			BC(2*CONECTIVITY(II,J))=1
		ENDDO    
    ENDDO
	READ(1,*)
    READ(1,*)!II	Atx	  Btx	 Ctx   
    DO I=1,NT1
		READ(1,*)II,BC_COEFFICIENTS(2*II-1,1),BC_COEFFICIENTS(2*II-1,2),BC_COEFFICIENTS(2*II-1,3)  
    ENDDO
	READ(1,*)
    READ(1,*)!II	Aty	  Bty	 Cty   
    DO I=1,NT2
		READ(1,*)II,BC_COEFFICIENTS(2*II,1),BC_COEFFICIENTS(2*II,2),BC_COEFFICIENTS(2*II,3) 
    ENDDO
    READ(1,*)
    READ(1,*)!CONTROLED PARAMETER (U=DISPLACEMENT, T=TRACTION, DU=CRACK_RELATIVE_DISPLACEMENT)   
    READ(1,*)!TYPE    NODE   DIRECTION      VALUE 
    READ(1,*)CONTROL_TYPE,CONTROLED_NODE,CONTROLED_DIRECTION,CONTROLED_VALUE 
!    
    CLOSE(1)     
!
    NDOF_U=0
    NDOF_DU=0
    DO I=1,N_NODES
        IF(COLLOCATION_TYPE(I).EQ."B")THEN
            IF(BC(2*I-1).EQ.1)THEN
                NDOF_U=NDOF_U+1  
            ENDIF
            IF(BC(2*I).EQ.1)THEN
                NDOF_U=NDOF_U+1  
            ENDIF  
        ELSE IF(COLLOCATION_TYPE(I).EQ."Ip")THEN
            NDOF_DU=NDOF_DU+2 
        ENDIF
    ENDDO
    NDOF_T=2*N_NODES-NDOF_U-2*NDOF_DU
!
    CONT=0
    IF(CONTROL_TYPE.EQ."U")THEN
        DO I=1,N_NODES
            IF(COLLOCATION_TYPE(I).EQ."B")THEN
                IF(BC(2*I-1).EQ.0)THEN
                    CONT=CONT+1
                ENDIF
                IF(I.EQ.CONTROLED_NODE)THEN
                    IF(CONTROLED_DIRECTION.EQ.1)THEN
                        CONTROLED_DOF=CONT    
                    ELSE IF(CONTROLED_DIRECTION.EQ.2)THEN
                        CONTROLED_DOF=CONT+1                 
                    ENDIF
                ENDIF
                IF(BC(2*I).EQ.0)THEN
                    CONT=CONT+1
                ENDIF
            ENDIF
        ENDDO
    ENDIF
    IF(CONTROL_TYPE.EQ."T")THEN
        CONT=NDOF_T
        DO I=1,N_NODES
            IF(COLLOCATION_TYPE(I).EQ."B")THEN
                IF(BC(2*I-1).EQ.1)THEN
                    CONT=CONT+1
                ENDIF
                IF(I.EQ.CONTROLED_NODE)THEN
                    IF(CONTROLED_DIRECTION.EQ.1)THEN
                        CONTROLED_DOF=CONT    
                    ELSE IF(CONTROLED_DIRECTION.EQ.2)THEN
                        CONTROLED_DOF=CONT+1                 
                    ENDIF
                ENDIF
                IF(BC(2*I).EQ.1)THEN
                    CONT=CONT+1
                ENDIF
            ENDIF
        ENDDO
    ENDIF
    IF(CONTROL_TYPE.EQ."DU")THEN
        CONT=NDOF_T+NDOF_U+NDOF_DU
        DO I=1,N_NODES
            IF(COLLOCATION_TYPE(I).EQ."Ip")THEN
                CONT=CONT+2
            ENDIF
            IF(I.EQ.CONTROLED_NODE)THEN
                IF(CONTROLED_DIRECTION.EQ.1)THEN
                    CONTROLED_DOF=CONT-1    
                ELSE IF(CONTROLED_DIRECTION.EQ.2)THEN
                    CONTROLED_DOF=CONT                 
                ENDIF
            ENDIF
        ENDDO         
    ENDIF


    ALLOCATE(VETOR_B_GDL_OFF_TOT(CONT*2),VETOR_B_GDL_ON_TOT(CONT*2))
    CONT_OFF = 0
    CONT_ON = 0
    DO IB=1,SIZE(VETOR_NODE_B)
        I = VETOR_NODE_B(IB)
        IF(BC(2*I-1).EQ.0)THEN
            CONT_OFF=CONT_OFF+1
            VETOR_B_GDL_OFF_TOT(CONT_OFF)=2*I-1
        ELSE
            CONT_ON=CONT_ON+1
            VETOR_B_GDL_ON_TOT(CONT_ON)=2*I-1
        ENDIF
        IF(BC(2*I).EQ.0)THEN
            CONT_OFF=CONT_OFF+1
            VETOR_B_GDL_OFF_TOT(CONT_OFF)=2*I
        ELSE
            CONT_ON=CONT_ON+1
            VETOR_B_GDL_ON_TOT(CONT_ON)=2*I
        ENDIF
    ENDDO

    ALLOCATE(VETOR_B_GDL_OFF(CONT_OFF),VETOR_B_GDL_ON(CONT_ON))
    VETOR_B_GDL_OFF(:) = VETOR_B_GDL_OFF_TOT(1:CONT_OFF)
    VETOR_B_GDL_ON(:) = VETOR_B_GDL_ON_TOT(1:CONT_ON)
    DEALLOCATE(VETOR_B_GDL_OFF_TOT,VETOR_B_GDL_ON_TOT)

    ALLOCATE(VET_Ip_ELEM_ORDER(SIZE(VETOR_NODE_Ip),3))

    DO I_AUX=1,SIZE(VETOR_NODE_Ip)
        I = VETOR_NODE_Ip(I_AUX)
        DO J=1,N_ELEM
            DO K=1,ORDER(J)+1
                IF(I.EQ.CONECTIVITY(J,K))THEN
                    VET_Ip_ELEM_ORDER(I_AUX,1) = I
                    VET_Ip_ELEM_ORDER(I_AUX,2) = J
                    VET_Ip_ELEM_ORDER(I_AUX,3) = K
                ENDIF
            ENDDO
        ENDDO
    ENDDO

    !WRITE(*,*)"VET_Ip_ELEM_ORDER = ",VET_Ip_ELEM_ORDER

    ! WRITE(*,*)"VETOR_B_GDL_ON = ",VETOR_B_GDL_ON
    ! WRITE(*,*)"VETOR_B_GDL_OFF = ",VETOR_B_GDL_OFF
!
 END SUBROUTINE INPUT_BOUNDARY_CONDITIONS    

! ========================================================================================!
![5]
SUBROUTINE INPUT_MATERIALS    
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!   
!
	IMPLICIT NONE 
!
    INTEGER::I,II
	OPEN(1,file='2_MATERIALS.txt',status='unknown')
!
    ALLOCATE(Young_Modulus(N_SUBREGIONS),Poisson(N_SUBREGIONS))
!        
	READ(1,*)!**************** MATERIALS DEFINITIONS *****************
	READ(1,*)!WARNING: RESPECT THIS FILE FORMAT.THE SPACE BETWEEN DATA
	READ(1,*)!AND THE NUMBER FORMAT
	READ(1,*)!
	READ(1,*)!INFORM THE MATERIAL PROPERTIES
	READ(1,*)!SUBREGION   Young_Modulus       Poisson 
    DO I=1,N_SUBREGIONS
		READ(1,*)II,Young_Modulus(II),Poisson(II)
	ENDDO
    READ(1,*)!
	READ(1,*)!INFORM THE COHESIVE PROPERTIES FOR MODE I
	READ(1,*)!Kn0          phin          smax                   
    READ(1,*)kn0,phin,smax
    READ(1,*)!
	READ(1,*)!INFORM THE COHESIVE PROPERTIES FOR MODE II
	READ(1,*)!Kto          phit          tmax             
    READ(1,*)kt0,phit,tmax
    ALLOCATE(Knn(N_NODES),Knt(N_NODES),Ktn(N_NODES),Ktt(N_NODES),Beta_MAX(N_NODES))
    Knn=0.D0
    Knt=0.D0
    Ktn=0.D0
    Ktt=0.D0
    Beta_MAX=0.D0
!    
    DO I=1,N_NODES
        IF(COLLOCATION_TYPE(I).NE."B")THEN
            Knn(I)=kn0
            Knt(I)=0.D0
            Ktn(I)=0.D0
            Ktt(I)=kt0
        ENDIF
    ENDDO
!    
    CLOSE(1)
!
! UPDATING THE ELASTIC CONSTANTS FOR PLANE STRESS (PT) CONDITION
!
    IF(PLANE_STATE.EQ."PT")THEN
        Young_Modulus=Young_Modulus*(1.D0-(Poisson/(1.D0+Poisson))**2.D0)
        Poisson=Poisson/(1.D0+Poisson)
    ENDIF
!
! COMPUTING COHESIVE PARAMETERS
!   
    deltan=2*phin/smax
    deltanc=smax/kn0
    Etan=(1.D0-deltanc/deltan)
!
    deltat=2*phit/tmax
    deltatc=tmax/kt0
    Etat=(1.D0-deltatc/deltat)
!
END SUBROUTINE INPUT_MATERIALS    

! ========================================================================================!
![6]
SUBROUTINE INPUT_POSTPROCESSING_MESH   
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!        
!
	IMPLICIT NONE 
! 
    INTEGER::I,J,II
    REAL(KIND=wp)::Z
!
	OPEN(1,file='5_POSTPROCESSING_MESH.txt',status='unknown')
!     	 
	READ(1,*)!**************** POSTPROCESSING MESH DEFINITIONS *****************
	READ(1,*)!WARNING: RESPECT THIS FILE FORMAT.THE SPACE BETWEEN DATA
	READ(1,*)!AND THE NUMBER FORMAT
	READ(1,*)!
	READ(1,*)!NUMBER OF ELEMENTS
	READ(1,*)N_POSPROCESSINGELEM
	READ(1,*)!
    READ(1,*)!NUMBER OF POINTS 
	READ(1,*)N_POINTS
	READ(1,*)!
    READ(1,*)!NUMBER OF BOUNDARY POINTS
	READ(1,*)NB
	READ(1,*)!
    NI=N_POINTS-NB
!
!	ALLOCATING SOME VARIABLES CONSIDERING THE INFORMED VALUES
!
	ALLOCATE(COORDBOUNDARYPOINTS(NB,2),COORDINTERNALPOINTS(NI,2),SUBREGION_POINTS(N_POINTS))
!     
	COORDBOUNDARYPOINTS=0.D0
	COORDINTERNALPOINTS=0
!
	READ(1,*)!POINTS COORDINATES
	READ(1,*)!		X		Y
!
	DO I=1,N_POINTS
        IF(I.LE.NB)THEN
			READ(1,*)II,COORDBOUNDARYPOINTS(I,1),COORDBOUNDARYPOINTS(I,2),Z,SUBREGION_POINTS(I)
        ELSE
			READ(1,*)II,COORDINTERNALPOINTS(I-NB,1),COORDINTERNALPOINTS(I-NB,2),Z,SUBREGION_POINTS(I)
        ENDIF
	ENDDO
	READ(1,*)!
!
	ALLOCATE(POSTPROCESSING_CONECTIVITY(N_POSPROCESSINGELEM,4))
	POSTPROCESSING_CONECTIVITY=0
!    
	READ(1,*)!INFORM THE CONECTIVITY FOR ALL POSTPROCESSING ELEMENTS
	READ(1,*)!NODES
!	
	DO I=1,N_POSPROCESSINGELEM
		READ(1,*)(POSTPROCESSING_CONECTIVITY(I,J),J=1,3)
	ENDDO
!
    CLOSE(1)    
!    
END SUBROUTINE INPUT_POSTPROCESSING_MESH

! ========================================================================================!
![7]
SUBROUTINE INPUT_DATA     
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!  
!
	IMPLICIT NONE 
!
	CALL INPUT_BOUNDARY_MESH    
!
	CALL INPUT_BOUNDARY_CONDITIONS
!
	CALL INPUT_MATERIALS
!
    CALL INPUT_POSTPROCESSING_MESH
!
END SUBROUTINE INPUT_DATA  
!
END MODULE MODULE_INPUT_DATA