! =====================================================================================!
!   MULTIDOMAIN BEM WITH GENERAL COHESIVE MODELS - PLANE STRAIN AND STRESS ANALYSIS    !   
!                                                                                      !
!	INSTITUTO TECNOLÓGICO DE AERONÁUTICA                                               !                                          
!	DIVISÃO DE ENGENHARIA CIVIL AERONÁUTICA                                            !                              
!   LABORATÓRIO DE MODELAGEM ESTRUTURAL                                                !                    
!	PRAÇA MARECHAL EDUARDO GOMES, 50. SÃO JOSÉ DOS CAMPOS-SP, BRAZIL.                  !
!   POSTAL CODE 12.228-900. PHONE/FAX +(55) 12 3947-6801                               !
!                                                                                      !                  
!   AUTHOR: SÉRGIO GUSTAVO FERREIRA CORDIRO 				                           !	       
!	cordeiro@ita.br                                                                    !
!   https://www.civil.ita.br/~cordeiro/	                                               !
!                                                                                      !
!	THIS COMPUTATIONAL CODE WAS DEVELOPED TO PERFORM LINEAR ELASTIC AND COHESIVE       !
!   FRACTURE ANALYSIS OF MULTIDOMAIN PROBLEMS IN PLANE STRAIN AND PLANE STRESS PROBLEMS!                                                     
!                                                                                      !
! =====================================================================================!
    PROGRAM MAIN_COHESIVE_DBEM_2D                                                           !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!   
    USE MODULE_DATA                                                                                   !
    USE MODULE_INPUT_DATA                                                              !
    USE MODULE_DATA_TREATMENT                                                       !
    USE MODULE_ANALYSIS                                                                !
    USE MODULE_POSTPROCESSING    
    USE omp_lib                                     !
!                                                                                      !                                            
    IMPLICIT NONE                                                                      !
	INTEGER::I            
    INTEGER::proc_num, thread_num                                                             !
	REAL(KIND=wp)::TI,TF                                                                      !                   
                              
    ! =====================================================================
    proc_num = omp_get_num_procs ( )
    thread_num = omp_get_max_threads ( )

    WRITE(*,*)'  The number of processors available = ', proc_num
    WRITE(*,*)'  The number of threads available    = ', thread_num
    ! =====================================================================
                                                              !
    CALL CPU_TIME(TI)                                                                  !
                                                                                       !
    CALL INPUT_DATA                                                                    !
!                                                                                      !
    CALL TREATMENT_OF_DATA                                                             !
!                                                                                      !
    CALL ANALYSIS                                                                      !
!                                                                                      !
    CALL POSTPROCESSING                                                                !
!                                                                                      !
    DO I=1,22                                                                          !
        WRITE(*,*)''                                                                   !
    ENDDO                                                                              !
!                                                                                      !
    WRITE(*,*)''                                                                       !
    WRITE(*,*)'SOLUTION IS DONE'                                                       !  
    WRITE(*,*)''                                                                       !
!                                                                                      !
    CALL CPU_TIME(TF)                                                                  !
!                                                                                      !
    WRITE(*,*)''                                                                       !
    WRITE(*,100)'COMPUTATIONAL WORK, SECONDS',(TF-TI)                                  !
    !PAUSE                                                                              !
!                                                                                      !
    100 FORMAT(A28,4x,E12.6)             
                                                  !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!        
    END PROGRAM MAIN_COHESIVE_DBEM_2D                                                       !
! =====================================================================================!





