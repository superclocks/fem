! P 712--715
C===============================================================================
C=============================== SUB: 1  =======================================
         SUBROUTINE ALLOCAT(M1,M2,M3,M4,M5,M8,M11,M12,N1,N2,N3,N4,N5,

     $  N6,N7,N8, N9, N10, N11,N12,N13,N14,N15,N16,N17,N18,N19,N20,N21,
     $    N22,N23,N24,N25,N26,N27 ,N31 ,N32,N33 ,N34 ,N35, N36 ,N37,N38)
C*******************************************************************************
C       INPUT BASIC PARAMETERS FROM FILE 'IN.DAT'
C              ALLOCATE DYANMICAL STORAGE SPACE
C*******************************************************************************
             IMPLICIT REAL * 8(A-H,O-Z)
               COMMON/COM1/MND, NUMEL, NUMPT,MBAND, NMATI
               COMMON/COM2/NF, NFSTR, MSOLV, MPROB, MTYPE, NVA
               COMMON/COMN/NFIX, NPC, GRAV
               COMMON/COM3/MND2,NUMPT2
C                   ! MND--MAXIMAL NODE NUMBER IN ALL ELEMENTS
C                   ! NUMEL--NUMBER OF GLOBAL ELEMENTS
C                   ! NUMPT--NUMBER OF GLOBAL NODES
C                   ! MBAND--HALF BANDWIDTH(INCLUDING DIAGONAL ELEMENT)
C###########从in.dat中读出第一和第二行数据#############  
                   
		  READ(5, * )
            READ(5, * )  MND,NUMEL,NUMPT,MBAND  
            READ(5, * )
            READ(5, * ) NFIX, NPc, MPROB, MSOLV
C######################################################
C                   ! NFIX--NUMBER OF NODES SUBJECTED TO CONSTRIANT
C                   ! NPC--NUMBER OF NODES SUBJECTED TO EQUIVALENT LOAD
C===============================================================================
C         MPROB= 1--PLANE STRESS PROBLEM,      MPROB=2--PLANE STRAIN PROBLEM
C         MPROB= 3--AXISYMMETRIC PROBLEM,      MPROB=4--MINDLIN PLATE PROBLEM
C-------------------------------------------------------------------------------
C         MSOLV= 1-STATIC ANALYSIS
C         MSOLV= 2-DYNAMIC RESPONSE ANALYSIS BY CENTRAL DIFFERENCE METHOD
C         MSOLV= 3-DYNAMIC RESPONSE ANALYSIS BY NEWMARK METHOD
C         MSOLV= 4--DYNAMIC CHARACTER ANALYSIS BY INVERSE ITERATION METHOD
C         MSOLV= 5-DYNAMIC CHARACTER ANALYSIS BY SUBSPACE ITERATION METHOD
C===============================================================================
C###########从in.dat中读出第三行数据######################
      READ(5, * )
      IF(MSOLV. NE. 4. OR. MSOLV. NE. 5) READ(5,*) NMATI,GRAV,MTYPE
      IF(MSOLV. EQ. 4. OR. MSOLV. EQ. 5) READ(5,*) NMATI,GRAV,MTYPE,NVA
C#########################################################
C               ! NMATI--KIND OF MATERIALS
C               ! GRAV--GRAVITY ACCELERATION
C               ! NVA--NUMBER OF EIGENVALUES
C-------------------------------------------------------------------------------
C         MTYPE--CONTROL KEY FOR OUTPUT RESULTS
C         MTYPE= O--OUTPUT RESULTS INCLUDE GLOBAL MATRIXES AND STRESS AT
C                GAUSS POINTS
C         MTYPE=l-OUTPUT RESULTS INCLUDE GLOBAL MATRIXES
C         MTYPE= 2-OUTPUT RESULTS INCLUDE STRESSES AT GAUSS POINTS
C-------------------------------------------------------------------------------
            IF(MPROB. EQ. 1. OR. MPROB. EQ. 2)    THEN
            NF= 2       !  NF--NUMBER OF NODAL FREEDOMS(DISPLACEMENT COMPONENTS)
            NFSTR=3     !  NFSTR--NUMBER OF STRESS COMPONENTS
            ENDIF
            IF(MPROB. EQ. 3)  THEN
            NF=2
            NFSTR = 4
            ENDIF
            IF(MPROB. EQ. 4)  THEN
            NF=3
            NFSTR = 5
            ENDIF
C*******************************************************************************
C      ALLOCATE STORAGE SPACE FOR INPUT DATA
C*******************************************************************************
            M1=1
            M2=M1
            M3=M2+NUMEL* 14   ! IELEM(NUMEL,14)--CODE OF MATERIAL AND NODES
		                    ! IELME为II单元的结点数  
            M4=M3+NFIX* 4     ! IFIXD(NFIX,4)-CODE OF NODES HAVING CONSTRIANT
							! 第II个约束号所约束的结点号
            M5=M4+NPC * 4     ! ILOAD(NPC,4)--CODE OF NODES HAVING LOAD
							! 第II个等效结点载荷作用的结点号
            N1= 1
            N2=N1+NMATI * 4   ! VAMATI(NMATI,4)--PARAMETERS OF MATERIALS
							
            N3=N2+NUMPT* 3   ! VCOOD(NUMPT,3)--GLOBAL NODE COORDINATES
            N4=N3+NFIX* 3     ! VFIXD(NFIX,3)-CONSTRIANED DISPLACEMENTS
            N5=N4+NPC * 3      ! VLOAD(NPC,3)-VALUES OF EQUIVALENT LOAD
C*******************************************************************************
C     ALLOCATE STORAGE SPACE FOR ELEMENTAL MATRIXES AND GLOBAL MATRIXES
C*******************************************************************************
            MND2= MND* NF              ! MND2--NUMBER OF FREEDOMS tN A ELEMENT
            NUMPT2 = NUMPT * NF        ! NUMPT2--NUMBER OF GLOBAL FREEDOMS
            M8=M5
            M11 = M8 + NUMPT           ! IADD(NUMPT)--USED FOR NODAL STRESS
            M12=M11+MND                !  IEL(MND)--NODE CODE IN A ELEMENT
            N6=N5
            N7=N6+NFSTR * MND2        ! VSG(NFSTR,MND2)-ELEMENT STRESS MATRIX AT
                                     !                INTEGRATION POINTS
            N8=N7+NFSTR* MND2		   ! VSN(NFSTR,MND2)--ELEMENT STRESS AT NODES
            N9 = N8 + NUMPT * NFSTR  ! SSN(NUMPT, NFSTR)-STRESS AT NODES
            N10 = N9+ 9 * NFSTR	 ! SSS(9,NFSTR)-STRESSES AT INTEGRATION POINTS
            N11=N10+4 * NFSTR		   ! VSS(4,NFSTR)-STRESSES AT DESIGNATED
                                     !                INTEGRATION POINTS

            N12=N11+NUMPT2 * MBAND   ! GMM(NUMPT2,MBAND)-GLOBAL MASS
C                                                         !               MATRIX
            N13=N12+NUMPT2 * MBAND   ! GKK(NUMPT2,MBAND)-GLOBAL STIFFNESS
C                                    !          MATRIX
            N14 = N13 + NUMPT2	   !  GP (NUMPT2)-GLOBAL LOAD VECTOR
            N15 = N14 + NUMPT2	   !  GU(NUMPT2)--GLOBAL DISPLACEMENT VECTOR
            N16=N15+MND *3		   !  VXY(MND,3)--NODE COORDINATE IN ELEMENT
            N17=N16+MND2 * MND2	   !  VMM(MND2 ,MND2)-ELEMENT MASS MATRIX
            N18=N17+MND2 * MND2	   !  VKK(MND2,MND2)--ELEMENT STIFFNESS MATRIX
            N19=N18+MND2			   !  VU(MND2)-NODE DISPLACEMENTS IN A ELEMENT
            N20=N19+4				   !  VMAE(4)-MATERIAL PARAMETERS IN A ELEMENT
C*******************************************************************************
C     ALLOCATE STORAGE SPACE FOR DYNAMIC RESPONSE ANALYSIS
C*******************************************************************************
            IF(MSOLV. EQ. 2. OR. MSOLV. EQ. 3) THEN
            N21=N20
            N22 = N21 + NUMPT2 * MBAND    ! DAMP(NUMPT2, MBAND)--DAMPLING MATRIX
            N23 = N22 + NUMPT2        ! U0 (NUMPT2)--INITIAL DISPLACEMENT VECTOR
            N24=N23+NUMPT2                 ! V0(NUMPT2)-INITIAL VELOCITY VECTOR
            N25=N24+NUMPT2               ! A(NUMPT2)-INITIAL ACCELERATION VECTOR
            N26 = N25 + NUMPT2 * MBAND     ! AW (NUMPT2, MBAND)--WORKING ARRAY
            N27 = N26+ NUMPT2              ! B(NUMPT2)--WORKING ARRAY
            ENDIF
C*******************************************************************************
C          ALLOCATE STORAGE SPACE FOR DYNAMIC CHARACTER ANALYSIS
C*******************************************************************************
            IF(MSOLV. EQ. 4. OR. MSOLV. EQ. 5)  THEN
            IF(NVA. EQ. 0)  THEN
            WRITE( * ,11)
11               FORMAT(' PLEASE  READ THE NUMBERS OF EIGENVALUE--NVA=')
            READ( * , * ) NVA
            ENDIF
            N31 = N20
          N32 = N31 + NUMPT2 * NVA  ! AA (NUMPT2, NVA)--INITIAL ITERATION VECTOR
            N33 = N32+NUMPT2 * NVA             ! BB(NUMPT2, NVA)--WORKING ARRAY
            N34=N33+NVA * NVA             ! GM(NVA,NVA)--MASS MATRIX IN SUBSPACE
            N35 =N34+NVA * NVA        ! GK(NVA,NVA)-STIFFNESS MATRIX IN SUBSPACE
            N36=N35+NVA * NVA         ! V(NVA,NVA)--EIGENVECTORS IN SUBSPACE
            N37 = N36 + NVA                       ! W1 ( NVA)-WORKING ARRAY
            N38= N37+NVA                     ! W2(NVA)--EIGENVALUES IN SUBSPACE
            ENDIF
            RETURN
            END