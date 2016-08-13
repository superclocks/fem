!!	P732--733

C  .........................  SUB: 3-1-6  ..........................
            SUBROUTINE ELEMENT_VD(MPROB,VD,VMAE)
C ****FORM  ELEMENT ELASTIC MATRIX [D'] ACCORDING TO TYPE OF PROBLEM
            IMPLICIT  REAL *8(A-H, O-Z)
            DIMENSION  VD(5,5) ,VMAE(4)
            E=VMAE(1)                  !   E-ELASTIC MODULUS
            V=VMAE(2)                  !   V-POSSION'S RATIO
C ******** CLEAR MEMORY FOR THE MATRIX: [D]
            DO 30 I=1,5
            DO 30 J=1,5
            VD(I,J)=0.0                 !      VD--ELASTIC MATRIX
30          CONTINUE
C ****COMPUTE [D] MATRIX FOR PLANE STRESS OR PLANE STRAIN PROBLEM
C ***
            MPROB=1 
            IF(MPROB.EQ.1.OR.MPROB.EQ.2)THEN
            IF(MPROB.EQ.2) E=E/(1-V*V)
            IF(MPROB.EQ.2) V=V/(1-V)
            D0=E/(1-V*V)
            VD(1,1)=D0                  ! MPROB= 1-PLANE STRESS
            VD(2,2)=D0
            VD(3,3)=D0*(1-V)/2           ! MPROB= 2-PLANE STRAIN
            VD(1,2)=D0*V
            VD(2,1)=D0*V
            ENDIF
C ***** COMPUTE [D] MATRIX FOR AXISYMMETRIC PROBLEM
C ***
	      IF(MPROB.EQ.3)THEN
            D0=E*(1-V)/(1+V)/(1-2*V)
            VD(1,1)=D0
            VD(2,2)=D0
            VD(3,3)=D0*(1-2*V)/2/(1-V)
            VD(4,4)=D0
            VD(2,1)=D0*V/(1-V)          ! MPROB = 3-AXISYMMETRIC
            VD(1,2)=D0*V/(1-V)
            VD(4,1)=D0*V/(1-V)
            VD(1,4)=D0*V/(1-V)
            VD(4,2)=D0*V/(1-V)
            VD(2,4)=D0*V/(1-V)
            ENDIF
C ********* COMPUTE [D] MATRIX FOR MINDLIN PLATE PROBLEM
            IF(MPROB.EQ.4)THEN
            TH=VMAE(4)
            D0=E*TH*TH*TH/12/(1-V*V)
            VD(1,1)=D0
            VD(2,2)=D0                            !  MPROB= 4--MINDLIN PLATE
            VD(3,3)=D0*(1-V)/2
            VD(1,2)=D0*V
            VD(2,1)=D0*V
            VD(4,4)=E/2/(1+V)*TH/(6.0/5.0)
            VD(5,5)=E/2/(1+V)*TH/(6.0/5.0)
            ENDIF
            RETURN
            END