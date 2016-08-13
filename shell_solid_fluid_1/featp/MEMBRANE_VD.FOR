
            SUBROUTINE MEMBRANE_VD(VD,VMAE) 
            IMPLICIT  REAL *8(A-H, O-Z)
            DIMENSION  VD(3,3) ,VMAE(4)
            E=VMAE(1)                  !   E-ELASTIC MODULUS
            V=VMAE(2)                  !   V-POSSION'S RATIO
	      
		  VD(:,:)=0.0
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
            RETURN
            END