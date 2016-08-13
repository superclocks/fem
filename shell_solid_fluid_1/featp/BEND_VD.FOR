       
	 
	 SUBROUTINE BEND_VD(VD,VMAE)
            IMPLICIT  REAL *8(A-H, O-Z)
            DIMENSION  VD(5,5) ,VMAE(4)
            E=VMAE(1)                  !   E-ELASTIC MODULUS
            V=VMAE(2)                  !   V-POSSION'S RATIO

            TH=VMAE(4)
            D0=E*TH*TH*TH/12/(1-V*V)
            VD(1,1)=D0
            VD(2,2)=D0                            !  MPROB= 4--MINDLIN PLATE
            VD(3,3)=D0*(1-V)/2
            VD(1,2)=D0*V
            VD(2,1)=D0*V
            VD(4,4)=E/2/(1+V)*TH/(6.0/5.0)
            VD(5,5)=E/2/(1+V)*TH/(6.0/5.0)
       
            RETURN
            END