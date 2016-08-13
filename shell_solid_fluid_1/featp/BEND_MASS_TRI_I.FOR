

	SUBROUTINE  BEND_MASS_TRI_I(VX,BM,VMAE,X,Y,Z)
	IMPLICIT  REAL* 8 (A-H, O-Z)
      COMMON/COM1/MND, NUMEL, NUMPT, MBAND, NMATI
      COMMON/COM2/NF, NFSTR, MSOLV, MPROB, MTYPE, NVA
      COMMON/COM3/MND2 , NUMPT2
      COMMON/ELEM/NODE, INTX, INTY
              
	DIMENSION  VX(3,3),BM(9,9),VMAE(4),AB(3),AC(3),BN(9)
     $ 
      
	AB(1)=VX(2,2)-VX(3,2)
	AB(2)=VX(3,2)-VX(1,2)
	AB(3)=VX(1,2)-VX(2,2)
	
	AC(1)=VX(3,1)-VX(2,1)
	AC(2)=VX(1,1)-VX(3,1)
	AC(3)=VX(2,1)-VX(1,1)
	
	BN(1)=X+X*X*Y+X*X*Z-X*Y*Y-X*Z*Z
	BN(2)=AB(2)*(X*X*Z+0.5*X*Y*Z)-AB(3)*(X*X*Y+0.5*X*Y*Z)
	BN(3)=AC(2)*(X*X*Z+0.5*X*Y*Z)-AC(3)*(X*X*Y+0.5*X*Y*Z) 
	BN(4)=Y+Y*Y*Z+X*Y*Y-Y*Z*Z-X*X*Y
      BN(5)=AB(3)*(X*Y*Y+0.5*X*Y*Z)-AB(1)*(Y*Y*Z+0.5*X*Y*Z)
	BN(6)=AC(3)*(X*Y*Y+0.5*X*Y*Z)-AC(1)*(Y*Y*Z+0.5*X*Y*Z)
	BN(7)=Z+X*Z*Z+Y*Z*Z-X*X*Z-Y*Y*Z
	BN(8)=AB(1)*(Y*Z*Z+0.5*X*Y*Z)-AB(2)*(X*Z*Z+0.5*X*Y*Z)
	BN(9)=AC(1)*(Y*Z*Z+0.5*X*Y*Z)-AC(2)*(X*Z*Z+0.5*X*Y*Z)

      BM(:,:)=0.0
	DO 10 I=1,9
	DO 10 J=1,9
	BM(I,J)=BM(I,J)+BN(I)*BN(J)
10    CONTINUE
	RETURN							   
	END
