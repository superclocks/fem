     
	SUBROUTINE TRI_TRANS(VX,TCOS)
      IMPLICIT  REAL* 8 (A-H, O-Z)
	DIMENSION VX(3,3),TCOS(3,3),AX(3),AY(3),AZ(3),A(3),AA(3,3),A1(3)

      AX(1)=VX(2,1)-VX(1,1) 
	AX(2)=VX(2,2)-VX(1,2) 
	AX(3)=VX(2,3)-VX(1,3)   !局部坐标下的X轴方向余弦(12矢量)

      A(1)=(VX(3,1)-VX(2,1)) 
	A(2)=(VX(3,2)-VX(2,2))  
	A(3)=(VX(3,3)-VX(2,3))   !局部坐标下的边23方向余弦

	AZ(1)=AX(2)*A(3)-AX(3)*A(2)
      AZ(2)=AX(3)*A(1)-AX(1)*A(3)
      AZ(3)=AX(1)*A(2)-AX(2)*A(1)   !局部坐标下的Z轴方向余弦(与面垂直)

	AY(1)=AX(3)*AZ(2)-AX(2)*AZ(3)
      AY(2)=AX(1)*AZ(3)-AX(3)*AZ(1)
	AY(3)=AX(2)*AZ(1)-AX(1)*AZ(2)  !局部坐标下的Y轴方向余弦

      XL=SQRT(AX(1)*AX(1)+AX(2)*AX(2)+AX(3)*AX(3))
	AX(1)=AX(1)/XL
      AX(2)=AX(2)/XL
	AX(3)=AX(3)/XL

	YL=SQRT(AY(1)*AY(1)+AY(2)*AY(2)+AY(3)*AY(3))
      AY(1)=AY(1)/YL
      AY(2)=AY(2)/YL
	AY(3)=AY(3)/YL
	
	ZL=SQRT(AZ(1)*AZ(1)+AZ(2)*AZ(2)+AZ(3)*AZ(3))
	AZ(1)=AZ(1)/ZL
      AZ(2)=AZ(2)/ZL
	AZ(3)=AZ(3)/ZL
      
	TCOS(:,1)=AX(:)
	TCOS(:,2)=AY(:)
	TCOS(:,3)=AZ(:)

	AA(:,1)=VX(:,1)-VX(1,1)
	AA(:,2)=VX(:,2)-VX(1,2)
      AA(:,3)=VX(:,3)-VX(1,3)

	DO 10 I=1,3
	DO 10 J=1,3
	VX(I,J)=0.0
      DO 10 K=1,3
	VX(I,J)=VX(I,J)+AA(I,K)*TCOS(K,J)
10    CONTINUE

	RETURN
	END