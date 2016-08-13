
      SUBROUTINE SOLID_TRANS(SXY,SCOS)
      IMPLICIT  REAL* 8 (A-H,O-Z)
	DIMENSION SXY(8,3),SX(8,3),SCOS(3,3),AX(3),AY(3),AZ(3),AA(3)

      AX(1)=SXY(2,1)-SXY(1,1) 
	AX(2)=SXY(2,2)-SXY(1,2)  
	AX(3)=SXY(2,3)-SXY(1,3)    !局部坐标下的X方向余弦

	AA(1)=SXY(3,1)-SXY(2,1) 
	AA(2)=SXY(3,2)-SXY(2,2)  
	AA(3)=SXY(3,3)-SXY(2,3)    !局部坐标下的工作坐标方向余弦,帮助形成Z方向余弦

	AZ(1)=AA(2)*AX(3)-AA(3)*AX(2)
      AZ(2)=AA(3)*AX(1)-AA(1)*AX(3)
      AZ(3)=AA(1)*AX(2)-AA(2)*AX(1)   !局部坐标下的Z轴方向余弦(与面垂直)

      AY(1)=AZ(2)*AX(3)-AZ(3)*AX(2)
	AY(2)=AZ(3)*AX(1)-AZ(1)*AX(3)
	AY(3)=AZ(1)*AX(2)-AZ(2)*AX(1)

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

      SCOS(:,1)=AX(:)
	SCOS(:,2)=AY(:)
	SCOS(:,3)=AZ(:)

      SX(:,1)=SXY(:,1)-SXY(1,1)   
	SX(:,2)=SXY(:,2)-SXY(1,2)  
      SX(:,3)=SXY(:,3)-SXY(1,3) 

	DO 10 I=1,8
	DO 10 J=1,3
	SXY(I,J)=0.0
      DO 10 K=1,3
	SXY(I,J)=SXY(I,J)+SX(I,K)*SCOS(K,J)
10    CONTINUE

	END