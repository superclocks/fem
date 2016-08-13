
	SUBROUTINE MEMBRANE_I(VX,X,Y,Z,AK,VMAE,CC)
	USE IMSL
	IMPLICIT REAL * 8(A-H,O-Z)
              
	DIMENSION VX(3,2),TR(3,3),BM(3,5),CC(5,9),AX(3),AY(3),AL(3),
     $          AXY(3),
     $          AA(3),AB(3),AC(3),AD(3),AE(3),VN(6),VDN(2,6),
     $          AH(2,9),ADH(4,9),B(3,9),SB(5,3),VMAE(4),VD(3,3),AK(5,5)
         
      TR(1,1)=VX(2,1)*VX(3,2)-VX(3,1)*VX(2,2)
	TR(1,2)=VX(2,2)-VX(3,2)
	TR(1,3)=VX(3,1)-VX(2,1)
	TR(2,1)=VX(3,1)*VX(1,2)-VX(1,1)*VX(3,2)
	TR(2,2)=VX(3,2)-VX(1,2)
	TR(2,3)=VX(1,1)-VX(3,1)
	TR(3,1)=VX(1,1)*VX(2,2)-VX(2,1)*VX(1,2)
	TR(3,2)=VX(1,2)-VX(2,2)
	TR(3,3)=VX(2,1)-VX(1,1)
	SA=VX(2,1)*VX(3,2)-VX(3,1)*VX(2,2)-VX(1,1)*VX(3,2)+
     $   VX(3,1)*VX(1,2)+VX(1,1)*VX(2,2)-VX(2,1)*VX(1,2)

	TR(:,:)=TR(:,:)/ABS(SA)
	TR(:,:)=.I.TR
	AXY(1)=TR(1,1)*(X)+TR(1,2)*Y+TR(1,3)*Z
      AXY(2)=TR(2,1)*(X)+TR(2,2)*Y+TR(2,3)*Z
	AXY(3)=TR(3,1)*(X)+TR(3,2)*Y+TR(3,3)*Z
      
	BM(:,:)=0.0
	BM(1,1)=1.0
	BM(1,2)=AXY(3)
	BM(2,3)=1.0
	BM(2,4)=AXY(2)
	BM(3,2)=-AXY(2)
	BM(3,4)=-AXY(3)
	BM(3,5)=1.0

      AX(1)=VX(1,1)-VX(2,1)
	AX(2)=VX(2,1)-VX(3,1)
	AX(3)=VX(3,1)-VX(1,1)

	AY(1)=VX(1,2)-VX(2,2)
	AY(2)=VX(2,2)-VX(3,2)
      AY(3)=VX(3,2)-VX(1,2)

	CC(:,:)=0.0

	CC(1,1)=AY(2)
	CC(1,3)=0.5*VX(1,2)*AY(2)
	CC(1,4)=AY(3)
	CC(1,6)=0.5*VX(2,2)*AY(3)
	CC(1,7)=AY(1)
	CC(1,9)=0.5*VX(3,2)*AY(1)

	CC(2,3)=-0.5*AY(2)
	CC(2,6)=-0.5*AY(3)
	CC(2,9)=-0.5*AY(1)

	CC(3,2)=-AX(2)
	CC(3,3)= 0.5*VX(1,1)*AX(2)
	CC(3,5)=-AX(3)
	CC(3,6)= 0.5*VX(2,1)*AX(3)
	CC(3,8)=-AX(1)
	CC(3,9)= 0.5*VX(3,1)*AX(1)

	CC(4,3)=-0.5*AX(2)
	CC(4,6)=-0.5*AX(3)
	CC(4,9)=-0.5*AX(1)

	CC(5,1)=-AX(2)
	CC(5,2)= AY(2)
	CC(5,3)= 0.5*(VX(1,1)*(-AY(2))+VX(1,2)*(-AX(2)))
	CC(5,4)=-AX(3)
	CC(5,5)= AY(3)
	CC(5,6)= 0.5*(VX(2,1)*(-AY(3))+VX(2,2)*(-AX(3)))
      CC(5,7)=-AX(1)
	CC(5,8)= AY(1)
	CC(5,9)= 0.5*(VX(3,1)*(-AY(1))+VX(3,2)*(-AX(1)))

	CC(:,:)=CC(:,:)/ABS(SA)
      CALL MEMBRANE_VD(VD,VMAE)
      
	SB(:,:)=0.0
	DO 10 I=1,5
	DO 10 J=1,3
	DO 10 KK=1,3
      SB(I,J)=SB(I,J)+BM(KK,I)*VD(KK,J)
10    CONTINUE
      
	AK(:,:)=0.0
      DO 11 I=1,5
	DO 11 J=1,5
	DO 11 KK=1,3
      AK(I,J)=AK(I,J)+SB(I,KK)*BM(KK,J)*ABS(SA)*VMAE(4)*0.5
11    CONTINUE
      RETURN
	END
