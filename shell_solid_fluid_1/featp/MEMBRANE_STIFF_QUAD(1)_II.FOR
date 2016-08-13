
	SUBROUTINE MEMBRANE_STIFF_QUAD1_II(VXY,QKK1,QMM1,VMAE,STA,QCOS)
	
      IMPLICIT   REAL * 8 (A-H,O-z)
      DIMENSION  IEL(4),VXY(4,3),VK(12,12),VM(12,12),AJ(2,2),AJJ(2,2),
     $        QKK1(12,12),QMM1(12,12),VDN(2,4),QCOS(3,3),
     $        VD(3,3),VMAE(4),AN(2,4),
     $        AA(12,12),AII(12),
     $        AA1(12,12)

      QKK1(:,:)=0.0
	QMM1(:,:)=0.0

	INTX=2
	INTY=2
	DO 40 I=1,INTX
	DO 40 J=1,INTY
	CALL GAUSS_INTEGRATION (INTX,INTY,I,J,X,Y,WXY)
      CALL MEMBRANE_STIFF_QUAD1_I(VXY,VK,VM,X,Y,VMAE,WXY,IEL)
	QKK1(:,:)=QKK1(:,:)+VK(:,:)
      QMM1(:,:)=QMM1(:,:)+VM(:,:)
40    CONTINUE

C******************加入弹簧类型的转动自由度****************
      AN(1,1)=-0.25
	AN(1,2)= 0.25
	AN(1,3)= 0.25
	AN(1,4)=-0.25
      AN(2,1)=-0.25
	AN(2,2)=-0.25
	AN(2,3)= 0.25
	AN(2,4)= 0.25

      DO 41 I=1,2
	DO 41 J=1,2
	AJ(I,J)=0.0
	DO 41 K=1,4
      AJ(I,J)=AJ(I,J)+AN(I,K)*VXY(K,J)
41    CONTINUE
      
	SJ0=AJ(1,1)*AJ(2,2)-AJ(1,2)*AJ(2,1)

      AJJ(1,1)= AJ(2,2)/SJ0
	AJJ(1,2)=-AJ(1,2)/SJ0
	AJJ(2,1)=-AJ(2,1)/SJ0
	AJJ(2,2)= AJ(1,1)/SJ0

	DO 42 I=1,2
	DO 42 J=1,4
	VDN(I,J)=0.0
	DO 42 K=1,2
      VDN(I,J)=VDN(I,J)+AJJ(I,K)*AN(K,J) !形函数对局部坐标的导数
42    CONTINUE

      AII(1)=-VDN(2,1)*0.5
	AII(2)=VDN(1,1)*0.5
	AII(3)=-0.25

      AII(4)=-VDN(2,2)*0.5
	AII(5)=VDN(1,2)*0.5
	AII(6)=-0.25

	AII(7)=-VDN(2,3)*0.5
	AII(8)=VDN(1,3)*0.5
	AII(9)=-0.25

	AII(10)=-VDN(2,4)*0.5
	AII(11)=VDN(1,4)*0.5
	AII(12)=-0.25

	
      DO 45 I=1,12
	DO 45 J=1,12
      AA(I,J)=AII(I)*AII(J)
45    CONTINUE
!      OPEN(111,FILE='AA.DAT',STATUS='UNKNOWN')
!	WRITE(111,112)((AA(I,J),J=1,12),I=1,12)
!112   FORMAT(12E18.8)
      G=VMAE(1)/((1.0+VMAE(2))*2.0)   !剪切模量
	
	S=0.5*(VXY(2,1)*VXY(3,2)-VXY(3,1)*VXY(2,2)+
     $       VXY(3,1)*VXY(1,2)-VXY(1,1)*VXY(3,2)+
     $	   VXY(1,1)*VXY(2,2)-VXY(2,1)*VXY(1,2)+
     $	   VXY(4,1)*VXY(1,2)-VXY(1,1)*VXY(4,2)+
     $	   VXY(1,1)*VXY(3,2)-VXY(3,1)*VXY(1,2)+
     $	   VXY(3,1)*VXY(4,2)-VXY(4,1)*VXY(3,2))
      AA1(:,:)=0.01*ABS(S)*G*AA(:,:)*VMAE(4)
      
 	QKK1(:,:)=QKK1(:,:)+AA1(:,:)  !在原刚度中加入弹簧转动刚度
C***************************************************************

C******************加入一个小的刚度防止数值不稳定***************
      X=G*VMAE(4)*ABS(S)*1E-8
      QKK1(3,3)=QKK1(3,3)+X
	QKK1(6,6)=QKK1(6,6)+X
	QKK1(9,9)=QKK1(9,9)+X
	QKK1(12,12)=QKK1(12,12)+X
C***************************************************************
      RETURN
	END