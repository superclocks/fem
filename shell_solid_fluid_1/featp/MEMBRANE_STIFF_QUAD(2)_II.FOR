
	SUBROUTINE MEMBRANE_STIFF_QUAD2_II(VXY,QKK1,QMM1,VMAE)
	USE IMSL
      IMPLICIT   REAL * 8 (A-H,O-z)
     	COMMON/COM1/MND, NUMEL, NUMPT, MBAND, NMATI
      COMMON/COM2/NF, NFSTR, MSOLV, MPROB, MTYPE, NVA
      COMMON/COM3/MND2 , NUMPT2
      COMMON/ELEM/NODE, INTX, INTY

      DIMENSION  IEL(NODE),VXY(4,3),VK(18,18),QKK1(12,12),QMM1(12,12),
     $        VSS(12,16),AKK(12,12),QM(16,16),akk1(12,12),
     $        VD(3,3),VMAE(4),QK(18,18),VT(16,12),VM(12,12),VP(12,12),
     $       AA(12,12),AAI(8,8),AA2(4,4),AII(12),AII2(4),
     $        AA3(12,12),AKK2(12,12),VX(8,3),VK0(12,16),VK00(12,12),
     $        AN(2,4),AJ(2,2),AJJ(2,2),VDN(2,4),AA1(12,12),EK(16,16),
     $        EK1(2,2),EK2(2,16),EK3(16,4),AAA(12,12),AII1(12)

      QK(:,:)=0.0
	QMM1(:,:)=0.0
      INTX=2
      INTy=2
	DO 40 I=1,INTX
	DO 40 J=1,INTY
	CALL GAUSS_INTEGRATION (INTX,INTY,I,J,X,Y,WXY)
      CALL MEMBRANE_STIFF_QUAD2_I(VXY,VK,VM,X,Y,VMAE,WXY)
	QK(:,:)=QK(:,:)+VK(:,:)
      QMM1(:,:)=QMM1(:,:)+VM(:,:)
40    CONTINUE

      EK1(:,:)=.I.QK(17:18,17:18)
!	X=0.0
!	Y=0.0
!      CALL MEMBRANE_STIFF_QUAD2_I(VXY,VK,VM,X,Y,VMAE,WXY)
	EK2(:,:)=QK(17:18,1:16)
	
	DO 70 I=1,16
	DO 70 J=1,2
	EK3(I,J)=0.0
	DO 70 K=1,2
	EK3(I,J)=EK3(I,J)+EK2(K,I)*EK1(K,J)
70    CONTINUE 

      DO 71 I=1,16
	DO 71 J=1,16
	EK(I,J)=0.0
	DO 71 K=1,2
      EK(I,J)=EK(I,J)+EK3(I,K)*EK2(K,J)
71    CONTINUE

      QK(1:16,1:16)=QK(1:16,1:16)-EK(:,:)

	VT(:,:)=0.0
	VT(1,1)=1.0
	VT(2,2)=1.0
	VT(3,4)=1.0
	VT(4,5)=1.0
	VT(5,7)=1.0
	VT(6,8)=1.0
	VT(7,10)=1.0
	VT(8,11)=1.0
	VT(9,1)=0.5
	VT(9,4)=0.5
	VT(9,3)=-0.125*(VXY(2,2)-VXY(1,2))
	VT(9,6)= 0.125*(VXY(2,2)-VXY(1,2))
	VT(10,2)=0.5
	VT(10,5)=0.5
	VT(10,3)= 0.125*(VXY(2,1)-VXY(1,1))
	VT(10,6)=-0.125*(VXY(2,1)-VXY(1,1))
	VT(11,4)=0.5
	VT(11,7)=0.5
      VT(11,6)=-0.125*(VXY(3,2)-VXY(2,2))
      VT(11,9)= 0.125*(VXY(3,2)-VXY(2,2))
      VT(12,5)=0.5 
	VT(12,8)=0.5
      VT(12,6)= 0.125*(VXY(3,1)-VXY(2,1))
      VT(12,9)=-0.125*(VXY(3,1)-VXY(2,1))
      VT(13,7)=0.5
	VT(13,10)=0.5
	VT(13,9)= -0.125*(VXY(4,2)-VXY(3,2))
      VT(13,12)= 0.125*(VXY(4,2)-VXY(3,2))
	VT(14,8)=0.5 
	VT(14,11)=0.5
	VT(14,9)= 0.125*(VXY(4,1)-VXY(3,1))
      VT(14,12)=-0.125*(VXY(4,1)-VXY(3,1))
	VT(15,1)=0.5
	VT(15,10)=0.5
      VT(15,3)= 0.125*(VXY(1,2)-VXY(4,2)) 
	VT(15,12)= -0.125*(VXY(1,2)-VXY(4,2))
      VT(16,2)=0.5
	VT(16,11)=0.5
      VT(16,3)=-0.125*(VXY(1,1)-VXY(4,1)) 
	VT(16,12)= 0.125*(VXY(1,1)-VXY(4,1))

	DO 45 I=1,12
	DO 45 J=1,16
      VK0(I,J)=0.0
	DO 45 K=1,16
	VK0(I,J)=VK0(I,J)+VT(K,I)*QK(K,J)
45    CONTINUE

      DO 46 I=1,12
	DO 46 J=1,12
	QKK1(I,J)=0.0
	DO 46 K=1,16
      QKK1(I,J)=QKK1(I,J)+VK0(I,K)*VT(K,J)
46    CONTINUE
      OPEN(111,FILE='QKK1.DAT',STATUS='UNKNOWN')
	WRITE(111,112)((QKK1(I,J),J=1,12),I=1,12)
112   FORMAT(12E18.8)
!      AKK(1:3,1:3)=QKK1(10:12,10:12)
!      AKK(1:3,4:12)=QKK1(10:12,1:9)
!      AKK(4:6,1:3)=QKK1(1:3,10:12)
!	AKK(4:6,4:12)=QKK1(1:3,1:9)
!	AKK(7:9,1:3)=QKK1(4:6,10:12)
!      AKK(7:9,4:12)=QKK1(4:6,1:9)
!	AKK(10:12,1:3)=QKK1(7:9,10:12)
!      AKK(10:12,4:12)=QKK1(7:9,1:9)

	

!      OPEN(50,FILE='mosti.dat',STATUS='OLD')
!      OPEN(51,FILE= 'AKK', STATUS= 'UNKNOWN')
!      READ(50,*)((AKK2(I,J),I=1,12),J=1,12)
!      AKK1(:,:)=AKK2(:,:)-AKK(:,:)
!      DO 42 I=1,12
!      WRITE(51,52)(AKK1(I,J),J=1,12)
!42    CONTINUE
!52    FORMAT(12E18.8)
C******************加入罚数,消除零能模式和沙漏问题****************
!===========================消除零能模式==========================
      AN(1,1)=-0.25
	AN(1,2)= 0.25
	AN(1,3)= 0.25
	AN(1,4)=-0.25
      AN(2,1)=-0.25
	AN(2,2)=-0.25
	AN(2,3)= 0.25
	AN(2,4)= 0.25

      DO 60 I=1,2
	DO 60 J=1,2
	AJ(I,J)=0.0
	DO 60 K=1,4
      AJ(I,J)=AJ(I,J)+AN(I,K)*VXY(K,J)
60    CONTINUE
      
	SJ0=AJ(1,1)*AJ(2,2)-AJ(1,2)*AJ(2,1)
      AJJ(1,1)= AJ(2,2)/SJ0
	AJJ(1,2)=-AJ(1,2)/SJ0
	AJJ(2,1)=-AJ(2,1)/SJ0
	AJJ(2,2)= AJ(1,1)/SJ0

	DO 61 I=1,2
	DO 61 J=1,4
	VDN(I,J)=0.0
	DO 61 K=1,2
      VDN(I,J)=VDN(I,J)+AJJ(I,K)*AN(K,J) !形函数对局部坐标的导数
61    CONTINUE

      VDN(:,:)=0.5*VDN(:,:)

      AII(1)=-VDN(2,1)
	AII(2)=VDN(1,1)
	AII(3)=-0.25

      AII(4)=-VDN(2,2)
	AII(5)=VDN(1,2)
	AII(6)=-0.25

	AII(7)=-VDN(2,3)
	AII(8)=VDN(1,3)
	AII(9)=-0.25

	AII(10)=-VDN(2,4)
	AII(11)=VDN(1,4)
	AII(12)=-0.25

      DO 62 I=1,12
	DO 62 J=1,12
      AA(I,J)=AII(I)*AII(J)
62    CONTINUE

      G=VMAE(1)/((1.0+VMAE(2))*2.0)   !剪切模量
	
	S=0.5*(VXY(2,1)*VXY(3,2)-VXY(3,1)*VXY(2,2)+
     $       VXY(3,1)*VXY(1,2)-VXY(1,1)*VXY(3,2)+
     $	   VXY(1,1)*VXY(2,2)-VXY(2,1)*VXY(1,2)+
     $	   VXY(4,1)*VXY(1,2)-VXY(1,1)*VXY(4,2)+
     $	   VXY(1,1)*VXY(3,2)-VXY(3,1)*VXY(1,2)+
     $	   VXY(3,1)*VXY(4,2)-VXY(4,1)*VXY(3,2))
      AA1(:,:)=(1E-6)*ABS(S)*G*AA(:,:)*VMAE(4)
      
!=================================================================
!============================消除沙漏=============================
      AII1(:)=0.0
	AII1(3)= 0.25
	AII1(6)=-0.25
	AII1(9)= 0.25
	AII1(12)=-0.25
	DO 63 I=1,12
	DO 63 J=1,12
      AA(I,J)=0.0
	AA(I,J)=AII1(I)*AII1(J) 
63    CONTINUE
      
	AAA(:,:)=(1E-3)*ABS(S)*G*VMAE(4)*AA(:,:)



!================================================================= 
C****************加入罚数,消除零能模式和沙漏问题结束*************
      QKK1=QKK1+AA1+AAA

      RETURN
	END