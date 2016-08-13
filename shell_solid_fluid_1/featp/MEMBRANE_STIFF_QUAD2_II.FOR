
	SUBROUTINE MEMBRANE_STIFF_QUAD(2)_II(VXY,AKK,VMAE,IEL)
	USE IMSL
      IMPLICIT   REAL * 8 (A-H,O-z)
     	COMMON/COM1/MND, NUMEL, NUMPT, MBAND, NMATI
      COMMON/COM2/NF, NFSTR, MSOLV, MPROB, MTYPE, NVA
      COMMON/COM3/MND2 , NUMPT2
      COMMON/ELEM/NODE, INTX, INTY

      DIMENSION  IEL(NODE),VXY(4,3),VK(16,16),
     $        VSS(12,16),AKK(12,12),QM(16,16),akk1(12,12),
     $        VD(3,3),VMAE(4),QK(16,16),VT(16,12),VM(16,16),VP(12,12),
     $        VPP(12,12),AA(12,12),AAI(8,8),AA2(4,4),AII(12),AII2(4),
     $        AA3(12,12),AKK2(12,12),VX(8,3),VK0(12,16),VK00(12,12)

      QK(:,:)=0.0
	QM(:,:)=0.0
	VPP(:,:)=0.0
!      INTX=3
!      INTy=3
	DO 40 I=1,INTX
	DO 40 J=1,INTY
	CALL GAUSS_INTEGRATION (INTX,INTY,I,J,X,Y,WXY)
      CALL MEMBRANE_STIFF_QUAD2_I(VXY,VK,VP,VM,X,Y,VMAE,WXY,IEL)
	QK(:,:)=QK(:,:)+VK(:,:)
      VPP(:,:)=VPP(:,:)+VP(:,:)
      QM(:,:)=QM(:,:)+VM(:,:)
40    CONTINUE

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
	VK00(I,J)=0.0
	DO 46 K=1,16
      VK00(I,J)=VK00(I,J)+VK0(I,K)*VT(K,J)
46    CONTINUE
      AKK(1:3,1:3)=VK00(10:12,10:12)
      AKK(1:3,4:12)=VK00(10:12,1:9)
      AKK(4:6,1:3)=VK00(1:3,10:12)
	AKK(4:6,4:12)=VK00(1:3,1:9)
	AKK(7:9,1:3)=VK00(4:6,10:12)
      AKK(7:9,4:12)=VK00(4:6,1:9)
	AKK(10:12,1:3)=VK00(7:9,10:12)
      AKK(10:12,4:12)=VK00(7:9,1:9)

	

      OPEN(50,FILE='mosti.dat',STATUS='OLD')
      OPEN(51,FILE= 'AKK', STATUS= 'UNKNOWN')
      READ(50,*)((AKK2(I,J),I=1,12),J=1,12)
      AKK1(:,:)=AKK2(:,:)-AKK(:,:)
      DO 42 I=1,12
      WRITE(51,52)(AKK1(I,J),J=1,12)
42    CONTINUE
52    FORMAT(12E18.8)
C******************加入罚数,消除零能模式和沙漏问题****************
      AII(1)=31.25
	AII(2)=-31.25
	AII(3)=-0.25

      AII(4)=31.25
	AII(5)=31.25
	AII(6)=-0.25

	AII(7)=-31.25
	AII(8)=31.25
	AII(9)=-0.25

	AII(10)=-31.25
	AII(11)=-31.25
	AII(12)=-0.25

	
      DO 41 I=1,12
	DO 41 J=1,12
      AA(I,J)=AII(I)*AII(J)
41    CONTINUE
      G=VMAE(1)/((1.0+VMAE(2))*2.0)
      AA3(:,:)=(6.4E-10)*G*AA(:,:)  !+VK00(:,:)
      
	AKK(1:3,1:3)=AA3(10:12,10:12)
      AKK(1:3,4:12)=AA3(10:12,1:9)
      AKK(4:6,1:3)=AA3(1:3,10:12)
	AKK(4:6,4:12)=AA3(1:3,1:9)
	AKK(7:9,1:3)=AA3(4:6,10:12)
      AKK(7:9,4:12)=AA3(4:6,1:9)
	AKK(10:12,1:3)=AA3(7:9,10:12)
      AKK(10:12,4:12)=AA3(7:9,1:9)

      
      
!	AII2(:)=(VXY(1,2)-VXY(2,2)+VXY(3,2)-VXY(4,2))*AAI(4,:)+
!     $        (VXY(2,1)-VXY(1,1)+VXY(4,1)-VXY(3,1))*AAI(4,:)



 
C****************加入罚数,消除零能模式和沙漏问题结束*************
      

      RETURN
	END