
      SUBROUTINE  ELEMENT_MATRIX_QUAD(VXY,QCOS,IEL,VMM1,VKK1,VMAE)
      IMPLICIT REAL*8 (A-H, O-Z)
	COMMON/KEYOPT3/KEYOPT_3
	COMMON/HIGH/H
      DIMENSION  VXY(4,3),VXYY(6,3),IEL(4),IELL(6),VMAE(4),VMM1(24,24),
     $          VKK1(24,24),VK1(24,24),VK2(24,24),VM1(24,24),
     $      VM2(24,24),DKK1(12, 12),QCOS(3,3),TR(24,24),VK(9,9),VX(3,3),
     $          VM(9,9),DMM1(12,12),QKK1(12,12),QMM1(12,12),STA(6,6)
C*********************************************************************************
C++++++++++++++++++++++Shell63单元的弯曲部分+++++++++++++++++++++++++++++++++++++++++++
C************由DKT单元组成的DKQ单元刚度矩阵及质量矩阵*****************************
      VXYY(1:4,:)=VXY(:,:)
	VXYY(5:6,:)=VXY(1:2,:)
	IELL(:)=[1,2,3,4,1,2]
	DKK1(:,:)=0.0
	DMM1(:,:)=0.0
      
      DO 40 II=1,4
      VX(:,:)=VXYY(II:II+2,:)
      CALL DKT_II(VX,VK,VMAE)
      CALL BEND_MASS_TRI_II(VX,VM,VMAE)
	DO 41 I1=1,3
	DO 41 J1=1,3
      I2=(IELL(II+I1-1)-1)*3+1
	J2=(IELL(II+J1-1)-1)*3+1
      DKK1(I2:I2+2,J2:J2+2)=DKK1(I2:I2+2,J2:J2+2)
     $            +VK((I1-1)*3+1:(I1-1)*3+3,(J1-1)*3+1:(J1-1)*3+3)
	DMM1(I2:I2+2,J2:J2+2)=DMM1(I2:I2+2,J2:J2+2)
     $            +VM((I1-1)*3+1:(I1-1)*3+3,(J1-1)*3+1:(J1-1)*3+3)
41    CONTINUE
40    CONTINUE
      DKK1(:,:)=DKK1(:,:)*0.5 !DKK1(:,:)为DKQ单元的单元刚度矩阵
	DMM1(:,:)=DMM1(:,:)*0.5 !DMM1(:,:)为DKQ单元的单元质量矩阵
C***********由DKT单元组成的DKQ单元刚度矩阵及质量结束*******
C**********************************************************************************

C++++++++++++++++++++++Shell63单元的膜部分+++++++++++++++++++++++++++++++++++++++++
!      KEYOPT_3=1    !控制参数
      SELECT CASE(KEYOPT_3)
C***************************带有额外形函数的四边形单元*****************************
      CASE(0)
	CALL MEMBRANE_STIFF_QUAD0_II(VXY,QKK1,QMM1,VMAE,STA,QCOS)
C**********************************************************************************
C********************带弹簧类型的转动自由度,不带额外形函数膜部分*******************
      CASE(1)
      CALL MEMBRANE_STIFF_QUAD1_II(VXY,QKK1,QMM1,VMAE,STA,QCOS)
						!QKK1(:,:)为四边形膜单元的刚度矩阵
						!QMM1(:,:)为四边形膜单元的质量矩阵
C***********************************膜部分结束*************************************
C***************************带Allman转角自由度类型的膜单元*************************
      CASE(2)
	CALL MEMBRANE_STIFF_QUAD2_II(VXY,QKK1,QMM1,VMAE)
C**********************************************************************************
      ENDSELECT
C********************将膜单元和弯曲单元矩阵组装成壳单元矩阵************************
      VK1(:,:)=0.0
	VM1(:,:)=0.0
      DO 42 I=1,4
	DO 42 J=1,4

      I2=(I-1)*6+1
	J2=(J-1)*6+1
      I3=(I-1)*3+1
	J3=(J-1)*3+1

	VK1(I2:I2+1,J2:J2+1)=QKK1(I3:I3+1,J3:J3+1)
	VK1(I2+5,J2:J2+1)=QKK1(I3+2,J3:J3+1)
      VK1(I2:I2+1,J2+5)=QKK1(I3:I3+1,J3+2)
	VK1(6*I,6*J)=QKK1(3*I,3*J)
	VK1(I2+2:I2+4,J2+2:J2+4)=DKK1(I3:I3+2,J3:J3+2)

	VM1(I2:I2+1,J2:J2+1)=QMM1(I3:I3+1,J3:J3+1)
	VM1(I2+5,J2:J2+1)=QMM1(I3+2,J3:J3+1)
      VM1(I2:I2+1,J2+5)=QMM1(I3:I3+1,J3+2)
	VM1(6*I,6*J)=QMM1(3*I,3*J)
	VM1(I2+2:I2+4,J2+2:J2+4)=DMM1(I3:I3+2,J3:J3+2)
42    CONTINUE
C*****************************当翘曲时加入调节矩阵********************************   
!	STA(:,:)=0.0
!	DO 60 I=1,6
!      STA(I,I)=1.0
!60    CONTINUE

!	DO 61 I=1,4
!	SELECT CASE(I)
!	CASE(1)
!      STA(1,5)=H
!	STA(2,4)=H
!	VK1(1:6,1:6)=VK1(1:6,1:6)+STA(:,:)
!	CASE(2)
!      STA(1,5)=-H
!	STA(2,4)=-H
!	VK1(7:12,7:12)=VK1(1:6,1:6)+STA(:,:)
!	CASE(3)
!      STA(1,5)=H
!	STA(2,4)=H
!	VK1(13:18,13:18)=VK1(1:6,1:6)+STA(:,:)
!	CASE(4)
 !     STA(1,5)=-H
!	STA(2,4)=-H
!	VK1(19:24,19:24)=VK1(1:6,1:6)+STA(:,:)
!	ENDSELECT
!61    CONTINUE
C*********************************************************************************
!     VMM1(:,:)=VM1(:,:)
!	VK1(:,:)=VK1(:,:)+STA(:,:)
C*******************膜单元和弯曲单元矩阵组装成壳单元矩阵结束********************
C**********************************************************************************
C****************************形成变换矩阵TR(24,24)****************************** 
      TR(:,:)=0.0
	DO 43 I=1,8
      II=(I-1)*3+1
      TR(II:II+2,II:II+2)=QCOS(:,:)
43    CONTINUE
C**********************************************************************************
C*******************将局部坐标下的单元矩阵变换到整体坐标系中***********************
      DO 44 I=1,24
	DO 44 J=1,24
	VK2(I,J)=0.0
	VM2(I,J)=0.0
	DO 44 K=1,24
	VK2(I,J)=VK2(I,J)+TR(I,K)*VK1(K,J)
	VM2(I,J)=VM2(I,J)+TR(I,K)*VM1(K,J)
44    CONTINUE

      DO 45 I=1,24
	DO 45 J=1,24
      VKK1(I,J)=0.0
      VMM1(I,J)=0.0
	DO 45 K=1,24
	VKK1(I,J)=VKK1(I,J)+VK2(I,K)*TR(J,K)
	VMM1(I,J)=VMM1(I,J)+VM2(I,K)*TR(J,K)
45    CONTINUE
C**********************************************************************************
     
      RETURN
      END