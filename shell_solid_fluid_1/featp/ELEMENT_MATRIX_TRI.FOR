
      SUBROUTINE ELEMENT_MATRIX_TRI(VX,TCOS,TKK2,TMM2,VMAE)
      IMPLICIT   REAL*8 (A-H,O-z)
	COMMON/KEYOPT3/KEYOPT_3
	DIMENSION TCOS(3,3),TKK1(9,9),TMM1(9,9),TMM2(18,18),TKK2(18,18),
     $          VK(9,9),VM(9,9),VX(3,3),TR(18,18),VK2(18,18),VM2(18,18),
     $         VKK1(18,18),VMM1(18,18) 
      
C*****************************************************************************
C++++++++++++++++++++++Shell63单元的弯曲部分+++++++++++++++++++++++
C*******************由DKT单元的刚度矩阵及质量矩阵**************	
	CALL DKT_II(VX,VK,VMAE)            !VK(9,9)为DKT单元刚度矩阵
      CALL BEND_MASS_TRI_II(VX,VM,VMAE)  !VM(9,9)为DKT单元质量矩阵
C*******************由DKT单元的刚度矩阵及质量结束**************
C*****************************************************************************

C++++++++++++++++++++++Shell63单元的膜部分+++++++++++++++++++++++
!      KEYOPT_3=0
      SELECT CASE(KEYOPT_3)
C********************加入弹簧类型的转动刚度和额外形函数********************
      CASE(0)
      CALL MEMBRANE_STIFF_TRI1_II(VX,TKK1,TMM1,VMAE)
						!TKK1(:,:)为四边形膜单元的刚度矩阵
						!TMM1(:,:)为四边形膜单元的质量矩阵
C*****************************************************************************
C*******************加入弹簧类型转动刚度，不加额外形函数**********************
	CASE(1)
      CALL MEMBRANE_STIFF_TRI1_II(VX,TKK1,TMM1,VMAE)
C*****************************************************************************
C**************************Allman类型转动自由度*******************************
	CASE(2)
	CALL MEMBRANE_STIFF_TRI2_II(VX,TKK1,TMM1,VMAE)
      ENDSELECT
C***********************************膜部分结束********************************

C********************将膜单元和弯曲单元矩阵组装成壳单元矩阵*******************
      TKK2(:,:)=0.0
      TMM2(:,:)=0.0
      DO 42 I=1,3
	DO 42 J=1,3

      I2=(I-1)*6+1
	J2=(J-1)*6+1
      I3=(I-1)*3+1
	J3=(J-1)*3+1

	TKK2(I2:I2+1,J2:J2+1)=TKK1(I3:I3+1,J3:J3+1)
	TKK2(I2+5,J2:J2+1)=TKK1(I3+2,J3:J3+1)
      TKK2(I2:I2+1,J2+5)=TKK1(I3:I3+1,J3+2)
	TKK2(6*I,6*J)=TKK1(3*I,3*J)
	TKK2(I2+2:I2+4,J2+2:J2+4)=VK(I3:I3+2,J3:J3+2)   
	                                     !TKK2(18,18)为Shell63的刚度矩阵

	TMM2(I2:I2+1,J2:J2+1)=TMM1(I3:I3+1,J3:J3+1)
	TMM2(I2+5,J2:J2+1)=TMM1(I3+2,J3:J3+1)
      TMM2(I2:I2+1,J2+5)=TMM1(I3:I3+1,J3+2)
	TMM2(6*I,6*J)=TMM1(3*I,3*J)
	TMM2(I2+2:I2+4,J2+2:J2+4)=VM(I3:I3+2,J3:J3+2)   
	                                     !TMM2(18,18)为Shell63的质量矩阵
42    CONTINUE
C*******************膜单元和弯曲单元矩阵组装成壳单元矩阵结束********************
C**********************************************************************************
C****************************形成变换矩阵TR(18,18)****************************** 
      TR(:,:)=0.0
	DO 43 I=1,6
      II=(I-1)*3+1
      TR(II:II+2,II:II+2)=TCOS(:,:)
43    CONTINUE
C**********************************************************************************
C*******************将局部坐标下的单元矩阵变换到整体坐标系中********************
      DO 44 I=1,18
	DO 44 J=1,18
	VK2(I,J)=0.0
	VM2(I,J)=0.0
	DO 44 K=1,18
	VK2(I,J)=VK2(I,J)+TR(I,K)*TKK2(K,J)
	VM2(I,J)=VM2(I,J)+TR(I,K)*TMM2(K,J)
44    CONTINUE

      DO 45 I=1,18
	DO 45 J=1,18
      VKK1(I,J)=0.0
      VMM1(I,J)=0.0
	DO 45 K=1,18
	VKK1(I,J)=VKK1(I,J)+VK2(I,K)*TR(J,K)
	VMM1(I,J)=VMM1(I,J)+VM2(I,K)*TR(J,K)
45    CONTINUE
      TKK2(:,:)=VKK1(:,:)
	TMM2(:,:)=VMM1(:,:)
C****************局部坐标下的单元矩阵变换到整体坐标系中程序结束*****************

	RETURN
	END