
      SUBROUTINE FLUID30_6_II(FXY,FKK,FMM,FSQQQ)
      IMPLICIT  REAL* 8 (A-H, O-Z)
	COMMON  DENS,SONC
 !     COMMON/COM1/MND, NUMEL, NUMPT, MBAND, NMATI
 !     COMMON/COM2/NF, NFSTR, MSOLV, MPROB, MTYPE, NVA
 !     COMMON/COM3/MND2 , NUMPT2
 !     COMMON/ELEM/NODE, INTX, INTY
      DIMENSION     FXY(8,3),FKK(8,8),FK(8,8),FMM(6,6),FM(6,6),
     $         FSQQ(8,24),FSQQQ(8,24),XN(3),YN(3),ZN(3),TRA(3,3),
     $         PXY(6,3),XY(6,3),FN(3),FKK1(8,8),FKK2(8,8)

	FKK(:,:)=0.0
      FMM(:,:)=0.0
	FSQQQ(:,:)=0.0
C***************判断最大面积面,目的为了对单元重新编码************

!     程序待添加

C****************************************************************


C***************************坐标转换*****************************
	XN(1)=FXY(2,1)-FXY(1,1)
	XN(2)=FXY(2,2)-FXY(1,2)
	XN(3)=FXY(2,3)-FXY(1,3)   !X轴向量

	FN(1)=FXY(1,1)-FXY(3,1)   
	FN(2)=FXY(1,2)-FXY(3,2)  !面IJK(L)或123(4)的边IK(L)或13(4)的向量
	FN(3)=FXY(1,3)-FXY(3,3)

	ZN(1)=FN(3)*XN(2)-FN(2)*XN(3)
	ZN(2)=FN(1)*XN(3)-FN(3)*XN(1)
	ZN(3)=FN(2)*XN(1)-FN(1)*XN(2) !Z轴向量

	YN(1)=ZN(2)*XN(3)-ZN(3)*XN(2)
	YN(2)=ZN(3)*XN(1)-ZN(1)*XN(3)
	YN(3)=ZN(1)*XN(2)-ZN(2)*XN(1) !Y轴向量
	
	AMOD=SQRT(XN(1)*XN(1)+XN(2)*XN(2)+XN(3)*XN(3))  !X轴向量的模
      AMOD1=SQRT(YN(1)*YN(1)+YN(2)*YN(2)+YN(3)*YN(3)) !Y轴向量的模
	AMOD2=SQRT(ZN(1)*ZN(1)+ZN(2)*ZN(2)+ZN(3)*ZN(3)) !Z轴向量的模

      XN(1)=XN(1)/AMOD 
	XN(2)=XN(2)/AMOD
	XN(3)=XN(3)/AMOD  !X轴单位向量

      YN(1)=YN(1)/AMOD1 
	YN(2)=YN(2)/AMOD1
	YN(3)=YN(3)/AMOD1  !Y轴单位向量

	ZN(1)=ZN(1)/AMOD2 
	ZN(2)=ZN(2)/AMOD2
	ZN(3)=ZN(3)/AMOD2  !Z轴单位向量

	TRA(:,1)=XN(:)
	TRA(:,2)=YN(:)
	TRA(:,3)=ZN(:)     !TRA为两个坐标的转换矩阵

	PXY(1,:)=FXY(1,:)-FXY(1,:)
	PXY(2,:)=FXY(2,:)-FXY(1,:)
	PXY(3,:)=FXY(3,:)-FXY(1,:)
!	PXY(4,:)=FXY(4,:)-FXY(1,:)
	PXY(4,:)=FXY(5,:)-FXY(1,:)
	PXY(5,:)=FXY(6,:)-FXY(1,:)
	PXY(6,:)=FXY(7,:)-FXY(1,:)
!	PXY(8,:)=FXY(8,:)-FXY(1,:)   !PXY为相对坐标矩阵

!	PXY=.T.PXY

      DO 30 I=1,6
	DO 30 J=1,3
	XY(I,J)=0.0
	DO 30 K=1,3
      XY(I,J)=XY(I,J)+PXY(I,K)*TRA(K,J)
30    CONTINUE
C*************************坐标转换结束***************************
C****************************************************************
C***********************计算六点高斯积分*************************
	DO 20 I=1,6
	SELECT CASE(I)
	CASE(1)
	AL1= 0.66666 66666 66666
	AL2= 0.16666 66666 66666
	R=   0.577350269189626
      CASE(2)
      AL1= 0.16666 66666 66666
	AL2= 0.66666 66666 66666
	R=   0.577350269189626
	CASE(3)
	AL1= 0.16666 66666 66666
	AL2= 0.16666 66666 66666
	R=   0.577350269189626
      CASE(4)
	AL1= 0.66666 66666 66666
	AL2= 0.16666 66666 66666
	R=  -0.577350269189626

	CASE(5)
      AL1= 0.16666 66666 66666
	AL2= 0.66666 66666 66666
	R=  -0.577350269189626
	CASE(6)
      AL1= 0.16666 66666 66666
	AL2= 0.16666 66666 66666
	R=  -0.577350269189626
	END SELECT
	 
	CALL FLUID30_6_I(XY,AL1,AL2,R,FK,FM,FSQQ)
      FKK(:,:)=FKK(:,:)+FK(:,:)
      FMM(:,:)=FMM(:,:)+FM(:,:)
	FSQQQ(:,:)=FSQQQ(:,:)+FSQQ(:,:)
C*****************计算高斯积分结束*******************
20    CONTINUE

!	DO 21 II=1,2
!	SELECT CASE(II)
!	CASE(1)
!      AL1= 1.0/3.0
!	AL2= 1.0/3.0
!	R=   0.577350269189626
!	CASE(2)
!	AL1= 1.0/3.0
!	AL2= 1.0/3.0
!	R=  -0.577350269189626
!	END SELECT
!	FKK(:,:)=FKK(:,:)+FK(:,:)*(-27.0/48.0)
 !     FMM(:,:)=FMM(:,:)+FM(:,:)
!	FSQQQ(:,:)=FSQQQ(:,:)+FSQQ(:,:)

!21    CONTINUE
      
!	DO 30 I=1,8
!	DO 30 J=1,8
	



      RETURN
	END