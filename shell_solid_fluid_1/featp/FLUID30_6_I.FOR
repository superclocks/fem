
	SUBROUTINE FLUID30_6_I(XY,AL1,AL2,R,FK,FM,FSQQ)
	USE IMSL
      IMPLICIT  REAL* 8 (A-H, O-Z)
	COMMON  DENS,SONC
 !     COMMON/COM1/MND, NUMEL, NUMPT, MBAND, NMATI
 !     COMMON/COM2/NF, NFSTR, MSOLV, MPROB, MTYPE, NVA
 !     COMMON/COM3/MND2 , NUMPT2
 !     COMMON/ELEM/NODE, INTX, INTY
      DIMENSION   FXY(6,3),VN(6),VNN(3,24),VDN(3,6),VJJ(3,3),VJ(3,3),
     $          VD(3,6),FK(6,6),FM(6,6),FN(3),FN1(3),FN2(3),FSQQ(8,24),
     $          FSQ(8,3),XY(6,3)
C***************计算Fluid30单元的刚度、质量矩阵***********************
!===============Fluid30六节点单元的形函数==================
	VN(1)=AL1*(1.0-R)*0.5
	VN(2)=AL2*(1.0-R)*0.5
	VN(3)=(1.0-AL1-AL2)*(1.0-R)*0.5
!	VN(4)=(1.0-AL1-AL2)*(1.0-R)*0.5
	VN(4)=AL1*(1.0+R)*0.5
	VN(5)=AL2*(1.0+R)*0.5
	VN(6)=(1.0-AL1-AL2)*(1.0+R)*0.5
!	VN(8)=(1.0-AL1-AL2)*(1.0+R)*0.5
!====================================================
!=========Fluid30六节点单元的形函数对自然坐标的导数========
	VDN(1,1)= (1.0-R)*0.5
      VDN(1,2)= 0.0
	VDN(1,3)=-(1.0-R)*0.5
!	VDN(1,4)=-(1.0-R)*0.5
	VDN(1,4)= (1.0+R)*0.5
	VDN(1,5)= 0.0
	VDN(1,6)=-(1.0+R)*0.5
!	VDN(1,8)=-(1.0+R)*0.5

      VDN(2,1)=0.0
      VDN(2,2)= (1.0-R)*0.5
	VDN(2,3)=-(1.0-R)*0.5
!	VDN(2,4)=-(1.0-R)*0.5
	VDN(2,4)=0.0
	VDN(2,5)= (1.0+R)*0.5
	VDN(2,6)=-(1.0+R)*0.5
!	VDN(2,8)=-(1.0+R)*0.5

	VDN(3,1)=-AL1*0.5
      VDN(3,2)=-AL2*0.5
	VDN(3,3)=-(1.0-AL1-AL2)*0.5
!	VDN(3,4)=-(1.0-AL1-AL2)*0.5
	VDN(3,4)= AL1*0.5
	VDN(3,5)= AL2*0.5
	VDN(3,6)= (1.0-AL1-AL2)*0.5
!	VDN(3,8)= (1.0-AL1-AL2)*0.5
!====================================================
!=============计算雅克比矩阵及雅克比行列式===========
	DO 10 I=1,3
	DO 10 J=1,3
	VJJ(I,J)=0.0
	DO 10 K=1,6
      VJJ(I,J)=VJJ(I,J)+VDN(I,K)*XY(K,J) !雅克比矩阵
10    CONTINUE

      SJ8=VJJ(1,1)*(VJJ(2,2)*VJJ(3,3)-VJJ(2,3)*VJJ(3,2))-
     $    VJJ(1,2)*(VJJ(2,1)*VJJ(3,3)-VJJ(2,3)*VJJ(3,1))+
     $    VJJ(1,3)*(VJJ(2,2)*VJJ(3,1)-VJJ(2,1)*VJJ(3,2))  !雅克比行列式
!====================================================
      VJ(:,:)=0.0
	VJ(:,:)=.I.VJJ  !计算雅克比矩阵的逆阵

      DO 11 I=1,3
	DO 11 J=1,6
	VD(I,J)=0.0
      DO 11 K=1,3
      VD(I,J)=VD(I,J)+VJ(I,K)*VDN(K,J)   !形函数对全局坐标的导数
11    CONTINUE
      
      DO 12 I=1,6
	DO 12 J=1,6
      FK(I,J)=0.0
	FM(I,J)=0.0
	FM(I,J)=FM(I,J)+(VN(I)*VN(J)*ABS(SJ8)/(1500*1500))/6.0 !(SONC*SONC)!质量矩阵
	DO 12 K=1,3
	FK(I,J)=FK(I,J)+VD(K,I)*VD(K,J)*ABS(SJ8)/6.0     !刚度矩阵
12    CONTINUE
C********************************************************************* 
C************************计算耦合矩阵******************
!====================计算耦合面的外法线向量============ 
  !   	FN1(1)=FXY(6,1)-FXY(5,1)
!	FN1(2)=FXY(6,2)-FXY(5,2)
!	FN1(3)=FXY(6,3)-FXY(5,3)

!	FN2(1)=FXY(7,1)-FXY(6,1)
!	FN2(2)=FXY(7,2)-FXY(6,2)
!	FN2(3)=FXY(7,3)-FXY(6,3)

!	FN(1)=FN1(2)*FN2(3)-FN2(2)*FN1(3)
!	FN(2)=FN1(3)*FN2(1)-FN1(1)*FN2(3)
!	FN(3)=FN1(1)*FN2(2)-FN1(2)*FN2(1)
!=======================================================
!===================形成结构的形函数矩阵================
      fsqq(:,:)=0.0
!	IF((X>0.AND.Y>0.AND.Z>0).OR.(X>0.AND.Y>0.AND.Z<0).OR.
 !    $ (X>0.AND.Y<0.AND.Z>0).OR.(X>0.AND.Y<0.AND.Z<0))THEN
!	X=1.0
!	VN(1)=(1.0-Y)*(1.0-Z)*0.25
!	VN(2)=(1.0+Y)*(1.0-Z)*0.25
!	VN(3)=(1.0+Y)*(1.0+Z)*0.25
!	VN(4)=(1.0-Y)*(1.0+Z)*0.25
!	VN(5)=(1.0-Y)*(1.0-Z)*0.25
!	VN(6)=(1.0+Y)*(1.0-Z)*0.25
!	VN(7)=(1.0+Y)*(1.0+Z)*0.25
!	VN(8)=(1.0-Y)*(1.0+Z)*0.25

!	VDN(1,1)=-(1.0-Y)*(1.0-Z)*0.125
!      VDN(1,2)=-(1.0+Y)*(1.0-Z)*0.125
!	VDN(1,3)=-(1.0+Y)*(1.0+Z)*0.125
!	VDN(1,4)=-(1.0-Y)*(1.0+Z)*0.125
!	VDN(1,5)= (1.0-Y)*(1.0-Z)*0.125
!	VDN(1,6)= (1.0+Y)*(1.0-Z)*0.125
!	VDN(1,7)= (1.0+Y)*(1.0+Z)*0.125
!	VDN(1,8)= (1.0-Y)*(1.0+Z)*0.125

 !     VDN(2,1)=-(1.0-X)*(1.0-Z)*0.125
 !     VDN(2,2)= (1.0-X)*(1.0-Z)*0.125
!	VDN(2,3)= (1.0-X)*(1.0+Z)*0.125
!	VDN(2,4)=-(1.0-X)*(1.0+Z)*0.125
!	VDN(2,5)=-(1.0+X)*(1.0-Z)*0.125
!	VDN(2,6)= (1.0+X)*(1.0-Z)*0.125
!	VDN(2,7)= (1.0+X)*(1.0+Z)*0.125
!	VDN(2,8)=-(1.0+X)*(1.0+Z)*0.125

!	VDN(3,1)=-(1.0-X)*(1.0-Y)*0.125
!      VDN(3,2)=-(1.0-X)*(1.0+Y)*0.125
!	VDN(3,3)= (1.0-X)*(1.0+Y)*0.125
!	VDN(3,4)= (1.0-X)*(1.0-Y)*0.125
!	VDN(3,5)=-(1.0+X)*(1.0-Y)*0.125
!	VDN(3,6)=-(1.0+X)*(1.0+Y)*0.125
!	VDN(3,7)= (1.0+X)*(1.0+Y)*0.125
!	VDN(3,8)= (1.0+X)*(1.0-Y)*0.125

 !     VNN(:,:)=0.0
!	DO 20 I=1,8
!	VNN(1,(I-1)*3+1)=VN(I)
!	VNN(2,(I-1)*3+2)=VN(I)
!	VNN(3,(I-1)*3+3)=VN(I)
!20    CONTINUE
!=======================================================
!=====================计算耦合矩阵======================
!	DO 21 I=1,8
!	DO 21 J=1,3
!      FSQ(I,J)=0.0
!	FSQ(I,J)=FSQ(I,J)+VN(I)*FN(J)
!21    CONTINUE

!      DO 22 I=1,8
!	DO 22 J=1,24
!	FSQQ(I,J)=0.0
!	DO 22 K=1,3
!	FSQQ(I,J)=FSQQ(I,J)+FSQ(I,K)*VNN(K,J)*0.25
!22    CONTINUE
!      ENDIF
!=======================================================
C*******************************************************    
	RETURN
     	END
