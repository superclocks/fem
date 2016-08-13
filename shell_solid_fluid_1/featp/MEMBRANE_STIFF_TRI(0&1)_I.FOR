
	SUBROUTINE MEMBRANE_STIFF_TRI1_I(VX,VK,VM,X,Y,Z,VMAE)
	IMPLICIT   REAL*8 (A-H,O-z)
      DIMENSION  IEL(4) ,VN(3),VNN(3,9),VDN(2,3),
     $          VX(3,3),VD(3,3),
     $         VJJ(2,2),VJ(2,2),VDD(2,3),
     $         VB(3,9),VSS(3,9),VK(9,9),
     $         VMAE(4),VM(9,9)
	VN(1)=X
	VN(2)=Y
	VN(3)=Z
	VDN(1,1)= 1.0
	VDN(1,2)= 0.0
	VDN(1,3)=-1.0
	VDN(2,1)= 0.0
	VDN(2,2)= 1.0
	VDN(2,3)=-1.0

	DO 10 I=1,2
	DO 10 J=1,2
      VJJ(I,J)=0.0
	DO 10 K=1,3
      VJJ(I,J)=VJJ(I,J)+VDN(I,K)*VX(K,J)
10    CONTINUE

      SJ=VJJ(1,1)*VJJ(2,2)-VJJ(1,2)*VJJ(2,1)  !雅克比行列式

	VJ(1,1)= VJJ(2,2)/SJ
	VJ(1,2)=-VJJ(1,2)/SJ
	VJ(2,1)=-VJJ(2,1)/SJ
	VJ(2,2)= VJJ(1,1)/SJ      !雅克比矩阵的逆阵

      DO 11 I=1,2
	DO 11 J=1,3
	VDD(I,J)=0.0
	DO 11 K=1,2
	VDD(I,J)=VDD(I,J)+VJ(I,K)*VDN(K,J)	!4节点形函数对全局坐标点的导数
11    CONTINUE

      VB(:,:)=0.0
	DO 12 I=1,3
	VB(1,(I-1)*3+1)=VDD(1,I)
      VB(2,(I-1)*3+2)=VDD(2,I)
      VB(3,(I-1)*3+1)=VDD(2,I)
	VB(3,(I-1)*3+2)=VDD(1,I)    !应变矩阵
12    CONTINUE

      CALL MEMBRANE_VD(VD,VMAE)  !调用弹性矩阵子函数
      DO 13 I=1,3
	DO 13 J=1,9   
	VSS(I,J)=0.0
	DO 13 K=1,3
      VSS(I,J)=VSS(I,J)+VD(I,K)*VB(K,J)
13    CONTINUE

      DO 14 I=1,9
	DO 14 J=1,9
	VK(I,J)=0.0
	DO 14 K=1,3
      VK(I,J)=VK(I,J)+VB(K,I)*VSS(K,J)*SJ*VMAE(4) !计算刚度矩阵
14    CONTINUE
	
	VNN(:,:)=0.0
	DO 15 I=1,3
	VNN(1,(I-1)*3+1)=VN(I)
      VNN(2,(I-1)*3+2)=VN(I)
!	VNN(3,(I-1)*3+3)=VN(I)
15    CONTINUE

      DO 16 I=1,9
	DO 16 J=1,9
	VM(I,J)=0.0
	DO 16 K=1,3
	VM(I,J)=VM(I,J)+VNN(K,I)*VNN(K,J)*SJ*VMAE(3)*VMAE(4) !计算质量矩阵
16    CONTINUE      

      RETURN
	END