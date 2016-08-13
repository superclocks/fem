      
	SUBROUTINE QUAD_COUP_I(TC,VXY,SN,X,Y)
      IMPLICIT  REAL*8 (A-H,O-Z)
      DIMENSION  VXY(4,3),TC(12,4),VN(4),VDN(2,4),VJJ(2,2),VNN(3,12),
     $           SN(3),FS(12)

      VN(1)= (1-X)*(1-Y)*0.25
      VN(2)= (1+X)*(1-Y)*0.25
      VN(3)= (1+X)*(1+Y)*0.25
      VN(4)= (1-X)*(1+Y)*0.25    !�Ľڵ㵥Ԫ�κ���
	
      VDN(1,1) =-(1-Y)*0.25
      VDN(1,2) = (1-Y)*0.25
      VDN(1,3) = (1+Y)*0.25
      VDN(1,4) =-(1+Y)*0.25


      VDN(2,1) =-(1-X)*0.25
      VDN(2,2) =-(1+X)*0.25
      VDN(2,3) = (1+X)*0.25
      VDN(2,4) = (1-X)*0.25      !�Ľڵ㵥Ԫ�κ�����ȫ������ĵ���

	DO 10 I=1,2
	DO 10 J=1,2
      VJJ(I,J)=0.0
	DO 10 K=1,4
      VJJ(I,J)=VJJ(I,J)+VDN(I,K)*VXY(K,J)   !�ſ˱Ⱦ���
10    CONTINUE

      SJ=VJJ(1,1)*VJJ(2,2)-VJJ(1,2)*VJJ(2,1)  !�ſ˱�����ʽ

      VNN(:,:)=0.0
	DO 20 I=1,4
	VNN(1,(I-1)*3+1)=VN(I)   !�κ�����ȫ������ĵ���
	VNN(2,(I-1)*3+2)=VN(I)
	VNN(3,(I-1)*3+3)=VN(I)
20    CONTINUE
!=======================================================
!=====================������Ͼ���======================
	DO 21 I=1,12
      FS(I)=0.0
	DO 21 J=1,3
	FS(I)=FS(I)+VNN(J,I)*SN(J)
21    CONTINUE

      DO 22 I=1,12
	DO 22 J=1,4
	TC(I,J)=0.0
	TC(I,J)=TC(I,J)+FS(I)*VN(J)*SJ
22    CONTINUE
    
	RETURN
	END