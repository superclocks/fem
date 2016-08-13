
	SUBROUTINE TRI_COUP_I(TC,VX,FN,X,Y,Z)
      IMPLICIT  REAL*8 (A-H,O-Z)
      DIMENSION  VX(3,3),TC(9,3),VN(3),VDN(2,3),VJJ(2,2),VNN(3,9),
     $           FN(3),FS(9)

	VN(1)=X  !�κ���
	VN(2)=Y
	VN(3)=Z
	
	VDN(1,1)= 1.0   !�κ����Ծֲ�����ĵ���
	VDN(1,2)= 0.0
	VDN(1,3)=-1.0
	VDN(2,1)= 0.0
	VDN(2,2)= 1.0
	VDN(2,3)=-1.0

	DO 10 I=1,2
	DO 10 J=1,2
      VJJ(I,J)=0.0
	DO 10 K=1,3
      VJJ(I,J)=VJJ(I,J)+VDN(I,K)*VX(K,J)   !�ſ˱Ⱦ���
10    CONTINUE

      SJ=VJJ(1,1)*VJJ(2,2)-VJJ(1,2)*VJJ(2,1)  !�ſ˱�����ʽ

      VNN(:,:)=0.0
	DO 20 I=1,3
	VNN(1,(I-1)*3+1)=VN(I)   !�κ�����ȫ������ĵ���
	VNN(2,(I-1)*3+2)=VN(I)
	VNN(3,(I-1)*3+3)=VN(I)
20    CONTINUE
!=======================================================
!=====================������Ͼ���======================
	DO 21 I=1,9
      FS(I)=0.0
	DO 21 J=1,3
	FS(I)=FS(I)+VNN(J,I)*FN(J)
21    CONTINUE

      DO 22 I=1,9
	DO 22 J=1,3
	TC(I,J)=0.0
	TC(I,J)=TC(I,J)+FS(I)*VN(J)*SJ/6.0
22    CONTINUE
    
	RETURN
	END