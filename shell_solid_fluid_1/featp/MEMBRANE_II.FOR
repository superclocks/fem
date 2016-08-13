

      SUBROUTINE MEMBRANE_II(VX,TKK1,TMM1,VMAE)
	IMPLICIT REAL * 8(A-H,O-Z)
 
      DIMENSION   VX(3,3),VBK(5,5),VK(9,5),VMAE(4),Q(9),SR(9,9),
     $             TKK1(9,9),TMM1(12,12),VM(12,12),VKK1(9,12),
     $    VMM1(9,12), AT(12,9),CC(5,9),AK(5,5)

      SA=0.5*(VX(2,1)*VX(3,2)-VX(3,1)*VX(2,2)-VX(1,1)*VX(3,2)+
     $   VX(3,1)*VX(1,2)+VX(1,1)*VX(2,2)-VX(2,1)*VX(1,2))
	 
	TH=VMAE(4)
	TKK1(:,:)=0.0
	TMM1(:,:)=0.0
	VBK(:,:)=0.0 
	B=1.0/3.0
	X=2.0/3.0
	Y=1.0/6.0
	Z=1.0/6.0
      CALL MEMBRANE_I(VX,X,Y,Z,AK,VMAE,CC)
      VBK(:,:)=VBK(:,:)+AK(:,:)*B
!	TMM1(:,:)=TMM1(:,:)+VM(:,:)*B*TH*VMAE(3)
      X=1.0/6.0
	Y=2.0/3.0
	Z=1.0/6.0
      CALL MEMBRANE_I(VX,X,Y,Z,AK,VMAE,CC)
      VBK(:,:)=VBK(:,:)+AK(:,:)*B
!      TMM1(:,:)=TMM1(:,:)+VM(:,:)*B*TH*VMAE(3)
      X=1.0/6.0
      Y=1.0/6.0
	Z=2.0/3.0
      CALL MEMBRANE_I(VX,X,Y,Z,AK,VMAE,CC)
      VBK(:,:)=VBK(:,:)+AK(:,:)*B

	DO 10 I=1,9
	DO 10 J=1,5
	VK(I,J)=0.0
      DO 10 K=1,5
      VK(I,J)=VK(I,J)+CC(K,I)*VBK(K,J)
10    CONTINUE

      DO 11 I=1,9
	DO 11 J=1,9
	DO 11 K=1,5
      TKK1(I,J)=TKK1(I,J)+VK(I,K)*CC(K,J)
11    CONTINUE

	Q(1)=-0.25*(VX(2,1)-VX(3,1))/ABS(SA)
      Q(2)=-0.25*(VX(2,2)-VX(3,2))/ABS(SA)
	Q(3)= 1.0/3.0

      Q(4)=-0.25*(VX(3,1)-VX(1,1))/ABS(SA)
      Q(5)=-0.25*(VX(3,2)-VX(1,2))/ABS(SA)
	Q(6)= 1.0/3.0

	Q(7)=-0.25*(VX(1,1)-VX(2,1))/ABS(SA)
      Q(8)=-0.25*(VX(1,2)-VX(2,2))/ABS(SA)
	Q(9)= 1.0/3.0

      DO 12 I=1,9
	DO 12 J=1,9
	SR(I,J)=0.0
      SR(I,J)=SR(I,J)+Q(I)*Q(J)
12    CONTINUE
      G=VMAE(1)/((1.0+VMAE(2))*2.0)

      TKK1(:,:)=TKK1(:,:)+G*ABS(SA)*(1E-2)*SR*TH










      RETURN
	END