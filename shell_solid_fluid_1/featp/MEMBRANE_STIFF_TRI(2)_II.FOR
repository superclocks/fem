
      !该子程序待优化
      SUBROUTINE MEMBRANE_STIFF_TRI2_II(VX,TKK1,TMM1,VMAE)
      IMPLICIT REAL * 8(A-H,O-Z)
 
      DIMENSION   VX(3,3),VBK(9,9),VK(9,9),VMAE(4),Q(9),SR(9,9),
     $             TKK1(9,9),TMM1(9,9),VM(9,9),VKK1(9,12),
     $    VMM1(9,12), AT(12,9)

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
      CALL MEMBRANE_STIFF_TRI2_I(VX,X,Y,Z,VK,VM,VMAE)
      VBK(:,:)=VBK(:,:)+VK(:,:)*B*TH*ABS(SA)
	TMM1(:,:)=TMM1(:,:)+VM(:,:)*B*TH*VMAE(3)*ABS(SA)
      X=1.0/6.0
	Y=2.0/3.0
	Z=1.0/6.0
      CALL MEMBRANE_STIFF_TRI2_I(VX,X,Y,Z,VK,VM,VMAE)
      VBK(:,:)=VBK(:,:)+VK(:,:)*B*TH*ABS(SA)
      TMM1(:,:)=TMM1(:,:)+VM(:,:)*B*TH*VMAE(3)*ABS(SA)
      X=1.0/6.0
	Y=1.0/6.0
	Z=2.0/3.0
      CALL MEMBRANE_STIFF_TRI2_I(VX,X,Y,Z,VK,VM,VMAE)
      VBK(:,:)=VBK(:,:)+VK(:,:)*B*TH*ABS(SA)
	TMM1(:,:)=TMM1(:,:)+VM(:,:)*B*TH*VMAE(3)*ABS(SA)


      DO 20 I=1,3
      TMM1(3*I,:)=0.0
	TMM1(:,I*3)=0.0
20    CONTINUE

!     	AT(:,:)=0.0
!	AT(1,1)=1.0
!	AT(2,2)=1.0
!	AT(3,4)=1.0
!	AT(4,5)=1.0
!	AT(5,7)=1.0
!	AT(6,8)=1.0

!	AT(7,1)= 0.5
!	AT(7,3)= 0.125*(VX(1,2)-VX(2,2))
!	AT(7,4)= 0.5
!	AT(7,6)=-0.125*(VX(1,2)-VX(2,2))

!	AT(8,2)= 0.5
!	AT(8,3)= 0.125*(VX(2,1)-VX(1,1))
!	AT(8,5)= 0.5
!	AT(8,6)=-0.125*(VX(2,1)-VX(1,1))

!	AT(9,4)= 0.5
!	AT(9,6)= 0.125*(VX(2,2)-VX(3,2))
!	AT(9,7)= 0.5
!	AT(9,9)=-0.125*(VX(2,2)-VX(3,2))

!	AT(10,5)= 0.5
!	AT(10,6)= 0.125*(VX(3,1)-VX(2,1))
!	AT(10,8)= 0.5
!	AT(10,9)=-0.125*(VX(3,1)-VX(2,1))

!	AT(11,1)= 0.5
!	AT(11,3)=-0.125*(VX(3,2)-VX(1,2))
!	AT(11,7)= 0.5
!	AT(11,9)= 0.125*(VX(3,2)-VX(1,2))
!
!	AT(12,2)= 0.5
!	AT(12,3)=-0.125*(VX(1,1)-VX(3,1))
 !     AT(12,8)= 0.5
!	AT(12,9)= 0.125*(VX(1,1)-VX(3,1))
!	
!	DO 15 I=1,9
!      DO 15 J=1,12
!	VKK1(I,J)=0.0
!	VMM1(I,J)=0.0
!	DO 15 K=1,12
 !     VKK1(I,J)=VKK1(I,J)+AT(K,I)*VBK(K,J)
!	VMM1(I,J)=VMM1(I,J)+AT(K,I)*TMM1(K,J)
!15    CONTINUE

!      DO 16 I=1,9
!      DO 16 J=1,9
!	TKK1(I,J)=0.0
!	TMM1(I,J)=0.0
!	DO 16 K=1,12
!      TKK1(I,J)=TKK1(I,J)+VKK1(I,K)*AT(K,J)
!	TMM1(I,J)=TMM1(I,J)+VMM1(I,K)*AT(K,J)
!16    CONTINUE
       

	Q(1)=-0.25*(VX(2,1)-VX(3,1))/ABS(SA)
      Q(2)=-0.25*(VX(2,2)-VX(3,2))/ABS(SA)
	Q(3)= 1.0/3.0

      Q(4)=-0.25*(VX(3,1)-VX(1,1))/ABS(SA)
      Q(5)=-0.25*(VX(3,2)-VX(1,2))/ABS(SA)
	Q(6)= 1.0/3.0

	Q(7)=-0.25*(VX(1,1)-VX(2,1))/ABS(SA)
      Q(8)=-0.25*(VX(1,2)-VX(2,2))/ABS(SA)
	Q(9)= 1.0/3.0

      DO 10 I=1,9
	DO 10 J=1,9
	SR(I,J)=0.0
      SR(I,J)=SR(I,J)+Q(I)*Q(J)
10    CONTINUE
      G=VMAE(1)/((1.0+VMAE(2))*2.0)
      TKK1(:,:)=VBK(:,:)+G*ABS(SA)*(1E-2)*SR*TH
	RETURN
	END
