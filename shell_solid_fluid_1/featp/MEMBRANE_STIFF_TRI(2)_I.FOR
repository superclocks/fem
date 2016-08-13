      
	SUBROUTINE MEMBRANE_STIFF_TRI2_I(VX,X,Y,Z,VK,VM,VMAE)
      IMPLICIT REAL * 8(A-H,O-Z)
      COMMON/COM1/SJ

      DIMENSION   VX(3,3),VXX(6,3),VK(9,9),VKK(12,12),SS(3,12),
     $            VN(6),VDN(2,6),VDD(2,6),AT(12,9),
     $            AJ(2,2),AJJ(2,2),VB(3,12),VD(3,3),VMAE(4),
     $            VKK1(9,12),VNN(2,12),VM(9,9),VMM(12,12),VMM1(9,12)              
	 
	VN(1)=(2.0*X-1.0)*X
	VN(2)=(2.0*Y-1.0)*Y
	VN(3)=(2.0*Z-1.0)*Z
	VN(4)=4.0*X*Y
	VN(5)=4.0*Y*Z
	VN(6)=4.0*Z*X

!	VDN(1,1)= 4.0*Y+4.0*Z-3.0
!	VDN(1,2)= 4.0*Y-1.0
!	VDN(1,3)= 0.0
!	VDN(1,4)=-8.0*Y-4.0*Z+4.0
!	VDN(1,5)= 4.0*Z
!	VDN(1,6)=-4.0*Z

!      VDN(2,1)= 4.0*Y+4.0*Z-3.0
!	VDN(2,2)=0.0
!	VDN(2,3)= 4.0*Z-1.0
!	VDN(2,4)=-4.0*Y
!	VDN(2,5)= 4.0*Y
!	VDN(2,6)=-8.0*Z-4.0*Y+4.0

	VDN(1,1)= 4.0*X-1.0
	VDN(1,2)= 0.0
	VDN(1,3)= 1.0-4.0*Z
	VDN(1,4)= 4.0*Y
	VDN(1,5)=-4.0*Y
	VDN(1,6)= 4.0*Z-4.0*X

      VDN(2,1)= 0.0
	VDN(2,2)= 4.0*Y-1.0
	VDN(2,3)= 1.0-4.0*Z
	VDN(2,4)= 4.0*X
	VDN(2,5)= 4.0*Z-4.0*Y
	VDN(2,6)=-4.0*X
	

      VXX(1:3,:)=VX(1:3,:)

	VXX(4,1)=0.5*(VX(1,1)+VX(2,1))
	VXX(5,1)=0.5*(VX(2,1)+VX(3,1))
	VXX(6,1)=0.5*(VX(3,1)+VX(1,1))

      VXX(4,2)=0.5*(VX(1,2)+VX(2,2))
	VXX(5,2)=0.5*(VX(2,2)+VX(3,2))
	VXX(6,2)=0.5*(VX(3,2)+VX(1,2))
	

	DO 10 I=1,2
	DO 10 J=1,2
	AJ(I,J)=0.0
	DO 10 K=1,6
	AJ(I,J)=AJ(I,J)+VDN(I,K)*VXX(K,J)
10    CONTINUE
      
	SJ=AJ(1,1)*AJ(2,2)-AJ(1,2)*AJ(2,1)

	AJJ(1,1)= AJ(2,2)/SJ
	AJJ(1,2)=-AJ(1,2)/SJ
	AJJ(2,1)=-AJ(2,1)/SJ
	AJJ(2,2)= AJ(1,1)/SJ

	DO 11 I=1,2
	DO 11 J=1,6
	VDD(I,J)=0.0
      DO 11 K=1,2
      VDD(I,J)=VDD(I,J)+AJJ(I,K)*VDN(K,J)
11    CONTINUE
      
	VB(:,:)=0.0
	DO 12 I=1,6
      VB(1,(I-1)*2+1)=VDD(1,I)
	VB(2,(I-1)*2+2)=VDD(2,I)
      VB(3,(I-1)*2+1)=VDD(2,I)
      VB(3,(I-1)*2+2)=VDD(1,I)
12    CONTINUE
       
	CALL MEMBRANE_VD(VD,VMAE)

	DO 13 I=1,3
	DO 13 J=1,12
      SS(I,J)=0.0
      DO 13 K=1,3
	SS(I,J)=SS(I,J)+VD(I,K)*VB(K,J)
13    CONTINUE

	DO 14 I=1,12
	DO 14 J=1,12
	VKK(I,J)=0.0
	DO 14 K=1,3
	VKK(I,J)=VKK(I,J)+VB(K,I)*SS(K,J)  !*SJ
14    CONTINUE

      DO 20 I=1,6
	VNN(1,(I-1)*2+1)=VN(I)
      VNN(2,(I-1)*2+2)=VN(I)
20    CONTINUE

      DO 21 I=1,12
	DO 21 J=1,12
	VMM(I,J)=0.0
	DO 21 K=1,2
      VMM(I,J)=VMM(I,J)+VNN(K,I)*VNN(K,J)  !*SJ
21    CONTINUE

!      VK(:,:)=VKK(:,:)   !*SJ
!	VM(:,:)=VMM(:,:)   !*SJ
      
	AT(:,:)=0.0
	AT(1,1)=1.0
	AT(2,2)=1.0
	AT(3,4)=1.0
	AT(4,5)=1.0
	AT(5,7)=1.0
	AT(6,8)=1.0

	AT(7,1)= 0.5
	AT(7,3)= 0.125*(VX(1,2)-VX(2,2))
	AT(7,4)= 0.5
	AT(7,6)=-0.125*(VX(1,2)-VX(2,2))

	AT(8,2)= 0.5
	AT(8,3)= 0.125*(VX(2,1)-VX(1,1))
	AT(8,5)= 0.5
	AT(8,6)=-0.125*(VX(2,1)-VX(1,1))

	AT(9,4)= 0.5
	AT(9,6)= 0.125*(VX(2,2)-VX(3,2))
	AT(9,7)= 0.5
	AT(9,9)=-0.125*(VX(2,2)-VX(3,2))

	AT(10,5)= 0.5
	AT(10,6)= 0.125*(VX(3,1)-VX(2,1))
	AT(10,8)= 0.5
	AT(10,9)=-0.125*(VX(3,1)-VX(2,1))

	AT(11,1)= 0.5
	AT(11,3)=-0.125*(VX(3,2)-VX(1,2))
	AT(11,7)= 0.5
	AT(11,9)= 0.125*(VX(3,2)-VX(1,2))

	AT(12,2)= 0.5
	AT(12,3)=-0.125*(VX(1,1)-VX(3,1))
      AT(12,8)= 0.5
	AT(12,9)= 0.125*(VX(1,1)-VX(3,1))
	
	DO 15 I=1,9
      DO 15 J=1,12
	VKK1(I,J)=0.0
	VMM1(I,J)=0.0
	DO 15 K=1,12
      VKK1(I,J)=VKK1(I,J)+AT(K,I)*VKK(K,J)
	VMM1(I,J)=VMM1(I,J)+AT(K,I)*VMM(K,J)
15    CONTINUE

      DO 16 I=1,9
      DO 16 J=1,9
	VK(I,J)=0.0
	VM(I,J)=0.0
	DO 16 K=1,12
      VK(I,J)=VK(I,J)+VKK1(I,K)*AT(K,J)
	VM(I,J)=VM(I,J)+VMM1(I,K)*AT(K,J)
16    CONTINUE

      RETURN
	END