      SUBROUTINE  TRI_COUP_II(VX,FN,TQ)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION VX(3,3),TQ(9,3),TC(9,3),FN(3),A12(3),A23(3)
      
!	A12(1)=VX(2,1)-VX(1,1)
!	A12(2)=VX(2,2)-VX(1,2)
!	A12(3)=VX(2,3)-VX(1,3)

!	A23(1)=VX(3,1)-VX(2,1)
!	A23(2)=VX(3,2)-VX(2,2)
!	A23(3)=VX(3,3)-VX(2,3)

!	FN(1)=A12(2)*A23(3)-A12(3)*A23(2)
!	FN(2)=A12(3)*A23(1)-A12(1)*A23(3)
!	FN(3)=A12(1)*A23(2)-A12(2)*A23(1)

!	ZL=SQRT(FN(1)*FN(1)+FN(2)*FN(2)+FN(3)*FN(3))

 !     FN(1)=FN(1)/ZL
  !    FN(2)=FN(2)/ZL
!	FN(3)=FN(3)/ZL

	TQ(:,:)=0.0
      X=2.0/3.0
	Y=1.0/6.0
	Z=1.0/6.0
	CALL TRI_COUP_I(TC,VX,FN,X,Y,Z)
	TQ(:,:)=TQ(:,:)+TC(:,:)
	X=1.0/6.0
	Y=2.0/3.0
	Z=1.0/6.0
	CALL TRI_COUP_I(TC,VX,FN,X,Y,Z)
	TQ(:,:)=TQ(:,:)+TC(:,:)
	X=1.0/6.0
	Y=1.0/6.0
	Z=2.0/3.0
	CALL TRI_COUP_I(TC,VX,FN,X,Y,Z)
	TQ(:,:)=TQ(:,:)+TC(:,:)
	RETURN
	END