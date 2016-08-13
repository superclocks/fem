      SUBROUTINE DKT_II(VX,VK,VMAE)
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION VX(3,3),AK(9,9),VK(9,9),VMAE(4)
 
      VK(:,:)=0.0
	B3=1.0/6.0
	Y=1.0/2.0
	Z=0.0
      CALL DKT_I(VX,Y,Z,AK,VMAE,SA)
	VK(:,:)=VK(:,:)+B3*AK(:,:)*ABS(SA)
      Y=1.0/2.0
	Z=1.0/2.0
      CALL DKT_I(VX,Y,Z,AK,VMAE,SA)
	VK(:,:)=VK(:,:)+B3*AK(:,:)*ABS(SA)
      Y=0.0
	Z=1.0/2.0
      CALL DKT_I(VX,Y,Z,AK,VMAE,SA)
	VK(:,:)=VK(:,:)+B3*AK(:,:)*ABS(SA)

!      VK1(1:3,1:3)=VK(1:3,1:3)
!      VK1(1:3,4:6)=VK(1:3,7:9)
!	VK1(1:3,7:9)=VK(1:3,4:6)
!      VK1(4:6,1:3)=VK(7:9,1:3)
!      VK1(4:6,4:6)=VK(7:9,7:9)!
!	VK1(4:6,7:9)=VK(7:9,4:6)
!      VK1(7:9,1:3)=VK(4:6,1:3)
!      VK1(7:9,4:6)=VK(4:6,7:9)
!	VK1(7:9,7:9)=VK(4:6,4:6)
!       OPEN(100,  FILE ='VK1.DAT', STATUS= 'UNKNOWN' )
!	WRITE(100,101)((VK1(I,J),I=1,9),J=1,9)
101   FORMAT(9E18.8)
	RETURN
	END


