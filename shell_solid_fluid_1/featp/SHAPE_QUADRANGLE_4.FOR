!!	P726--727

C  .......................  SUB: 3-1-2  .............................
      SUBROUTINE SHAPE_QUADRANGLE_4(NODE,VXY,X, Y, IEL, VN, VDN)
C     COMPUTE SHAPE FUNCTION OF QUADRANGLE ELEMENT AT INTEGRATION POINT
      IMPLICIT REAL * 8 (A-H,O-z)
      DIMENSION  IEL(NODE) ,VN(9) ,VDN(3,9),VXY(4,3)
C ***
C *******CLEAR THE MEMORY FOR SHAPE FUNCTION AND ITS DERIVATIVES
C ***
      DO 10 I=1,9
         VN(I) =0.0                          ! VN-SHAPE FUNCTION
      DO 10 J=1,3
         VDN (J, I)= 0. 0                   !  VDN--LOCAL DERIVATIVE OF VN
10    CONTINUE
C ********* SET FUNCTION VALUE FOR QUADRANGLE ELEMENT OF 4 NODES
C ***
      VN(1)= (1+(VXY(1,1)/ABS(VXY(1,1)))*X)*
     $       (1+(VXY(1,2)/ABS(VXY(1,2)))*Y)*0.25
      VN(2)= (1+(VXY(2,1)/ABS(VXY(2,1)))*X)*
     $       (1+(VXY(2,2)/ABS(VXY(2,2)))*Y)*0.25
      VN(3)= (1+(VXY(3,1)/ABS(VXY(3,1)))*X)*
     $	   (1+(VXY(3,2)/ABS(VXY(3,2)))*Y)*0.25
      VN(4)= (1+(VXY(4,1)/ABS(VXY(4,1)))*X)*
     $	   (1+(VXY(4,2)/ABS(VXY(4,2)))*Y)*0.25
     
      VDN(1,1) =(VXY(1,1)/ABS(VXY(1,1)))*
     $	      (1+(VXY(1,2)/ABS(VXY(1,2)))*Y)*0.25
      VDN(1,2) =(VXY(2,1)/ABS(VXY(2,1)))*
     $	      (1+(VXY(2,2)/ABS(VXY(2,2)))*Y)*0.25
      VDN(1,3) =(VXY(3,1)/ABS(VXY(3,1)))*
     $          (1+(VXY(3,2)/ABS(VXY(3,2)))*Y)*0.25
      VDN(1,4) =(VXY(4,1)/ABS(VXY(4,1)))*
     $          (1+(VXY(4,2)/ABS(VXY(4,2)))*Y)*0.25
      VDN(2,1) =(VXY(1,2)/ABS(VXY(1,2)))*
     $          (1+(VXY(1,1)/ABS(VXY(1,1)))*X)*0.25
      VDN(2,2) =(VXY(2,2)/ABS(VXY(2,2)))*
     $          (1+(VXY(2,1)/ABS(VXY(2,1)))*X)*0.25
      VDN(2,3) =(VXY(3,2)/ABS(VXY(3,2)))*
     $          (1+(VXY(3,1)/ABS(VXY(3,1)))*X)*0.25
      VDN(2,4) =(VXY(4,2)/ABS(VXY(4,2)))*
     $          (1+(VXY(4,1)/ABS(VXY(4,1)))*X)*0.25
      
	RETURN
      END