!!	P733--734

C  ........................  SUB: 3-1-7  .............................
            SUBROUTINE ELEMENT_JACOBI(MND,VXY,VDN,SJ,VD0,VDN1)
C***********************************************************************
C   GET THE DETERMINANT OF JACOBI MATRIX AND CARTESIAN DERIVATIVES
C***********************************************************************
            IMPLICIT  REAL *8(A-H,O-Z)
            DIMENSION  VXY(MND,3),VXY1(8,3),VDN(3,9),VD0(3,9),VDN1(2,4)
            DIMENSION  VJJ(2,2),VJ1(2,2)
C ***
C ******** FORM JACOBI MATRIX
C ***			
            VXY1(:,:)=0.0
	      VXY1(1:4,:)=VXY(1:4,:)
	      VXY1(5,1)=0.5*(VXY(1,1)+VXY(2,1))
		  VXY1(6,1)=0.5*(VXY(2,1)+VXY(3,1))
		  VXY1(7,1)=0.5*(VXY(3,1)+VXY(4,1))
		  VXY1(8,1)=0.5*(VXY(4,1)+VXY(1,1))
		  
		  VXY1(5,2)=0.5*(VXY(1,2)+VXY(2,2))
		  VXY1(6,2)=0.5*(VXY(2,2)+VXY(3,2))
		  VXY1(7,2)=0.5*(VXY(3,2)+VXY(4,2))
		  VXY1(8,2)=0.5*(VXY(4,2)+VXY(1,2))

		  DO 11 II=1,2
            DO 11 JJ=1,2
            VJJ(II,JJ)=0.0                         !  VJJ-JACOBI MATRIX [J]    
		  DO 11 KK=1,MND+4
            VJJ(II,JJ) =VJJ(II,JJ)+VDN(II,KK)*VXY1(KK,JJ)
11          CONTINUE
C ***
C *********FORM THE DETERMINANT OF JACOBI MATRIX
            SJ=VJJ(1,1)*VJJ(2,2)-VJJ(1,2)*VJJ(2,1)     ! SJ-|J|
C ***
C ******** COMPUTE THE INVERSE OF JACOBI MATRIX
C ***
            VJ1(1,1)=+VJJ(2,2)/SJ
            VJ1(1,2)=-VJJ(1,2)/SJ       ! VJ1--INVERSE OF [J]
            VJ1(2,1)=-VJJ(2,1)/SJ
            VJ1(2,2)=+VJJ(1,1)/SJ
C ***
C ************ COMPUTE THE CARTESIAN DERIVATIVES OF SHAPE FUNCTION
            DO 51 II=1,2
            DO 51 JJ=1,9
            VD0(II,JJ) =0.0             !VD0-GLOBAL DERIVATIVE OF VN
            DO 51 KK=1,2
            VD0 (II ,JJ)=VD0(II,JJ)+VJ1(II,KK)*VDN(KK ,JJ)
51          CONTINUE
            RETURN
            END