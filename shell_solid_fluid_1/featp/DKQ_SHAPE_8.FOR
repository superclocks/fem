
            SUBROUTINE     DKQ_SHAPE_8  (X, Y,VN, VDN)
C     COMPUTE SHAPE FUNCTION OF QUADRANGLE ELEMENT AT INTEGRATION POINT
            IMPLICIT   REAL * 8 (A-H,O-z)
            DIMENSION  VN(9) ,VDN(3,9)
C ***
C *******CLEAR THE MEMORY FOR SHAPE FUNCTION AND ITS DERIVATIVES
C ***
            DO 10 I=1,9
            VN(I) =0.0                          ! VN-SHAPE FUNCTION
            DO 10 J=1,3
            VDN (J, I)= 0. 0                   !  VDN--LOCAL DERIVATIVE OF VN
10         CONTINUE
C ********* SET FUNCTION VALUE FOR QUADRANGLE ELEMENT OF 4 NODES
C ***
            VN(1)= -0.25*((1-X)*(1-Y)*(1+X+Y))
            VN(2)= -0.25*((1+X)*(1-Y)*(1-X+Y))
            VN(3)= -0.25*((1+X)*(1+Y)*(1-X-Y))
            VN(4)= -0.25*((1-X)*(1+Y)*(1+X-Y))
		  VN(5)= 0.5*(1-X*X)*(1-Y)
		  VN(6)= 0.5*(1+X)*(1-Y*Y)
		  VN(7)= 0.5*(1-X*X)*(1+Y)
		  VN(8)= 0.5*(1-X)*(1-Y*Y)

            VDN(1,1) = 0.25*(2*X+Y)*(1-Y)
            VDN(1,2) = 0.25*(2*X-Y)*(1-Y)
            VDN(1,3) = 0.25*(2*X+Y)*(1+Y)
            VDN(1,4) = 0.25*(2*X-Y)*(1+Y)
            VDN(1,5) = -X*(1-Y)
            VDN(1,6) =  0.5*(1-Y*Y)
            VDN(1,7) = -X*(1+Y)
            VDN(1,8) = -0.5*(1-Y*Y)
	      
		  VDN(2,1) = 0.25*(2*Y+X)*(1-X)
            VDN(2,2) = 0.25*(2*Y-X)*(1+X)
            VDN(2,3) = 0.25*(2*Y+X)*(1+X)
            VDN(2,4) = 0.25*(2*Y-X)*(1-X)
            VDN(2,5) = -0.5*(1-X*X)
            VDN(2,6) = -Y*(1+X)
            VDN(2,7) =  0.5*(1-X*X)
            VDN(2,8) = -Y*(1-X)
 !         
            RETURN

            END