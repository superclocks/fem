!!	P725--726

C  .........................  SUB: 3-1-1  .............................
            SUBROUTINE    GAUSS_INTEGRATION (INTX,INTY,I,J,X,Y,WXY)
C                  GET THE INFORMATION OF GAUSS INTEGRATION POINT
            IMPLICIT  REAL * 8 (A-H, O-Z)
            DIMENSION     GXY(3,3) ,WG(3,3)
C ********* GAUSS INTEGRATION CONSTANTS FOR 1, 2 AND 3 POINTS
            GXY(1,1) = 0.0
            WG(1,1)= 2.0
            GXY(1,2) =-0. 577350269189626
            GXY(2,2)=  0. 577350269189626
            WG(1,2) = 1.0
            WG(2,2)= 1.0
            GXY(1,3) =-0. 774596669241483
	      GXY(2,3)= 0.0
            GXY(3,3)= 0.774596669241483
            WG(1,3)= 0.555555555555556
            WG(2,3)= 0.888888888888889
            WG(3,3)= 0.555555555555556
C *** GET PARAMETERS OF INTEGRATION POINT
            X= GXY(I, INTX)
            Y= GXY(J , INTY)
            WXY= WG(I, INTX) * WG(J, INTY)
            RETURN
            END