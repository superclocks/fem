!!	P734--736

C  ....................  SUB: 3-1-8  ..............................
            SUBROUTINE    ELEMENT_VB (MPROB, MND, VXY,VN,VD0,VDD,VB,SR)
C          FORM ELEMENT STRAIN MATRIX [BI  ACCORDING TO TYPE OF PROBLEM
            IMPLICIT  REAL *8 (A-H, O-Z)
            DIMENSION    VXY(MND,2),VN(9),VD0(3,9) ,VB(5,27),VX(4),
     $    VY(4),VL(4),VC(4),VS(4),VDD(3,9)
                                                          !  VXY-NODAL COORDINATES OF ELEMENT
                                                          !  VN--SHAPE FUNCTION
                                                          !  VD0--GLOBAL DERIVATIVE OF VN
                                                          !  VB-ELEMENTAL STRAIN MATRIX
C ***
C *******CLEAR MEMORY FOR THE MATRIX [BI
            VX(1)=VXY(2,1)-VXY(1,1)
		  VX(2)=VXY(3,1)-VXY(2,1)
		  VX(3)=VXY(4,1)-VXY(3,1)
	      VX(4)=VXY(1,1)-VXY(4,1)
		
	      VY(1)=VXY(2,2)-VXY(1,2)
		  VY(2)=VXY(3,2)-VXY(2,2)
		  VY(3)=VXY(4,2)-VXY(3,2)
	      VY(4)=VXY(1,2)-VXY(4,2)

	VL(1)=SQRT(VX(1)*VX(1)+VY(1)*VY(1))
      VL(2)=SQRT(VX(2)*VX(2)+VY(2)*VY(2))
	VL(3)=SQRT(VX(3)*VX(3)+VY(3)*VY(3))
	VL(4)=SQRT(VX(4)*VX(4)+VY(4)*VY(4))

	VC(1)=VY(1)/VL(1)
      VC(2)=VY(2)/VL(2)
	VC(3)=VY(3)/VL(3)
	VC(4)=VY(4)/VL(4)

	VS(1)=-VX(1)/VL(1)
      VS(2)=-VX(2)/VL(2)
	VS(3)=-VX(3)/VL(3)
	VS(4)=-VX(4)/VL(4)
C ***
            MPROB=1
            DO 30 II=1,27
            DO 30 JJ=1,5
            VB(JJ, II) = 0.0
30          CONTINUE
C***
C *** COMPUTE [B] MATRIX FOR PLANE STRESS OR PLANE STRAIN PROBLEM
             
             IF(MPROB. EQ. 1. OR. MPROB. EQ. 2)THEN
             SR= 1.0
             DO 40 II=1,4
             VB(1,(II-1) * 3+1)=VDD(1,II)   ! MPROB=1--PLANE STRESS
             VB(2,(II-1) * 3+2)=VDD(2,II)
             VB(3,(II-1) * 3+1)=VDD(2,II)   ! MPROB=2--PLANE STRAIN
             VB(3,(II-1) * 3+2)=VDD(1,II)

40           CONTINUE
             VB(1,3)=VL(4)*(-VC(4))*VD0(1,8)-VL(1)*VC(1)*VD0(1,5)
		   VB(2,3)=VL(4)*(-VS(4))*VD0(2,8)-VL(1)*VS(1)*VD0(2,5)
		   VB(3,3)=VL(4)*(-VC(4))*VD0(2,8)-VL(1)*VC(1)*VD0(2,5)+
     $	           VL(4)*(-VS(4))*VD0(1,8)-VL(1)*VS(1)*VD0(1,5)
		   VB(1,6)=VL(1)*(-VC(1))*VD0(1,5)-VL(2)*VC(2)*VD0(1,6)
	       VB(2,6)=VL(1)*(-VS(1))*VD0(2,5)-VL(2)*VS(2)*VD0(2,6)
		   VB(3,6)=VL(1)*(-VC(1))*VD0(2,5)-VL(2)*VC(2)*VD0(2,6)+
     $			   VL(1)*(-VS(1))*VD0(1,5)-VL(2)*VS(2)*VD0(1,6)
		   VB(1,9)=VL(2)*(-VC(2))*VD0(1,6)-VL(3)*VC(3)*VD0(1,7)
	       VB(2,9)=VL(2)*(-VS(2))*VD0(2,6)-VL(3)*VS(3)*VD0(2,7)
	       VB(3,9)=VL(2)*(-VC(2))*VD0(2,6)-VL(3)*VC(3)*VD0(2,7)+
     $               VL(2)*(-VS(2))*VD0(1,6)-VL(3)*VS(3)*VD0(1,7)
		   VB(1,12)=VL(3)*(-VC(3))*VD0(1,7)-VL(4)*VC(4)*VD0(1,8)
             VB(2,12)=VL(3)*(-VS(3))*VD0(2,7)-VL(4)*VS(4)*VD0(2,8)
             VB(3,12)=VL(3)*(-VC(3))*VD0(2,7)-VL(4)*VC(4)*VD0(2,8)+
     $                VL(3)*(-VS(3))*VD0(1,7)-VL(4)*VS(4)*VD0(1,8)

             ENDIF

C ******
C ****** COMPUTE [B] MATRIX FOR AXISYMMETRIC PROBLEM
C ***
             IF(MPROB. EQ. 3)  THEN
             SR=0. 0
             DO 45 II=1,MND
             SR=SR+VXY(II,1) * VN(II)
45          CONTINUE                ! MPROB= 3--AXISYMMETRIC
             DO 50 II=1,9
             VB(1,(II-1) * 2+1) =VD0(1,II)
             VB(2,(II-1) * 2+2) =VD0(2,II)
             VB(3,(II-1) * 2+1) =VD0(2,II)
             VB(3,(II-1) * 2+2) =VD0(1,II)
            IF(ABS(SR). LT. 1E-20) VB(4,(II-1) *2+1)=0.0
            IF(ABS(SR). GE. 1E-20) VB(4,(II-1) *2+1)=VN(II)/SR
50           CONTINUE
            SR=SR* 2 * 3.14159265
            ENDIF
C ******* COMPUTE [B] MATRIX FOR THE MINDLIN PLATE PROBLEM
C*********
            IF(MPROB. EQ. 4)  THEN
            SR  =  1.0                              !  MPROB=4--MINDLIN PLATE
            DO 70 II=1,9
            VB(1, (II-1) * 3+1) =-VD0(1,II)
            VB(2, (II-1) * 3+2) =-VD0(2,II)
            VB(3, (II-1) * 3+1) =-VD0(2,II)
            VB(3, (II-1) * 3+2) =-VD0(1,II)
            VB(4, (II-1) * 3+1) =-VN(II)
            VB(4, (II-1) * 3+3) = VD0(1,II)            	  
		  VB(5, (II-1) * 3+2) =-VN(II)
            VB(5, (II-1) * 3+3) = VD0(2,II)
70        CONTINUE
            ENDIF
            RETURN
            END