
      SUBROUTINE MEMBRANE_STIFF_QUAD2_I(VXY,VK,VM,X,Y,VMAE,WXY)
      IMPLICIT   REAL * 8 (A-H,O-z)


      DIMENSION  VN(9),VN0(4),VDN(2,10),VDN0(2,4),
     $          VXY(4,3),VX(8,3),VD0(2,4),
     $         VJJ(2,2),VJJ0(2,2),VJ(2,2),VJ0(2,2),VDD(2,10),VDD0(2,4),
     $         VB(3,18),VSS(3,18),VK(18,18),
     $         VD(3,3),VMAE(4),VDN1(2,4),VM(12,12),VNN(3,12),
     $         VXX(4),VL(4),VC(4),VS(4),BI(12),VJ1(2,2),AB(3,2)
       
	AB(1,1)=0.25*(-VXY(1,1)+VXY(2,1)+VXY(3,1)-VXY(4,1))
	AB(1,2)=0.25*(-VXY(1,2)+VXY(2,2)+VXY(3,2)-VXY(4,2))

	AB(2,1)=0.25*(VXY(1,1)-VXY(2,1)+VXY(3,1)-VXY(4,1))
	AB(2,2)=0.25*(VXY(1,2)-VXY(2,2)+VXY(3,2)-VXY(4,2))

	AB(3,1)=0.25*(-VXY(1,1)-VXY(2,1)+VXY(3,1)+VXY(4,1))
	AB(3,2)=0.25*(-VXY(1,2)-VXY(2,2)+VXY(3,2)+VXY(4,2))

	AJ0=AB(1,1)*AB(3,2)-AB(3,1)*AB(1,2)
	A=AB(1,1)*AB(2,2)-AB(2,1)*AB(1,2)
	B=AB(2,1)*AB(3,2)-AB(3,1)*AB(2,2)


C*******计算4节点的形函数、导数、雅克比行列式及形函数对全局坐标的导数*******
      VN0(:)=0.0
	VDN0(:,:)=0.0
      VN0(1)= (1-X)*(1-Y)*0.25
      VN0(2)= (1+X)*(1-Y)*0.25
      VN0(3)= (1+X)*(1+Y)*0.25
      VN0(4)= (1-X)*(1+Y)*0.25
     
      VDN0(1,1) =-(1-Y)*0.25
      VDN0(1,2) = (1-Y)*0.25
      VDN0(1,3) = (1+Y)*0.25
      VDN0(1,4) =-(1+Y)*0.25


      VDN0(2,1) =-(1-X)*0.25
      VDN0(2,2) =-(1+X)*0.25
      VDN0(2,3) = (1+X)*0.25
      VDN0(2,4) = (1-X)*0.25

      DO 20 I=1,2
	DO 20 J=1,2
	VJJ0(I,J)=0.0
	DO 20 K=1,4
      VJJ0(I,J)=VJJ0(I,J)+VDN0(I,K)*VXY(K,J)
20    CONTINUE
	
	SJ0=VJJ0(1,1)*VJJ0(2,2)-VJJ0(1,2)*VJJ0(2,1)

	VJ0(1,1)= VJJ0(2,2)/SJ0
	VJ0(1,2)=-VJJ0(1,2)/SJ0
	VJ0(2,1)=-VJJ0(2,1)/SJ0
	VJ0(2,2)= VJJ0(1,1)/SJ0

	DO 21 I=1,2
	DO 21 J=1,4
	VDD0(I,J)=0.0
	DO 21 K=1,2
	VDD0(I,J)=VDD0(I,J)+VJ0(I,K)*VDN0(K,J)	!4节点形函数对全局坐标点的导数
21    CONTINUE

C****************************************************************************
C************计算8节点的形函数、导数、雅克比行列式及形函数对全局坐标的导数********
      VN(:)=0.0
	VDN(:,:)=0.0
	VX(:,:)=0.0
      VN(1)= (1-X)*(1-Y)*(-X-Y-1.0)*0.25
      VN(2)= (1+X)*(1-Y)*(X-Y-1.0)*0.25
      VN(3)= (1+X)*(1+Y)*(X+Y-1.0)*0.25
      VN(4)= (1-X)*(1+Y)*(-X+Y-1.0)*0.25
      VDN(1,1) = (1-Y)*(2.0*X+Y)*0.25 
      VDN(1,2) = (1-Y)*(2.0*X-Y)*0.25
      VDN(1,3) = (1+Y)*(2.0*X+Y)*0.25
      VDN(1,4) = (1+Y)*(2.0*X-Y)*0.25

      VDN(2,1) = (1-X)*(2.0*Y+X)*0.25
      VDN(2,2) = (1+X)*(2.0*Y-X)*0.25
      VDN(2,3) = (1+X)*(2.0*Y+X)*0.25
      VDN(2,4) = (1-X)*(2.0*Y-X)*0.25

      VN(5)= (1-X*X)*(1-Y)*0.5
      VN(6)= (1-Y*Y)*(1+X)*0.5
      VN(7)= (1-X*X)*(1+Y)*0.5
      VN(8)= (1-Y*Y)*(1-X)*0.5 

      VDN(1,5)= (-2*X)*(1-Y)*0.5
      VDN(1,6)= (1-Y*Y)*(+1)*0.5
      VDN(1,7)= (-2*X)*(1+Y)*0.5
      VDN(1,8)= (1-Y*Y)*(-1)*0.5
   
     
      VDN(2,5)= (1-X*X)*(-1)*0.5
      VDN(2,6)= (-2*Y)*(1+X)*0.5
      VDN(2,7)= (1-X*X)*(+1)*0.5
      VDN(2,8)= (-2*Y)*(1-X)*0.5
      
	VDN(1,9)=-2.0*X
!     VDN(1,10)=-2.0*X
!	VDN(2,9)=-2.0*Y
	VDN(2,10)=-2.0*Y
	
	VX(1:4,:)=VXY(1:4,:)
	VX(5,1)=(VXY(1,1)+VXY(2,1))*0.5
	VX(6,1)=(VXY(3,1)+VXY(2,1))*0.5
	VX(7,1)=(VXY(4,1)+VXY(3,1))*0.5
      VX(8,1)=(VXY(1,1)+VXY(4,1))*0.5

	VX(5,2)=(VXY(1,2)+VXY(2,2))*0.5
	VX(6,2)=(VXY(3,2)+VXY(2,2))*0.5
	VX(7,2)=(VXY(4,2)+VXY(3,2))*0.5
      VX(8,2)=(VXY(1,2)+VXY(4,2))*0.5

	DO 11 I=1,2
	DO 11 J=1,2
	VJJ(I,J)=0.0
	DO 11 K=1,8
      VJJ(I,J)=VJJ(I,J)+VDN(I,K)*VX(K,J)
11    CONTINUE

      SJ=VJJ(1,1)*VJJ(2,2)-VJJ(1,2)*VJJ(2,1)
	VJ(1,1)= VJJ(2,2)/SJ
	VJ(1,2)=-VJJ(1,2)/SJ
	VJ(2,1)=-VJJ(2,1)/SJ
	VJ(2,2)= VJJ(1,1)/SJ

	SJ1=AB(1,1)*AB(3,2)-AB(1,2)*AB(3,1)
	VJ1(1,1)= AB(3,2)/SJ
	VJ1(1,2)=-AB(1,2)/SJ
	VJ1(2,1)=-AB(3,1)/SJ
	VJ1(2,2)= AB(1,1)/SJ     !零点雅克比矩阵逆阵


	DO 12 I=1,2
	DO 12 J=1,10
	VDD(I,J)=0.0
	DO 12 K=1,2
	VDD(I,J)=VDD(I,J)+VJ(I,K)*VDN(K,J)	!8节点形函数对全局坐标点的导数
12    CONTINUE

      DO 30 I=1,2
	DO 30 J=9,10
	VDD(I,J)=0.0
	DO 30 K=1,2
      VDD(I,J)=VDD(I,J)+VJ1(I,K)*VDN(K,J)   !零点雅克比矩阵
30    CONTINUE
C****************************************************************************


      VB(:,:)=0.0
	DO 13 I=1,8
	VB(1,(I-1)*2+1)=VDD(1,I)
      VB(2,(I-1)*2+2)=VDD(2,I)
      VB(3,(I-1)*2+1)=VDD(2,I)
	VB(3,(I-1)*2+2)=VDD(1,I)
13    CONTINUE

      AL12=0.5*SQRT((VXY(2,1)-VXY(1,1))*(VXY(2,1)-VXY(1,1))+
     $	         (VXY(2,2)-VXY(1,2))*(VXY(2,2)-VXY(1,2)))
	AL23=0.5*SQRT((VXY(3,1)-VXY(2,1))*(VXY(3,1)-VXY(2,1))+
     $             (VXY(3,2)-VXY(2,2))*(VXY(3,2)-VXY(2,2)))
	AL34=0.5*SQRT((VXY(4,1)-VXY(3,1))*(VXY(4,1)-VXY(3,1))+
     $	         (VXY(4,2)-VXY(3,2))*(VXY(4,2)-VXY(3,2)))
      AL41=0.5*SQRT((VXY(1,1)-VXY(4,1))*(VXY(1,1)-VXY(4,1))+
     $	         (VXY(1,2)-VXY(4,2))*(VXY(1,2)-VXY(4,2)))
      
!	IF(X.LT.0.0.AND.Y.LT.0.0)THEN
!	VB(1,17)=-2.0*X/AL12
!	VB(2,18)=-2.0*Y/AL34
!	ELSEIF(X.LT.0.0.AND.Y.GT.0.0)THEN
	VB(1,17)=VDD(1,9)
!      VB(1,18)=VDD(1,10)

!	VB(2,17)=VDD(2,9)
	VB(2,18)=VDD(2,10)

	VB(3,17)=VDD(2,9)
	VB(3,18)=VDD(1,10)

!	ELSEIF(X.GT.0.0.AND.Y.LT.0.0)THEN
 !     VB(1,17)=-2.0*X/AL12
!	VB(2,18)=-2.0*Y/AL23
!	ELSEIF(X.GT.0.0.AND.Y.GT.0.0)THEN
!      VB(1,17)=-2.0*X/AL34
!	VB(2,18)=-2.0*Y/AL23
!	ENDIF
      
	CALL MEMBRANE_VD(VD,VMAE)
      DO 14 I=1,3
	DO 14 J=1,18   !12
	VSS(I,J)=0.0
	DO 14 K=1,3
      VSS(I,J)=VSS(I,J)+VD(I,K)*VB(K,J)
14    CONTINUE

      DO 15 I=1,18   !12
	DO 15 J=1,18    !12
	VK(I,J)=0.0
	DO 15 K=1,3
      VK(I,J)=VK(I,J)+VB(K,I)*VSS(K,J)*WXY*SJ*VMAE(4)
15    CONTINUE
      
      
      VNN(:,:)=0.0
	DO 17 I=1,4
      VNN(1,(I-1)*3+1)=VN0(I)
      VNN(2,(I-1)*3+2)=VN0(I)
17    CONTINUE
	VM(:,:)=0.0
	DO 16 I=1,12
	DO 16 J=1,12
      DO 16 K=1,3
      VM(I,J)=VM(I,J)+VNN(K,I)*VNN(K,J)*VMAE(3)*SJ*VMAE(4)
16    CONTINUE
      

    	RETURN
	END
