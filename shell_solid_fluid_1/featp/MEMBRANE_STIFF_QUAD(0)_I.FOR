	
	
      SUBROUTINE MEMBRANE_STIFF_QUAD0_I(VXY,VK,VM,X,Y,VMAE,WXY,IEL)

      IMPLICIT   REAL * 8 (A-H,O-z)
      DIMENSION  IEL(4) ,VN(4),VNN(3,12),VDN(2,6),VDNN(2,2),VJ1(2,2),
     $          VXY(4,3),VD(3,3),AB(3,2),
     $         VJJ(2,2),VJ(2,2),VDD(2,6),VJJ1(2,2),
     $         VB(3,12),VSS(3,12),VK(12,12),
     $         VMAE(4),VM(12,12)
    
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
      VN(:)=0.0
	VDN(:,:)=0.0
      VN(1)= (1-X)*(1-Y)*0.25
      VN(2)= (1+X)*(1-Y)*0.25
      VN(3)= (1+X)*(1+Y)*0.25
      VN(4)= (1-X)*(1+Y)*0.25    !四节点单元形函数
      
      VDN(1,1) =-(1-Y)*0.25
      VDN(1,2) = (1-Y)*0.25
      VDN(1,3) = (1+Y)*0.25
      VDN(1,4) =-(1+Y)*0.25


      VDN(2,1) =-(1-X)*0.25
      VDN(2,2) =-(1+X)*0.25
      VDN(2,3) = (1+X)*0.25
      VDN(2,4) = (1-X)*0.25 !四节点单元形函数对全局坐标的导数

!	IF(AJ0.EQ.0)THEN       
      VDN(1,5) =-2.0*X
      VDN(2,6) =-2.0*Y 
!      ELSE
!      VDN(1,5)=-2.0*X+2.0*A/(3.0*AJ0)
!      VDN(1,6)=-2.0*A/(3.0*AJ0)
!	VDN(2,5)=-2.0*B/(3.0*AJ0)
!	VDN(2,6)=-2.0*Y+2.0*B/(3.0*AJ0)  !WU的函数
 
 !     VDN(1,5)=-2.0*(X+(B/AJ0)*X*Y)*(1.0+(B/AJ0)*Y)
 !     VDN(1,6)=-2.0*(Y-(A/AJ0)*X*Y)*(-(A/AJ0)*Y)
!	VDN(2,5)=-2.0*(X+(B/AJ0)*X*Y)*((B/AJ0)*X)
!	VDN(2,6)=-2.0*(Y-(A/AJ0)*X*Y)*(1.0-(A/AJ0)*X)   !Wilson的改进函数

!	ENDIF


      DO 20 I=1,2
	DO 20 J=1,2
	VJJ(I,J)=0.0
	DO 20 K=1,4
      VJJ(I,J)=VJJ(I,J)+VDN(I,K)*VXY(K,J)  !雅克比变换
20    CONTINUE
	
	SJ=VJJ(1,1)*VJJ(2,2)-VJJ(1,2)*VJJ(2,1)  !雅克比行列式

	VJ(1,1)= VJJ(2,2)/SJ
	VJ(1,2)=-VJJ(1,2)/SJ
	VJ(2,1)=-VJJ(2,1)/SJ
	VJ(2,2)= VJJ(1,1)/SJ                 !雅克比矩阵的逆阵

	VJ1(1,1)= AB(3,2)/SJ
	VJ1(1,2)=-AB(1,2)/SJ
	VJ1(2,1)=-AB(3,1)/SJ
	VJ1(2,2)= AB(1,1)/SJ     !零点雅克比矩阵逆阵,求逆犯规法


	DO 21 I=1,2
	DO 21 J=1,4
	VDD(I,J)=0.0
	DO 21 K=1,2
	VDD(I,J)=VDD(I,J)+VJ(I,K)*VDN(K,J)	!4节点形函数对全局坐标点的导数
21    CONTINUE
       
	DO 22 I=1,2
	DO 22 J=5,6
	VDD(I,J)=0.0
	DO 22 K=1,2
      VDD(I,J)=VDD(I,J)+VJ1(I,K)*VDN(K,J)
22    CONTINUE

      VB(:,:)=0.0
	DO 13 I=1,4
	VB(1,(I-1)*2+1)=VDD(1,I)
      VB(2,(I-1)*2+2)=VDD(2,I)
      VB(3,(I-1)*2+1)=VDD(2,I)
	VB(3,(I-1)*2+2)=VDD(1,I)    !应变矩阵
13    CONTINUE
      
!	VB(1,9)= -2.0*X/AL12
!	VB(2,12)=-2.0*Y/AL41
!	VB(3,10)=-2.0*Y/AL23
!	VB(3,11)=-2.0*X/AL34

!      VJJ1(1,1)=0.25*(-VXY(1,1)+VXY(2,1)+VXY(3,1)-VXY(4,1))
!      VJJ1(1,2)=0.25*(-VXY(1,2)+VXY(2,2)+VXY(3,2)-VXY(4,2))
!      VJJ1(2,1)=0.25*(-VXY(1,1)-VXY(2,1)+VXY(3,1)+VXY(4,1))
!      VJJ1(2,2)=0.25*(-VXY(1,2)-VXY(2,2)+VXY(3,2)+VXY(4,2))

!**************Wilson的方法******************
!      VB(1,9)= -2.0*X/VJJ1(1,1)
!	VB(2,12)=-2.0*Y/VJJ1(2,2)
!	VB(3,10)=-2.0*Y/VJJ1(2,2)
!	VB(3,11)=-2.0*X/VJJ1(1,1)
!      IF(VJJ1(2,1).NE.0.0)THEN
!	VB(1,10)=-2.0*Y/VJJ1(2,1)
!	VB(3,12)=-2.0*Y/VJJ1(2,1)
!      ENDIF

!	IF(VJJ1(1,2).NE.0.0)THEN
!	VB(2,11)=-2.0*X/VJJ1(1,2)
 !     VB(3,9)=-2.0*X/VJJ1(1,2)
 !     ENDIF
!********************************************
!****************WU的方法********************
      VB(1,9)=  VDD(1,5)
	VB(1,11)= VDD(1,6)
	VB(2,10)= VDD(2,5)
	VB(2,12)= VDD(2,6)
      VB(3,9)= VDD(2,5)
	VB(3,10)= VDD(1,5)
	VB(3,11)= VDD(2,6)
	VB(3,12)= VDD(1,6)
!------------------------------------------
!      IF(VJJ(2,1).NE.0)THEN
!	VB(1,9)=  VDN(1,5)/VJJ(1,1)+VDN(2,5)/VJJ(2,1)
!	VB(1,10)= VDN(1,6)/VJJ(1,1)+VDN(2,6)/VJJ(2,1)
!	VB(3,11)= VDN(1,5)/VJJ(1,1)+VDN(2,5)/VJJ(2,1)
!	VB(3,12)= VDN(1,6)/VJJ(1,1)+VDN(2,6)/VJJ(2,1)
!	ELSE
!      VB(1,9)=  VDN(1,5)/VJJ(1,1)
!	VB(1,10)= VDN(1,6)/VJJ(1,1)
!	VB(3,11)= VDN(1,5)/VJJ(1,1)
!	VB(3,12)= VDN(1,6)/VJJ(1,1)

!	ENDIF
      
!	IF(VJJ(1,2).NE.0)THEN
!	VB(2,11)= VDN(1,5)/VJJ(1,2)+VDN(2,5)/VJJ(2,2)
!	VB(2,12)= VDD(1,6)/VJJ(1,2)+VDN(2,6)/VJJ(2,2)
!      VB(3,9)= VDN(1,5)/VJJ(1,2)+VDN(2,5)/VJJ(2,2)
!	VB(3,10)= VDD(1,6)/VJJ(1,2)+VDN(2,6)/VJJ(2,2)
!      ELSE 
!	VB(2,11)=VDN(2,5)/VJJ(2,2)
!	VB(2,12)=VDN(2,6)/VJJ(2,2)
!      VB(3,9)= VDN(2,5)/VJJ(2,2)
!	VB(3,10)=VDN(2,6)/VJJ(2,2)

!	ENDIF


!********************************************
      
	CALL MEMBRANE_VD(VD,VMAE)  !调用弹性矩阵子函数
      DO 14 I=1,3
	DO 14 J=1,12   
	VSS(I,J)=0.0
	DO 14 K=1,3
      VSS(I,J)=VSS(I,J)+VD(I,K)*VB(K,J)
14    CONTINUE

      DO 15 I=1,12
	DO 15 J=1,12
	VK(I,J)=0.0
	DO 15 K=1,3
      VK(I,J)=VK(I,J)+VB(K,I)*VSS(K,J)*WXY*SJ*VMAE(4) !计算刚度矩阵
15    CONTINUE
     	
	VNN(:,:)=0.0
	DO 16 I=1,4
	VNN(1,(I-1)*3+1)=VN(I)
      VNN(2,(I-1)*3+2)=VN(I)
!	VNN(3,(I-1)*3+3)=VN(I)
16    CONTINUE

      DO 17 I=1,12
	DO 17 J=1,12
	VM(I,J)=0.0
	DO 17 K=1,3
	VM(I,J)=VM(I,J)+VNN(K,I)*VNN(K,J)*SJ*VMAE(3)*VMAE(4) !计算质量矩阵
17    CONTINUE      
	
      


    	RETURN
	END
