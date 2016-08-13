
      SUBROUTINE MEMBRANE_STIFF_QUAD(2)_I(VXY,VX,VK,VP,VM,X,Y,VMAE,WXY,IEL)
      IMPLICIT   REAL * 8 (A-H,O-z)

	COMMON/COM1/MND, NUMEL, NUMPT, MBAND, NMATI
      COMMON/COM2/NF, NFSTR, MSOLV, MPROB, MTYPE, NVA
      COMMON/COM3/MND2 , NUMPT2
      COMMON/ELEM/NODE, INTX, INTY

      DIMENSION  IEL(NODE) ,VN(9),VN0(4),VDN(3,9),VDN0(2,4),
     $          VXY(4,3),VX(8,3),VD0(2,4),
     $         VJJ(2,2),VJJ0(2,2),VJ(2,2),VJ0(2,2),VDD(2,8),VDD0(2,4),
     $         VB(3,16),VSS(3,16),VK(16,16),
     $         VD(3,3),VMAE(4),VDN1(2,4),VM(16,16),VNN(2,16),
     $         VXX(4),VY(4),VL(4),VC(4),VS(4),BI(12),VP(12,12)



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

!	VDN1(1:2,1:4)=VDN(1:2,1:4)

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

	DO 12 I=1,2
	DO 12 J=1,8
	VDD(I,J)=0.0
	DO 12 K=1,2
	VDD(I,J)=VDD(I,J)+VJ(I,K)*VDN(K,J)	!8节点形函数对全局坐标点的导数
12    CONTINUE
C****************************************************************************


      VB(:,:)=0.0
	DO 13 I=1,8
	VB(1,(I-1)*2+1)=VDD(1,I)
      VB(2,(I-1)*2+2)=VDD(2,I)
      VB(3,(I-1)*2+1)=VDD(2,I)
	VB(3,(I-1)*2+2)=VDD(1,I)
13    CONTINUE

!	DO 13 I=1,4
!	VB(1,(I-1)*3+1)=VDD0(1,I)
!      VB(2,(I-1)*3+2)=VDD0(2,I)
!      VB(3,(I-1)*3+1)=VDD0(2,I)
!	VB(3,(I-1)*3+2)=VDD0(1,I)
!13    CONTINUE


!      VXX(1)=VXY(2,1)-VXY(1,1)
!	VXX(2)=VXY(3,1)-VXY(2,1)
!	VXX(3)=VXY(4,1)-VXY(3,1)
!	VXX(4)=VXY(1,1)-VXY(4,1)
		
!	VY(1)=VXY(2,2)-VXY(1,2)
!      VY(2)=VXY(3,2)-VXY(2,2)
!      VY(3)=VXY(4,2)-VXY(3,2)
!	VY(4)=VXY(1,2)-VXY(4,2)

!	VL(1)=SQRT(VXX(1)*VXX(1)+VY(1)*VY(1))
!      VL(2)=SQRT(VXX(2)*VXX(2)+VY(2)*VY(2))
!	VL(3)=SQRT(VXX(3)*VXX(3)+VY(3)*VY(3))
!	VL(4)=SQRT(VXX(4)*VXX(4)+VY(4)*VY(4))

!	VC(1)=VY(1)/VL(1)
!      VC(2)=VY(2)/VL(2)
!	VC(3)=VY(3)/VL(3)
!	VC(4)=VY(4)/VL(4)

!	VS(1)=-VXX(1)/VL(1)
!      VS(2)=-VXX(2)/VL(2)
!	VS(3)=-VXX(3)/VL(3) 
!	VS(4)=-VXX(4)/VL(4)

!      VB(1,3)=(VY(4)*VDD(1,8)-VY(1)*VDD(1,5))/8
!      VB(2,3)=(-VXX(4)*VDD(2,8)+VXX(1)*VDD(2,5))/8
!      VB(3,3)=(VY(4)*VDD(2,8)-VY(1)*VDD(2,5)+
!     $	    -VXX(4)*VDD(1,8)+VXX(1)*VDD(1,5))/8
!      VB(1,6)=(VY(1)*VDD(1,5)-VY(2)*VDD(1,6))/8
!	VB(2,6)=(-VXX(1)*VDD(2,5)+VXX(2)*VDD(2,6))/8
!      VB(3,6)=(VY(1)*VDD(2,5)-VY(2)*VDD(2,6)+
!     $        -VXX(1)*VDD(1,5)+VXX(2)*VDD(1,6))/8
!      VB(1,9)=(VY(2)*VDD(1,6)-VY(3)*VDD(1,7))/8
!	VB(2,9)=(-VXX(2)*VDD(2,6)+VXX(3)*VDD(2,7))/8
!	VB(3,9)=(VY(2)*VDD(2,6)-VY(3)*VDD(2,7)+
!     $        -VXX(2)*VDD(1,6)+VXX(3)*VDD(1,7))/8
!      VB(1,12)=(VY(3)*VDD(1,7)-VY(4)*VDD(1,8))/8
!      VB(2,12)=(-VXX(3)*VDD(2,7)+VXX(4)*VDD(2,8))/8
!      VB(3,12)=(VY(3)*VDD(2,7)-VY(4)*VDD(2,8)+
!     $         -VXX(3)*VDD(1,7)+VXX(4)*VDD(1,8))/8

!	 BI(:)=0.0
!      BI(1)=-0.5*VDD0(2,1)
!	BI(2)= 0.5*VDD0(1,1)
!	BI(3)= -0.0625*(VL(4)*(-VC(4))*VDD(2,8)-VL(1)*VC(1)*VDD(2,5))
!     $    +0.0625*(VL(4)*(-VS(4))*VDD(1,8)-VL(1)*VS(1)*VDD(1,5))-VN0(1)
	
!	BI(4)=-0.5*VDD0(2,2)
!	BI(5)= 0.5*VDD0(1,2)
!	BI(6)=-0.0625*(VL(1)*(-VC(1))*VDD(2,5)-VL(2)*VC(2)*VDD(2,6))
!     $    +0.0625*(VL(1)*(-VS(1))*VDD(1,5)-VL(2)*VS(2)*VDD(1,6))-VN0(2)

!	BI(7)=-0.5*VDD0(2,3)
!      BI(8)= 0.5*VDD0(1,3)
!	BI(9)=-0.0625*(VL(2)*(-VC(2))*VDD(2,6)-VL(3)*VC(3)*VDD(2,7))
!     $    +0.0625*(VL(2)*(-VS(2))*VDD(1,6)-VL(3)*VS(3)*VDD(1,7))-VN0(3)
	 
!	BI(10)=-0.5*VDD0(2,4)
!	BI(11)= 0.5*VDD0(1,4)
 !     BI(12)= -0.0625*(VL(3)*(-VC(3))*VDD(2,7)-VL(4)*VC(4)*VDD(2,8))
!     $    +0.0625*(VL(3)*(-VS(3))*VDD(1,7)-VL(4)*VS(4)*VDD(1,8))-VN0(4)

	DO 30 I=1,12
      DO 30 J=1,12
	VP(I,J)=0.0
	VP(I,J)=VP(I,J)+BI(I)*BI(J)*SJ
30    CONTINUE
      
	CALL MEMBRANE_VD(MPROB,VD,VMAE)
      DO 14 I=1,3
	DO 14 J=1,16   !12
	VSS(I,J)=0.0
	DO 14 K=1,3
      VSS(I,J)=VSS(I,J)+VD(I,K)*VB(K,J)
14    CONTINUE

      DO 15 I=1,16   !12
	DO 15 J=1,16    !12
	VK(I,J)=0.0
	DO 15 K=1,3
      VK(I,J)=VK(I,J)+VB(K,I)*VSS(K,J)*WXY*SJ*VMAE(4)
15    CONTINUE
      
      
      VNN(:,:)=0.0
	DO 17 I=1,8
      VNN(1,(I-1)*2+1)=VN(I)
      VNN(2,(I-1)*2+2)=VN(I)
17    CONTINUE
	VM(:,:)=0.0
	DO 16 I=1,16
	DO 16 J=1,16
      DO 16 K=1,2
      VM(I,J)=VM(I,J)+VNN(K,I)*VNN(K,J)*VMAE(3)*SJ*0.1
16    CONTINUE
      

    	RETURN
	END
