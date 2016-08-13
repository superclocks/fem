	
	SUBROUTINE DKQ(VXY,X,Y,VN,VDN,B,A)


	IMPLICIT REAL * 8(A-H,O-Z)
               COMMON/COM1/MND, NUMEL, NUMPT,MBAND, NMATI
               COMMON/COM2/NF, NFSTR, MSOLV, MPROB, MTYPE, NVA
               COMMON/COMN/NFIX, NPC, GRAV
               COMMON/COM3/MND2,NUMPT2
C                   ! MND--MAXIMAL NODE NUMBER IN ALL ELEMENTS
C                   ! NUMEL--NUMBER OF GLOBAL ELEMENTS
C                   ! NUMPT--NUMBER OF GLOBAL NODES
C                   ! MBAND--HALF BANDWIDTH(INCLUDING DIAGONAL ELEMENT)
	DIMENSION VXY(NUMPT,3),AX(4),AY(4),AL(4),AC(4),AS(4),AJ(2,2),
     $          AJJ(2,2),AA(4),AB(4),ACK(4),AD(4),AE(4),AH(2,12),
     $          VN(9),VDN(3,9),B(3,12),ADH(4,12)
	ADH(:,:)=0.0
	AX(1)=VXY(1,1)-VXY(2,1)
	AX(2)=VXY(2,1)-VXY(3,1)
	AX(3)=VXY(3,1)-VXY(4,1)
	AX(4)=VXY(4,1)-VXY(1,1)

	AY(1)=VXY(1,2)-VXY(2,2)
	AY(2)=VXY(2,2)-VXY(3,2)
	AY(3)=VXY(3,2)-VXY(4,2)
	AY(4)=VXY(4,2)-VXY(1,2)

 !     AX(3)=VXY(1,1)-VXY(2,1)
!	AX(4)=VXY(2,1)-VXY(3,1)
!	AX(1)=VXY(3,1)-VXY(4,1)
!	AX(2)=VXY(4,1)-VXY(1,1)

!	AY(3)=VXY(1,2)-VXY(2,2)
!	AY(4)=VXY(2,2)-VXY(3,2)
!	AY(1)=VXY(3,2)-VXY(4,2)
!	AY(2)=VXY(4,2)-VXY(1,2)

	AL(1)=SQRT(AX(1)*AX(1)+AY(1)*AY(1))
	AL(2)=SQRT(AX(2)*AX(2)+AY(2)*AY(2))
	AL(3)=SQRT(AX(3)*AX(3)+AY(3)*AY(3))
	AL(4)=SQRT(AX(4)*AX(4)+AY(4)*AY(4))

!	AS(1)=AX(1)/AL(1)
!	AS(2)=AX(2)/AL(2)
!	AS(3)=AX(3)/AL(3)
!	AS(4)=AX(4)/AL(4)

!	AC(1)=-AY(1)/AL(1)
!	AC(2)=-AY(2)/AL(2)
!	AC(3)=-AY(3)/AL(3)
!	AC(4)=-AY(4)/AL(4)

	AJ(1,1)=0.25*(-AX(1)+AX(3)+Y*(AX(1)+AX(3)))
	AJ(1,2)=0.25*(-AY(1)+AY(3)+Y*(AY(1)+AY(3)))
	AJ(2,1)=0.25*(-AX(2)+AX(4)+X*(AX(1)+AX(3)))
	AJ(2,2)=0.25*(-AY(2)+AY(4)+X*(AY(1)+AY(3)))
	A=0.125*((VXY(4,2)-VXY(2,2))*(VXY(3,1)-VXY(1,1))-
     $         (VXY(3,2)-VXY(1,2))*(VXY(4,1)-VXY(2,1)))+
     $   0.125*X*((VXY(3,2)-VXY(4,2))*(VXY(2,1)-VXY(1,1))-
     $         (VXY(2,2)-VXY(1,2))*(VXY(3,1)-VXY(4,1)))+
     $   0.125*Y*((VXY(4,2)-VXY(1,2))*(VXY(3,1)-VXY(2,1))-
     $         (VXY(3,2)-VXY(2,2))*(VXY(4,1)-VXY(1,1)))
!	A=0.25*0.25*(AJ(1,1)*AJ(2,2)-AJ(1,2)*AJ(2,1))
	AJJ(1,1)= AJ(2,2)/A
	AJJ(1,2)=-AJ(1,2)/A
	AJJ(2,1)=-AJ(2,1)/A
	AJJ(2,2)= AJ(2,2)/A

	AA(1)=-AX(1)/(AL(1)*AL(1))
	AA(2)=-AX(2)/(AL(2)*AL(2))
	AA(3)=-AX(3)/(AL(3)*AL(3))
	AA(4)=-AX(4)/(AL(4)*AL(4))

	AB(1)=(0.75*AX(1)*AY(1))/(AL(1)*AL(1))
	AB(2)=(0.75*AX(2)*AY(2))/(AL(2)*AL(2))
	AB(3)=(0.75*AX(3)*AY(3))/(AL(3)*AL(3))
	AB(4)=(0.75*AX(4)*AY(4))/(AL(4)*AL(4))

	ACK(1)=(0.25*AX(1)*AX(1)-0.5*AY(1)*AY(1))/(AL(1)*AL(1))
	ACK(2)=(0.25*AX(2)*AX(2)-0.5*AY(2)*AY(2))/(AL(2)*AL(2))
	ACK(3)=(0.25*AX(3)*AX(3)-0.5*AY(3)*AY(3))/(AL(3)*AL(3))
	ACK(4)=(0.25*AX(4)*AX(4)-0.5*AY(4)*AY(4))/(AL(4)*AL(4))

	AD(1)=-AY(1)/(AL(1)*AL(1))
	AD(2)=-AY(2)/(AL(2)*AL(2))
	AD(3)=-AY(3)/(AL(3)*AL(3))
	AD(4)=-AY(4)/(AL(4)*AL(4))

	AE(1)=(-0.5*AX(1)*AX(1)+0.25*AY(1)*AY(1))/(AL(1)*AL(1))
	AE(2)=(-0.5*AX(2)*AX(2)+0.25*AY(2)*AY(2))/(AL(2)*AL(2))
	AE(3)=(-0.5*AX(3)*AX(3)+0.25*AY(3)*AY(3))/(AL(3)*AL(3))
	AE(4)=(-0.5*AX(4)*AX(4)+0.25*AY(4)*AY(4))/(AL(4)*AL(4))

!      CALL SHAPE_QUADRANGLE_8(NODE,X,Y,IEL,VN,VDN)

	AH(1,1)=1.5*(AA(1)*VN(5)-AA(4)*VN(8))
	AH(1,2)=AB(1)*VN(5)+AB(4)*VN(8)
	AH(1,3)=VN(1)-ACK(1)*VN(5)-ACK(4)*VN(8)
	AH(1,4)=1.5*(AA(2)*VN(6)-AA(1)*VN(5))
	AH(1,5)=AB(2)*VN(6)+AB(1)*VN(5)
	AH(1,6)=VN(2)-ACK(2)*VN(6)-ACK(1)*VN(5)
	AH(1,7)=1.5*(AA(3)*VN(7)-AA(2)*VN(6))
	AH(1,8)=AB(3)*VN(7)+AB(2)*VN(6)
	AH(1,9)=VN(3)-ACK(3)*VN(7)-ACK(2)*VN(6)
	AH(1,10)=1.5*(AA(4)*VN(8)-AA(3)*VN(7))
	AH(1,11)=AB(4)*VN(8)+AB(3)*VN(7)
	AH(1,12)=VN(4)-ACK(4)*VN(8)-ACK(3)*VN(7)

	AH(2,1)=1.5*(AD(1)*VN(5)-AD(4)*VN(8))
	AH(2,2)=-VN(1)+AE(1)*VN(5)+AE(4)*VN(8)
	AH(2,3)=-AB(1)*VN(5)-AB(4)*VN(8)
	AH(2,4)=1.5*(AD(2)*VN(6)-AD(1)*VN(5))
	AH(2,5)=-VN(2)+AE(2)*VN(6)+AE(1)*VN(5)
	AH(2,6)=-AB(2)*VN(6)-AB(1)*VN(5)
	AH(2,7)=1.5*(AD(3)*VN(7)-AD(2)*VN(6))
	AH(2,8)=-VN(3)+AE(3)*VN(7)+AE(2)*VN(6)
	AH(2,9)=-AB(3)*VN(7)-AB(2)*VN(6)
	AH(2,10)=1.5*(AD(4)*VN(8)-AD(3)*VN(7))
	AH(2,11)=-VN(4)+AE(4)*VN(8)+AE(3)*VN(7)
	AH(2,12)=-AB(4)*VN(8)-AB(3)*VN(7)

	ADH(1,1)=1.5*(AA(1)*VDN(1,5)-AA(4)*VDN(1,8))
	ADH(1,2)=AB(1)*VDN(1,5)+AB(4)*VDN(1,8)
	ADH(1,3)=VDN(1,1)-ACK(1)*VDN(1,5)-ACK(4)*VDN(1,8)
	ADH(1,4)=1.5*(AA(2)*VDN(1,6)-AA(1)*VDN(1,5))
	ADH(1,5)=AB(2)*VDN(1,6)+AB(1)*VDN(1,5)
	ADH(1,6)=VDN(1,2)-ACK(2)*VDN(1,6)-ACK(1)*VDN(1,5)
	ADH(1,7)=1.5*(AA(3)*VDN(1,7)-AA(2)*VDN(1,6))
	ADH(1,8)=AB(3)*VDN(1,7)+AB(2)*VDN(1,6)
	ADH(1,9)=VDN(1,3)-ACK(3)*VDN(1,7)-ACK(2)*VDN(1,6)
	ADH(1,10)=1.5*(AA(4)*VDN(1,8)-AA(3)*VDN(1,7))
	ADH(1,11)=AB(4)*VDN(1,8)+AB(3)*VDN(1,7)
	ADH(1,12)=VDN(1,4)-ACK(4)*VDN(1,8)-ACK(3)*VDN(1,7)

	ADH(2,1)=1.5*(AA(1)*VDN(2,5)-AA(4)*VDN(2,8))
	ADH(2,2)=AB(1)*VDN(2,5)+AB(4)*VDN(2,8)
	ADH(2,3)=VDN(2,1)-ACK(1)*VDN(2,5)-ACK(4)*VDN(2,8)
	ADH(2,4)=1.5*(AA(2)*VDN(2,6)-AA(1)*VDN(2,5))
	ADH(2,5)=AB(2)*VDN(2,6)+AB(1)*VDN(2,5)
	ADH(2,6)=VDN(2,2)-ACK(2)*VDN(2,6)-ACK(1)*VDN(2,5)
	ADH(2,7)=1.5*(AA(3)*VDN(2,7)-AA(2)*VDN(2,6))
	ADH(2,8)=AB(3)*VDN(2,7)+AB(2)*VDN(2,6)
	ADH(2,9)=VDN(2,3)-ACK(3)*VDN(2,7)-ACK(2)*VDN(2,6)
	ADH(2,10)=1.5*(AA(4)*VDN(2,8)-AA(3)*VDN(2,7))
	ADH(2,11)=AB(4)*VDN(2,8)+AB(3)*VDN(2,7)
	ADH(2,12)=VDN(2,4)-ACK(4)*VDN(2,8)-ACK(3)*VDN(2,7)

	ADH(3,1)=1.5*(AD(1)*VDN(1,5)-AD(4)*VDN(1,8))
	ADH(3,2)=-VDN(1,1)+AE(1)*VDN(1,5)+AE(4)*VDN(1,8)
	ADH(3,3)=-AB(1)*VDN(1,5)-AB(4)*VDN(1,8)
	ADH(3,4)=1.5*(AD(2)*VDN(1,6)-AD(1)*VDN(1,5))
	ADH(3,5)=-VDN(1,2)+AE(2)*VDN(1,6)+AE(1)*VDN(1,5)
	ADH(3,6)=-AB(2)*VDN(1,6)-AB(1)*VDN(1,5)
	ADH(3,7)=1.5*(AD(3)*VDN(1,7)-AD(2)*VDN(1,6))
	ADH(3,8)=-VDN(1,3)+AE(3)*VDN(1,7)+AE(2)*VDN(1,6)
	ADH(3,9)=-AB(3)*VDN(1,7)-AB(2)*VDN(1,6)
	ADH(3,10)=1.5*(AD(4)*VDN(1,8)-AD(3)*VDN(1,7))
	ADH(3,11)=-VDN(1,4)+AE(4)*VDN(1,8)+AE(3)*VDN(1,7)
	ADH(3,12)=-AB(4)*VDN(1,8)-AB(3)*VDN(1,7)

	ADH(4,1)=1.5*(AD(1)*VDN(2,5)-AD(4)*VDN(2,8))
	ADH(4,2)=-VDN(2,1)+AE(1)*VDN(2,5)+AE(4)*VDN(2,8)
	ADH(4,3)=-AB(1)*VDN(2,5)-AB(4)*VDN(2,8)
	ADH(4,4)=1.5*(AD(2)*VDN(2,6)-AD(1)*VDN(2,5))
	ADH(4,5)=-VDN(2,2)+AE(2)*VDN(2,6)+AE(1)*VDN(2,5)
	ADH(4,6)=-AB(2)*VDN(2,6)-AB(1)*VDN(2,5)
	ADH(4,7)=1.5*(AD(3)*VDN(2,7)-AD(2)*VDN(2,6))
	ADH(4,8)=-VDN(2,3)+AE(3)*VDN(2,7)+AE(2)*VDN(2,6)
	ADH(4,9)=-AB(3)*VDN(2,7)-AB(2)*VDN(2,6)
	ADH(4,10)=1.5*(AD(4)*VDN(2,8)-AD(3)*VDN(2,7))
	ADH(4,11)=-VDN(2,4)+AE(4)*VDN(2,8)+AE(3)*VDN(2,7)
	ADH(4,12)=-AB(4)*VDN(2,8)-AB(3)*VDN(2,7)
	
	B(:,:)=0.0
	B(1,:)=AJJ(1,1)*ADH(1,:)+AJJ(1,2)*ADH(2,:)
	B(2,:)=AJJ(2,1)*ADH(3,:)+AJJ(2,2)*ADH(4,:)
	B(3,:)=AJJ(1,1)*ADH(3,:)+AJJ(1,2)*ADH(4,:)+
     $       AJJ(2,1)*ADH(1,:)+AJJ(2,2)*ADH(2,:)


	RETURN
	END