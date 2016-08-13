 
 
      SUBROUTINE DKT_I(VX,Y,Z,AK,VMAE,SA)
	IMPLICIT REAL * 8(A-H,O-Z)
               COMMON/COM1/MND, NUMEL, NUMPT,MBAND, NMATI
               COMMON/COM2/NF, NFSTR, MSOLV, MPROB, MTYPE, NVA
               COMMON/COMN/NFIX, NPC, GRAV
               COMMON/COM3/MND2,NUMPT2
C                   ! MND--MAXIMAL NODE NUMBER IN ALL ELEMENTS
C                   ! NUMEL--NUMBER OF GLOBAL ELEMENTS
C                   ! NUMPT--NUMBER OF GLOBAL NODES
C                   ! MBAND--HALF BANDWIDTH(INCLUDING DIAGONAL ELEMENT)
	DIMENSION VX(3,3),AX(3),AY(3),AL(3),
     $          AA(3),AB(3),AC(3),AD(3),AE(3),VN(6),VDN(2,6),
     $          AH(2,9),ADH(4,9),B(3,9),SB(9,3),VMAE(4),VD(5,5),AK(9,9)

                                !    (4)I--------------I(3)
	VN(1)=2*(1-Y-Z)*(0.5-Y-Z) !       I              I 
	VN(2)=Y*(2*Y-1)           !       I              I
	VN(3)=Z*(2*Z-1)           !       I              I
	VN(4)=4*Y*Z				  !       I              I
	VN(5)=4*Z*(1-Y-Z)		  !       I              I
	VN(6)=4*Y*(1-Y-Z)         !    (1)I--------------I(2) 六节点单元形函数
                                !
      VDN(1,1)=-3+4*Y+4*Z
	VDN(1,2)=4*Y-1
	VDN(1,3)=0
	VDN(1,4)=4*Z
	VDN(1,5)=-4*Z
	VDN(1,6)=4-8*Y-4*Z

	VDN(2,1)=-3+4*Y+4*Z
	VDN(2,2)=0
	VDN(2,3)=4*Z-1
	VDN(2,4)=4*Y
	VDN(2,5)=4-4*Y-8*Z
	VDN(2,6)=-4*Y

	AX(1)=VX(2,1)-VX(3,1)   !
	AX(2)=VX(3,1)-VX(1,1)   !
	AX(3)=VX(1,1)-VX(2,1)   !               /I(3)
	                        !              / I
							!             /  I
	AY(1)=VX(2,2)-VX(3,2)   !         (5)/   I(4)
	AY(2)=VX(3,2)-VX(1,2)   !           /    I
	AY(3)=VX(1,2)-VX(2,2)   !          /     I
	                        !      (1)/------I(2)

	AL(1)=AX(1)*AX(1)+AY(1)*AY(1)
	AL(2)=AX(2)*AX(2)+AY(2)*AY(2)
	AL(3)=AX(3)*AX(3)+AY(3)*AY(3)

	AA(1)=-AX(1)/AL(1)
	AA(2)=-AX(2)/AL(2)
	AA(3)=-AX(3)/AL(3)

	AB(1)=(0.75*AX(1)*AY(1))/AL(1)
      AB(2)=(0.75*AX(2)*AY(2))/AL(2)
      AB(3)=(0.75*AX(3)*AY(3))/AL(3)

	AC(1)=(0.25*AX(1)*AX(1)-0.5*AY(1)*AY(1))/AL(1)
	AC(2)=(0.25*AX(2)*AX(2)-0.5*AY(2)*AY(2))/AL(2)
	AC(3)=(0.25*AX(3)*AX(3)-0.5*AY(3)*AY(3))/AL(3)

	AD(1)=-AY(1)/AL(1)
	AD(2)=-AY(2)/AL(2)
	AD(3)=-AY(3)/AL(3)

	AE(1)=(0.25*AY(1)*AY(1)-0.5*AX(1)*AX(1))/AL(1)
	AE(2)=(0.25*AY(2)*AY(2)-0.5*AX(2)*AX(2))/AL(2)
	AE(3)=(0.25*AY(3)*AY(3)-0.5*AX(3)*AX(3))/AL(3)

!	AH(1,1)=1.5*(AA(3)*VN(6)-AA(2)*VN(5))
!	AH(1,2)=AB(2)*VN(5)+AB(3)*VN(6)
!	AH(1,3)=VN(1)-AC(2)*VN(5)-AC(3)*VN(6)
!	AH(1,4)=1.5*(AA(1)*VN(4)-AA(3)*VN(6))
!	AH(1,5)=AB(3)*VN(6)+AB(1)*VN(4)
!	AH(1,6)=VN(2)-AC(3)*VN(6)-AC(1)*VN(4)
!	AH(1,7)=1.5*(AA(2)*VN(5)-AA(1)*VN(4))
!	AH(1,8)=AB(1)*VN(4)+AB(2)*VN(5)
!	AH(1,9)=VN(3)-AC(1)*VN(4)-AC(2)*VN(5)

!	AH(2,1)=1.5*(AD(3)*VN(6)-AD(2)*VN(5))
!	AH(2,2)=-VN(1)+AE(2)*VN(5)+AE(3)*VN(6)
!	AH(2,3)=-AB(2)*VN(5)-AB(3)*VN(6)
!	AH(2,4)=1.5*(AD(1)*VN(4)-AD(3)*VN(6))
!	AH(2,5)=-VN(2)+AE(3)*VN(6)+AE(1)*VN(4)
!	AH(2,6)=-AB(3)*VN(6)-AB(1)*VN(4)
!	AH(2,7)=1.5*(AD(2)*VN(5)-AD(1)*VN(4))
!	AH(2,8)=-VN(3)+AE(1)*VN(4)+AE(2)*VN(5)
!	AH(2,9)=-AB(1)*VN(4)-AB(2)*VN(5)

	ADH(1,1)=1.5*(AA(3)*VDN(1,6)-AA(2)*VDN(1,5))
	ADH(1,2)=AB(2)*VDN(1,5)+AB(3)*VDN(1,6)
	ADH(1,3)=VDN(1,1)-AC(2)*VDN(1,5)-AC(3)*VDN(1,6)
	ADH(1,4)=1.5*(AA(1)*VDN(1,4)-AA(3)*VDN(1,6))
	ADH(1,5)=AB(3)*VDN(1,6)+AB(1)*VDN(1,4)
	ADH(1,6)=VDN(1,2)-AC(3)*VDN(1,6)-AC(1)*VDN(1,4)
	ADH(1,7)=1.5*(AA(2)*VDN(1,5)-AA(1)*VDN(1,4))
	ADH(1,8)=AB(1)*VDN(1,4)+AB(2)*VDN(1,5)
	ADH(1,9)=VDN(1,3)-AC(1)*VDN(1,4)-AC(2)*VDN(1,5)

	ADH(2,1)=1.5*(AA(3)*VDN(2,6)-AA(2)*VDN(2,5))
	ADH(2,2)=AB(2)*VDN(2,5)+AB(3)*VDN(2,6)
	ADH(2,3)=VDN(2,1)-AC(2)*VDN(2,5)-AC(3)*VDN(2,6)
	ADH(2,4)=1.5*(AA(1)*VDN(2,4)-AA(3)*VDN(2,6))
	ADH(2,5)=AB(3)*VDN(2,6)+AB(1)*VDN(2,4)
	ADH(2,6)=VDN(2,2)-AC(3)*VDN(2,6)-AC(1)*VDN(2,4)
	ADH(2,7)=1.5*(AA(2)*VDN(2,5)-AA(1)*VDN(2,4))
	ADH(2,8)=AB(1)*VDN(2,4)+AB(2)*VDN(2,5)
	ADH(2,9)=VDN(2,3)-AC(1)*VDN(2,4)-AC(2)*VDN(2,5)

	ADH(3,1)=1.5*(AD(3)*VDN(1,6)-AD(2)*VDN(1,5))
	ADH(3,2)=-VDN(1,1)+AE(2)*VDN(1,5)+AE(3)*VDN(1,6)
	ADH(3,3)=-AB(2)*VDN(1,5)-AB(3)*VDN(1,6)
	ADH(3,4)=1.5*(AD(1)*VDN(1,4)-AD(3)*VDN(1,6))
	ADH(3,5)=-VDN(1,2)+AE(3)*VDN(1,6)+AE(1)*VDN(1,4)
	ADH(3,6)=-AB(3)*VDN(1,6)-AB(1)*VDN(1,4)
	ADH(3,7)=1.5*(AD(2)*VDN(1,5)-AD(1)*VDN(1,4))
	ADH(3,8)=-VDN(1,3)+AE(1)*VDN(1,4)+AE(2)*VDN(1,5)
	ADH(3,9)=-AB(1)*VDN(1,4)-AB(2)*VDN(1,5)

      ADH(4,1)=1.5*(AD(3)*VDN(2,6)-AD(2)*VDN(2,5))
	ADH(4,2)=-VDN(2,1)+AE(2)*VDN(2,5)+AE(3)*VDN(2,6)
	ADH(4,3)=-AB(2)*VDN(2,5)-AB(3)*VDN(2,6)
	ADH(4,4)=1.5*(AD(1)*VDN(2,4)-AD(3)*VDN(2,6))
	ADH(4,5)=-VDN(2,2)+AE(3)*VDN(2,6)+AE(1)*VDN(2,4)
	ADH(4,6)=-AB(3)*VDN(2,6)-AB(1)*VDN(2,4)
	ADH(4,7)=1.5*(AD(2)*VDN(2,5)-AD(1)*VDN(2,4))
	ADH(4,8)=-VDN(2,3)+AE(1)*VDN(2,4)+AE(2)*VDN(2,5)
	ADH(4,9)=-AB(1)*VDN(2,4)-AB(2)*VDN(2,5)

      B(:,:)=0.0
	B(1,:)=AY(2)*ADH(1,:)+AY(3)*ADH(2,:)
	B(2,:)=-AX(2)*ADH(3,:)-AX(3)*ADH(4,:)
	B(3,:)=-AX(2)*ADH(1,:)-AX(3)*ADH(2,:)
     $       +AY(2)*ADH(3,:)+AY(3)*ADH(4,:)

	SA=(-AY(1)*AX(2)+AY(2)*AX(1))
	B(:,:)=B(:,:)/ABS(SA)
	CALL BEND_VD(VD,VMAE)
	SB(:,:)=0.0
	DO 10 I=1,9
	DO 10 J=1,3
	DO 10 KK=1,3
      SB(I,J)=SB(I,J)+B(KK,I)*VD(KK,J)
10    CONTINUE
      
	AK(:,:)=0.0
      DO 11 I=1,9
	DO 11 J=1,9
	DO 11 KK=1,3
      AK(I,J)=AK(I,J)+SB(I,KK)*B(KK,J)
11    CONTINUE

	RETURN
	END
