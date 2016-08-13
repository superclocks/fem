
      SUBROUTINE SOLID0_II(SXY,VMAE,GSK,GSM)
	USE IMSL
	IMPLICIT  REAL*8 (A-H,O-Z)
      DIMENSION   SXY(8,3),GSK(24,24),GSM(24,24),SK(33,33),
     $            SM(24,24),VMAE(4),GSKK(33,33),AK(9,9),AKK(24,9),
     $            AKKK(24,24)


      GSK(:,:)=0.0
	GSM(:,:)=0.0

      DO 10 I=1,8
      SELECT CASE(I)
	CASE(1)
      X= 0.57735 02691 89626
	Y= 0.57735 02691 89626
	Z=-0.57735 02691 89626
      CALL SOLID0_I(SXY,X,Y,Z,VMAE,SK,SM)
      GSKK(:,:)=GSKK(:,:)+SK(:,:)
	GSM(:,:)=GSM(:,:)+SM(:,:)

      CASE(2)
      X=-0.57735 02691 89626
	Y= 0.57735 02691 89626
	Z=-0.57735 02691 89626
	CALL SOLID0_I(SXY,X,Y,Z,VMAE,SK,SM)
      GSKK(:,:)=GSKK(:,:)+SK(:,:)
	GSM(:,:)=GSM(:,:)+SM(:,:)

	CASE(3)
      X=-0.57735 02691 89626
	Y=-0.57735 02691 89626
	Z=-0.57735 02691 89626
	CALL SOLID0_I(SXY,X,Y,Z,VMAE,SK,SM)
      GSKK(:,:)=GSKK(:,:)+SK(:,:)
	GSM(:,:)=GSM(:,:)+SM(:,:)

	CASE(4)
	X= 0.57735 02691 89626
	Y=-0.57735 02691 89626
	Z=-0.57735 02691 89626
	CALL SOLID0_I(SXY,X,Y,Z,VMAE,SK,SM)
      GSKK(:,:)=GSKK(:,:)+SK(:,:)
	GSM(:,:)=GSM(:,:)+SM(:,:)

	CASE(5)
      X= 0.57735 02691 89626
	Y= 0.57735 02691 89626
	Z= 0.57735 02691 89626
	CALL SOLID0_I(SXY,X,Y,Z,VMAE,SK,SM)
      GSKK(:,:)=GSKK(:,:)+SK(:,:)
	GSM(:,:)=GSM(:,:)+SM(:,:)

	CASE(6)
      X=-0.57735 02691 89626
	Y= 0.57735 02691 89626
	Z= 0.57735 02691 89626
	CALL SOLID0_I(SXY,X,Y,Z,VMAE,SK,SM)
      GSKK(:,:)=GSKK(:,:)+SK(:,:)
	GSM(:,:)=GSM(:,:)+SM(:,:)

	CASE(7)
      X=-0.57735 02691 89626
	Y=-0.57735 02691 89626
	Z= 0.57735 02691 89626
	CALL SOLID0_I(SXY,X,Y,Z,VMAE,SK,SM)
      GSKK(:,:)=GSKK(:,:)+SK(:,:)
	GSM(:,:)=GSM(:,:)+SM(:,:)

	CASE(8)
      X= 0.57735 02691 89626
	Y=-0.57735 02691 89626
	Z= 0.57735 02691 89626
	CALL SOLID0_I(SXY,X,Y,Z,VMAE,SK,SM)
      GSKK(:,:)=GSKK(:,:)+SK(:,:)
	GSM(:,:)=GSM(:,:)+SM(:,:)
	END SELECT
10    CONTINUE
      
	AK(:,:)=.I.GSKK(25:33,25:33)

      AKK(:,:)=0.0
	DO 11 I=1,24
	DO 11 J=1,9
	DO 11 K=1,9
      AKK(I,J)=AKK(I,J)+GSKK(I,K+24)*AK(K,J)
11    CONTINUE

      DO 12 I=1,24
	DO 12 J=1,24
	DO 12 K=1,9
      AKKK(I,J)=AKKK(I,J)+AKK(I,K)*GSKK(K+24,J)
12    CONTINUE

      GSK(:,:)=GSKK(1:24,1:24)-AKKK(:,:)



	END