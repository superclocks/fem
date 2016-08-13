
      SUBROUTINE SOLID1_II(SXY,VMAE,GSK,GSM,KEYOPT2)
	IMPLICIT  REAL*8 (A-H,O-Z)
      DIMENSION   SXY(8,3),GSK(24,24),GSM(24,24),SK(24,24),
     $            SM(24,24),VMAE(4)


      GSK(:,:)=0.0
	GSM(:,:)=0.0
      KEYOPT2=0
      SELECT CASE(KEYOPT2)
	CASE(0)   !八点高斯积分
      DO 10 I=1,8
      SELECT CASE(I)
	CASE(1)
      X= 0.57735 02691 89626
	Y= 0.57735 02691 89626
	Z=-0.57735 02691 89626
      CALL SOLID1_I(SXY,X,Y,Z,VMAE,SK,SM)
      GSK(:,:)=GSK(:,:)+SK(:,:)
	GSM(:,:)=GSM(:,:)+SM(:,:)

      CASE(2)
      X=-0.57735 02691 89626
	Y= 0.57735 02691 89626
	Z=-0.57735 02691 89626
	CALL SOLID1_I(SXY,X,Y,Z,VMAE,SK,SM)
      GSK(:,:)=GSK(:,:)+SK(:,:)
	GSM(:,:)=GSM(:,:)+SM(:,:)

	CASE(3)
      X=-0.57735 02691 89626
	Y=-0.57735 02691 89626
	Z=-0.57735 02691 89626
	CALL SOLID1_I(SXY,X,Y,Z,VMAE,SK,SM)
      GSK(:,:)=GSK(:,:)+SK(:,:)
	GSM(:,:)=GSM(:,:)+SM(:,:)

	CASE(4)
	X= 0.57735 02691 89626
	Y=-0.57735 02691 89626
	Z=-0.57735 02691 89626
	CALL SOLID1_I(SXY,X,Y,Z,VMAE,SK,SM)
      GSK(:,:)=GSK(:,:)+SK(:,:)
	GSM(:,:)=GSM(:,:)+SM(:,:)

	CASE(5)
      X= 0.57735 02691 89626
	Y= 0.57735 02691 89626
	Z= 0.57735 02691 89626
	CALL SOLID1_I(SXY,X,Y,Z,VMAE,SK,SM)
      GSK(:,:)=GSK(:,:)+SK(:,:)
	GSM(:,:)=GSM(:,:)+SM(:,:)

	CASE(6)
      X=-0.57735 02691 89626
	Y= 0.57735 02691 89626
	Z= 0.57735 02691 89626
	CALL SOLID1_I(SXY,X,Y,Z,VMAE,SK,SM)
      GSK(:,:)=GSK(:,:)+SK(:,:)
	GSM(:,:)=GSM(:,:)+SM(:,:)

	CASE(7)
      X=-0.57735 02691 89626
	Y=-0.57735 02691 89626
	Z= 0.57735 02691 89626
	CALL SOLID1_I(SXY,X,Y,Z,VMAE,SK,SM)
      GSK(:,:)=GSK(:,:)+SK(:,:)
	GSM(:,:)=GSM(:,:)+SM(:,:)

	CASE(8)
      X= 0.57735 02691 89626
	Y=-0.57735 02691 89626
	Z= 0.57735 02691 89626
	CALL SOLID1_I(SXY,X,Y,Z,VMAE,SK,SM)
      GSK(:,:)=GSK(:,:)+SK(:,:)
	GSM(:,:)=GSM(:,:)+SM(:,:)

	END SELECT
10    CONTINUE

	CASE(1)   !一点高斯积分
      X= 0.0
	Y= 0.0
	Z= 0.0
	CALL SOLID1_I(SXY,X,Y,Z,VMAE,SK,SM)
      GSK(:,:)=GSK(:,:)+SK(:,:)*2.0
	GSM(:,:)=GSM(:,:)+SM(:,:)*2.0
      END SELECT
	END