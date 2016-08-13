      SUBROUTINE CHECK_SHAPE(MR,NUM,ISO)
      IMPLICIT   REAL * 8(A-H,O-Z)
	DIMENSION MR(8),MRR(4)

	SELECT CASE(NUM)
	CASE(1)
 	MRR(1)=MR(2)
	MRR(2)=MR(1)
	MRR(3)=MR(4)
	MRR(4)=MR(3)
	CASE(2)
    	MRR(1)=MR(1)
	MRR(2)=MR(2)
	MRR(3)=MR(6)
	MRR(4)=MR(5)
	CASE(3)
     	MRR(1)=MR(2)
	MRR(2)=MR(3)
	MRR(3)=MR(7)
	MRR(4)=MR(6)
	CASE(4)
    	MRR(1)=MR(3)
	MRR(2)=MR(4)
	MRR(3)=MR(8)
	MRR(4)=MR(7)
	CASE(5)
     	MRR(1)=MR(4)
	MRR(2)=MR(1)
	MRR(3)=MR(5)
	MRR(4)=MR(8)
	CASE(6)
    	MRR(1)=MR(5)
	MRR(2)=MR(6)
	MRR(3)=MR(7)
	MRR(4)=MR(8)
	ENDSELECT

	IF(MRR(1).EQ.MRR(2).OR.MRR(1).EQ.MRR(3).OR.MRR(1).EQ.MRR(4).OR.
     $   MRR(2).EQ.MRR(3).OR.MRR(2).EQ.MRR(4).OR.MRR(3).EQ.MRR(4))THEN
      ISO=3
	ELSE
      ISO=4
	ENDIF



	RETURN
	END