

      SUBROUTINE REORDER_SOLID_NODE_NUM(SXY,VNODE,NODE)
      IMPLICIT   REAL * 8(A-H,O-Z)
      COMMON/NUMOFF/NODENUM     !节点总数
      DIMENSION  SXY(8,3),NODE(8),NODEE(8),XY(4,3),VX(3,2),
     $           VXX(3),XYY(4,3),VY(3,2),VYY(3),NN(7),VNODE(NODENUM,3)
      
	DO 10 I=1,3
	SELECT CASE(I)
      CASE(1)
      XY(1,:)=SXY(2,:)
	XY(2,:)=SXY(1,:)
	XY(3,:)=SXY(4,:)
	XY(4,:)=SXY(3,:)
	VX(1,1)=XY(2,1)-XY(1,1)
	VX(2,1)=XY(2,2)-XY(1,2)
	VX(3,1)=XY(2,3)-XY(1,3)
	VX(1,2)=XY(3,1)-XY(2,1)
	VX(2,2)=XY(3,2)-XY(2,2)
	VX(3,2)=XY(3,3)-XY(2,3)
      
	VXX(1)=VX(2,1)*VX(3,2)-VX(3,1)*VX(2,2)
	VXX(2)=VX(3,1)*VX(1,2)-VX(1,1)*VX(3,2)
	VXX(3)=VX(1,1)*VX(2,2)-VX(2,1)*VX(1,2)
      
	XYY(1,:)=SXY(5,:)
	XYY(2,:)=SXY(6,:)
	XYY(3,:)=SXY(7,:)
	XYY(4,:)=SXY(8,:)
	VY(1,1)=XYY(2,1)-XYY(1,1)
	VY(2,1)=XYY(2,2)-XYY(1,2)
	VY(3,1)=XYY(2,3)-XYY(1,3)
	VY(1,2)=XYY(3,1)-XYY(2,1)
	VY(2,2)=XYY(3,2)-XYY(2,2)
	VY(3,2)=XYY(3,3)-XYY(2,3)
	VYY(1)=VY(2,1)*VY(3,2)-VY(3,1)*VY(2,2)
	VYY(2)=VY(3,1)*VY(1,2)-VY(1,1)*VY(3,2)
	VYY(3)=VY(1,1)*VY(2,2)-VY(2,1)*VY(1,2)


	COS=(VXX(1)*VYY(1)+VXX(2)*VYY(2)+VXX(3)*VYY(3))/
     $    (SQRT(VXX(1)*VXX(1)+VXX(2)*VXX(2)+VXX(3)*VXX(3))*
     $     SQRT(VYY(1)*VYY(1)+VYY(2)*VYY(2)+VYY(3)*VYY(3)))

	IF(ABS(ABS(COS)-1.0).LE.1E-6)THEN
      NODEE(:)=NODE(:)
	NODE(1)=NODEE(2)
	NODE(2)=NODEE(1)
	NODE(3)=NODEE(4)
	NODE(4)=NODEE(3)
	NODE(5)=NODEE(6)
	NODE(6)=NODEE(5)
	NODE(7)=NODEE(8)
	NODE(8)=NODEE(7)

	SXY(1,:)=VNODE(NODE(1),:)
	SXY(2,:)=VNODE(NODE(2),:)
	SXY(3,:)=VNODE(NODE(3),:)
	SXY(4,:)=VNODE(NODE(4),:)
	SXY(5,:)=VNODE(NODE(5),:)
	SXY(6,:)=VNODE(NODE(6),:)
	SXY(7,:)=VNODE(NODE(7),:)
	SXY(8,:)=VNODE(NODE(8),:)
	EXIT
	ENDIF

	
	CASE(2)
	XY(1,:)=SXY(1,:)
	XY(2,:)=SXY(2,:)
	XY(3,:)=SXY(6,:)
	XY(4,:)=SXY(5,:)
	VX(1,1)=XY(2,1)-XY(1,1)
	VX(2,1)=XY(2,2)-XY(1,2)
	VX(3,1)=XY(2,3)-XY(1,3)
	VX(1,2)=XY(3,1)-XY(2,1)
	VX(2,2)=XY(3,2)-XY(2,2)
	VX(3,2)=XY(3,3)-XY(2,3)
      
	VXX(1)=VX(2,1)*VX(3,2)-VX(3,1)*VX(2,2)
	VXX(2)=VX(3,1)*VX(1,2)-VX(1,1)*VX(3,2)
	VXX(3)=VX(1,1)*VX(2,2)-VX(2,1)*VX(1,2)
      
	XYY(1,:)=SXY(3,:)
	XYY(2,:)=SXY(4,:)
	XYY(3,:)=SXY(8,:)
	XYY(4,:)=SXY(7,:)
		VY(1,1)=XYY(2,1)-XYY(1,1)
	VY(2,1)=XYY(2,2)-XYY(1,2)
	VY(3,1)=XYY(2,3)-XYY(1,3)
	VY(1,2)=XYY(3,1)-XYY(2,1)
	VY(2,2)=XYY(3,2)-XYY(2,2)
	VY(3,2)=XYY(3,3)-XYY(2,3)
	VYY(1)=VY(2,1)*VY(3,2)-VY(3,1)*VY(2,2)
	VYY(2)=VY(3,1)*VY(1,2)-VY(1,1)*VY(3,2)
	VYY(3)=VY(1,1)*VY(2,2)-VY(2,1)*VY(1,2)


	COS=(VXX(1)*VYY(1)+VXX(2)*VYY(2)+VXX(3)*VYY(3))/
     $    (SQRT(VXX(1)*VXX(1)+VXX(2)*VXX(2)+VXX(3)*VXX(3))*
     $     SQRT(VYY(1)*VYY(1)+VYY(2)*VYY(2)+VYY(3)*VYY(3)))

	IF(ABS(ABS(COS)-1.0).LE.1E-6)THEN
      NODEE(:)=NODE(:)
	NODE(1)=NODEE(1)
	NODE(2)=NODEE(2)
	NODE(3)=NODEE(6)
	NODE(4)=NODEE(5)
	NODE(5)=NODEE(4)
	NODE(6)=NODEE(3)
	NODE(7)=NODEE(7)
	NODE(8)=NODEE(8)


	SXY(1,:)=VNODE(NODE(1),:)
	SXY(2,:)=VNODE(NODE(2),:)
	SXY(3,:)=VNODE(NODE(3),:)
	SXY(4,:)=VNODE(NODE(4),:)
	SXY(5,:)=VNODE(NODE(5),:)
	SXY(6,:)=VNODE(NODE(6),:)
	SXY(7,:)=VNODE(NODE(7),:)
	SXY(8,:)=VNODE(NODE(8),:)
	EXIT
	ENDIF

      
	CASE(3)
	XY(1,:)=SXY(2,:)
	XY(2,:)=SXY(3,:)
	XY(3,:)=SXY(7,:)
	XY(4,:)=SXY(6,:)
	VX(1,1)=XY(2,1)-XY(1,1)
	VX(2,1)=XY(2,2)-XY(1,2)
	VX(3,1)=XY(2,3)-XY(1,3)
	VX(1,2)=XY(3,1)-XY(2,1)
	VX(2,2)=XY(3,2)-XY(2,2)
	VX(3,2)=XY(3,3)-XY(2,3)
      
	VXX(1)=VX(2,1)*VX(3,2)-VX(3,1)*VX(2,2)
	VXX(2)=VX(3,1)*VX(1,2)-VX(1,1)*VX(3,2)
	VXX(3)=VX(1,1)*VX(2,2)-VX(2,1)*VX(1,2)
      
	XYY(1,:)=SXY(4,:)
	XYY(2,:)=SXY(1,:)
	XYY(3,:)=SXY(5,:)
	XYY(4,:)=SXY(8,:)
      VY(1,1)=XYY(2,1)-XYY(1,1)
	VY(2,1)=XYY(2,2)-XYY(1,2)
	VY(3,1)=XYY(2,3)-XYY(1,3)
	VY(1,2)=XYY(3,1)-XYY(2,1)
	VY(2,2)=XYY(3,2)-XYY(2,2)
	VY(3,2)=XYY(3,3)-XYY(2,3)
	VYY(1)=VY(2,1)*VY(3,2)-VY(3,1)*VY(2,2)
	VYY(2)=VY(3,1)*VY(1,2)-VY(1,1)*VY(3,2)
	VYY(3)=VY(1,1)*VY(2,2)-VY(2,1)*VY(1,2)


	COS=(VXX(1)*VYY(1)+VXX(2)*VYY(2)+VXX(3)*VYY(3))/
     $    (SQRT(VXX(1)*VXX(1)+VXX(2)*VXX(2)+VXX(3)*VXX(3))*
     $     SQRT(VYY(1)*VYY(1)+VYY(2)*VYY(2)+VYY(3)*VYY(3)))

	IF(ABS(ABS(COS)-1.0).LE.1E-6)THEN
      NODEE(:)=NODE(:)
	NODE(1)=NODEE(2)
	NODE(2)=NODEE(3)
	NODE(3)=NODEE(7)
	NODE(4)=NODEE(6)
	NODE(5)=NODEE(1)
	NODE(6)=NODEE(4)
	NODE(7)=NODEE(8)
	NODE(8)=NODEE(5)

	SXY(1,:)=VNODE(NODE(1),:)
	SXY(2,:)=VNODE(NODE(2),:)
	SXY(3,:)=VNODE(NODE(3),:)
	SXY(4,:)=VNODE(NODE(4),:)
	SXY(5,:)=VNODE(NODE(5),:)
	SXY(6,:)=VNODE(NODE(6),:)
	SXY(7,:)=VNODE(NODE(7),:)
	SXY(8,:)=VNODE(NODE(8),:)
	EXIT
	ENDIF
	ENDSELECT
10    CONTINUE

	RETURN
	END