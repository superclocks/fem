
      SUBROUTINE ASSEM(VNODE,IELEM,IELES,IELEF,VMAE,FMAE,NN,SKK,SMM,
     $                 THIC,N0,MARK,M,N,NDIM) 
      IMPLICIT   REAL * 8(A-H,O-Z)
	COMMON/COM/NNN
      DIMENSION  NN(7),VMAE(4),FMAE(2),VNODE(NN(1),3),IELEM(NN(3),4),
     $           IELES(NN(4),8),FXY(8,3),FKK(8,8),FMM(8,8),THIC(NN(6)),
     $           VKK1(24,24),VMM1(24,24),IEL(4),VXY(4,3),VX(3,3),
     $           TCOS(3,3),QCOS(3,3),SKK(NDIM,NDIM),MARK(NNN,2),
     $           SMM(NDIM,NDIM),AK(24,24),AM(24,24),N0(NN(2)),
     $           TMM2(18,18),TKK2(18,18),TQ(9,3),TQQ(3,9),FN(3),
     $           FSKK(24,24),QQ(12,4),QQQ(4,12),SN(3),MR(8)

      DIMENSION  SXY(8,3),SCOS(3,3),GSK(24,24),GSM(24,24),NODE(8),
     $          IELEF(NN(5),8)

	SKK(:,:)=0.0
	SMM(:,:)=0.0
C====================������弰��Ͼ������ӵ����������======================
C======================���������ο��嵥Ԫ������Ͼ���==================
      DO 10 II=1,NN(3)
	VMAE(4)=THIC(N0(II))       !��ÿ����Ԫʵ������ֵ
      IF(IELEM(II,3).EQ.IELEM(II,4))THEN
	VX(1,:)=VNODE(IELEM(II,1),:)
	VX(2,:)=VNODE(IELEM(II,2),:)
	VX(3,:)=VNODE(IELEM(II,3),:)  !��ÿ�����嵥Ԫ�ڵ㸳����ֵ
	CALL TRI_TRANS(VX,TCOS)  !����������任

C-----------------------------------------------------------------------	
C---------------------------���������οǵ�Ԫ------------------------------
      CALL ELEMENT_MATRIX_TRI(VX,TCOS,TKK2,TMM2,VMAE) !����ṹ��Ԫ����
C-----------------------------------------------------------------------
C-----------------------�ѵ�Ԫ�����������նȾ���----------------------
      DO 11 I=1,3
	DO 11 J=1,3

	I2=(IELEM(II,I)-1)*6+1  !
      J2=(IELEM(II,J)-1)*6+1
      I3=(I-1)*6+1
	J3=(J-1)*6+1
	SKK(I2:I2+5,J2:J2+5)=SKK(I2:I2+5,J2:J2+5)+TKK2(I3:I3+5,J3:J3+5)
      SMM(I2:I2+5,J2:J2+5)=SMM(I2:I2+5,J2:J2+5)+TMM2(I3:I3+5,J3:J3+5)
11    CONTINUE
C----------------------------------------------------------------------   
C======================================================================
C======================�����ı��ο��嵥Ԫ����==========================
	ELSE
      VXY(1,:)=VNODE(IELEM(II,1),:)
	VXY(2,:)=VNODE(IELEM(II,2),:)
	VXY(3,:)=VNODE(IELEM(II,3),:)
	VXY(4,:)=VNODE(IELEM(II,4),:)   !��ÿ�����嵥Ԫ�ڵ㸳����ֵ

	CALL  QUAD_TRANS(VXY,QCOS)

C----------------------------------------------------------------------- 
C--------------------------�����ı��οǵ�Ԫ-----------------------------   
	CALL  ELEMENT_MATRIX_QUAD(VXY,QCOS,IEL,VMM1,VKK1,VMAE)
C-----------------------------------------------------------------------
C------------------�ѵ�Ԫ�նȾ�����뵽����նȾ�����-------------------
	DO 12 I=1,4
	DO 12 J=1,4
      I2=(IELEM(II,I)-1)*6+1  !
      J2=(IELEM(II,J)-1)*6+1
      I3=(I-1)*6+1
	J3=(J-1)*6+1    
	SKK(I2:I2+5,J2:J2+5)=SKK(I2:I2+5,J2:J2+5)+VKK1(I3:I3+5,J3:J3+5)
      SMM(I2:I2+5,J2:J2+5)=SMM(I2:I2+5,J2:J2+5)+VMM1(I3:I3+5,J3:J3+5)            
12    CONTINUE
C-----------------------------------------------------------------------
      ENDIF
10    CONTINUE
C===============================================================================
C=========================����ʵ�嵥Ԫ�������ӵ����������====================
      DO 15 II=1,NN(4)
      
      DO 16 I=1,8
      SXY(I,:)=VNODE(IELES(II,I),:)
16    CONTINUE
      NODE(:)=IELES(II,:)

!-----------�����嵥Ԫ�ڵ������ӳ���������Ԫ��ʽ��ִ�д��ӳ���--------
	IF(NODE(3).NE.NODE(4).AND.NODE(7).NE.NODE(8))THEN
      CALL REORDER_SOLID_NODE_NUM(SXY,VNODE,NODE)   !ʵ�嵥Ԫ�����ӳ���
      ENDIF
!-----------------------------------------------------------------------
      CALL ELEMENT_MATRIX_SOLID(SXY,SCOS,GSK,GSM,VMAE)

      DO 17 I=1,8
	DO 17 J=1,8
      IF(NODE(I).LE.M.AND.NODE(J).LE.M)THEN
	I2=(NODE(I)-1)*6+1     !
      J2=(NODE(J)-1)*6+1
	I3=(I-1)*3+1
	J3=(J-1)*3+1
	SKK(I2:I2+2,J2:J2+2)=SKK(I2:I2+2,J2:J2+2)+GSK(I3:I3+2,J3:J3+2)
      SMM(I2:I2+2,J2:J2+2)=SMM(I2:I2+2,J2:J2+2)+GSM(I3:I3+2,J3:J3+2)

      ELSEIF(NODE(I).GT.M.AND.NODE(J).LE.M)THEN
      I2=M*6+(NODE(I)-M-1)*3+1    !
      J2=(NODE(J)-1)*6+1
	I3=(I-1)*3+1
	J3=(J-1)*3+1
	SKK(I2:I2+2,J2:J2+2)=SKK(I2:I2+2,J2:J2+2)+GSK(I3:I3+2,J3:J3+2)
      SMM(I2:I2+2,J2:J2+2)=SMM(I2:I2+2,J2:J2+2)+GSM(I3:I3+2,J3:J3+2)

	ELSEIF(NODE(I).LE.M.AND.NODE(J).GT.M)THEN
      I2=(NODE(I)-1)*6+1    !
      J2=M*6+(NODE(J)-M-1)*3+1
	I3=(I-1)*3+1
	J3=(J-1)*3+1
	SKK(I2:I2+2,J2:J2+2)=SKK(I2:I2+2,J2:J2+2)+GSK(I3:I3+2,J3:J3+2)
      SMM(I2:I2+2,J2:J2+2)=SMM(I2:I2+2,J2:J2+2)+GSM(I3:I3+2,J3:J3+2)

	ELSE
      I2=M*6+(NODE(I)-M-1)*3+1    !
      J2=M*6+(NODE(J)-M-1)*3+1
	I3=(I-1)*3+1
	J3=(J-1)*3+1
	SKK(I2:I2+2,J2:J2+2)=SKK(I2:I2+2,J2:J2+2)+GSK(I3:I3+2,J3:J3+2)
      SMM(I2:I2+2,J2:J2+2)=SMM(I2:I2+2,J2:J2+2)+GSM(I3:I3+2,J3:J3+2)

	ENDIF

17    CONTINUE

15    CONTINUE
C===============================================================================
C======================�������嵥Ԫ�ն�,�������󲢼��뵽���������==============
      IF(N.GE.M)THEN
      IDIM=M*6+(N-M)*3
	ELSE
      IDIM=M*6
	ENDIF

      DO 20 I=1,NN(5)
	DO 21 II=1,8
      FXY(II,:)=VNODE(IELEF(I,II),:)        !�����嵥Ԫÿ���ڵ㸳ֵ
21    CONTINUE
      NODE(:)=IELEF(I,:)
	IF(NODE(3).NE.NODE(4).AND.NODE(7).NE.NODE(8))THEN
      CALL REORDER_SOLID_NODE_NUM(FXY,VNODE,NODE)
      ENDIF
      CALL FLUID30_8_II(FXY,FKK,FMM)
	DO 22 II=1,8
      DO 22 JJ=1,8
       
	SKK(IDIM+NODE(II),IDIM+NODE(JJ))=
     $     SKK(IDIM+NODE(II),IDIM+NODE(JJ))+FKK(II,JJ) 
      SMM(IDIM+NODE(II),IDIM+NODE(JJ))=
     $     SMM(IDIM+NODE(II),IDIM+NODE(JJ))+FMM(II,JJ) 
22    CONTINUE
20    CONTINUE
C===============================================================================
C============================������Ͼ��󲢼��뵽���������=====================
      DO 30 I=1,NNN

      DO 31 II=1,8
      FXY(II,:)=VNODE(IELEF(MARK(I,1)-NN(3)-NN(4),II),:)
	                                 !�����嵥Ԫÿ���ڵ㸳ֵ
31    CONTINUE

      MR(:)=IELEF(MARK(I,1)-NN(3)-NN(4),:)   !������������йص���ڵ����ҳ�
      NUM=MARK(I,2)              !�������

      CALL CHECK_SHAPE(MR,NUM,ISO)  !�ж��������״
C------------------------�����������������Ͼ��󲢼������������--------------- 
      IF(ISO.EQ.3)THEN
      CALL CACUL_OUT_NORMAL_TRI(FXY,FN,NUM,MR,VX)
      CALL TRI_COUP_II(VX,FN,TQ)

      TQ(:,:)=(-TQ(:,:))
	TQQ(:,:)=-FMAE(1)*TRANSPOSE(TQ(:,:)) !���������е���Ͼ���
	DO 33 I1=1,3
	DO 33 J1=1,3

	IF(MR(I1).LE.M)THEN

      I2=(MR(I1)-1)*6+1  !
	IF(N.GE.M)THEN
      J2=MR(J1)+M*6+(N-M)*3
	ELSE
      J2=MR(J1)+M*6
	ENDIF 
      SKK(I2:I2+2,J2)=SKK(I2:I2+2,J2)+TQ((I1-1)*3+1:(I1-1)*3+3,J1)
      SMM(J2,I2:I2+2)=SMM(J2,I2:I2+2)+TQQ(I1,(J1-1)*3+1:(J1-1)*3+3)     
	ELSE
      I2=M*6+(MR(I1)-M-1)*3+1  !

      IF(N.GE.M)THEN
      J2=MR(J1)+M*6+(N-M)*3
	ELSE
      J2=MR(J1)+M*6
	ENDIF 
      SKK(I2:I2+2,J2)=SKK(I2:I2+2,J2)+TQ((I1-1)*3+1:(I1-1)*3+3,J1)
      SMM(J2,I2:I2+2)=SMM(J2,I2:I2+2)+TQQ(I1,(J1-1)*3+1:(J1-1)*3+3)     
	
	ENDIF          
      

!      I2=(MR(I1)-1)*7+1  !���������ʱ
!      J2=(MR(J1)-1)*7+1
      
!      SKK(I2:I2+2,J2+6)=SKK(I2:I2+2,J2+6)+TQ((I1-1)*3+1:(I1-1)*3+3,J1)

!      SMM(I2+6,J2:J2+2)=SMM(I2+6,J2:J2+2)+TQQ(I1,(J1-1)*3+1:(J1-1)*3+3)  
33    CONTINUE
C------------------------------------------------------------------------------
C-----------------------�����ı����������Ͼ��󲢼������������---------------
	ELSE
      CALL CACUL_OUT_NORMAL_QUAD(FXY,SN,NUM,MR,VXY)

     	CALL QUAD_COUP_II(VXY,SN,QQ)  !������Ͼ���,QQΪ�ı�����Ͼ���
      QQ(:,:)=(-QQ(:,:))
	QQQ(:,:)=-FMAE(1)*TRANSPOSE(QQ(:,:)) !���������е���Ͼ���	
	DO 32 I1=1,4
	DO 32 J1=1,4
	IF(MR(I1).LE.M)THEN

      I2=(MR(I1)-1)*6+1  !
	IF(N.GE.M)THEN
      J2=MR(J1)+M*6+(N-M)*3
	ELSE
      J2=MR(J1)+M*6
	ENDIF
!      SKK(I2:I2+2,J2)=SKK(I2:I2+2,J2)+QQ((I1-1)*3+1:(I1-1)*3+3,J1)
!      SMM(I2,J2:J2+2)=SMM(I2,J2:J2+2)+QQQ(I1,(J1-1)*3+1:(J1-1)*3+3)     
      SKK(I2:I2+2,J2)=SKK(I2:I2+2,J2)+QQ((I1-1)*3+1:(I1-1)*3+3,J1)
      SMM(J2,I2:I2+2)=SMM(J2,I2:I2+2)+QQQ(I1,(J1-1)*3+1:(J1-1)*3+3)     
	ELSE
      I2=M*6+(MR(I1)-M-1)*3+1  !

      IF(N.GE.M)THEN
      J2=MR(J1)+M*6+(N-M)*3
	ELSE
      J2=MR(J1)+M*6
	ENDIF
	     
!      SKK(I2:I2+2,J2)=SKK(I2:I2+2,J2)+QQ((I1-1)*3+1:(I1-1)*3+3,J1)
!      SMM(I2,J2:J2+2)=SMM(I2,J2:J2+2)+QQQ(I1,(J1-1)*3+1:(J1-1)*3+3)     
      SKK(I2:I2+2,J2)=SKK(I2:I2+2,J2)+QQ((I1-1)*3+1:(I1-1)*3+3,J1)
      SMM(J2,I2:I2+2)=SMM(J2,I2:I2+2)+QQQ(I1,(J1-1)*3+1:(J1-1)*3+3)     
	
	ENDIF          
32    CONTINUE
      ENDIF
C-------------------------------------------------------------------------------
30    CONTINUE
C===============================================================================

      OPEN(5,FILE='SKK.DAT',STATUS='UNKNOWN')
      OPEN(6,FILE='SMM.DAT',STATUS='UNKNOWN')
!	OPEN(8,FILE ='aSKK.DAT', STATUS='OLD')
!      OPEN(9,FILE ='aSmm.DAT', STATUS='old')
!      OPEN(7,FILE ='p.DAT', STATUS='unknown')
      
 
	WRITE(5,100)((SKK(I,J),J=1,NDIM),I=1,NDIM)
      WRITE(6,100)((SMM(I,J),J=1,NDIM),I=1,NDIM)
!      WRITE(7,101)(p(I),i=1,3360)
100   FORMAT(2716E18.8)
!101   FORMAT(E18.8)

      RETURN
      END