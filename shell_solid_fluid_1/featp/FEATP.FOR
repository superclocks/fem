CCC
ccc
      PROGRAM  main

      IMPLICIT   REAL * 8(A-H ,O-Z)
	COMMON/KEYOPT3/KEYOPT_3
      COMMON/COM/NNN
	COMMON/NUMOFF/NODENUM
 !     COMMON/STR/VMAE(1),VMAE(2),VMAE(3)  !�ṹ������ģ�������ɱȣ��ܶ�
	PARAMETER(PI=3.1415926)                  
	COMMON  DENS,SONC
      DIMENSION  NN(8),VXY(4,3),FKK(8,8),FMM(8,8),FSQQQ(8,24),
     $         GKK(15,15),GMM(15,15),VMAE(4),FMAE(2),UX(50)
	REAL*8,POINTER::P(:,:)

      CHARACTER(LEN=80)::S,S1(100)
	ALLOCATABLE VNODE(:,:),IELEM(:,:),IELEF(:,:),IELES(:,:),SKK(:,:),
     $            SMM(:,:),THIC(:),N0(:),A(:,:),MAR(:,:),MARK(:,:)

	ALLOCATABLE MPTR(:),IND(:),VAL(:),RHS(:),NCP(:,:),NET(:,:)
C************************����Ԫ��Ϣ�ļ�*******************************
!	
!      CALL ASSEM_SOLID45()    !ʵ�嵥Ԫ���Խӿڳ���
!      CALL ASSEM_FLUID30()    !���嵥Ԫ���Խӿڳ���
	OPEN(200,  FILE= 'ux.dat', STATUS= 'UNKNOWN')
      OPEN(201,  FILE= 'uy.dat', STATUS= 'UNKNOWN')
      OPEN(202,  FILE= 'uz.dat', STATUS= 'UNKNOWN')
      OPEN(203,  FILE= 'rotx.dat', STATUS= 'UNKNOWN')
	OPEN(204, FILE= 'roty.dat', STATUS= 'UNKNOWN')
	OPEN(205, FILE= 'rotz.dat', STATUS= 'UNKNOWN')
      OPEN(206, FILE= 'pre.dat', STATUS= 'UNKNOWN')

      OPEN(100,FILE ='double_cylindrical2_5.dat', STATUS= 'OLD',
     $              IOSTAT=ERROR)
      KEYOPT_3=0   !���Ʋ��������ܶ�Ӧ��Ansys��keyopt3
     
      DO 10 I=1,8
      READ(100,*)S,S,NN(I) !NN(1)Ϊ�ڵ�����,NN(2)Ϊ��Ԫ����
	                     !NN(3)Ϊ���嵥Ԫ��,NN(4)Ϊʵ�嵥Ԫ��
	                     !NN(5)Ϊ���嵥Ԫ��,NN(6)Ϊ���ϲ�������    
						 !NN(7)Ϊʵ��������,NN(8)Ϊ��Ԫ��������
10    CONTINUE
      NODENUM=NN(1)
    
      ALLOCATE(MAR((NN(3)+NN(4))*2+2,2))   !��̬����������־����,��ʵ�ʵ�Ҫ��
      ALLOCATE(THIC(NN(7)))          !��̬�������ʵ��������
      ALLOCATE(VNODE(NN(1),3))       !��̬����ڵ���������ڴ�
	ALLOCATE(IELEM(NN(3),4))       !��̬����ǵ�Ԫ��ž����ڴ�
	ALLOCATE(IELES(NN(4),8))       !��̬����ʵ�嵥Ԫ��ž����ڴ�
      ALLOCATE(IELEF(NN(5),8))       !��̬�������嵥Ԫ��ž����ڴ�
      ALLOCATE(NET(NN(8),2))         !��̬�ֵ�Ԫ��������

      ALLOCATE(N0(NN(2)))            !��̬���乤�������ڴ�
	ALLOCATE(NCP(NN(7),3))
      THIC(:)=0.0
	VNODE(:,:)=0.0
      IELEM(:,:)=0.0
	IELES(:,:)=0.0
	IELEF(:,:)=0.0
	NET(:,:)=0.0
      N0(:)=0.0

      DO 11 II=1,NN(8)
      READ(100,*)S,NET(II,1),NET(II,2) !����Ԫ����
11    CONTINUE

      DO 12 II=1,NN(7)
      READ(100,*)I0,I0,THIC(II) !���ṹ���
12    CONTINUE
      
	DO 13 II=1,NN(1)
 	READ(100,*)I0,I1,I2,(VNODE(II,JJ),JJ=1,3) !VNODEΪ�ڵ��������
13    CONTINUE

      DO 14 II=1,NN(3)
	READ(100,*)I0,I0,N0(II),I0,I0,I0,I0,I0,I0,I0,I0,
     $       (IELEM(II,JJ),JJ=1,4)
14    CONTINUE                        !IELEMΪÿ���ǵ�Ԫ�ڵ����
      
	DO 15 II=1,NN(4)      
	READ(100,*)I0,I0,N0(II+NN(3)),I0,I0,I0,I0,I0,I0,I0,I0,
     $       (IELES(II,JJ),JJ=1,8)
15    CONTINUE	               !IELEFΪÿ�����嵥Ԫ�ڵ����

	DO 16 II=1,NN(5)      
	READ(100,*)I0,I0,N0(II+NN(3)+NN(4)),I0,I0,I0,I0,I0,I0,I0,I0,
     $       (IELEF(II,JJ),JJ=1,8)
16    CONTINUE	               !IELEFΪÿ�����嵥Ԫ�ڵ����
!      WRITE(*,*)N0(1537),N0(4744)                     

      READ(100,*)S,S,S,S,S,S,VMAE(1)  !���ṹ����ģ��
	READ(100,*)S,S,S,S,S,S,VMAE(2)  !���ṹ���ɱ�
      READ(100,*)S,S,S,S,S,S,VMAE(3)  !���ṹ�ܶ�
	READ(100,*)S,S,S,S,S,S,FMAE(1)  !�������ܶ�
      READ(100,*)S,S,S,S,S,S,FMAE(2)  !����������
      
!	DO 16 II=1,3
!      READ(100,*)S  !,MARK(II,1),MARK(II,2)
!17    CONTINUE
!      WRITE(*,*)MARK(737,1)
!      DO 16 II=1,NN(7)
!      READ(100,*)S,S,S,NCP(II,1),S,NCP(II,2),NCP(II,3)
!16    CONTINUE


      II=1
      DO WHILE(.TRUE.)
      READ(100,*,IOSTAT=ERROR)S,MAR(II,1),MAR(II,2)
	II=II+1
	IF(ERROR/=0)EXIT
	ENDDO
	NNN=II-2

	ALLOCATE(MARK(NNN,2))    !����ȷ����Ҫ��������йصı߽������ڴ�
	MARK(:,:)=MAR(1:NNN,:)
	DEALLOCATE(MAR)          !�ͷ������ڴ�
	DENS=FMAE(1)
	SONC=FMAE(2)

      M=MAXVAL(IELEM(:,:))         !�ǵ�Ԫ�����ڵ��
	N=MAXVAL(IELES(:,:))         !ʵ�嵥Ԫ�����ڵ��
      MAX=MAXVAL(IELEF(:,:))       !���嵥Ԫ�����ڵ��
      MIN=MINVAL(IELEF(:,:))       !���嵥Ԫ����С�ڵ��
      IF(N.GE.M)THEN
      NDIM=M*6+(N-M)*3+MAX-MIN+1   !�������ά��
	ELSE
      NDIM=M*6+MAX-MIN+1            !�������ά��
	ENDIF
	
     
	ALLOCATE(SKK(NDIM,NDIM)) !��̬��������նȾ����ڴ�
	                               
	ALLOCATE(SMM(NDIM,NDIM)) !��̬�����������������ڴ�
	                               
C************************************************************************

!======================�����������ӳ����ͷ������ڴ�==================
      START=SECNDS(0.0)
      CALL ASSEM(VNODE,IELEM,IELES,IELEF,VMAE,FMAE,NN,SKK,SMM,
     $          THIC,N0,MARK,M,N,NDIM)	
      ASSEM_TIME=SECNDS(START)
      WRITE(*,*)'assembly time',ASSEM_TIME
	DEALLOCATE(THIC)
      DEALLOCATE(VNODE)
	DEALLOCATE(IELEM)
      DEALLOCATE(IELES)
	DEALLOCATE(IELEF)  !
      DEALLOCATE(N0)
C************************************************************************

      NS=0 ! 10
	NST=0  !  !ʩ�ӱ߽������Ľڵ���
	IF(NS.NE.0.0.AND.NST.NE.0.0)THEN
	SKK((NS-1)*7+1:NST*7,:)=0.0
	SKK(:,(NS-1)*7+1:NST*7)=0.0 
	SMM((NS-1)*7+1:NST*7,:)=0.0
	SMM(:,(NS-1)*7+1:NST*7)=0.0 
      
	DO 24 I=(NS-1)*7+1,NST*7
	SKK(I,I)=1.0  !-AK(:,:)
	SMM(I,I)=1.0  !-AM(:,:)
24    CONTINUE
      ENDIF

	NS=0 !396
	NST=0 !396  !ʩ�ӱ߽������Ľڵ���
	IF(NS.NE.0.0.AND.NST.NE.0.0)THEN
	SKK((NS-1)*7+1:NST*7,:)=0.0
	SKK(:,(NS-1)*7+1:NST*7)=0.0 
	SMM((NS-1)*7+1:NST*7,:)=0.0
	SMM(:,(NS-1)*7+1:NST*7)=0.0 
      
	DO 25 I=(NS-1)*7+1,NST*7
	SKK(I,I)=1.0  !-AK(:,:)
	SMM(I,I)=1.0  !-AM(:,:)
25    CONTINUE
      ENDIF
!========================================================================

C*********************************���������*****************************
!+++++++++++++++++++++++++++++���������������+++++++++++++++++++++++++
      NUMBER=1       !г��Ӧ������
	FB=499.95         !��ʼƵ��
	FO=499.95         !����Ƶ��
 !     NDIM=M*6+N
	
	ALLOCATE(A(NDIM,NDIM))
	DO 60 II=1,NUMBER
      
	START=SECNDS(0.0)
      F=FB+((FO-FB)/NUMBER)*II
	WRITE(*,*)"substep",II,"begine "," frequency=",F
	A(:,:)=SKK(:,:)-4*PI*PI*F*F*SMM(:,:)
!	DEALLOCATE(SMM)
!      ALLOCATE(A(NDIM,NDIM))
!      A=SKK
      
      NNS=COUNT(A(:,:).NE.0)  !�������Ԫ�ظ���
	
      ALLOCATE(MPTR(NDIM+1)) !NPTR�ĳ�����N+1,PTR[I]��¼��i�е�һ������Ԫ
	             !��λ��,���һ��Ԫ��PTR[N��1]=NN+1��������I�еĳ��Ⱦ���PTR[I+1]-PTR[I]��
	ALLOCATE(IND(NNS))  !IND���м���ÿ������Ԫ��λ�� 
	ALLOCATE(VAL(NNS))  !VAL���м���ÿ������Ԫ��ֵ
	IF(II.EQ.1)THEN
	ALLOCATE(RHS(NDIM))  !�غ�����
      ENDIF
	
	K=1

      DO 50 J=1,NDIM
	DO 51 I=1,NDIM
	IF(A(I,J).NE.0.0)THEN
      IND(K)=I
	VAL(K)=A(I,J)      
	K=K+1
	ENDIF
51    CONTINUE

      IF(J.EQ.1)THEN
	MPTR(J)=1
	ELSE
	MPTR(J)=MPTR(J-1)+COUNT(A(:,J-1).NE.0.0)
	ENDIF
50    CONTINUE

      MPTR(NDIM+1)=NNS+1
     	RHS(:)=0.0
      RHS(1)=10.0
 !     RHS(3131)=1.0
!	RHS(15)=10.0
  !    RHS(20)=10.0
      PRE_MANAGE=SECNDS(START)
	WRITE(*,*)"pre_manage time",PRE_MANAGE
!++++++++++++++++++++++++++++������ӿ��ӳ���+++++++++++++++++++++++++++++++
      CALL SOLVER(NDIM,MPTR,IND,VAL,RHS,NNS,M_TYPE)	
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !    UX(II)=RHS(4477)
      
!	WRITE(200,600)RHS(541)
!	WRITE(201,600)RHS(542)
!	WRITE(202,600)RHS(543)
!	WRITE(203,600)RHS(544)
!	WRITE(204,600)RHS(545)
!	WRITE(205,600)RHS(546)
!      write(*,*)rhs(2950)

      DEALLOCATE(MPTR) !NPTR�ĳ�����N+1,PTR[I]��¼��i�е�һ������Ԫ
	             !��λ��,���һ��Ԫ��PTR[N��1]=NN+1��������I�еĳ��Ⱦ���PTR[I+1]-PTR[I]��
	DEALLOCATE(IND)  !IND���м���ÿ������Ԫ��λ�� 
	DEALLOCATE(VAL)  !VAL���м���ÿ������Ԫ��ֵ

60    CONTINUE
      
	DO 200 I=1,ndim  !NN(1)
	IF(I.LE.M)THEN
     	WRITE(200,600)RHS(1+(I-1)*6)
	WRITE(201,600)RHS(2+(I-1)*6)
	WRITE(202,600)RHS(3+(I-1)*6)
	WRITE(203,600)RHS(4+(I-1)*6)
	WRITE(204,600)RHS(5+(I-1)*6)
	WRITE(205,600)RHS(6+(I-1)*6)
	ELSEIF(I.GT.M.AND.I.LE.N)THEN
      WRITE(200,600)RHS(6*M+1+(I-M-1)*3)
	WRITE(201,600)RHS(6*M+2+(I-M-1)*3)
	WRITE(202,600)RHS(6*M+3+(I-M-1)*3)
	ENDIF
200   CONTINUE
      
      WRITE(206,600)(RHS(60+I),I=1,11)
      
600   FORMAT(E18.8)

      DEALLOCATE(SKK) !�ͷŸնȾ����ڴ�
	DEALLOCATE(SMM) !�ͷ����������ڴ�
      DEALLOCATE(A)   !�ͷŹ��������ڴ�
C***************************************************************************
C***************************************************************************
      STOP
      END