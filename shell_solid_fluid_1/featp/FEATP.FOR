CCC
ccc
      PROGRAM  main

      IMPLICIT   REAL * 8(A-H ,O-Z)
	COMMON/KEYOPT3/KEYOPT_3
      COMMON/COM/NNN
	COMMON/NUMOFF/NODENUM
 !     COMMON/STR/VMAE(1),VMAE(2),VMAE(3)  !结构的杨氏模量，泊松比，密度
	PARAMETER(PI=3.1415926)                  
	COMMON  DENS,SONC
      DIMENSION  NN(8),VXY(4,3),FKK(8,8),FMM(8,8),FSQQQ(8,24),
     $         GKK(15,15),GMM(15,15),VMAE(4),FMAE(2),UX(50)
	REAL*8,POINTER::P(:,:)

      CHARACTER(LEN=80)::S,S1(100)
	ALLOCATABLE VNODE(:,:),IELEM(:,:),IELEF(:,:),IELES(:,:),SKK(:,:),
     $            SMM(:,:),THIC(:),N0(:),A(:,:),MAR(:,:),MARK(:,:)

	ALLOCATABLE MPTR(:),IND(:),VAL(:),RHS(:),NCP(:,:),NET(:,:)
C************************读单元信息文件*******************************
!	
!      CALL ASSEM_SOLID45()    !实体单元测试接口程序
!      CALL ASSEM_FLUID30()    !流体单元测试接口程序
	OPEN(200,  FILE= 'ux.dat', STATUS= 'UNKNOWN')
      OPEN(201,  FILE= 'uy.dat', STATUS= 'UNKNOWN')
      OPEN(202,  FILE= 'uz.dat', STATUS= 'UNKNOWN')
      OPEN(203,  FILE= 'rotx.dat', STATUS= 'UNKNOWN')
	OPEN(204, FILE= 'roty.dat', STATUS= 'UNKNOWN')
	OPEN(205, FILE= 'rotz.dat', STATUS= 'UNKNOWN')
      OPEN(206, FILE= 'pre.dat', STATUS= 'UNKNOWN')

      OPEN(100,FILE ='double_cylindrical2_5.dat', STATUS= 'OLD',
     $              IOSTAT=ERROR)
      KEYOPT_3=0   !控制参数，功能对应于Ansys中keyopt3
     
      DO 10 I=1,8
      READ(100,*)S,S,NN(I) !NN(1)为节点总数,NN(2)为单元总数
	                     !NN(3)为壳体单元数,NN(4)为实体单元数
	                     !NN(5)为流体单元数,NN(6)为材料参数种类    
						 !NN(7)为实常数种类,NN(8)为单元类型种类
10    CONTINUE
      NODENUM=NN(1)
    
      ALLOCATE(MAR((NN(3)+NN(4))*2+2,2))   !动态分配耦合面标志数组,比实际的要大
      ALLOCATE(THIC(NN(7)))          !动态分配材料实常数数组
      ALLOCATE(VNODE(NN(1),3))       !动态分配节点坐标矩阵内存
	ALLOCATE(IELEM(NN(3),4))       !动态分配壳单元编号矩阵内存
	ALLOCATE(IELES(NN(4),8))       !动态分配实体单元编号矩阵内存
      ALLOCATE(IELEF(NN(5),8))       !动态分配流体单元编号矩阵内存
      ALLOCATE(NET(NN(8),2))         !动态分单元种类数组

      ALLOCATE(N0(NN(2)))            !动态分配工作数组内存
	ALLOCATE(NCP(NN(7),3))
      THIC(:)=0.0
	VNODE(:,:)=0.0
      IELEM(:,:)=0.0
	IELES(:,:)=0.0
	IELEF(:,:)=0.0
	NET(:,:)=0.0
      N0(:)=0.0

      DO 11 II=1,NN(8)
      READ(100,*)S,NET(II,1),NET(II,2) !读单元种类
11    CONTINUE

      DO 12 II=1,NN(7)
      READ(100,*)I0,I0,THIC(II) !读结构厚度
12    CONTINUE
      
	DO 13 II=1,NN(1)
 	READ(100,*)I0,I1,I2,(VNODE(II,JJ),JJ=1,3) !VNODE为节点坐标矩阵
13    CONTINUE

      DO 14 II=1,NN(3)
	READ(100,*)I0,I0,N0(II),I0,I0,I0,I0,I0,I0,I0,I0,
     $       (IELEM(II,JJ),JJ=1,4)
14    CONTINUE                        !IELEM为每个壳单元节点编码
      
	DO 15 II=1,NN(4)      
	READ(100,*)I0,I0,N0(II+NN(3)),I0,I0,I0,I0,I0,I0,I0,I0,
     $       (IELES(II,JJ),JJ=1,8)
15    CONTINUE	               !IELEF为每个流体单元节点编码

	DO 16 II=1,NN(5)      
	READ(100,*)I0,I0,N0(II+NN(3)+NN(4)),I0,I0,I0,I0,I0,I0,I0,I0,
     $       (IELEF(II,JJ),JJ=1,8)
16    CONTINUE	               !IELEF为每个流体单元节点编码
!      WRITE(*,*)N0(1537),N0(4744)                     

      READ(100,*)S,S,S,S,S,S,VMAE(1)  !读结构杨氏模量
	READ(100,*)S,S,S,S,S,S,VMAE(2)  !读结构泊松比
      READ(100,*)S,S,S,S,S,S,VMAE(3)  !读结构密度
	READ(100,*)S,S,S,S,S,S,FMAE(1)  !读流体密度
      READ(100,*)S,S,S,S,S,S,FMAE(2)  !读流体声速
      
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

	ALLOCATE(MARK(NNN,2))    !分配确切需要的与耦合有关的边界条件内存
	MARK(:,:)=MAR(1:NNN,:)
	DEALLOCATE(MAR)          !释放无用内存
	DENS=FMAE(1)
	SONC=FMAE(2)

      M=MAXVAL(IELEM(:,:))         !壳单元中最大节点号
	N=MAXVAL(IELES(:,:))         !实体单元中最大节点号
      MAX=MAXVAL(IELEF(:,:))       !流体单元中最大节点号
      MIN=MINVAL(IELEF(:,:))       !流体单元中最小节点号
      IF(N.GE.M)THEN
      NDIM=M*6+(N-M)*3+MAX-MIN+1   !整体矩阵维数
	ELSE
      NDIM=M*6+MAX-MIN+1            !整体矩阵维数
	ENDIF
	
     
	ALLOCATE(SKK(NDIM,NDIM)) !动态分配整体刚度矩阵内存
	                               
	ALLOCATE(SMM(NDIM,NDIM)) !动态分配整体质量矩阵内存
	                               
C************************************************************************

!======================求解整体矩阵子程序并释放无用内存==================
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
	NST=0  !  !施加边界条件的节点数
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
	NST=0 !396  !施加边界条件的节点数
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

C*********************************调用求解器*****************************
!+++++++++++++++++++++++++++++对整体矩阵排序处理+++++++++++++++++++++++++
      NUMBER=1       !谐响应步长数
	FB=499.95         !起始频率
	FO=499.95         !截至频率
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
      
      NNS=COUNT(A(:,:).NE.0)  !计算非零元素个数
	
      ALLOCATE(MPTR(NDIM+1)) !NPTR的长度是N+1,PTR[I]记录第i列第一个非零元
	             !的位置,最后一个元素PTR[N＋1]=NN+1。这样第I列的长度就是PTR[I+1]-PTR[I]。
	ALLOCATE(IND(NNS))  !IND按列计算每个非零元素位置 
	ALLOCATE(VAL(NNS))  !VAL按列计算每个非零元素值
	IF(II.EQ.1)THEN
	ALLOCATE(RHS(NDIM))  !载荷向量
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
!++++++++++++++++++++++++++++求解器接口子程序+++++++++++++++++++++++++++++++
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

      DEALLOCATE(MPTR) !NPTR的长度是N+1,PTR[I]记录第i列第一个非零元
	             !的位置,最后一个元素PTR[N＋1]=NN+1。这样第I列的长度就是PTR[I+1]-PTR[I]。
	DEALLOCATE(IND)  !IND按列计算每个非零元素位置 
	DEALLOCATE(VAL)  !VAL按列计算每个非零元素值

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

      DEALLOCATE(SKK) !释放刚度矩阵内存
	DEALLOCATE(SMM) !释放质量矩阵内存
      DEALLOCATE(A)   !释放工作数组内存
C***************************************************************************
C***************************************************************************
      STOP
      END