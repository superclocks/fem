      
	SUBROUTINE ASSEM_SOLID45()
      
      IMPLICIT REAL * 8(A-H,O-Z)

      DIMENSION  NN(5),SXY(8,3),VMAE(4),SCOS(3,3),
     $         GSK(24,24),GSM(24,24)
      CHARACTER(LEN=80)::S,S1(100)
	ALLOCATABLE FNODE(:,:),IELEM(:,:),GFKK(:,:),GFMM(:,:),N0(:,:)
C************************读单元信息文件*************************
      OPEN(500,FILE ='solid1.dat', STATUS= 'OLD')

      DO 11 I=1,5
      READ(500,*)S,S,NN(I) !NN(1)为节点总数,NN(2)为单元总数
11    CONTINUE

      ALLOCATE(FNODE(NN(1),3))    !动态分配节点坐标矩阵内存
	ALLOCATE(IELEM(NN(2),8))     !动态分配单元编号矩阵内存
	ALLOCATE(GFKK(NN(1)*3,NN(1)*3)) !动态分配整体刚度矩阵内存
	ALLOCATE(GFMM(NN(1)*3,NN(1)*3)) !动态分配整体质量矩阵内存
	ALLOCATE(N0(NN(2),10))      !动态分配工作数组内存

	GFKK(:,:)=0.0   !对整体刚度矩阵清零
	GFMM(:,:)=0.0   !对整体质量矩阵清零
 
	READ(500,*)(I,I1,I2,(FNODE(II,JJ),JJ=1,3),II=1,NN(1)) 
                       !FNODE为节点坐标矩阵
	READ(500,*)((N0(II,JJ),JJ=1,10),I,(IELEM(II,JJ),JJ=1,8),II=1,NN(2))
	                                   !IELEM为每个单元节点编码
	VMAE(:)=0.0
	READ(500,*)S,S,S,S,S,S,VMAE(1)  !杨氏模量
	READ(500,*)S,S,S,S,S,S,VMAE(2)  !泊松比
      READ(500,*)S,S,S,S,S,S,VMAE(3)  !密度
      
      DEALLOCATE(N0)   !释放N0
C*********************读单元信息文件结束**********************
C====================================================================
C*******************单元循环并由单刚集成总刚******************
      DO 50 IE=1,NN(2)   !NN(2)为单元总数
      
	DO 51 I=1,8

      SXY(I,:)=FNODE(IELEM(IE,I),:)
51    CONTINUE
      

      KEYOPT2=0
      CALL ELEMENT_MATRIX_SOLID(SXY,SCOS,GSK,GSM,VMAE)
C*****************由单元刚度、质量矩阵集成整体刚度质量矩阵**************
      
	DO 52 I=1,8
      DO 52 J=1,8
      II=(IELEM(IE,I)-1)*3+1
      JJ=(IELEM(IE,J)-1)*3+1
	I1=(I-1)*3+1
	J1=(J-1)*3+1
	GFKK(II:II+2,JJ:JJ+2)=GFKK(II:II+2,JJ:JJ+2)+GSK(I1:I1+2,J1:J1+2)
	GFMM(II:II+2,JJ:JJ+2)=GFMM(II:II+2,JJ:JJ+2)+GSM(I1:I1+2,J1:J1+2)
52    CONTINUE
50    CONTINUE
C*****************************矩阵集成过程结束**************************
C=========================================================================
C************************把单元矩阵和整体矩阵写入文件*******************
!	OPEN(101,FILE='FKK.DAT',STATUS= 'UNKNOWN')
!      OPEN(102,FILE='FMM.DAT',STATUS= 'UNKNOWN')!
!	OPEN(103,FILE='GFKK.DAT',STATUS= 'UNKNOWN')
!	OPEN(104,FILE='GFMM.DAT',STATUS= 'UNKNOWN')
!	WRITE(101,110)(II,(FKK(II,JJ),JJ=1,8),II=1,8)
!      WRITE(102,110)(II,(FMM(II,JJ),JJ=1,8),II=1,8)
!	WRITE(103,111)(II,(GFKK(II,JJ),JJ=1,NN(1)),II=1,NN(1))
!	WRITE(104,111)(II,(GFMM(II,JJ),JJ=1,NN(1)),II=1,NN(1))
!110    FORMAT(1I,8E18.8)
!111    FORMAT(1I,4E18.8)
	END