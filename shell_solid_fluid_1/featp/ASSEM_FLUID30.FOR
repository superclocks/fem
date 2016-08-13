
      SUBROUTINE ASSEM_FLUID30()
	IMPLICIT REAL * 8(A-H,O-Z)
	COMMON  DENS,SONC
      DIMENSION  NN(5),FXY(8,3),FKK(8,8),FMM(8,8),FSQQQ(8,24),
     $         GKK(15,15),GMM(15,15)
      CHARACTER(LEN=80)::S,S1(100)
	ALLOCATABLE FNODE(:,:),IELEM(:,:),GFKK(:,:),GFMM(:,:),N0(:,:)
C************************读单元信息文件*************************
      OPEN(100,FILE ='FILE3n.DAT', STATUS= 'OLD')
      
	DO 10 I=1,7
      READ(100,*)
10    CONTINUE

      DO 11 I=1,5
      READ(100,*)S,S,NN(I) !NN(1)为节点总数,NN(2)为单元总数
11    CONTINUE

      DO 12 I=1,31
	READ(100,*)s1(i)
12    CONTINUE  
      ALLOCATE(FNODE(NN(1),3))    !动态分配节点坐标矩阵内存
	ALLOCATE(IELEM(NN(2),8))     !动态分配单元编号矩阵内存
	ALLOCATE(GFKK(NN(1),NN(1))) !动态分配整体刚度矩阵内存
	ALLOCATE(GFMM(NN(1),NN(1))) !动态分配整体质量矩阵内存
	ALLOCATE(N0(NN(2),10))      !动态分配工作数组内存

	GFKK(:,:)=0.0   !对整体刚度矩阵清零
	GFMM(:,:)=0.0   !对整体质量矩阵清零
 
	READ(100,*)(I,I1,I2,(FNODE(II,JJ),JJ=1,3),II=1,NN(1)) 
                       !FNODE为节点坐标矩阵
	DO 13 I=1,4
	READ(100,*)
13    CONTINUE

	READ(100,*)((N0(II,JJ),JJ=1,10),I,(IELEM(II,JJ),JJ=1,8),II=1,NN(2))
	                                   !IELEM为每个单元节点编码
	READ(100,*)
	READ(100,*)
	READ(100,*)
	READ(100,*)S,S,S,S,S,S,DENS   !读流体密度
      READ(100,*)
	READ(100,*)S,S,S,S,S,S,SONC  !读流体声速
      
      DEALLOCATE(N0)   !释放N0
C*********************读单元信息文件结束**********************
C====================================================================
C*******************单元循环并由单刚集成总刚******************
      DO 50 IE=1,NN(2)   !NN(2)为单元总数
      
	DO 51 I=1,8
!	DO 51 J=1,3
      FXY(I,:)=FNODE(IELEM(IE,I),:)
51    CONTINUE
      
!	IF(IELEM(IE,3).EQ.IELEM(IE,4).AND.IELEM(IE,5).EQ.IELEM(IE,6).EQ.
!     $    IELEM(IE,7).EQ.IELEM(IE,8))THEN
     	                                     !调用该函数形成4节点单元矩阵
!	ELSEIF(IELEM(IE,3).EQ.IELEM(IE,4).AND.IELEM(IE,7).EQ.IELEM(IE,8))
 !    $	THEN
 !     CALL FLUID30_6_II(FXY,FKK,FMM,FSQQQ) !调用该函数形成6节点单元矩阵
!	ELSE
      CALL FLUID30_8_II(FXY,FKK,FMM) !调用该函数形成8节点单元矩阵
!	ENDIF


C*****************由单元刚度、质量矩阵集成整体刚度质量矩阵**************
      
	DO 52 I=1,8
      DO 52 J=1,8
      II=IELEM(IE,I)
      JJ=IELEM(IE,J)
	GFKK(II,JJ)=GFKK(II,JJ)+FKK(I,J)
	GFMM(II,JJ)=GFMM(II,JJ)+FMM(I,J)
52    CONTINUE
50    CONTINUE
C*****************************矩阵集成过程结束**************************
C=========================================================================
C************************把单元矩阵和整体矩阵写入文件*******************
	OPEN(101,FILE='FKK.DAT',STATUS= 'UNKNOWN')
      OPEN(102,FILE='FMM.DAT',STATUS= 'UNKNOWN')
	OPEN(103,FILE='GFKK.DAT',STATUS= 'UNKNOWN')
	OPEN(104,FILE='GFMM.DAT',STATUS= 'UNKNOWN')
	WRITE(101,110)(II,(FKK(II,JJ),JJ=1,8),II=1,8)
      WRITE(102,110)(II,(FMM(II,JJ),JJ=1,8),II=1,8)
	WRITE(103,111)(II,(GFKK(II,JJ),JJ=1,NN(1)),II=1,NN(1))
	WRITE(104,111)(II,(GFMM(II,JJ),JJ=1,NN(1)),II=1,NN(1))
110    FORMAT(1I,8E18.8)
111    FORMAT(1I,4E18.8)
C******************************矩阵写入文件结束*************************
!      GKK(:,:)=GFKK(:,:)
!	GMM(:,:)=GFMM(:,:)
!	WRITE(*,*)GFKK(38,35),GFKK(38,36),GFKK(38,37)

	RETURN
	END
      