
      SUBROUTINE ASSEM_FLUID30()
	IMPLICIT REAL * 8(A-H,O-Z)
	COMMON  DENS,SONC
      DIMENSION  NN(5),FXY(8,3),FKK(8,8),FMM(8,8),FSQQQ(8,24),
     $         GKK(15,15),GMM(15,15)
      CHARACTER(LEN=80)::S,S1(100)
	ALLOCATABLE FNODE(:,:),IELEM(:,:),GFKK(:,:),GFMM(:,:),N0(:,:)
C************************����Ԫ��Ϣ�ļ�*************************
      OPEN(100,FILE ='FILE3n.DAT', STATUS= 'OLD')
      
	DO 10 I=1,7
      READ(100,*)
10    CONTINUE

      DO 11 I=1,5
      READ(100,*)S,S,NN(I) !NN(1)Ϊ�ڵ�����,NN(2)Ϊ��Ԫ����
11    CONTINUE

      DO 12 I=1,31
	READ(100,*)s1(i)
12    CONTINUE  
      ALLOCATE(FNODE(NN(1),3))    !��̬����ڵ���������ڴ�
	ALLOCATE(IELEM(NN(2),8))     !��̬���䵥Ԫ��ž����ڴ�
	ALLOCATE(GFKK(NN(1),NN(1))) !��̬��������նȾ����ڴ�
	ALLOCATE(GFMM(NN(1),NN(1))) !��̬�����������������ڴ�
	ALLOCATE(N0(NN(2),10))      !��̬���乤�������ڴ�

	GFKK(:,:)=0.0   !������նȾ�������
	GFMM(:,:)=0.0   !������������������
 
	READ(100,*)(I,I1,I2,(FNODE(II,JJ),JJ=1,3),II=1,NN(1)) 
                       !FNODEΪ�ڵ��������
	DO 13 I=1,4
	READ(100,*)
13    CONTINUE

	READ(100,*)((N0(II,JJ),JJ=1,10),I,(IELEM(II,JJ),JJ=1,8),II=1,NN(2))
	                                   !IELEMΪÿ����Ԫ�ڵ����
	READ(100,*)
	READ(100,*)
	READ(100,*)
	READ(100,*)S,S,S,S,S,S,DENS   !�������ܶ�
      READ(100,*)
	READ(100,*)S,S,S,S,S,S,SONC  !����������
      
      DEALLOCATE(N0)   !�ͷ�N0
C*********************����Ԫ��Ϣ�ļ�����**********************
C====================================================================
C*******************��Ԫѭ�����ɵ��ռ����ܸ�******************
      DO 50 IE=1,NN(2)   !NN(2)Ϊ��Ԫ����
      
	DO 51 I=1,8
!	DO 51 J=1,3
      FXY(I,:)=FNODE(IELEM(IE,I),:)
51    CONTINUE
      
!	IF(IELEM(IE,3).EQ.IELEM(IE,4).AND.IELEM(IE,5).EQ.IELEM(IE,6).EQ.
!     $    IELEM(IE,7).EQ.IELEM(IE,8))THEN
     	                                     !���øú����γ�4�ڵ㵥Ԫ����
!	ELSEIF(IELEM(IE,3).EQ.IELEM(IE,4).AND.IELEM(IE,7).EQ.IELEM(IE,8))
 !    $	THEN
 !     CALL FLUID30_6_II(FXY,FKK,FMM,FSQQQ) !���øú����γ�6�ڵ㵥Ԫ����
!	ELSE
      CALL FLUID30_8_II(FXY,FKK,FMM) !���øú����γ�8�ڵ㵥Ԫ����
!	ENDIF


C*****************�ɵ�Ԫ�նȡ��������󼯳�����ն���������**************
      
	DO 52 I=1,8
      DO 52 J=1,8
      II=IELEM(IE,I)
      JJ=IELEM(IE,J)
	GFKK(II,JJ)=GFKK(II,JJ)+FKK(I,J)
	GFMM(II,JJ)=GFMM(II,JJ)+FMM(I,J)
52    CONTINUE
50    CONTINUE
C*****************************���󼯳ɹ��̽���**************************
C=========================================================================
C************************�ѵ�Ԫ������������д���ļ�*******************
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
C******************************����д���ļ�����*************************
!      GKK(:,:)=GFKK(:,:)
!	GMM(:,:)=GFMM(:,:)
!	WRITE(*,*)GFKK(38,35),GFKK(38,36),GFKK(38,37)

	RETURN
	END
      