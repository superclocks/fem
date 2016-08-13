!!	 P744--745
c..........  SUB: 6  .............................
            SUBROUTINE DYNAM (GMM, GKK, GP, DAMP, U0, V0, A, AW, B)
c*****************************************************************
C      SOLVE THE DYNAMIC PROBLEM BY  CENTRAL DIFFERENCE  METHOD
C              OR NEWMARK METHOD
C    CALL SUBROUTINES: DECOMPOS,BACKSUBS,CENTER OR NEWMARK
c*****************************************************************
            IMPLICIT REAL * 8(A-H,O-Z)
              COMMON/COM1/MND, NUMEL, NUMPT, MBAND, NMATI
              COMMON/COM2/NF, NFSTR, MSOLV, MPROB, MTYPE, NVA
              COMMON/COM3/MND2, NUMPT2
              COMMON/DYN/XMUV, OMEGA, CC1,CC2, TT, DT, ALPHA, DELTA
         DIMENSION  GMM(NUMPT2, MBAND) ,GKK(NUMPT2 ,MBAND) ,GP(NUMPT2) ,
     $ A( N UMPT2) , U0 ( NUMPT2 ), V0 (NU MPT2) , AW ( NU MPT2, MBAND),
     $                   B( NU MPT2) , DAMP ( N UMPT2, MBAND)

c ***
C ********  COMPUTE DAMPLING MATRIX (DAMP)
c ***
            DO 50 I=1,NUMPT2
            DO 50 J = 1, MBAND
            DAMP(I,J) =CC1 * GMM(I,J)+CC2 * GKK(I,J)
50      CONTINUE
c ***
C *********** FORM INITIAL DISPLACEMENT (U0) AND VELOCITY (V0)
C ***
            DO  60  I=1,NUMPT2
            U0(I) =0.0
            V0(I) =0.0
60      CONTINUE
          IF(XMUV. EQ. 1. OR. XMUV. EQ. 3) THEN
            READ(5, * )
            READ(5, * ) (U0(I),I=1,NUMPT2)
          ENDIF
          IF(XMUV. EQ. 2. OR.X MUV. EQ. 3) THEN
            READ(5, * )
            READ(5, * ) (V0(I),I=1,NUMPT2)
          ENDIF
            WRITE(6, * )' INITIAL DISPLACEMENTS--U0:'
            WRITE(6,62) (U0(I) ,I= 1,NUMPT2)
            WRITE(6, * ) ' INITIAL VELOCITIES--V0: '
            WRITE(6,62) (V0(I) ,I=1 ,NUMPT2)
62      FORMAT(2X,4E15.6)
c *********************************************************** 
C                 COMPUTE  INITIAL  ACCELERATION--A(NUMPT2)
c *********************************************************** 
            DO  71  I=1,NUMPT2
71          A(I) = GP(I)
            DO 70 I=1,NUMPT2
            IK = I-MBAND+ 1
            IF(IK. LT. 1) IK=1
            DO 70 J=IK,I
            A(I) =A(I)-GKK(J,I-J+1) * U0(J)-DAMP(J,I-J+1) * V0(J)
70      CONTINUE
            DO 72 I=1,NUMPT2
            IK = I+ MBAND-1
            IF(IK. GT. NUMPT2) IK= NUMPT2
            DO 72 J=I+1 ,IK
            A(I) =A(I)-GKK(I,J-I+1) * U0(J)-DAMP(I,J-I+1) * V0(J)
72      CONTINUE
            DO  73  I=1,NUMPT2
            DO 73 J=1,MBAND
            AW(I,J) =GMM(I,J)
73      CONTINUE
              CALL DECOMPOS (NUMPT2, MBAND, AW)
              CALL BACKSUBS  (NUMPT2, MBAND, AW, A)
C ******
C ********* BY USING CENTRAL DIFFERENCE METHOD WHEN MSOLV= g
C ***
           IF(MSOLV. EQ. 2)  CALL CENTER  (GMM,GKK,DAMP,GP,U0,V0,A,AW,B)
C ***
C *******  BY USING NEWMARK METHOD WHEN MSOLV=3
C ***
          IF(MSOLV. EQ. 3)  CALL NEWMARK  (GMM,GKK,DAMP,GP,U0,V0,A,AW,B)
            RETURN
            END