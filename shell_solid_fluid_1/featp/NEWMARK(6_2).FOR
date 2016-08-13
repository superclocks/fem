!!	P747--749

C  ......................  SUB: 6-2  ..................................
            SUBROUTINE  NEWMARK  (GMM, GK K, DAMP, GP, U0, V0, A, AW, B)
C **********************************************************************
C           SOLVE DYNAMIC RESPONSE BY NEWMARK METHOD
C             CALL SUBROUTINES: DECOMPOS, BACKSUBS
C **********************************************************************
            IMPLICIT REAL * 8 (A-H, O-Z)
              COMMON/COM1/MND, NUMEL, NUMPT, MBAND, NMATI
              COMMON/COM2/NF, NFSTR, MSOLV, MPROB, MTYPE, NVA
              COMMON/COM3/MND2 , NUMPT2
            COMMON/DYN/XMUV, OMEGA, CC1, CC2, TT, DT, ALP HA, DELTA
            DIMENSION    GMM(NUMPT2, MBAND) , GKK(NUMPT2,MBAND) ,
     $      AW (NUMPT2 , MBAND) , DAMP ( NUMPT2, MBAND) , GP (NUMPT2) ,
     $                    U0 (NUMPT2), V0 (NUMPT2), A(NUMPT2),B(NUMPT2)
              CLOSE(31, STATUS= 'DELETE')
            OPEN(31, FILE= 'OUT_NMK', STATUS= 'UNKNOWN')
            WRITE(31, * )' DYNAMIC RESPONSE RESULT BY NEWMARK METHOD'
            WRITE(31, * ) ' TOTAL TIME= ', TT
            WRITE( * ,101)
101      FORMAT(/6X,' #  SOLVE DYNAMIC RESPONSE BY NEWMARK METHOD  # #')
            WRITE( * , 102)
102      FORMAT(/6X, '% OUTPUT DYNAMIC DISPLACEMENT IN FILE<OUT_NMK>')
C ***
C ************** INITIAL COMPUTATIONS
C ***
              C0= 1.0/(ALPHA * DT * DT)
              C1 = DELTA/(ALPH A * DT)
              C2=1. 0/(ALPHA * DT)
              C3=0.5/ALPHA-1. 0
              C4 = DELTA/ALPHA-1.0
              C5=DT/2.0 * (DELTA/ALPHA-2.0)
              C6=DT * (1.0-DELTA)
              C7 = DELTA * DT
C **********************************************************************
C         COMPUTE  K'=K+C0 * M+C1 * C
C **********************************************************************
            DO 40 I=1 ,NUMPT2
            B(I) = GP(I)
            DO 40 J=1 ,MBAND
            AW(I,J) =GKK(I,J)
40     GKK(I,J) =GKK(I,J) +C0 * GMM(I,J)+C1 * DAMP(I,J)
C **********************************************************************
C         TRIANGLE DECOMPOSITION    OF THE MATRIX.-  [GKK]
C **********************************************************************
              CALL DECOMPOS (NUMPT2,MBAND,GKK)
C ***
C************** COMPUTATIONS FOR EACH TIME STEP
C ***
            DO 300 Y=0,TT,DT
            IF(OMEGA. GT. 0.0) AP=SIN(OMEGA * Y)
            IF(OMEGA. LT. 0.0) AP=COS(OMEGA * Y)
            DO 41 I= 1, NUMPT2
            GP(I) =B(I)
            IF(OMEGA. NE. 0.0) GP(I)=GP(I)*AP
41          CONTINUE
            DO 50 I= 1,NUMPT2
            IK = I-MBAND+ 1
            IF(IK. LT. 1) IK=1
            DO 50 J=IK,I
            GP(I) =GP(I)+GMM(J ,I-J+ 1)*(C0*U0(J)+C2*V0(J)+C3 * A(J))+
     $                     DAMP(J,I-J+1)*(C1 * U0(J)+C4* V0(J)+C5*A(J))
50      CONTINUE
            DO 60 I=1,NUMPT2
            IK = I+ MBAND-1
            IF(IK. GT. NUMPT2) IK=NUMPT2
            DO 60 J=I+1,IK
      GP(I) = GP(I) +GMM(I,J-I+1)*(C0* U0(J) +C2 * V0(J)+C3*A(J)) +
     $              DAMP(I,J-I+I)*(C1 * U0(J) +C4 * V0(J) +C5*A(J))
60      CONTINUE
              CALL BACKSUBS (NUMPT2, MBAND, GKK, GP)
            NSTEP = Y/DT+ 1
            WRITE(31,61) NSTEP,DT,Y+DT
            IF(MPROB. EQ. 4) THEN
            WRITE(31,73)
            WRITE(31,74) (GP(I) ,I=1 ,NUMPT2)
            ELSE
           WRITE(31,71)
           WRITE( 31,72) (GP(I), I= 1 ,NUMPT2)
           ENDIF
           DO 70 I= 1, NUMPT2
           A(I) =C0 * (GP(I)-U0(I))-C2*V0(I)-C3*A(I)
           V0(I) =V0(I)+C6*A(I)+C7*A(I)
           U0(I) =GP(I)
70      CONTINUE
300      CONTINUE
61       FORMAT( 2X, 'NO. OF STEP = ', I5,3X, 'STEP LENGTH =', E15.6,3X,
     $                  'AT TIME= ', E15.6)
71      FORMAT(2X, 'DISPLACEMENT: '/2 (13X, 'X-', 13X, 'Y-'))
72      FORMAT(4X,4E16.8)
73      FORMAT(2X, 'DISPLACEMENT: '/11X 'THETA-X', 11X, 'THETA-Y',
     $                   12X, 'W-Z')
74              FORMAT(6X,3E16.8)
              CLOSE ( 31 )
            RETURN
            END