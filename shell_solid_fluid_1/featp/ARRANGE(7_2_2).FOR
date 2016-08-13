!!	P761--762

C  ========================  SUB: 7-2-2  ===================================
      SUBROUTINE ARRANGE(W ,V, N)
C  *************************************************************************
C         ARRANGE  EIGENVALUE INTO ORDER
C  *************************************************************************
           IMPLICIT  REAL* 8(A-H,O-Z)
           DIMENSION    W(N),V(N,N)
           DO 13 I=1,N-1
           K=I
           P=W(I)
           DO 11 J=I+1,N
           IF(W(J).LT.P)  THEN
           K=J
           P=W(J)
           ENDIF
11    CONTINUE
           IF(K.NE.I)  THEN
           W(K)=W(I)
           W(I) =P
           DO 12 J=1,N
           P=V(J,I)
           V(J,I)=V(J,K)
           V(J,K)=P
12    CONTINUE
           ENDIF
13    CONTINUE
           RETURN
           END
