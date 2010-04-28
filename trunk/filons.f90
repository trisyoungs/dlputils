	program filons 
	implicit none
	real*8 :: gr(

	! Load in test rdf
	

C    *******************************************************************
C    ** FOURIER SINE TRANSFORM BY FILON'S METHOD                      **
C    **                                                               **
C    ** A SPATIAL CORRELATION FUNCTION, H(R), IS TRANSFORMED TO       **
C    ** HHAT(K) IN RECIPROCAL SPACE.                                  **
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** FILON, PROC ROY SOC EDIN, 49 38, 1928.                        **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** REAL    KVEC                THE WAVENUMBER                    **
C    ** REAL    RMAX                MAXIMUM DIST IN CORREL. FUNCTION  **
C    ** REAL    ALPHA, BETA, GAMMA  FILON PARAMETERS                  **
C    ** REAL    H(NMAX)             THE CORRELATION FUNCTION          **
C    ** REAL    HHAT(NMAX)          THE 3-D TRANSFORM                 **
C    ** REAL    DR                  INTERVAL BETWEEN POINTS IN H      **
C    ** REAL    DK                  INTERVAL BETWEEN POINTS IN HHAT   **
C    ** INTEGER NMAX                NO. OF INTERVALS                  **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THE ROUTINE REQUIRES THAT THE NUMBER OF INTERVALS, NMAX, IS   **
C    ** EVEN AND CHECKS FOR THIS CONDITION. THE FIRST VALUE OF H(R)   **
C    ** IS AT R=0. THE MAXIMUM R FOR THE CORRELATION FUNCTION IS      **
C    ** RMAX=DR*NMAX. FOR AN ACCURATE TRANSFORM H(RMAX)=0.            **
C    *******************************************************************

       SUBROUTINE FILONS ( DR, DK, NMAX, H, HHAT )

        INTEGER     NMAX
        REAL        DR, DK, H(0:NMAX), HHAT(0:NMAX)

        REAL        RMAX, K, THETA, SINTH, COSTH
        REAL        SINSQ, COSSQ, THSQ, THCUB, ALPHA, BETA, GAMMA
        REAL        SE, SO, FOURPI, R
        INTEGER     IR, IK

C    *******************************************************************

C    ** CHECKS NMAX IS EVEN **

        IF ( MOD ( NMAX, 2 ) .NE. 0 ) THEN

           STOP ' NMAX SHOULD BE EVEN '

        ENDIF

        FOURPI = 16.0 * ATAN ( 1.0 )
        RMAX   = REAL ( NMAX ) * DR

C    ** LOOP OVER K **

        DO 30 IK = 0, NMAX

           K  = REAL ( IK ) * DK
           THETA = K * DR

C       ** CALCULATE THE FILON PARAMETERS **

           SINTH = SIN ( THETA )
           COSTH = COS ( THETA )
           SINSQ = SINTH * SINTH
           COSSQ = COSTH * COSTH
           THSQ  = THETA * THETA
           THCUB = THSQ * THETA

           IF ( THETA. EQ. 0.0 ) THEN

              ALPHA = 0.0
              BETA  = 2.0 / 3.0
              GAMMA = 4.0 / 3.0

            ELSE

              ALPHA = ( 1.0 / THCUB )
     :               * ( THSQ + THETA * SINTH * COSTH - 2.0 * SINSQ )
              BETA  = ( 2.0 / THCUB )
     :               * ( THETA * ( 1.0 + COSSQ ) -2.0 * SINTH * COSTH )
              GAMMA = ( 4.0 / THCUB ) * ( SINTH - THETA * COSTH )

           ENDIF

C       ** THE INTEGRAND IS H(R) * R FOR THE 3-D TRANSFORM **

C       ** DO THE SUM OVER THE EVEN ORDINATES **

           SE = 0.0

           DO 10 IR = 0, NMAX, 2

              R = REAL ( IR ) * DR
              SE = SE + H(IR) * R * SIN ( K * R )

10         CONTINUE

C       ** SUBTRACT HALF THE FIRST AND LAST TERMS **
C       ** HERE THE FIRST TERM IS ZERO            **

           SE = SE - 0.5 * ( H(NMAX) * RMAX * SIN ( K * RMAX ) )

C       ** DO THE SUM OVER THE ODD ORDINATES **

           SO = 0.0

           DO 20 IR = 1, NMAX - 1, 2

              R = REAL ( IR ) * DR
              SO = SO + H(IR) * R * SIN ( K * R )

20         CONTINUE

           HHAT(IK) = ( - ALPHA * H(NMAX) * RMAX * COS ( K * RMAX)
     :                 + BETA * SE + GAMMA * SO ) * DR

C       ** INCLUDE NORMALISING FACTOR **

           HHAT(IK) = FOURPI * HHAT(IK) / K

30      CONTINUE

        RETURN
        END


