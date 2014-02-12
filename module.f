!=========================================================================================
!=========================================================================================
! Contains subroutines: JACOBI, EIGSRT
!          function : THRJ 

       module f77_subroutine
         contains
C******* DIAGONALIZE THE REAL SYMMETRIC MATRIX *******************
C
C            ..... good for small matrices .....
C
      SUBROUTINE JACOBI(A,N,NP,D,V,NROT)
! A is the matrix that you want to diagonalize
! N, NP is the size of A
! D is the column matrix of eigenvalues
! V is the eigenvector matrix whose columns are normalized eigenvectors of A
! NROT is the number of rotations has been performed during the diagonalization process (out)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=250)
      DIMENSION A(NP,NP),D(NP),V(NP,NP),B(NMAX),Z(NMAX)
      DO 12 IP=1,N
        DO 11 IQ=1,N
          V(IP,IQ)=0.
11      CONTINUE
        V(IP,IP)=1.
12    CONTINUE
      DO 13 IP=1,N
        B(IP)=A(IP,IP)
        D(IP)=B(IP)
        Z(IP)=0.
13    CONTINUE
      NROT=0
      DO 24 I=1,50
        SM=0.
        DO 15 IP=1,N-1
          DO 14 IQ=IP+1,N
            SM=SM+ABS(A(IP,IQ))
14        CONTINUE
15      CONTINUE
        IF(SM.EQ.0.)RETURN
        IF(I.LT.4)THEN
          TRESH=0.2*SM/N**2
        ELSE
          TRESH=0.
        ENDIF
        DO 22 IP=1,N-1
          DO 21 IQ=IP+1,N
            G=100.*ABS(A(IP,IQ))
            IF((I.GT.4).AND.(ABS(D(IP))+G.EQ.ABS(D(IP)))
     *         .AND.(ABS(D(IQ))+G.EQ.ABS(D(IQ))))THEN
              A(IP,IQ)=0.
            ELSE IF(ABS(A(IP,IQ)).GT.TRESH)THEN
              H=D(IQ)-D(IP)
              IF(ABS(H)+G.EQ.ABS(H))THEN
                T=A(IP,IQ)/H
              ELSE
                THETA=0.5*H/A(IP,IQ)
                T=1./(ABS(THETA)+SQRT(1.+THETA**2))
                IF(THETA.LT.0.)T=-T
              ENDIF
              C=1./SQRT(1+T**2)
              S=T*C
              TAU=S/(1.+C)
              H=T*A(IP,IQ)
              Z(IP)=Z(IP)-H
              Z(IQ)=Z(IQ)+H
              D(IP)=D(IP)-H
              D(IQ)=D(IQ)+H
              A(IP,IQ)=0.
              DO 16 J=1,IP-1
                G=A(J,IP)
                H=A(J,IQ)
                A(J,IP)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
16            CONTINUE
              DO 17 J=IP+1,IQ-1
                G=A(IP,J)
                H=A(J,IQ)
                A(IP,J)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
17            CONTINUE
              DO 18 J=IQ+1,N
                G=A(IP,J)
                H=A(IQ,J)
                A(IP,J)=G-S*(H+G*TAU)
                A(IQ,J)=H+S*(G-H*TAU)
18            CONTINUE
              DO 19 J=1,N
                G=V(J,IP)
                H=V(J,IQ)
                V(J,IP)=G-S*(H+G*TAU)
                V(J,IQ)=H+S*(G-H*TAU)
19            CONTINUE
              NROT=NROT+1
            ENDIF
21        CONTINUE
22      CONTINUE
        DO 23 IP=1,N
          B(IP)=B(IP)+Z(IP)
          D(IP)=B(IP)
          Z(IP)=0.
23      CONTINUE
24    CONTINUE
      PAUSE '50 iterations should never happen'
      RETURN
      END SUBROUTINE
C
C**** SORT THE EIGENVALUES AND VECTORS OF THE DIAGONALIZED MATRIX *********
C
      SUBROUTINE EIGSRT(D,V,N,NP)
C
C  Modified by R. Krems to sort the eigenvalues and 
C  eigenvectors in the ascending rather then descending 
C  order
C               April 20, 2000, Goteborg
C
C
! D is the column matrix of eigenvalues
! V is the eigenvector matrix whose columns are eigenvectors
! N, NP is the size of V
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION D(NP),V(NP,NP)
      DO 13 I=1,N-1
        K=I
        P=D(I)
        DO 11 J=I+1,N
          IF(D(J).LE.P)THEN
            K=J
            P=D(J)
          ENDIF
11      CONTINUE
        IF(K.NE.I)THEN
          D(K)=D(I)
          D(I)=P
          DO 12 J=1,N
            P=V(J,I)
            V(J,I)=V(J,K)
            V(J,K)=P
12        CONTINUE
        ENDIF
13    CONTINUE
      RETURN
      END SUBROUTINE
   
C**********************************************************************
C
      DOUBLE PRECISION FUNCTION THRJ(F1,F2,F3,G1,G2,G3)
C  the input variables should be of double precision types ----added by Ping
!  calculate the value of three-j symbol
!     / F1     F2     F3 \
!    |                    |
!     \ G1     G2     G3 /
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     SMALL CHANGES 31 JUL 95 (SG)
      SAVE MUNG,X,Y
      PARAMETER (MXIX=302)
      DIMENSION X(MXIX),Y(MXIX)
      DATA MUNG/0/
      IF (MUNG.EQ.21) GO TO 69
      MUNG = 21
      X(1) = 0.D0
      DO 100 I = 1, MXIX-1
      A = I
      X(I+1) = LOG(A) +X(I)
      Y(I+1) = LOG(A)
  100 CONTINUE
   69 IF(F1-ABS(G1)) 1,13,13
   13 IF(F2-ABS(G2))1,14,14
   14 IF(F3-ABS(G3))1,15,15
   15 SUM=F1+F2+F3
      NSUM=SUM+.001D0
      IF(SUM-NSUM)2,2,1
    1 THRJ=0.D0
      RETURN
    2 IF(ABS(G1+G2+G3)-1.D-08)3,3,1
    3 IF(F1+F2-F3)1,4,4
    4 IF(F1+F3-F2)1,5,5
    5 IF(F2+F3-F1)1,6,6
    6 J1=2.D0*F3+2.001D0
      J2=F1+F2-F3+1.001D0
      J3=F1-F2+F3+1.001D0
      J4=-F1+F2+F3+1.001D0
      J5=F1+F2+F3+2.001D0
      J6=F1+G1+1.001D0
      J7=F1-G1+1.001D0
      J8=F2+G2+1.001D0
      J9=F2-G2+1.001D0
      J10=F3+G3+1.001D0
      J11=F3-G3+1.001D0
      IF(J5.GT.MXIX) THEN
        WRITE(6,601) J5,MXIX
  601   FORMAT(' *** DIMENSION ERROR IN THRJ - INDEX.GT.MXIX',2I5)
        STOP
      ENDIF
      R=0.5D0*(Y(J1)+X(J2)+X(J3)+X(J4)-X(J5)
     1+X(J6)+X(J7)+X(J8)+X(J9)+X(J10)+X(J11))
      SUM=0.D0
      F=-1
      KZ=-1
    7 KZ=KZ+1
      F=-F
      J1=KZ+1
      J2=F1+F2-F3-KZ+1.001D0
      IF(J2)20,20,8
    8 J3=F1-G1-KZ+1.001D0
      IF(J3)20,20,9
    9 J4=F2+G2-KZ+1.001D0
      IF(J4)20,20,10
   10 J5=F3-F2+G1+KZ+1.001D0
      IF(J5)7,7,11
   11 J6=F3-F1-G2+KZ+1.001D0
      IF(J6)7,7,12
   12 JMAX=MAX(J1,J2,J3,J4,J5,J6)
      IF(JMAX.GT.MXIX) THEN
        WRITE(6,601) JMAX,MXIX
        STOP
      ENDIF
      S=-(X(J1)+X(J2)+X(J3)+X(J4)+X(J5)+X(J6))
      SUM=SUM+F*EXP(R+S)
      GO TO 7
   20 INT=ABS(F1-F2-G3)+0.0001D0
      VAL=((-1.D0)**INT)*SUM/SQRT(2.D0*F3+1.D0)
      IF(ABS(VAL).LE.1.D-6) VAL=0.D0
      THRJ=VAL
      RETURN
      END FUNCTION
C
!***********************************************************************
! input A = a real symmetric matrix to be diagonalized
! due to symmetry of A, you don't need to specify every elements of A
! in this subroutine, you can only store the uptriangle elemnts in A
! ouput A = orthonormal eigenvector matrix
! N ----- dimension of A
! W ----- eigenvalue array
! if the subroutine fails, it will generate a error message: lapack_eig.err
      SUBROUTINE lapack_eig(A, N, W)
      INTEGER          N
      INTEGER          LDA
      INTEGER          LWMAX
      PARAMETER        ( LWMAX = 100000000  ) ! LWMAX >= 1 + 6*N + 2*N**2 
! don't delete () in the above line
!
!     .. Local Scalars ..
      INTEGER          INFO, LWORK, LIWORK
!
!     .. Local Arrays ..
      INTEGER          IWORK( LWMAX )
      DOUBLE PRECISION :: A( N, N )  
      double precision :: W( N ), WORK( LWMAX ) 

!
!     .. External Subroutines ..
      EXTERNAL         DSYEVD

!
!     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
!
      LDA=N
!
!     Query the optimal workspace.
!
      LWORK = -1
      LIWORK = -1
      CALL DSYEVD( 'Vectors', 'Upper', N, A, LDA, W, WORK, LWORK,
     $             IWORK, LIWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      LIWORK = MIN( LWMAX, IWORK( 1 ) )
!
!     Solve eigenproblem.
!
      CALL DSYEVD( 'Vectors', 'Upper', N, A, LDA, W, WORK, LWORK,
     $             IWORK, LIWORK, INFO )
!
!     Check for convergence.
!
!      open(unit=13,file="lapack_eig.err")
      IF( INFO.GT.0 ) THEN
         WRITE(13,*)'The algorithm failed to compute eigenvalues.'
!         close(13)
         STOP
      END IF
      end subroutine lapack_eig 

!**********************************************************************
! slove A*X=B
! dimension of A : N by N
! dimension of B : N by 1
! in output, B is the solution X 
      SUBROUTINE lapack_solve(A, N, B)
      INTEGER          N, NRHS
      PARAMETER        ( NRHS = 1 ) ! don't delete the ()
      INTEGER          LDA, LDB
!
!     .. Local Scalars ..
      INTEGER          INFO
!
!     .. Local Arrays ..
      INTEGER          IPIV( N )
      DOUBLE PRECISION A( N, N ), B( N, NRHS )

!     .. External Subroutines ..
      EXTERNAL         DGESV

!
      LDA = N
      LDB = N
!
!     Solve the equations A*X = B.
!
      CALL DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
!     Check for the exact singularity.
!
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The diagonal element of the triangular factor of A,'
         WRITE(*,*)'U(',INFO,',',INFO,') is zero, so that'
         WRITE(*,*)'A is singular; the solution could not be computed.'
         STOP
      END IF
      END subroutine lapack_solve 

!***********************************************************************
! this subroutine has a low requirement for memory
! thus it is recommended for matrix with dimension > 5000X5000 
! if encounter very large matrix, try lapack_eig2
! input A = a real symmetric matrix to be diagonalized
! due to symmetry of A, you don't need to specify every elements of A
! in this subroutine, you can only store the uptriangle elemnts in A
! note part of A will be destoryed on exit
!---------------------------------------------------------------------
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of A contains the
*          upper triangular part of the matrix A.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of A contains
*          the lower triangular part of the matrix A.
*          On exit, the lower triangle (if UPLO='L') or the upper
*          triangle (if UPLO='U') of A, including the diagonal, is
*          destroyed.
!----------------------------------------------------------------------
! ouput Z = orthonormal eigenvector matrix
! N ----- dimension of A
! W ----- eigenvalue array from small value to large value
! if the subroutine fails, it will generate a error message: lapack_eig.err
      SUBROUTINE lapack_eig2(A, N, W, Z)
      INTEGER          N
      INTEGER          LDA, LDZ
      INTEGER          LWMAX
      PARAMETER        ( LWMAX = 10000000  )  
! don't delete () in the above line
!
!     .. Local Scalars ..
      INTEGER          INFO, LWORK, LIWORK, IL, IU, M
      DOUBLE PRECISION ABSTOL, VL, VU
!
!     .. Local Arrays ..
      INTEGER          ISUPPZ( 2*N ), IWORK( LWMAX )
      DOUBLE PRECISION :: A( N, N )  
      double precision :: W( N ), WORK( LWMAX ), Z( N, N )

!
!     .. External Subroutines ..
      EXTERNAL         DSYEVR

!
!     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
!
      LDA=N
      LDZ=N
!     Negative ABSTOL means using the default value
      ABSTOL = -1.0
!
!     Query the optimal workspace.
!
      LWORK = -1
      LIWORK = -1
      CALL DSYEVR( 'Vectors', 'A', 'Upper', N, A, LDA, VL, VU, IL,
     $             IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK,
     $             LIWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      LIWORK = MIN( LWMAX, IWORK( 1 ) )
!
!     Solve eigenproblem.
!
      CALL DSYEVR( 'Vectors', 'A', 'Upper', N, A, LDA, VL, VU, IL,
     $             IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK,
     $             LIWORK, INFO )
!
!     Check for convergence.
!
      open(unit=13,file="lapack_eig.err")
      IF( INFO.GT.0 ) THEN
         WRITE(13,*)'The algorithm failed to compute eigenvalues.'
         close(13)
         STOP
      END IF
      end subroutine lapack_eig2

      end module f77_subroutine
!========================================================================================
!========================================================================================






!=========================================================================================
!=========================================================================================
      module unit_conversion
        implicit none
        contains       
!******************Convert Debye to atomic units*************************
      SUBROUTINE Dipole_in_atomic_unit(Dipole)
      	IMPLICIT NONE
      	DOUBLE PRECISION :: Dipole
      	Dipole=Dipole*0.3934302014076827
      END SUBROUTINE

!*****************Convert meter to atomic unit****************************
      SUBROUTINE Length_in_Bohr(Lattice_Constant)
      	IMPLICIT NONE
      	DOUBLE PRECISION :: lattice_constant
      	Lattice_Constant=Lattice_Constant*1.889726133921252D10
      END SUBROUTINE	

!*************Convert from energy atomic unit  to kHz**********************
      SUBROUTINE Energy_in_kHz(energy_in_atomic_unit)
      	IMPLICIT NONE
      	DOUBLE PRECISION :: energy_in_atomic_unit
        energy_in_atomic_unit=energy_in_atomic_unit*6.57968392072144D12
      END SUBROUTINE     

!*************Convert from energy atomic unit  to kHz**********************
      SUBROUTINE Energy_in_MHz(energy_in_atomic_unit)
      	IMPLICIT NONE
      	DOUBLE PRECISION :: energy_in_atomic_unit
        energy_in_atomic_unit=energy_in_atomic_unit*6.57968392072144D9
      END SUBROUTINE 
 
!*****************Convert meter to atomic unit****************************
      SUBROUTINE Field_in_atomic_unit(electric_field)
      	IMPLICIT NONE
      	DOUBLE PRECISION :: electric_field
      	electric_field=electric_field*1.944690567144141D-12
      END SUBROUTINE	
 

!*****************Convert Hz to atomic unit****************************
      SUBROUTINE Hz_to_atomic_unit(RotationalConstant)
      	IMPLICIT NONE
      	DOUBLE PRECISION :: RotationalConstant
      	RotationalConstant=RotationalConstant*2.418884324306202D-17
      END SUBROUTINE	

  
      end module unit_conversion
!==========================================================================================
!==========================================================================================





!==========================================================================================
!==========================================================================================
! contains subroutine New_EigState and Matrix_for_new_EigStates

      module in_electric_field
	    implicit none
        contains
!**********find new rotational eigenstates (in terms of old ones) in an electric field****
      SUBROUTINE New_EigStates(NewEigValue, coefficient_array,
     &                      Dipole, electric_field,
     &                      RotationalConstant,NumberState,
     &                      N, M) 
!  NewEigstate = Sum_i{ coefficient_array(i) * oldEigstate(i) }
!  NumberState is the number of basis sets (how many OldEigstates you want to include to
!                                           express the NewEignstates)	 
! |N,M> is the rotational states
	use f77_subroutine	      
      	IMPLICIT NONE
      	DOUBLE PRECISION :: Dipole, electric_field, RotationalConstant      	
      	INTEGER :: N, M
        DOUBLE PRECISION, ALLOCATABLE :: TotalHamiltonian(:,:)
        DOUBLE PRECISION, ALLOCATABLE :: EigValue(:)
        double precision :: NewEigValue
        DOUBLE PRECISION :: coefficient_array(NumberState)
      	INTEGER, INTENT(IN) ::  NumberState
!assign initial values for the matrix, otherwise it may cause errors in some system
        coefficient_array = 0.d0  
        ALLOCATE( TotalHamiltonian(NumberState,NumberState) ) 
       	ALLOCATE( EigValue(NumberState) )
!     .. Form the matrix ..		
       	call Matrix_for_new_EigStates(TotalHamiltonian,NumberState,N,M,
     &                  Dipole,electric_field,RotationalConstant)

!      .... solve the eigenvalue and eigenvector ...	
        call lapack_eig(TotalHamiltonian, NumberState, EigValue)

!      assign values for output
        NewEigValue=EigValue(N-ABS(M)+1)
! because the number of old states that couple with |N,M> are the following:
! |ABS(M),ABS(M)>, |ABS(M)+1,ABS(M)>, |ABS(M)+2,ABS(M)>, ..., |N,M>, |N+1, ABS(M)>,...
      	coefficient_array(1:NumberState)=TotalHamiltonian(1:NumberState, N-ABS(M)+1)
       		
       	DEALLOCATE(TotalHamiltonian,EigValue)
      END SUBROUTINE New_EigStates 



 	
!***********form matrix for calculating the new rotational eigenstates***************     
      SUBROUTINE Matrix_for_new_EigStates(MatrixName,MatrixDIM,N,M,
     &                      Dipole, electric_field,
     &                      RotationalConstant) !NOTE: MatrixDIM > N
! MatriName: Hamiltonian in matrix form 
! MatrixDIM: the size of MatrixName (MatrixDIM X MatrixDIM)
! |N,M> represent the rotational states	 
	use f77_subroutine
        IMPLICIT NONE
        DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        INTEGER :: I, J, K, L, MatrixDIM, N, M,  NStart, NEndWith
        DOUBLE PRECISION,allocatable :: InteractionHamiltonian(:,:) !the dimension shouldn't exceed 100
        DOUBLE PRECISION,allocatable :: MoleculeHamiltonian(:,:)  ! default: index starts from 1 to 100
        DOUBLE PRECISION :: MatrixName(1:MatrixDIM,1:MatrixDIM)
        DOUBLE PRECISION :: Dipole, electric_field 
        DOUBLE PRECISION :: RotationalConstant

	NStart=ABS(M)
	NEndWith=ABS(M) + MatrixDIM -1
	    
	allocate( InteractionHamiltonian(NStart : NEndWith, NStart : NEndWith) )
	allocate( MoleculeHamiltonian(NStart : NEndWith, NStart : NEndWith) )
		
	    IF ( N < ABS(M) ) THEN 
	    	WRITE(*,*) "Wrong! N should be larger than or equal to M."
	    ELSE IF (MatrixDIM <= N) THEN
	           WRITE(*,*) "Wrong! Matrix Dimension should be larger than N"
	    ELSE       
            DO J=NStart,NEndWith  
!           form the matrix for interaction Hamiltonian
!           for the expression for the matrix elements, see Chapter 2 of <<Cold molecule>>
                  DO I=NStart,J
                       InteractionHamiltonian(I,J)=(-1)*Dipole*(-1)**(M)
     &                 *electric_field*DSQRT((2.d0*I+1.d0)*(2.d0*J+1.d0))               
     &                 *THRJ(DFLOAT(I),1.D0,DFLOAT(J),-1.d0*DFLOAT(M),0.D0,DFLOAT(M))
     &                 *THRJ(DFLOAT(I),1.D0,DFLOAT(J),0.D0,0.D0,0.D0)
	                   InteractionHamiltonian(J,I)=InteractionHamiltonian(I,J) ! since the matrix is symmetric
                  END DO
            END DO
 
            DO L=ABS(M),NEndWith 
                  DO K=ABS(M),NEndWith
                    IF (L==K) THEN 
                     MoleculeHamiltonian(L,K)=2*Pi*RotationalConstant*(K+1)*K
! Energy = h*f = atomic_unit(h)*atomic_unit(f) = 2*Pi*atomic_unit(f)
                    ELSE 
	              MoleculeHamiltonian(L,K)=0.D0   
                    END IF
                  END DO
            END DO

	    MatrixName(1:MatrixDIM, 1:MatrixDIM)
     &	    = InteractionHamiltonian(NStart:NEndWith, NStart:NEndWith) 
     &        + MoleculeHamiltonian(NStart:NEndWith, NStart:NEndWith)
	   END IF
      
	  deallocate(InteractionHamiltonian, MoleculeHamiltonian)
	  	  
	  RETURN
      END SUBROUTINE Matrix_for_new_EigStates 
   
      end module in_electric_field
!========================================================================================
!========================================================================================










!=========================================================================================
!=========================================================================================
      module sums
        implicit none
        contains
!*******************calculate the half_sum1 of Cos[mka]/m^3*****************
! k --- wave vector, a --- lattice constant
! number_of_molecule = Number of molecules
      FUNCTION half_sum1(number_of_molecule, ka)
      	IMPLICIT NONE
      	DOUBLE PRECISION :: half_sum1,ka 
      	INTEGER :: number_of_molecule, I
      	half_sum1=0.D0
      	DO I=1, number_of_molecule/2 
           !do not use DFLOAT(I**3), may exceed the limit for integer 
      		half_sum1=half_sum1 + ( DCOS(I*ka) )/( DFLOAT(I)**3 ) 
      	END DO
      	RETURN
      END FUNCTION half_sum1
!*******************calculate the half_sum1 of 1/m^3*****************
      FUNCTION half_sum2(number_of_molecule)
      	IMPLICIT NONE
      	DOUBLE PRECISION :: half_sum2
      	INTEGER :: number_of_molecule, I
      	half_sum2=0.D0
      	DO I=1, number_of_molecule/2
      		half_sum2=half_sum2 + 1.d0/(DFLOAT(I)**3) 
                !do not use DFLOAT(I**3), may exceed the limit for integer
      	END DO
      	RETURN
      END FUNCTION half_sum2
!******************************************************************************

!*******************calculate the half_sum1 of Cos[mka]Cos[mk'a]/m^3*****************
! k --- wave vector, a --- lattice constant
! number_of_molecule = Number of molecules
      FUNCTION half_sum3(number_of_molecule, ka,ka_)
      	IMPLICIT NONE
      	DOUBLE PRECISION :: half_sum3,ka,ka_ 
      	INTEGER :: number_of_molecule, I
      	half_sum3=0.D0
      	DO I=1, number_of_molecule/2 
           !do not use DFLOAT(I**3), may exceed the limit for integer 
      		half_sum3=half_sum3 + ( DCOS(I*ka)*DCOS(I*ka_) )/(DFLOAT(I)**3) 
      	END DO
      	RETURN
      END FUNCTION half_sum3

! ************************************************************************************
      FUNCTION half_sum4(number_of_molecule, ka1,ka2,ka4)
      	IMPLICIT NONE
      	DOUBLE PRECISION :: half_sum4,ka1,ka2,ka4 
      	INTEGER :: number_of_molecule, I
      	half_sum4=0.D0
      	DO I=1, number_of_molecule/2 
           !do not use DFLOAT(I**3), may exceed the limit for integer 
      		half_sum4=half_sum4 + (    DCOS( I*(ka1-ka4) )  + DCOS( I*(ka2-ka4) )    )
     &                 /(DFLOAT(I)**3) 
      	END DO
      	RETURN
      END FUNCTION half_sum4

      END MODULE sums
!=========================================================================================
!=========================================================================================












!==========================================================================================
!==========================================================================================
      module dipole_dipole
        implicit none
        contains		
! please note that function in the same module knows each other
! so there is no need to specify them as external when they are being used		
!***************************************************************************      
      FUNCTION AminusBminus(NA,MA,NB,MB,NA_,MA_,NB_,MB_)
	    USE f77_subroutine
      	IMPLICIT NONE
      	DOUBLE PRECISION :: AminusBminus
      	INTEGER :: NA,MA,NB,MB,NA_,MA_,NB_,MB_
      	AminusBminus=THRJ(1.D0,1.D0,2.D0,-1.D0,-1.D0,2.D0)*
     & THRJ(DFLOAT(NA),1.D0,DFLOAT(NA_),-DFLOAT(MA),-1.D0,DFLOAT(MA_))*
     & THRJ(DFLOAT(NB),1.D0,DFLOAT(NB_),-DFLOAT(MB),-1.D0,DFLOAT(MB_))
        RETURN
      END  FUNCTION AminusBminus
 

!****************************************************************************      	      
      FUNCTION AplusBplus(NA,MA,NB,MB,NA_,MA_,NB_,MB_)
	    USE f77_subroutine
      	IMPLICIT NONE
      	DOUBLE PRECISION :: AplusBplus
      	INTEGER :: NA,MA,NB,MB,NA_,MA_,NB_,MB_
      	AplusBplus=THRJ(1.D0,1.D0,2.D0,1.D0,1.D0,-2.D0)*
      	
     & THRJ(DFLOAT(NA),1.D0,DFLOAT(NA_),-DFLOAT(MA),1.D0,DFLOAT(MA_))*
     
     & THRJ(DFLOAT(NB),1.D0,DFLOAT(NB_),-DFLOAT(MB),1.D0,DFLOAT(MB_))
     
        RETURN
      END  FUNCTION  AplusBplus    

!****************************************************************************      
      FUNCTION AplusBminus(NA,MA,NB,MB,NA_,MA_,NB_,MB_)
	    USE f77_subroutine
      	IMPLICIT NONE
      	DOUBLE PRECISION :: AplusBminus
      	INTEGER :: NA,MA,NB,MB,NA_,MA_,NB_,MB_
      	AplusBminus=THRJ(1.D0,1.D0,2.D0,1.D0,-1.D0,0.D0)*
      	
     & THRJ(DFLOAT(NA),1.D0,DFLOAT(NA_),-DFLOAT(MA),1.D0,DFLOAT(MA_))*
     
     & THRJ(DFLOAT(NB),1.D0,DFLOAT(NB_),-DFLOAT(MB),-1.D0,DFLOAT(MB_))
     
        RETURN
      END  FUNCTION AplusBminus     

!****************************************************************************      
      FUNCTION AminusBplus(NA,MA,NB,MB,NA_,MA_,NB_,MB_)
	    USE f77_subroutine
      	IMPLICIT NONE
      	DOUBLE PRECISION :: AminusBplus
      	INTEGER :: NA,MA,NB,MB,NA_,MA_,NB_,MB_
      	AminusBplus=THRJ(1.D0,1.D0,2.D0,-1.D0,1.D0,0.D0)*
      	
     & THRJ(DFLOAT(NA),1.D0,DFLOAT(NA_),-DFLOAT(MA),-1.D0,DFLOAT(MA_))*
     
     & THRJ(DFLOAT(NB),1.D0,DFLOAT(NB_),-DFLOAT(MB),1.D0,DFLOAT(MB_))
     
        RETURN
      END   FUNCTION  AminusBplus   


!****************************************************************************      
      FUNCTION A0B0(NA,MA,NB,MB,NA_,MA_,NB_,MB_)
	    USE f77_subroutine
      	IMPLICIT NONE
      	DOUBLE PRECISION :: A0B0
      	INTEGER :: NA,MA,NB,MB,NA_,MA_,NB_,MB_
      	A0B0=THRJ(1.D0,1.D0,2.D0,0.D0,0.D0,0.D0)*
      	
     & THRJ(DFLOAT(NA),1.D0,DFLOAT(NA_),-DFLOAT(MA),0.D0,DFLOAT(MA_))*
     
     & THRJ(DFLOAT(NB),1.D0,DFLOAT(NB_),-DFLOAT(MB),0.D0,DFLOAT(MB_))
        RETURN
      END   FUNCTION  A0B0   

!****************************************************************************      
      FUNCTION A0Bminus(NA,MA,NB,MB,NA_,MA_,NB_,MB_)

Cf2py intent(in) NA
Cf2py intent(in) MA
Cf2py intent(in) NB
Cf2py intent(in) MB
Cf2py intent(in) NA_
Cf2py intent(in) MA_
Cf2py intent(in) NB_
Cf2py intent(in) MB_
        use f77_subroutine
	    IMPLICIT NONE
      	DOUBLE PRECISION :: A0Bminus
      	INTEGER :: NA,MA,NB,MB,NA_,MA_,NB_,MB_
!        DOUBLE PRECISION, external :: thrj
      	A0Bminus=THRJ(1.D0,1.D0,2.D0,0.D0,-1.D0,1.D0)*
      	
     & THRJ(DFLOAT(NA),1.D0,DFLOAT(NA_),-DFLOAT(MA),0.D0,DFLOAT(MA_))*
     
     & THRJ(DFLOAT(NB),1.D0,DFLOAT(NB_),-DFLOAT(MB),-1.D0,DFLOAT(MB_))
     
        RETURN
      END  FUNCTION  A0Bminus    

!****************************************************************************      
      FUNCTION A0Bplus(NA,MA,NB,MB,NA_,MA_,NB_,MB_)

Cf2py intent(in) NA
Cf2py intent(in) MA
Cf2py intent(in) NB
Cf2py intent(in) MB
Cf2py intent(in) NA_
Cf2py intent(in) MA_
Cf2py intent(in) NB_
Cf2py intent(in) MB_
        use f77_subroutine
	    IMPLICIT NONE
      	DOUBLE PRECISION :: A0Bplus
      	INTEGER :: NA,MA,NB,MB,NA_,MA_,NB_,MB_
!        DOUBLE PRECISION, external :: thrj
      	A0Bplus=THRJ(1.D0,1.D0,2.D0,0.D0,1.D0,-1.D0)*
      	
     & THRJ(DFLOAT(NA),1.D0,DFLOAT(NA_),-DFLOAT(MA),0.D0,DFLOAT(MA_))*
     
     & THRJ(DFLOAT(NB),1.D0,DFLOAT(NB_),-DFLOAT(MB),1.D0,DFLOAT(MB_))
     
        RETURN
      END  FUNCTION A0Bplus     
      
!****************************************************************************      
      FUNCTION AminusB0(NA,MA,NB,MB,NA_,MA_,NB_,MB_)

Cf2py intent(in) NA
Cf2py intent(in) MA
Cf2py intent(in) NB
Cf2py intent(in) MB
Cf2py intent(in) NA_
Cf2py intent(in) MA_
Cf2py intent(in) NB_
Cf2py intent(in) MB_
        use f77_subroutine
	    IMPLICIT NONE
      	DOUBLE PRECISION :: AminusB0
      	INTEGER :: NA,MA,NB,MB,NA_,MA_,NB_,MB_
!        DOUBLE PRECISION, external :: thrj
      	AminusB0=THRJ(1.D0,1.D0,2.D0,-1.D0,0.D0,1.D0)*
      	
     & THRJ(DFLOAT(NA),1.D0,DFLOAT(NA_),-DFLOAT(MA),-1.D0,DFLOAT(MA_))*
     
     & THRJ(DFLOAT(NB),1.D0,DFLOAT(NB_),-DFLOAT(MB),0.D0,DFLOAT(MB_))
     
        RETURN
      END  FUNCTION AminusB0

!****************************************************************************      
      FUNCTION AplusB0(NA,MA,NB,MB,NA_,MA_,NB_,MB_)

Cf2py intent(in) NA
Cf2py intent(in) MA
Cf2py intent(in) NB
Cf2py intent(in) MB
Cf2py intent(in) NA_
Cf2py intent(in) MA_
Cf2py intent(in) NB_
Cf2py intent(in) MB_
        use f77_subroutine
        IMPLICIT NONE
      	DOUBLE PRECISION :: AplusB0
      	INTEGER :: NA,MA,NB,MB,NA_,MA_,NB_,MB_
!        DOUBLE PRECISION, external :: thrj
      	AplusB0=THRJ(1.D0,1.D0,2.D0,1.D0,0.D0,-1.D0)*
      	
     & THRJ(DFLOAT(NA),1.D0,DFLOAT(NA_),-DFLOAT(MA),1.D0,DFLOAT(MA_))*
     
     & THRJ(DFLOAT(NB),1.D0,DFLOAT(NB_),-DFLOAT(MB),0.D0,DFLOAT(MB_))
     
        RETURN
      END  FUNCTION AplusB0
! all the above A+/-/0 B+/-/0 functions has been checked with my notes


      
      
!*******************************dipole-dipole interaction***********************************
! the electric field is perpendicular to the intermolecular axis
! calculate the dipole-dipole interaction between two molecules: 
!  < NA,MA|<NB,MB| V |NA',MA'>|NB',MB'> 
! the rotational states are eigenstates of the Hamiltonian when there is no external field    
      FUNCTION DDInteraction(Dipole,lattice_constant,NA,MA,NB,MB,NA_,MA_,NB_,MB_)
	    use f77_subroutine
      	IMPLICIT NONE
        DOUBLE PRECISION :: DDInteraction
!      	DOUBLE PRECISION, EXTERNAL :: AminusBminus, AplusBplus, 
!     & 	                              AplusBminus, AminusBplus, A0B0
        INTEGER :: NA,MA,NB,MB,NA_,MA_,NB_,MB_

        DOUBLE PRECISION :: Dipole, Lattice_Constant                             

! see my notes for the expression        
        DDInteraction=-3.D0*DSQRT(5.D0)*( (-1)**((-1)*(MA+MB)) )
        
     & *( (Dipole*Dipole) / (Lattice_Constant)**(3) ) 
     
     & *DSQRT( (2.D0*NA+1.D0)*(2.D0*NA_+1.D0)*(2.D0*NB+1.D0)*(2.D0*NB_+1.D0) )
     
     & *THRJ(DFLOAT(NA), 1.D0,DFLOAT(NA_),0.D0,0.D0,0.D0)
     
     & *THRJ(DFLOAT(NB), 1.D0,DFLOAT(NB_),0.D0,0.D0,0.D0)
     
     & *( 0.5*AplusBplus(NA,MA,NB,MB,NA_,MA_,NB_,MB_) 
     
     &    + 0.5*AminusBminus(NA,MA,NB,MB,NA_,MA_,NB_,MB_)
     
     &    - DSQRT(1.D0/6.D0)*AplusBminus(NA,MA,NB,MB,NA_,MA_,NB_,MB_)
     
     &    - DSQRT(1.D0/6.D0)*AminusBplus(NA,MA,NB,MB,NA_,MA_,NB_,MB_)
     
     &    - DSQRT(1.D0/6.D0)*A0B0(NA,MA,NB,MB,NA_,MA_,NB_,MB_) )
      
        RETURN
      END  FUNCTION DDInteraction

!***************************calculate matrix elements in electric field*********************
! the electric field is perpendicular to the intermolecular axis
! dipole-dipole interaction between two molecules in electric field
! <NA1,MA1|<NB1,MB1| V |NA2,MA2>|NB2,MB2>
! the states are dressed states, each is a linear combination of the original rotational states
      FUNCTION DDInteract_in_Field(Dipole,electric_field,RotationalConstant,Lattice_Constant,
     &                             NA1,MA1,NB1,MB1,NA2,MA2,NB2,MB2)  !<B-ground-A-excited|V|A-ground-B-excited>
	    USE in_electric_field ! including New_EigStates 
      	IMPLICIT NONE
	INTEGER, PARAMETER :: number_state=50
      	INTEGER :: NA1, MA1, NB1, MB1, NA2,MA2,NB2,MB2
      	INTEGER :: I, J, K, L
      	INTEGER :: final_state_number1, final_state_number2
      	INTEGER :: final_state_number3, final_state_number4
      	      	
      	DOUBLE PRECISION :: Dipole, electric_field, RotationalConstant, Lattice_Constant
      	DOUBLE PRECISION :: DDInteract_in_Field
!      	DOUBLE PRECISION, EXTERNAL :: DDInteraction
        DOUBLE PRECISION :: coefficient_array1(number_state)
        DOUBLE PRECISION :: coefficient_array2(number_state)
        DOUBLE PRECISION :: coefficient_array3(number_state)
        DOUBLE PRECISION :: coefficient_array4(number_state)
        DOUBLE PRECISION :: NewEigValue1,NewEigValue2,NewEigValue3,NewEigValue4
		
		final_state_number1=number_state
		final_state_number2=number_state
		final_state_number3=number_state
		final_state_number4=number_state
		
		
        CALL New_EigStates(NewEigValue1,coefficient_array1,
     &                      Dipole, electric_field,
     &                      RotationalConstant,final_state_number1,
     &                      NA1, MA1)
     
!        Write(40,*) coefficient_array1(1:final_state_number1)
        
        CALL New_EigStates(NewEigValue2,coefficient_array2,
     &                      Dipole, electric_field,
     &                      RotationalConstant,final_state_number2,
     &                      NB1, MB1)
     
        CALL New_EigStates(NewEigValue3,coefficient_array3,
     &                      Dipole, electric_field,
     &                      RotationalConstant,final_state_number3,
     &                      NA2, MA2)
    
        
        CALL New_EigStates(NewEigValue4,coefficient_array4,
     &                      Dipole, electric_field,
     &                      RotationalConstant,final_state_number4,
     &                      NB2, MB2)     
        DDInteract_in_Field=0.D0
! the coefficient_array has been normalized
        DO I=1, final_state_number1
        	DO J=1, final_state_number2
        		DO K=1, final_state_number3
        			DO L=1, final_state_number4
        				DDInteract_in_Field=DDInteract_in_Field + coefficient_array1(I)*coefficient_array2(J)
     &   	                                            *coefficient_array3(K)*coefficient_array4(L)	
     &   					            *DDInteraction(Dipole, Lattice_Constant,ABS(MA1)+I-1,MA1,ABS(MB1)+J-1,MB1,    
     &                                                                     ABS(MA2)+K-1,MA2,ABS(MB2)+L-1,MB2)
                    END DO
                 END DO   
            END DO
        END DO
        
        RETURN
      END  FUNCTION DDInteract_in_Field
!************************************************************************************************** 

!*******************************dipole-dipole interaction with theta and phi************************     
! assume the molecular axis is along the space-fixed x-axis
! theta is the angle between the direction of electric field and its projection in x-y plane
! phi is the angle between the projection of electric field and the x axis
! for the one dimensional case, we can assume phi = 0, for the case where electric field
! is perpendicular to intermoleculear axis, theta = pi/2; for the case where electric field 
! is parallel to intermolecular axis, theta = 0

      FUNCTION DDInteraction_angle(Dipole,lattice_constant,
     &                    NA,MA,NB,MB,NA_,MA_,NB_,MB_,theta,phi)

Cf2py intent(in) Dipole
Cf2py intent(in) lattice_constant
Cf2py intent(in) NA
Cf2py intent(in) MA
Cf2py intent(in) NB
Cf2py intent(in) MB
Cf2py intent(in) NA_
Cf2py intent(in) MA_
Cf2py intent(in) NB_
Cf2py intent(in) MB_
Cf2py intent(in) theta
Cf2py intent(in) phi
        use f77_subroutine
        IMPLICIT NONE
        DOUBLE complex :: DDInteraction_angle
! in fact, DDInteraction_angle is always real if we only consider M=0
! for M= +1/-1, I am not sure
!      	DOUBLE PRECISION, EXTERNAL :: AminusBminus, AplusBplus, 
!     & 	                              AplusBminus, AminusBplus, A0B0,
!     &                                AplusB0,AminusB0,A0Bplus,A0Bminus
        INTEGER :: NA,MA,NB,MB,NA_,MA_,NB_,MB_

        DOUBLE PRECISION :: Dipole, Lattice_Constant,prefactor  
        double precision :: theta, phi 
        double complex :: imag_unit
        
        imag_unit = (0.D0,1.D0)                          
       
        prefactor=-DSQRT(5.D0)*( (-1)**(MA+MB) )        
     & *( (Dipole*Dipole) / (2.D0*Lattice_Constant**3) )      
     & *DSQRT( (2.D0*NA+1.D0)*(2.D0*NA_+1.D0)*(2.D0*NB+1.D0)*(2.D0*NB_+1.D0) )     
     & *THRJ(DFLOAT(NA), 1.D0,DFLOAT(NA_),0.D0,0.D0,0.D0)     
     & *THRJ(DFLOAT(NB), 1.D0,DFLOAT(NB_),0.D0,0.D0,0.D0) 

        DDInteraction_angle = prefactor    
     & *(   3.D0* (DSIN(theta)**2) *exp(-2.D0*phi*imag_unit)*AminusBminus(NA,MA,NB,MB,NA_,MA_,NB_,MB_)      
     &    -6.D0*dsin(theta)*dcos(theta)*exp(-phi*imag_unit)*
     &         ( A0Bminus(NA,MA,NB,MB,NA_,MA_,NB_,MB_) + AminusB0(NA,MA,NB,MB,NA_,MA_,NB_,MB_) )    
     &    + DSQRT(6.D0)* (3*dcos(theta)**2 -1) *(  AplusBminus(NA,MA,NB,MB,NA_,MA_,NB_,MB_)     
     &         + AminusBplus(NA,MA,NB,MB,NA_,MA_,NB_,MB_) +  A0B0(NA,MA,NB,MB,NA_,MA_,NB_,MB_)  )  
     &    + 6.D0*dsin(theta)*dcos(theta)*exp(phi*imag_unit)*
     &         ( A0Bplus(NA,MA,NB,MB,NA_,MA_,NB_,MB_) + AplusB0(NA,MA,NB,MB,NA_,MA_,NB_,MB_) )
     &    + 3.D0*(dsin(theta)**2)*exp(2.D0*phi*imag_unit)*AplusBplus(NA,MA,NB,MB,NA_,MA_,NB_,MB_) 
     &     )

        RETURN
      END  FUNCTION DDInteraction_angle

!*******************************dipole-dipole interaction with theta************************ 
! This is for one-dimensional case, phi = 0
! for the special case where MA=0, MB=0, MA'=0, and MB'=0    
      FUNCTION DDInteraction_angle2(Dipole,lattice_constant,
     &                    NA,MA,NB,MB,NA_,MA_,NB_,MB_,theta)

Cf2py intent(in) Dipole
Cf2py intent(in) lattice_constant
Cf2py intent(in) NA
Cf2py intent(in) MA
Cf2py intent(in) NB
Cf2py intent(in) MB
Cf2py intent(in) NA_
Cf2py intent(in) MA_
Cf2py intent(in) NB_
Cf2py intent(in) MB_
Cf2py intent(in) theta
        use f77_subroutine
        use in_electric_field ! including New_EigStates 
        IMPLICIT NONE
        DOUBLE precision :: DDInteraction_angle2
!      	DOUBLE PRECISION, EXTERNAL :: AminusBminus, AplusBplus, 
!     & 	                              AplusBminus, AminusBplus, A0B0,
!     &                                AplusB0,AminusB0,A0Bplus,A0Bminus
        INTEGER :: NA,MA,NB,MB,NA_,MA_,NB_,MB_

        DOUBLE PRECISION :: Dipole, Lattice_Constant,prefactor  
        double precision :: theta 
                         
       
        prefactor=-DSQRT(5.D0)*( (-1)**(MA+MB) )        
     & *( (Dipole*Dipole) / (2.D0*Lattice_Constant**3) )      
     & *DSQRT( (2.D0*NA+1.D0)*(2.D0*NA_+1.D0)*(2.D0*NB+1.D0)*(2.D0*NB_+1.D0) )     
     & *THRJ(DFLOAT(NA), 1.D0,DFLOAT(NA_),0.D0,0.D0,0.D0)     
     & *THRJ(DFLOAT(NB), 1.D0,DFLOAT(NB_),0.D0,0.D0,0.D0) 

        DDInteraction_angle2 = prefactor *       
     &     DSQRT(6.D0)* (3*dcos(theta)**2 -1) *A0B0(NA,MA,NB,MB,NA_,MA_,NB_,MB_) 
   
        RETURN
! an alternative expression is
! DDInteraction_angle2 = (1-3*DCOS(theta)**2)*Dipole*Dipole/(Lattice_Constant**3)
!                        *DSQRT((2*NA + 1)*(2*NA_+1)*(2*NB+1)*(2*NB_+1))
!                        *THRJ(NA, 1, NA_, 0, 0, 0)**2
!                        *THRJ(NB, 1, NB_, 0, 0, 0)**2
      END  FUNCTION DDInteraction_angle2






!**********calculate matrix elements in electric field with theta and phi****************
! require two angles
      FUNCTION DDInteract_in_Field_angle(Dipole,electric_field,rot_const,Lattice_Constant,
     &                             NA1,MA1, NB1,MB1, NA2,MA2, NB2,MB2, theta,phi)  
!<B-ground-A-excited|V|A-ground-B-excited>

Cf2py intent(in) Dipole
Cf2py intent(in) electric_field
Cf2py intent(in) rot_const
Cf2py intent(in) Lattice_Constant
Cf2py intent(in) NA1
Cf2py intent(in) MA1
Cf2py intent(in) NB1
Cf2py intent(in) MB1
Cf2py intent(in) NA2
Cf2py intent(in) MA2
Cf2py intent(in) NB2
Cf2py intent(in) MB2
Cf2py intent(in) theta
Cf2py intent(in) phi
        use in_electric_field ! including New_EigStates 
        IMPLICIT NONE
        INTEGER, PARAMETER :: number_state= 50 !100 !30
      	INTEGER :: NA1, MA1,  NB1, MB1,  NA2,MA2,  NB2,MB2
      	INTEGER :: I, J, K, L
      	INTEGER :: final_state_number1, final_state_number2
      	INTEGER :: final_state_number3, final_state_number4      	      	
      	DOUBLE PRECISION :: Dipole, electric_field, rot_const 
        double precision :: Lattice_Constant
      	DOUBLE complex :: DDInteract_in_Field_angle ! for M=0, it is real
        DOUBLE PRECISION :: coefficient_array1(number_state)
        DOUBLE PRECISION :: coefficient_array2(number_state)
        DOUBLE PRECISION :: coefficient_array3(number_state)
        DOUBLE PRECISION :: coefficient_array4(number_state)
        DOUBLE PRECISION :: NewEigValue1,NewEigValue2,NewEigValue3,
     &                       NewEigValue4
        double precision :: theta, phi
!      	DOUBLE complex, EXTERNAL :: DDInteraction_angle
!        external :: New_EigStates 
		
        final_state_number1=number_state
        final_state_number2=number_state
        final_state_number3=number_state
        final_state_number4=number_state

		
        CALL New_EigStates(NewEigValue1,coefficient_array1,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number1,
     &                      NA1, MA1)
     
        
        CALL New_EigStates(NewEigValue2,coefficient_array2,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number2,
     &                      NB1, MB1)
     
        CALL New_EigStates(NewEigValue3,coefficient_array3,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number3,
     &                      NA2, MA2)
    
        
        CALL New_EigStates(NewEigValue4,coefficient_array4,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number4,
     &                      NB2, MB2) 
    
        DDInteract_in_Field_angle=(0.D0,0.D0)
        DO I=1, final_state_number1
          DO J=1, final_state_number2
            DO K=1, final_state_number3
              DO L=1, final_state_number4
        		DDInteract_in_Field_angle=DDInteract_in_Field_angle + 
     &               coefficient_array1(I)*coefficient_array2(J)
     &   	        *coefficient_array3(K)*coefficient_array4(L)	
     &              *DDInteraction_angle( Dipole, Lattice_Constant,
     &                         ABS(MA1)+I-1,MA1,ABS(MB1)+J-1,MB1,    
     &                         ABS(MA2)+K-1,MA2,ABS(MB2)+L-1,MB2,theta,phi )
                END DO
              END DO   
            END DO
          END DO
        
        RETURN
      END  FUNCTION DDInteract_in_Field_angle
!*******************************************************************************************


!***************************calculate matrix elements in electric field*********************
! this is for one-dimensional case, phi = 0
! for the special case where MA1=0, MB1=0, MA2=0, and MB2=0 
! only require one angle
      FUNCTION DDInteract_in_Field_angle2(Dipole,electric_field,rot_const,Lattice_Constant,
     &                             NA1,MA1,NB1,MB1,NA2,MA2,NB2,MB2,theta)  
!<B-ground-A-excited|V|A-ground-B-excited>

Cf2py intent(in) Dipole
Cf2py intent(in) electric_field
Cf2py intent(in) rot_const
Cf2py intent(in) Lattice_Constant
Cf2py intent(in) NA1
Cf2py intent(in) MA1
Cf2py intent(in) NB1
Cf2py intent(in) MB1
Cf2py intent(in) NA2
Cf2py intent(in) MA2
Cf2py intent(in) NB2
Cf2py intent(in) MB2
Cf2py intent(in) theta
Cf2py intent(in) phi
        use in_electric_field
        IMPLICIT NONE
        INTEGER, PARAMETER :: number_state= 50 !100 !30
      	INTEGER :: NA1, MA1, NB1, MB1, NA2,MA2,NB2,MB2
      	INTEGER :: I, J, K, L
      	INTEGER :: final_state_number1, final_state_number2
      	INTEGER :: final_state_number3, final_state_number4      	      	
      	DOUBLE PRECISION :: Dipole, electric_field, rot_const 
        double precision :: Lattice_Constant
      	DOUBLE precision :: DDInteract_in_Field_angle2
        DOUBLE PRECISION :: coefficient_array1(number_state)
        DOUBLE PRECISION :: coefficient_array2(number_state)
        DOUBLE PRECISION :: coefficient_array3(number_state)
        DOUBLE PRECISION :: coefficient_array4(number_state)
        DOUBLE PRECISION :: NewEigValue1,NewEigValue2,NewEigValue3,
     &                       NewEigValue4
        double precision :: theta
!      	DOUBLE precision, EXTERNAL :: DDInteraction_angle2
!        external :: New_EigStates 
		
        final_state_number1=number_state
        final_state_number2=number_state
        final_state_number3=number_state
        final_state_number4=number_state
		
		
        CALL New_EigStates(NewEigValue1,coefficient_array1,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number1,
     &                      NA1, MA1)
     
        
        CALL New_EigStates(NewEigValue2,coefficient_array2,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number2,
     &                      NB1, MB1)
     
        CALL New_EigStates(NewEigValue3,coefficient_array3,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number3,
     &                      NA2, MA2)
    
        
        CALL New_EigStates(NewEigValue4,coefficient_array4,
     &                      Dipole, electric_field,
     &                      rot_const,final_state_number4,
     &                      NB2, MB2) 
    
        DDInteract_in_Field_angle2=0.D0
        DO I=1, final_state_number1
          DO J=1, final_state_number2
            DO K=1, final_state_number3
              DO L=1, final_state_number4
        		DDInteract_in_Field_angle2=DDInteract_in_Field_angle2 + 
     &               coefficient_array1(I)*coefficient_array2(J)
     &   	        *coefficient_array3(K)*coefficient_array4(L)	
     &              *DDInteraction_angle2( Dipole, Lattice_Constant,
     &                         ABS(MA1)+I-1,MA1,ABS(MB1)+J-1,MB1,    
     &                         ABS(MA2)+K-1,MA2,ABS(MB2)+L-1,MB2,theta)
                END DO
              END DO   
            END DO
          END DO
        
        RETURN
      END  FUNCTION DDInteract_in_Field_angle2
!*******************************************************************************************     
      end module dipole_dipole 
!==========================================================================================
!==========================================================================================














!=========================================================================================
!=========================================================================================
      module exciton
        implicit none
        contains
!*************************************************************************
! the input electric_field is in atomic unit
! the output of exciton10_energy is kHz
! this is the dispersion energy (Excitonic energy vs wavevector)
! exciton10_energy is an array of excitonic energies for different ka (from 0 to pi)
! exciton10_energy is the output 
! only for theta = pi/2
      Subroutine Exciton10_dispersion(electric_field,number_of_molecule,exciton10_energy)
        use unit_conversion
        use dipole_dipole
        use f77_subroutine
        use sums
      	IMPLICIT NONE
      	DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        double precision :: interaction_energy
        ! ka =0, 2pi/N, 2*(2pi/N), 3*(2pi/N), ..., (N/2)*(2pi/N)
        ! so the dimension of the matrix is N/2 + 1
        double precision :: exciton10_energy(number_of_molecule/2 + 1)  
      	INTEGER ::  I 
      	INTEGER :: number_of_molecule
      	DOUBLE PRECISION :: RotationalConstant,electric_field
        double precision :: summation(number_of_molecule/2 + 1)
      	DOUBLE PRECISION :: Dipole, Lattice_Constant
        double precision :: ka(number_of_molecule/2 + 1)
      	RotationalConstant=11.7998D9/2.D0  !unit=Hz 1GHz=10^9 Hz
      	Dipole=5.529D0  !Debye
        Lattice_Constant=4.D-7  !400nm       
        
        CALL Dipole_in_atomic_unit(Dipole)
        CALL Length_in_Bohr(Lattice_Constant)
        CALL Hz_to_atomic_unit(RotationalConstant)

        interaction_energy=DDInteract_in_Field(Dipole,electric_field,
     &                 RotationalConstant,Lattice_Constant,1,0,0,0,0,0,1,0)                
        
        call Energy_in_kHz(interaction_energy) !convert atomic unit to kHz 
        do i=1, number_of_molecule/2 + 1
          ka(i) = (i-1)*pi/(number_of_molecule/2) 
          summation(i)=2*half_sum1(number_of_molecule, ka(i)) 
        end do

        exciton10_energy=interaction_energy*summation 
              
      END SUBROUTINE Exciton10_dispersion


!*************************************************************************
! the electric_field is in atomic unit
! the output exciton20_energy is in atomic unit
! exciton20_energy is the dispersion energy (purely excitonic energy)
! exciton20_energy is an array of excitonic energies for different ka (from 0 to pi)
! only for theta = pi/2 
      Subroutine Exciton20_dispersion(electric_field,number_of_molecule,exciton20_energy)
        use unit_conversion
        use dipole_dipole
        use f77_subroutine
        use sums
      	IMPLICIT NONE
      	DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        double precision :: interaction_energy
        ! ka =0, 2pi/N, 2*(2pi/N), 3*(2pi/N), ..., (N/2)*(2pi/N)
        double precision :: exciton20_energy(number_of_molecule/2 + 1)  
      	INTEGER ::  I 
      	INTEGER :: number_of_molecule
      	DOUBLE PRECISION :: RotationalConstant,electric_field
        double precision :: summation(number_of_molecule/2 + 1)
      	DOUBLE PRECISION :: Dipole, Lattice_Constant
        double precision :: ka(number_of_molecule/2 + 1)
      	RotationalConstant=11.7998D9/2.D0  !unit=Hz 1GHz=10^9 Hz
      	Dipole=5.529D0 !Debye
        Lattice_Constant=4.D-7 !400nm        
        
        CALL Dipole_in_atomic_unit(Dipole)
        CALL Length_in_Bohr(Lattice_Constant)
        CALL Hz_to_atomic_unit(RotationalConstant)

        
        interaction_energy=DDInteract_in_Field(Dipole,electric_field,
     &                 RotationalConstant,Lattice_Constant,2,0,0,0,0,0,2,0)                
        call Energy_in_kHz(interaction_energy) !convert atomic unit to kHz
        
        do i=1, number_of_molecule/2 + 1
           ka(i) = (i-1)*pi/(number_of_molecule/2) 
           summation(i)=2*half_sum1(number_of_molecule, ka(i)) 
        end do

       exciton20_energy=interaction_energy*summation 
              
      END SUBROUTINE Exciton20_dispersion

!*****************************calculate the matrix element <f|V|i> as a function of k***********************
!unit of M_vs_k is: kHz  Electric_field :: atomic unit
! only for theta = pi/2
      Subroutine M_vs_k(electric_field, Number_of_molecule,m_k_array) 
        use dipole_dipole
        use unit_conversion
        use sums
      	IMPLICIT NONE
      	DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        ! ka =0, 2pi/N, 2*(2pi/N), 3*(2pi/N), ..., (N/2)*(2pi/N)
      	DOUBLE PRECISION :: m_k_array(Number_of_molecule/2 + 1)
      	DOUBLE PRECISION :: ka(Number_of_molecule/2 + 1)
        double precision :: summation(Number_of_molecule/2 + 1)
      	double precision :: matrix_element
      	INTEGER :: Number_of_molecule,i
      	DOUBLE PRECISION :: RotationalConstant,electric_field
      	DOUBLE PRECISION :: Dipole, Lattice_Constant

      	RotationalConstant=11.7998D9/2.D0  !unit=Hz 1GHz=10^9 Hz
      	Dipole=5.529D0 !Debye
        Lattice_Constant=4.D-7 !400nm
                   
        CALL Dipole_in_atomic_unit(Dipole)
        CALL Length_in_Bohr(Lattice_Constant)
        CALL Hz_to_atomic_unit(RotationalConstant)
     
      	matrix_element = DDInteract_in_Field(Dipole,electric_field,RotationalConstant,Lattice_Constant,
     &                             2,0,0,0,1,0,1,0)
        CALL Energy_in_kHz(matrix_element)

        do i=1, number_of_molecule/2 + 1
           ka(i) = (i-1)*pi/(number_of_molecule/2)
           summation(i)=2*half_sum1(number_of_molecule, ka(i)) 
        end do

        m_k_array=matrix_element*summation
        return            
      END subroutine M_vs_k
    
!***************************calculate E_{k1} + E_{k2} - E_{k3} excluding the excitonic energy*******************
! the energy difference between dressed states: (E(|N=1,M=0>) + D1) + (E(|N=1,M=0>) + D1) - (E(|N=2,M=0>) + D2)
! where D1 and D2 are gas-condensed matter shift for |10> and |20> respectively
! D1 = sum_{m\=n}[<10|Vnm|10> - <00|Vnm|00> ] ; D2 = sum_{m\=n}[<20|Vnm|20> - <00|Vnm|00> ] 
! the return energy is in the unit of kHz
      Function delta_energy_exclude(electric_field, number_of_molecule) !unit of electric field: atomic unit 
        use unit_conversion
        implicit none
        double precision :: electric_field, delta_energy_exclude
        integer :: number_of_molecule
! E1G = E(|10>) - E(|00>); E2G = E(|20>) - E(|10>) here |00>, |10>, |20> are dressed states
        delta_energy_exclude= 2*( E1G(electric_field) + GPSA(electric_field,Number_of_molecule, 1) )  
     &                       - E2G(electric_field) - GPSA(electric_field,Number_of_molecule, 2)  
!       here delta_energy is in atomic unit
        CALL Energy_in_kHz(delta_energy_exclude)
        return ! return energy in kHz
      end function delta_energy_exclude

!***************************calculate E_{k1} + E_{k2} - E_{k3} excluding the excitonic energy*******************
! the energy difference between dressed states: (E(|N=1,M=0>) + D1) + (E(|N=1,M=0>) + D1) - (E(|N=2,M=0>) + D2)
! where D1 and D2 are gas-condensed matter shift for |10> and |20> respectively
! D1 = sum_{m\=n}[<10|Vnm|10> - <00|Vnm|00> ] ; D2 = sum_{m\=n}[<20|Vnm|20> - <00|Vnm|00> ] 
! the return energy is in the unit of kHz
! angle between electric field an the intermolecular axis
      Function delta_energy_ex_angle(electric_field, number_of_molecule, theta) !unit of electric field: atomic unit 
        use unit_conversion
        implicit none
        double precision :: electric_field, delta_energy_ex_angle
        double precision :: theta
        integer :: number_of_molecule
! E1G = E(|10>) - E(|00>); E2G = E(|20>) - E(|10>) here |00>, |10>, |20> are dressed states
        delta_energy_ex_angle= 2*( E1G(electric_field) + GPSA_angle(electric_field,Number_of_molecule,theta, 1) )  
     &                       - E2G(electric_field) - GPSA_angle(electric_field,Number_of_molecule,theta, 2)  
!       here delta_energy is in atomic unit
        CALL Energy_in_kHz(delta_energy_ex_angle)
        return ! return energy in kHz
      end function delta_energy_ex_angle



!*******energy difference between dressed states |10> and |00> at some electric field************
      FUNCTION E1G(electric_field) !unit of electric field: atomic units
        use in_electric_field
        IMPLICIT NONE
        DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        integer, parameter :: Number_state = 50
        DOUBLE PRECISION :: Dipole, RotationalConstant
        DOUBLE PRECISION :: electric_field, E1G
        DOUBLE PRECISION :: NewEigValue00,NewEigValue10
        DOUBLE PRECISION :: coefficient_array00(Number_state)
        DOUBLE PRECISION :: coefficient_array10(Number_state)
        INTEGER :: final_state_number00, final_state_number10
      
        final_state_number00 = Number_state 
        final_state_number10 = Number_state 
        Dipole=5.529D0*0.393430307D0 !convert to atomic unit
        RotationalConstant=(11.7998D9/2.D0)*2.418884324306202D-17 !convert to atomic unit
        CALL New_EigStates(NewEigValue00,coefficient_array00,
     &                      Dipole, electric_field,
     &                      RotationalConstant,final_state_number00,
     &                      0, 0)
        CALL New_EigStates(NewEigValue10,coefficient_array10,
     &                      Dipole, electric_field,
     &                      RotationalConstant,final_state_number10,
     &                      1, 0)
        E1G=NewEigValue10 - NewEigValue00
        return
      END function E1G

!*******energy difference between dressed states |20> and |00> at certain electric field************
      FUNCTION E2G(electric_field) !unit of electric field: atomic unit
        use in_electric_field
        IMPLICIT NONE
        DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        integer, parameter :: Number_state = 50
        DOUBLE PRECISION :: Dipole, RotationalConstant
        DOUBLE PRECISION :: electric_field, E2G
        DOUBLE PRECISION :: NewEigValue00,NewEigValue20
        DOUBLE PRECISION :: coefficient_array00(Number_state)
        DOUBLE PRECISION :: coefficient_array20(Number_state)
        INTEGER :: final_state_number00, final_state_number20

        final_state_number00= Number_state 
        final_state_number20 = Number_state 
        Dipole=5.529D0*0.393430307D0 !convert to atomic unit
        RotationalConstant=(11.7998D9/2.D0)*2.418884324306202D-17 !convert to atomic unit
        CALL New_EigStates(NewEigValue00,coefficient_array00,
     &                      Dipole, electric_field,
     &                      RotationalConstant,final_state_number00,
     &                      0, 0)
        CALL New_EigStates(NewEigValue20,coefficient_array20,
     &                      Dipole, electric_field,
     &                      RotationalConstant,final_state_number20,
     &                      2, 0)
        E2G=NewEigValue20 - NewEigValue00
        return
      END FUNCTION E2G
!************************************************************************************** 
  
!***********************gas-condensed matter shift for N=A exciton, A can be 1 or 2*****************
! D1 = sum_{m\=n}[<10|Vnm|10> - <00|Vnm|00> ] ; D2 = sum_{m\=n}[<20|Vnm|20> - <00|Vnm|00> ] 
      Function GPSA(electric_field,Number_of_molecule, A)
        use sums
        use dipole_dipole
        use unit_conversion
        IMPLICIT NONE
        integer :: Number_of_molecule 
        integer :: A
        double precision :: dipole, electric_field,RotationalConstant,Lattice_Constant
        double precision :: GPSA
      	RotationalConstant=11.7998D9/2.D0  !unit=Hz 1GHz=10^9 Hz
      	Dipole=5.529D0 !Debye
        Lattice_Constant=4.D-7 !400 nm
           
        !***************unit conversion*************        
        CALL Dipole_in_atomic_unit(Dipole)
        CALL Length_in_Bohr(Lattice_Constant)
        CALL Hz_to_atomic_unit(RotationalConstant)
        
        GPSA=( DDInteract_in_Field(Dipole,electric_field,RotationalConstant,Lattice_Constant,
     &                             A,0,0,0,A,0,0,0)
     & - DDInteract_in_Field(Dipole,electric_field,RotationalConstant,Lattice_Constant,
     &                             0,0,0,0,0,0,0,0) )*2*half_sum2(Number_of_molecule)

        return ! unit of GPSA: atomic unit
      end FUNCTION GPSA

!***********************gas-condensed matter shift for N=A exciton, A can be 1 or 2*****************
! D1 = sum_{m\=n}[<10|Vnm|10> - <00|Vnm|00> ] ; D2 = sum_{m\=n}[<20|Vnm|20> - <00|Vnm|00> ] 
! angle between electric field and the intermolecular axis
      Function GPSA_angle(electric_field,Number_of_molecule,theta, A)
        use sums
        use dipole_dipole
        use unit_conversion
        IMPLICIT NONE
        integer :: Number_of_molecule 
        integer :: A
        double precision :: dipole, electric_field,RotationalConstant,Lattice_Constant
        double precision :: GPSA_angle
        double precision :: theta
      	RotationalConstant=11.7998D9/2.D0  !unit=Hz 1GHz=10^9 Hz
      	Dipole=5.529D0 !Debye
        Lattice_Constant=4.D-7 !400 nm
           
        !***************unit conversion*************        
        CALL Dipole_in_atomic_unit(Dipole)
        CALL Length_in_Bohr(Lattice_Constant)
        CALL Hz_to_atomic_unit(RotationalConstant)
        
        GPSA_angle=( DDInteract_in_Field_angle2(Dipole,electric_field,RotationalConstant,Lattice_Constant,
     &                             A,0,0,0,A,0,0,0, theta)
     & - DDInteract_in_Field_angle2(Dipole,electric_field,RotationalConstant,Lattice_Constant,
     &                             0,0,0,0,0,0,0,0, theta) )*2*half_sum2(Number_of_molecule)

        return ! unit of GPSA: atomic unit
      end FUNCTION GPSA_angle


!*************************************************************************
! the unit of exciton10_energy is kHz
! calculate the sum of dispersion energies for two exciton with wavevectors k, -k
! exciton10_energy is an array of pure excitonic energies for different ka
! the dimension of exciton10_exciton energy is pair_number (input)
! normally, the pair_number = Number_of_molecule/2 + 1
! this subroutine is for the one-dimensional case and we are assuming phi = 0
! this subroutine is for M=0 excitons
      Subroutine Twoexciton_energy_angle(electric_field, num_mole,theta, 
     &                                   pair_number,k1_array, k2_array, twoexciton_energy_array)
        use unit_conversion
        use sums
        use dipole_dipole
        IMPLICIT NONE
      	DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        integer :: pair_number
        integer :: k1_array(pair_number), k2_array(pair_number)
        double precision :: interaction_energy
        double precision :: twoexciton_energy_array(pair_number) 
        double precision :: theta  
      	INTEGER ::  I 
      	INTEGER :: num_mole
      	DOUBLE PRECISION :: rot_const,electric_field
        double precision :: summation(pair_number)
      	DOUBLE PRECISION :: Dipole, Lattice_Constant
        double precision :: ka1(pair_number), ka2(pair_number)

      	rot_const=11.7998D9/2.D0  !unit=Hz 1GHz=10^9 Hz
      	Dipole=5.529D0 ! Debye
        Lattice_Constant=4.D-7  ! 400nm       
        
        CALL Dipole_in_atomic_unit(Dipole)
        CALL Length_in_Bohr(Lattice_Constant)
        CALL Hz_to_atomic_unit(rot_const)

        interaction_energy =  DDInteract_in_Field_angle2(Dipole,electric_field,rot_const,Lattice_Constant,
     &                             1,0, 0,0, 0,0, 1,0,theta)    ! for M=0 excitons            
        call Energy_in_kHz(interaction_energy) !convert atomic unit to kHz 
        write(*,*) "J =", interaction_energy

          open(unit=40,file='J_1001.dat')
          write(40,*) "According to the notation in Vektaris' paper on biexciton"
          write(40,*) "M =", interaction_energy, "kHz" 
          ! M here is the sum of the exciton propagation matrix elements <10|Vnm|01>
          close(40)

        do i=1, pair_number
          ka1(i) = k1_array(i)*pi/(num_mole/2)
          ka2(i) = k2_array(i)*pi/(num_mole/2) 
          summation(i)=2*( half_sum1(num_mole, ka1(i)) + half_sum1(num_mole, ka2(i)) ) 
! using nearest-neighbor-approximation
!          summation(i)=2*( half_sum1(2, ka1(i)) + half_sum1(2, ka2(i)) ) 
        end do


        twoexciton_energy_array= interaction_energy*summation
        ! th sum of pure excitonic energies for two noninteracting excitons 
              
      END SUBROUTINE Twoexciton_energy_angle 


!****************************************************************************************
! give the excitonic energy of a N=1 exciton with wavevector k1
      Function E_10(electric_field,number_of_molecule,theta,k1n)
        use unit_conversion
        use dipole_dipole
        use f77_subroutine
        use sums
      	IMPLICIT NONE
      	DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        double precision, save :: interaction_energy
        double precision :: E_10
        double precision :: theta  
        integer :: k1n
        double precision :: k1
      	INTEGER ::  I 
      	INTEGER :: number_of_molecule
      	DOUBLE PRECISION :: RotationalConstant,electric_field
      	DOUBLE PRECISION :: Dipole, Lattice_Constant
        integer, save :: icounter
        data icounter /0/

        if (icounter==0) then
          icounter = 1
      	  RotationalConstant=11.7998D9/2.D0  !unit=Hz 1GHz=10^9 Hz
      	  Dipole=5.529D0  !Debye
          Lattice_Constant=4.D-7  !400nm       
        
          CALL Dipole_in_atomic_unit(Dipole)
          CALL Length_in_Bohr(Lattice_Constant)
          CALL Hz_to_atomic_unit(RotationalConstant)

          interaction_energy=DDInteract_in_Field_angle2(Dipole,electric_field,
     &                 RotationalConstant,Lattice_Constant,1,0,0,0,0,0,1,0,theta)                
        
          call Energy_in_kHz(interaction_energy) !convert atomic unit to kHz 
        endif
        k1 = k1n*pi/(number_of_molecule/2)
        E_10=interaction_energy*2*half_sum1(number_of_molecule, k1) 
        ! using nearest-neighbor-approximation
!        E_10=interaction_energy*2*half_sum1(2, k1)
        return      
      END function E_10


!****************************************************************************************
! give the excitonic energy of a N=2 exciton with wavevector k2
      Function E_20(electric_field,number_of_molecule,theta, k2n)
        use unit_conversion
        use dipole_dipole
        use f77_subroutine
        use sums
      	IMPLICIT NONE
      	DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
        double precision, save :: interaction_energy2
        double precision :: E_20 
        double precision :: theta
        integer :: k2n
        double precision :: k2 
      	INTEGER :: number_of_molecule
      	DOUBLE PRECISION :: RotationalConstant,electric_field
      	DOUBLE PRECISION :: Dipole, Lattice_Constant
        integer, save :: icounter2
        data icounter2 /0/

        if (icounter2==0) then
          icounter2 = 1
      	  RotationalConstant=11.7998D9/2.D0  !unit=Hz 1GHz=10^9 Hz
      	  Dipole=5.529D0  !Debye
          Lattice_Constant=4.D-7  !400nm       
        
          CALL Dipole_in_atomic_unit(Dipole)
          CALL Length_in_Bohr(Lattice_Constant)
          CALL Hz_to_atomic_unit(RotationalConstant)

          interaction_energy2=DDInteract_in_Field_angle2(Dipole,electric_field,
     &                 RotationalConstant,Lattice_Constant,2,0,0,0,0,0,2,0,theta)                
        
          call Energy_in_kHz(interaction_energy2) !convert atomic unit to kHz 
        endif
        k2 = k2n*pi/(number_of_molecule/2)
        E_20=interaction_energy2*2*half_sum1(number_of_molecule, k2) 
        return      
      END function E_20



!******************calculate the matrix element of the transition matrix***********************
! k1n, k1n_ are the wavevectors of N=1 excitons in the unit of 2*pi/N
! k2n are the wavevector of N=2 exciton in the unit of 2*pi/N
! k1n, k1n_, and k2n are integers in range [-N/2, N/2]
!unit of transition_mat is: kHz  Electric_field :: atomic unit
      Function transition_mat(electric_field, Number_of_molecule, theta, k1n, k1n_, k2n) 
        use dipole_dipole
        use unit_conversion
        use sums
      	IMPLICIT NONE
      	DOUBLE PRECISION, PARAMETER :: Pi=3.141592653589793115997963468d0
      	DOUBLE PRECISION :: transition_mat
        double precision :: theta
        integer :: k1n, k1n_, k2n
        double precision :: k1, k1_, k2
      	double precision, save :: matrix_element
      	INTEGER :: Number_of_molecule,i
      	DOUBLE PRECISION :: RotationalConstant,electric_field
      	DOUBLE PRECISION :: Dipole, Lattice_Constant
        integer, save :: ntime 
        data ntime /0/
        
        if (ntime ==0) then
          ntime =1
      	  RotationalConstant=11.7998D9/2.D0  !unit=Hz 1GHz=10^9 Hz
      	  Dipole=5.529D0 !Debye
          Lattice_Constant=4.D-7 !400nm
                   
          CALL Dipole_in_atomic_unit(Dipole)
          CALL Length_in_Bohr(Lattice_Constant)
          CALL Hz_to_atomic_unit(RotationalConstant)
     
      	  matrix_element = DDInteract_in_Field_angle2(Dipole,electric_field,RotationalConstant,Lattice_Constant,
     &                             2,0,0,0,1,0,1,0,theta)
          CALL Energy_in_kHz(matrix_element)
          open(unit=33,file='tran_mat.dat')
          write(33,*) "Transition matrix (two nearest neighbors):", 2*matrix_element
          close(33)
        endif

!        if ( (k1n + k1n_==k2n).or.( k1n + k1n_==k2n + 2*(Number_of_molecule/2)) 
!     &     .or.(k1n + k1n_==k2n - 2*(Number_of_molecule/2)) ) then
          k1 = k1n*pi/(Number_of_molecule/2)
          k1_ = k1n_*pi/(Number_of_molecule/2)
          k2 = k2n*pi/(Number_of_molecule/2)
          !--------------------------------------------------------------------------------------------------!
!          if ( (k1n==k1n_).or.(k1n==k1n_+2*(Number_of_molecule/2))                                           !
!     &        .or.(k1n==k1n_-2*(Number_of_molecule/2)) ) then                                                !
!            transition_mat =( 2*half_sum1(number_of_molecule, k1) + 2*half_sum1(Number_of_molecule,k1_) )    !
!     &                       *matrix_element/sqrt(number_of_molecule*2.D0)                                   !
!          else                                                                                               !
!            transition_mat =( 2*half_sum1(number_of_molecule, k1) + 2*half_sum1(number_of_molecule, k1_) )   !
!     &                       *matrix_element/sqrt(number_of_molecule*1.D0)                                   !
!          endif                                                                                              !
          !--------------------------------------------------------------------------------------------------!

        ! using a different normalization scheme for the basis set
          transition_mat = ( 2*half_sum1(number_of_molecule, k1) + 2*half_sum1(number_of_molecule, k1_) )   
     &                       *matrix_element/sqrt(number_of_molecule*1.D0) 
!        else
!          transition_mat = 0.D0
!        endif

        return            
      END Function transition_mat

      end module exciton
!=======================================================================================
!======================================================================================= 



!=======================================================================================
!======================================================================================= 
!------------------------------------------------------------------------
! usage of this module: in the main program, first give the value for K_n and Number_of_molecule,
! then call subroutine two_exciton_pairs, after that call subroutine assign_values 
! up to this point, the pair_number, k1_array and k2_array will be obtained
! if you want to do the same thing with another value of K_n, you have to 
! call deallocation before repeating the above procedure
      module find_exciton_pairs 
        implicit none
        integer :: K_n 
        ! the value for K_n comes from outside, K_n = k1 + k2
        ! K_n is the quantum number for biexciton
        integer :: pair_number
        integer :: Number_of_molecule ! the value comes from outside
        integer, allocatable :: k1_array(:)
        integer, allocatable :: k2_array(:) 
        integer :: tmp1(10000) ! 10000 should be large enough ! to save values for k1_array
        integer :: tmp2(10000)

        contains
        ! find all possible pairs of excitons given K
        subroutine two_exciton_pairs
          implicit none
          integer :: k1_n, k2_n ! k1_n >= k2_n
          integer :: i, a, b
          tmp1 = 0
          tmp2 = 0
          pair_number = 0

          if (K_n /= 0) then
!------------------------------------------------------------------------------------------
            ! find out all exciton pairs for which k1_n > k2_n
            do k1_n = -Number_of_molecule/2, Number_of_molecule/2
              do k2_n = -Number_of_molecule/2, k1_n
                if ( ( (k1_n + k2_n==K_n).or.(k1_n + k2_n==K_n+2*(Number_of_molecule/2))
     &             .or.(k1_n + k2_n==K_n-2*(Number_of_molecule/2)) ).and. (k1_n > k2_n) ) then
                  pair_number = pair_number + 1
                  tmp1( pair_number) = k1_n
                  tmp2( pair_number) = k2_n 
                endif         
              end do
            end do
            ! find out the exciton pairs for which k1_n = k2_n
            if ( MOD(K_n,2) == 0) then
              pair_number = pair_number + 1
              tmp1(pair_number) = K_n/2
              tmp2(pair_number) = K_n/2
              !---------------------------------------------------------
              if ( ( K_n < 0 ).and.(K_n /=Number_of_molecule/2).and.(K_n /=-Number_of_molecule/2 ) ) then
                pair_number = pair_number + 1
                tmp1(pair_number) = ( K_n + 2*(Number_of_molecule/2) )/2
                tmp2(pair_number) = ( K_n + 2*(Number_of_molecule/2) )/2 
              elseif ( ( K_n > 0 ).and.(K_n /=Number_of_molecule/2).and.(K_n /=-Number_of_molecule/2 ) ) then
                pair_number = pair_number + 1
                tmp1(pair_number) = ( K_n - 2*(Number_of_molecule/2) )/2
                tmp2(pair_number) = ( K_n - 2*(Number_of_molecule/2) )/2                  
              endif
              !---------------------------------------------------------
            else
              write(*,*) "k1 and k2 cannot be equal in this case"
            endif
!-------------------------------------------------------------------------------------------
          else
!-------------------------------------------------------------------------------------------
            do k1_n = 0, Number_of_molecule/2
              pair_number = pair_number + 1
              tmp1( pair_number) = k1_n
              tmp2( pair_number) = -k1_n
            end do
!            pair_number = pair_number + 1
!            tmp1(pair_number) = 0
!            tmp2(pair_number) = 0            
!-------------------------------------------------------------------------------------------
          endif
              
        end subroutine two_exciton_pairs
        
        

        ! assign the values to k1_array
        subroutine assign_values
          implicit none
          integer :: L
          allocate(k1_array(pair_number))
          allocate(k2_array(pair_number))
          k1_array(1:pair_number) = tmp1(1:pair_number)
          k2_array(1:pair_number) = tmp2(1:pair_number)
          open(unit=21,file='exciton_pairs.dat', position='append')
          write(21,*) "K_n=", K_n
          do L =1, pair_number
            write(21,*) k1_array(L), k2_array(L),  k1_array(L) + k2_array(L)
          end do
          write(21,*) "#***************************************"
          write(21,*) "#***************************************"
          close(21)
        end subroutine assign_values

        !deallocate k1_array if you want to find the possible exciton pairs for other K_n
        subroutine deallocation
          implicit none
          deallocate(k1_array)
          deallocate(k2_array)
        end subroutine deallocation

      end module find_exciton_pairs
!------------------------------------------------------------------------
!=======================================================================================
!======================================================================================= 



           





! The following three subprograms can be used to find at what electric field 
! the transition between N=2 exciton and N=1 exciton can be expected
!*****************find the root of func(x) in range (x1,x2)*************************
      FUNCTION zriddr(func,x1,x2,xacc)
      INTEGER MAXIT
      double precision :: zriddr,x1,x2,xacc,func,UNUSED
      PARAMETER (MAXIT=6000,UNUSED=-1.11E30)
      EXTERNAL func
CU    USES func
      INTEGER j
      double precision :: fh,fl,fm,fnew,s,xh,xl,xm,xnew
      fl=func(x1)
      fh=func(x2)
      if((fl.gt.0..and.fh.lt.0.).or.(fl.lt.0..and.fh.gt.0.))then
        xl=x1
        xh=x2
        zriddr=UNUSED
        do 11 j=1,MAXIT
          xm=0.5*(xl+xh)
          fm=func(xm)
          s=sqrt(fm**2-fl*fh)
          if(s.eq.0.)return
          xnew=xm+(xm-xl)*(sign(1.,fl-fh)*fm/s)
          if (abs(xnew-zriddr).le.xacc) return
          zriddr=xnew
          fnew=func(zriddr)
          if (fnew.eq.0.) return
          if(sign(fm,fnew).ne.fm) then
            xl=xm
            fl=fm
            xh=zriddr
            fh=fnew
          else if(sign(fl,fnew).ne.fl) then
            xh=zriddr
            fh=fnew
          else if(sign(fh,fnew).ne.fh) then
            xl=zriddr
            fl=fnew
          else
            pause 'never get here in zriddr'
          endif
          if(abs(xh-xl).le.xacc) return
11      continue
        pause 'zriddr exceed maximum iterations'
      else if (fl.eq.0.) then
        zriddr=x1
      else if (fh.eq.0.) then
        zriddr=x2
      else
        pause 'root must be bracketed in zriddr'
      endif
      return
      END function zriddr 


!***********************************************************************************
! calculate the energy difference between two N=1 excitons with k=0 
! and one N=2 exciton with k=0
! the electric_field is in atomic unit
      function func(electric_field)
        use exciton
        implicit none
        double precision :: func
        double precision :: electric_field
        integer, parameter :: Number_of_molecule = 100
        double precision :: exciton10_energy(Number_of_molecule/2 + 1)
        double precision :: exciton20_energy(Number_of_molecule/2 + 1)       

        call Exciton10_dispersion(electric_field,number_of_molecule,exciton10_energy)
        call Exciton20_dispersion(electric_field,number_of_molecule,exciton20_energy) 
        func = delta_energy_exclude(electric_field, number_of_molecule)  
     &                      + 2*exciton10_energy(1)
     &                      - exciton20_energy(1)
        return   
      end function func 
!*************************************************************************************
      subroutine find_electric_field(start_value, end_value, electric_field)
! start_value, end_value, and electric_field are all in atomic units
        implicit none
        double precision :: electric_field
        double precision :: start_value, end_value 
        double precision :: xacc 
        double precision, external :: func
        double precision, external :: zriddr
        xacc = 1.D-22   ! determine the accuracy of the solution     
        electric_field = zriddr(func,start_value,end_value,xacc)        
      end subroutine find_electric_field






!******Form the evolution Hamiltonian matrix for the transition between N=1 exciton and N=2 exciton****
      subroutine form_evolution_matrix(electric_field, evolution_matrix) 
        use exciton
        use find_exciton_pairs 
        ! the module find_exciton_pairs will pass the following variables into this subroutine:
        ! Number_of_molecule, k1_array, k2_array     
        implicit none
        double precision, parameter :: Pi=3.141592653589793115997963468d0  
        double precision :: electric_field, record
        double precision, parameter :: theta=pi/2 ! E is perpendicular to intermolecular axis
        double precision :: evolution_matrix(Number_of_molecule/2 + 2, Number_of_molecule/2 + 2)
        ! the basis sets are |N=2 exciton with k=0>, |two N=1 excitons with k=-k=0*2pi/N>
        ! |two N=1 excitons with k=-k=1*2pi/N>,|two N=1 excitons with k=-k=2*2pi/N>, ...
        ! |two N=1 excitons with k=-k=(N/2)*2pi/N>, so the dimension of evolution_matrix is
        ! 1 + 1 + N/2  
        integer :: i, j, K
        integer :: ka1n, ka2n, ka3n, ka4n
        double precision :: delta_energy(Number_of_molecule/2 + 1)
        double precision :: exciton10_energy(Number_of_molecule/2 + 1) 
        double precision :: exciton20_energy(Number_of_molecule/2 + 1)
        double precision :: m_k_array(Number_of_molecule/2 + 1)
        double precision, external :: dyn_interact
 
        evolution_matrix = 0.D0 

        CALL Exciton10_dispersion(electric_field,Number_of_molecule,exciton10_energy) 
        CALL Exciton20_dispersion(electric_field,Number_of_molecule,exciton20_energy)
        CALL M_vs_k(electric_field, Number_of_molecule,m_k_array)

        record = delta_energy_exclude(electric_field, number_of_molecule) 

        do i=1, Number_of_molecule/2 + 1  
          delta_energy(i) = 2*exciton10_energy(i) - exciton20_energy(1) + record
          ! assuming the energy of N=2 exciton is zero
        end do

        do j=2, Number_of_molecule/2 + 2 
          evolution_matrix(j,j) = delta_energy(j-1) 
        end do
        
        do K=2,Number_of_molecule/2 + 2
           evolution_matrix(1, K) = 2*m_k_array(K-1)/sqrt(Number_of_molecule*1.D0) 
!           evolution_matrix(K, 1) = evolution_matrix(1, K) 
        end do

        evolution_matrix(1,2) = evolution_matrix(1,2)/sqrt(2.D0)  ! k=0 case (kronecker delta)
        evolution_matrix(1, Number_of_molecule/2+2) =evolution_matrix(1, Number_of_molecule/2+2)/sqrt(2.D0) 
!        evolution_matrix(2,1) = evolution_matrix(2,1)/sqrt(2.D0)  ! k=0 case 
        ! The corresponding matrix element is: 
        !  [M(k) + M(-k)]/[ Sqrt( N*(1 + KroneckerDelta_{k,-k}) ) ]
        do i=1, pair_number
           ka1n = k1_array(i)
           ka2n = k2_array(i)
          do j = i, pair_number
            ka3n = k1_array(j)
            ka4n = k2_array(j)
            evolution_matrix(i+1,j+1) = evolution_matrix(i+1,j+1) 
     &                  + dyn_interact(Number_of_molecule,electric_field,theta, ka1n,ka2n,ka3n,ka4n)
          end do
        end do

        write(*,*) "The shape of evolution_matrix:", shape(evolution_matrix)
        do i=1, Number_of_molecule/2+2
             write(90,* ) evolution_matrix(i,:)
             write(90,*) " "
        end do
        evolution_matrix = 2000*Pi*evolution_matrix

      end subroutine form_evolution_matrix



!******Form the evolution Hamiltonian matrix for the transition between N=1 exciton and N=2 exciton****
      subroutine form_evolution_matrix_K(electric_field, theta, evolution_matrix) 
        use exciton
        use find_exciton_pairs 
        ! the module find_exciton_pairs will pass the following variables into this subroutine:
        ! Number_of_molecule, k1_array, k2_array     
        implicit none
        double precision, parameter :: Pi=3.141592653589793115997963468d0  
        double precision :: electric_field, record
        double precision :: theta 
        double precision :: evolution_matrix(Number_of_molecule/2 + 2, Number_of_molecule/2 + 2)
        ! the basis sets are |N=2 exciton with k=0>, |two N=1 excitons with k=-k=0*2pi/N>
        ! |two N=1 excitons with k=-k=1*2pi/N>,|two N=1 excitons with k=-k=2*2pi/N>, ...
        ! |two N=1 excitons with k=-k=(N/2)*2pi/N>, so the dimension of evolution_matrix is
        ! 1 + 1 + N/2  
        integer :: i, j, K
        integer :: ka1n, ka2n, ka3n, ka4n
        double precision :: delta_energy(Number_of_molecule/2 + 1)
        double precision, external :: dyn_interact
 
        evolution_matrix = 0.D0 

        record = delta_energy_ex_angle(electric_field, number_of_molecule,theta) 
        record = record - E_20(electric_field, Number_of_molecule,theta, 2*K_n)

        do i=1, Number_of_molecule/2 + 1 
          ka1n = k1_array(i)
          ka2n = k2_array(i) 
          delta_energy(i) = E_10(electric_field,Number_of_molecule,theta,ka1n) 
     &                    + E_10(electric_field,Number_of_molecule,theta,ka2n) + record
          ! K_n = (ka1n + ka2n )/2, then the wavevector of N=2 exciton is 2*K_n
        end do

        do j=2, Number_of_molecule/2 + 2 
          evolution_matrix(j,j) = delta_energy(j-1) 
        end do
        ! transition between two types of excitons
        do K=1,Number_of_molecule/2 + 1
          ka1n = k1_array(K)
          ka2n = k2_array(K)            
          evolution_matrix(1, K+1) = transition_mat(electric_field, Number_of_molecule,theta, ka1n, ka2n, 2*K_n) 
        end do
        ! consider dynamical interaction
        do i=1, pair_number
           ka1n = k1_array(i)
           ka2n = k2_array(i)
          do j = i, pair_number
            ka3n = k1_array(j)
            ka4n = k2_array(j)
            evolution_matrix(i+1,j+1) = evolution_matrix(i+1,j+1) 
     &                  + dyn_interact(Number_of_molecule,electric_field,theta, ka1n,ka2n,ka3n,ka4n)
          end do
        end do

        write(*,*) "The shape of evolution_matrix:", shape(evolution_matrix)
        do i=1, Number_of_molecule/2+2
             write(90,* ) evolution_matrix(i,:)
             write(90,*) " "
        end do
        evolution_matrix = 2000*Pi*evolution_matrix

      end subroutine form_evolution_matrix_K



!******Form the evolution Hamiltonian matrix for the transition between N=1 exciton and N=2 exciton****
!
      subroutine form_evolution_matrix_K2(electric_field, theta, em_n1, evolution_matrix) 
        use exciton
        use find_exciton_pairs 
        ! the module find_exciton_pairs will pass the following variables into this subroutine:
        ! Number_of_molecule, k1_array, k2_array     
        implicit none
        double precision, parameter :: Pi=3.141592653589793115997963468d0  
        double precision :: electric_field 
        double precision :: theta 
        double precision :: em_n1(pair_number, pair_number) ! evolution matrix for N=1 exciton pairs
        double precision :: evolution_matrix(pair_number+1, pair_number+1) 
        integer :: i, j, K
        integer :: ka1n, ka2n, ka3n, ka4n
 
        evolution_matrix = 0.D0 
        ! define the energy of dressed state |N=2,M=0> as zero
        evolution_matrix(1,1) = E_20(electric_field, Number_of_molecule,theta, K_n) ! K_n is the wavevector of N=2 exciton
        open(unit=73, file='N2_dispersion.dat', position='append')
        write(73,*) K_n*pi/(Number_of_molecule/2), evolution_matrix(1,1)
        close(73)
        evolution_matrix(2:pair_number+1, 2:pair_number+1) = em_n1(1:pair_number, 1:pair_number)

        ! transition between two types of excitons
        do K=1,pair_number
          ka1n = k1_array(K)
          ka2n = k2_array(K)            
          evolution_matrix(1, K+1) = transition_mat(electric_field, Number_of_molecule,theta, ka1n, ka2n, K_n) 
        end do

        evolution_matrix = 2000*Pi*evolution_matrix

      end subroutine form_evolution_matrix_K2



!******Form the evolution Hamiltonian matrix for the transition between N=1 exciton and N=2 exciton****
! including kinematic interaction
      subroutine form_evolution_matrix_kinematic(electric_field, theta, 
     &                        xx_energy, eig_vector_full_n1, evolution_matrix)
        use mkl95_precision, only: WP => DP
        use mkl95_blas, only: DOT 
        use exciton
        use find_exciton_pairs 
        ! the module find_exciton_pairs will pass the following variables into this subroutine:
        ! Number_of_molecule, k1_array, k2_array     
        implicit none
        double precision, parameter :: Pi=3.141592653589793115997963468d0  
        double precision :: electric_field 
        double precision :: theta 
        double precision :: eig_vector_full_n1(pair_number, pair_number) ! evolution matrix for N=1 exciton pairs
        double precision :: evolution_matrix(pair_number, pair_number) 
        double precision :: coefficient_array(pair_number)
        double precision :: trans_array(pair_number)
        double precision :: xx_energy(pair_number-1)
        double precision, save :: E_gap
        integer :: i, j, K
        integer :: ka1n, ka2n, ka3n, ka4n
        integer :: calltimes
        data calltimes /0/

        if (calltimes == 0 ) then
          calltimes = 1
          ! the energy difference: 2*|N=1,M=0> - |N=2,M=0>
          ! define the energy of dressed state |N=2,M=0> as zero
          E_gap = delta_energy_ex_angle(electric_field, number_of_molecule,theta)
        endif

        evolution_matrix = 0.D0 
        ! define the energy of dressed state |N=2,M=0> as zero
        evolution_matrix(1,1) = E_20(electric_field, Number_of_molecule,theta, K_n) ! K_n is the wavevector of N=2 exciton
        open(unit=73, file='N2_dispersion.dat', position='append')
        write(73,*) K_n*pi/(Number_of_molecule/2), evolution_matrix(1,1)
        close(73)
        
        ! calculat the transition elements in the basis set of |k1, k2>
        ! then they will be converted into the new basise set which include the biexciton state
        do i=1,pair_number
          ka1n = k1_array(i)
          ka2n = k2_array(i)            
          trans_array(i) = transition_mat(electric_field, Number_of_molecule,theta, ka1n, ka2n, K_n) 
        end do
        
        do i = 1, pair_number-1
          coefficient_array(1:pair_number) = eig_vector_full_n1(1:pair_number, i) 
          evolution_matrix(1, i+1) = DOT(coefficient_array, trans_array)
          evolution_matrix(i+1, 1) = evolution_matrix(1, i+1)
        end do

        do i = 2, pair_number
          evolution_matrix(i,i) = xx_energy(i-1) 
        end do

!        evolution_matrix = 2000*Pi*evolution_matrix
        evolution_matrix = 1000*evolution_matrix
      end subroutine form_evolution_matrix_kinematic










! form the full basis set
! N must be odd
      subroutine basis_set(N, basisN1, basisN2)
        implicit none
        integer :: N
        integer :: basisN2(1:2*(N/2)) ! the elements of the array are k2n of the N=2 exciton
        integer :: basisN1(2*(N/2)+1 : N*(N-1)/2 + 2*(N/2), 2) 
        ! we assume k1n is associated with exciton 1, k1n_ with exciton 2
        integer :: i, j, k, m
        
        
        do i=1, 2*(N/2)
           basisN2(i) = -N/2 + (i-1)
           write(99,*) i, basisN2(i) 
        end do
           write(99,*) "============="
        j = 2*(N/2) +1 
        do k=-N/2, N/2-1 !-N/2 and N/2 are equivalent
          ! assuming k1n <= k1n_
          do m=k, N/2-1
!            write(*,*) j
            basisN1(j,1) = k ! k1n
            basisN1(j,2) = m ! k1n_
            write(99,*) j, basisN1(j,1), basisN1(j,2)  
            j = j + 1
          end do
        end do 

      end subroutine basis_set



!******Form the evolution Hamiltonian matrix for the transition between N=1 exciton and N=2 exciton****
! this version includes the full basis sets: all k states for N=2, all pairs for N=1
      subroutine form_evolution_matrix2(electric_field, N, theta,evolution_matrix) 
        use exciton     
        implicit none
        double precision, parameter :: Pi=3.141592653589793115997963468d0  
        double precision :: electric_field, record, theta
        integer :: N !number_of_molecule
        double precision :: evolution_matrix(N*(N-1)/2+2*(N/2), N*(N-1)/2+2*(N/2) )
        integer ::  basisN2(1:2*(N/2))! the elements of the array are k2n of the N=2 exciton
        integer ::  basisN1(2*(N/2)+1 : N*(N-1)/2 + 2*(N/2), 2)
        ! we assume k1n is associated with exciton 1, k1n_ with exciton 2
        integer :: i, j, K
        integer :: lk1n, lk1n_, rk1n, rk1n_, k2n ! l---left <|; r----right |>
        double precision, external :: dyn_interact !, E_10, E_20, transition_mat  are in module exciton
        external :: basis_set

        write(*,*) "0"

        write(*,*) "01"
        evolution_matrix = 0.D0 
        write(*,*) "02"
        record = delta_energy_ex_angle(electric_field, N,theta) 
        write(*,*) "03"
        ! form the basis set
        call basis_set(N, basisN1, basisN2) 
        write(*,*) "04"
! diagonal terms without considering the dynamical interactions
!-------------------------------------------------------------------------------
        do i=1, 2*(N/2)
              k2n = basisN2(i)
              evolution_matrix(i,i) = E_20(electric_field,N,theta,k2n) 
!              write(90,*) i, i, evolution_matrix(i,i) 
        end do

        do i=2*(N/2)+1, N*(N-1)/2 + 2*(N/2)
              lk1n = basisN1(i,1)
              lk1n_ = basisN1(i,2)
              evolution_matrix(i,i) = E_10(electric_field, N, theta,lk1n)
     &                           + E_10(electric_field, N, theta,lk1n_) 
     &                           + record
        end do
!-------------------------------------------------------------------------------
        write(*,*) "05"
! consider the dynamical interactions
        do i=2*(N/2)+1, N*(N-1)/2 + 2*(N/2) 
          do j=i, N*(N-1)/2 + 2*(N/2)
          lk1n = basisN1(i,1)
          lk1n_ = basisN1(i,2)
          rk1n = basisN1(j,1)
          rk1n_ = basisN1(j,2) 
          evolution_matrix(i,j) = evolution_matrix(i,j)  
     &     + dyn_interact(N,electric_field,theta, lk1n,lk1n_,rk1n,rk1n_) 
!          write(91,*)  lk1n, lk1n_,"|",rk1n,rk1n_,"|", evolution_matrix(i,j)                        
          end do
        end do

        write(*,*) "06"
! consider the transiton between N=2 and N=1 excitons
        do i=1, 2*(N/2)
          do j= 2*(N/2)+1, N*(N-1)/2 + 2*(N/2)
             k2n = basisN2(i)
             rk1n = basisN1(j,1)
             rk1n_ = basisN1(j,2) 
             evolution_matrix(i,j)=transition_mat(electric_field, N, theta,rk1n, rk1n_, k2n)
!             write(90,*) k2n, rk1n, rk1n_, evolution_matrix(i,j)
          end do
        end do
        write(*,*) "07"
        evolution_matrix = 2000*Pi*evolution_matrix
!        write(*,*) "The shape of evolution_matrix:", shape(evolution_matrix)
!        do i=1, 2*(N/2)
!          do j=2*(N/2)+1 ,N*(N-1)/2 + 2*(N/2) 
!             write(90,"(1000(f14.1))" ) basisN2(i), basisN1(j,1), basisN1(j,2), evolution_matrix(i,j)
!             write(90,*) " "
!          end do
!        end do

        write(*,*) "08"
      end subroutine form_evolution_matrix2


! need to write an interface in the calling program, otherwise it may lead to errors
      subroutine probability_dis2(time, eig_vector_matrix, eig_value, initial_coefficients,
     &                            coefficients_squared_k, K_n)
! ndim is the dimension of initial_coefficients, it is a even integer
! i*hbar*dot(C) = HC, i*hbar*Ut*dot(C) = Ut*H*C = Ut*H*U*Ut*C where Ut*H*U = D (diagonalized matrix)
! assume A = Ut*C, then i*hbar*dot(A) = D*A, A(t) = exp(D*t/(i*hbar))*A(t=0)
! convert A back to C, we have C(t) = U*exp(D*t/(i*hbar))*Ut*C(t=0) 
! Note U is the eigenvectors matrix whose columns are eigenvectors of H, it is the output of the lapack subroutine 
        implicit none
        double precision, parameter :: Pi=3.141592653589793115997963468d0
        integer, save :: ndim2 
        integer :: i
        integer :: K_n
        integer, save :: last_K_n
        double precision ::  time
        double precision :: coefficients_squared_k( : ) ! assumed-shape, require explicit interface in calling program
        double precision, allocatable :: coefficients_real( : )
        double precision, allocatable :: coefficients_imaginary( : )
        double precision :: eig_vector_matrix( :, : )
!        double precision, allocatable :: eig_vector_matrix_t(:,:)     
        double precision :: eig_value( : )
        double precision :: initial_coefficients( : )
        double precision, allocatable :: exponential_Dt_real(:,:)
        double precision, allocatable :: exponential_Dt_imaginary(:,:)
        double precision, allocatable, save :: Utc(:)
        double precision :: sum_coefficient
        integer, save :: counter
        integer, save :: ireg
        data counter /1/
        data ireg /1/
        data last_K_n /100000000000/
        external :: dgemv ! from blas, parallel version, should be faster than intrinsic matmul
        double precision, allocatable :: intermediate(:) 
        
        if ( (last_K_n /= K_n).and.(ireg/=1) ) then
          counter = 1
          deallocate(Utc)
        endif        

        if (counter == 1) then
          counter = 2
          ndim2 = Size(initial_coefficients)
          allocate(Utc(ndim2))
          call dgemv('T',ndim2,ndim2,1.D0,eig_vector_matrix,ndim2,initial_coefficients,1,0.D0,Utc,1)
        endif
         
        allocate(coefficients_real(ndim2))
        allocate(coefficients_imaginary(ndim2))
        allocate(exponential_Dt_real(ndim2,ndim2))
        allocate(exponential_Dt_imaginary(ndim2,ndim2))
        allocate(intermediate(ndim2))

        exponential_Dt_real = 0.D0
        exponential_Dt_imaginary = 0.D0
        ! set the initial values because the following loop only assign the diagonal elements
        do i=1, ndim2
          exponential_Dt_real(i,i) = COS(eig_value(i)*time)
          exponential_Dt_imaginary(i,i) = SIN(-eig_value(i)*time)
        end do

!        coefficients_real = MATMUL(  eig_vector_matrix, MATMUL( exponential_Dt_real, Utc )  )
        ! intermediate = exponential_Dt_real*Utc
        call dgemv('N',ndim2,ndim2,1.D0,exponential_Dt_real,ndim2,Utc,1,0.D0,intermediate,1 )
        ! coefficients_real=eig_vector_matrix*intermediate
        call dgemv('N',ndim2,ndim2,1.D0,eig_vector_matrix,ndim2,intermediate,1,0.D0,coefficients_real,1 )
!        coefficients_imaginary = MATMUL(  eig_vector_matrix, MATMUL( exponential_Dt_imaginary, Utc )  )
        call dgemv('N',ndim2,ndim2,1.D0,exponential_Dt_imaginary,ndim2,Utc,1,0.D0,intermediate,1 )
        call dgemv('N',ndim2,ndim2,1.D0,eig_vector_matrix,ndim2,intermediate,1,0.D0,coefficients_imaginary,1 )

        deallocate(exponential_Dt_real )
        deallocate(exponential_Dt_imaginary)
        deallocate(intermediate)  
        coefficients_squared_k = coefficients_real**2 + coefficients_imaginary**2
        deallocate(coefficients_real)
        deallocate(coefficients_imaginary)
        sum_coefficient = SUM(coefficients_squared_k)
        ! normalization
        coefficients_squared_k = coefficients_squared_k/sum_coefficient 
        ! save the value for K_n
        last_K_n = K_n
        ireg = ireg + 1
      end subroutine probability_dis2




      subroutine probability_dis(time, eig_vector_matrix, eig_value, initial_coefficients, 
     &                           coefficients_squared, coefficients_squared_x)
! i*hbar*dot(C) = HC, i*hbar*Ut*dot(C) = Ut*H*C = Ut*H*U*Ut*C where Ut*H*U = D (diagonalized matrix)
! assume A = Ut*C, then i*hbar*dot(A) = D*A, A(t) = exp(D*t/(i*hbar))*A(t=0)
! convert A back to C, we have C(t) = U*exp(D*t/(i*hbar))*Ut*C(t=0) 
! Note U is the eigenvectors matrix whose columns are eigenvectors of H, it is the output of the lapack subroutine 
        implicit none
        double precision, parameter :: Pi=3.141592653589793115997963468d0
        integer, parameter :: Number_of_molecule = 51
        integer :: i, j, K
        double precision ::  time
        double precision :: coefficients_squared(Number_of_molecule/2+2)
        double precision :: coefficients_real(Number_of_molecule/2+2)
        double precision :: coefficients_imaginary(Number_of_molecule/2+2)
        double precision :: eig_vector_matrix(Number_of_molecule/2+2, Number_of_molecule/2 + 2)
        double precision :: eig_vector_matrix_t(Number_of_molecule/2+2, Number_of_molecule/2 + 2)     
        double precision :: eig_value(Number_of_molecule/2+2)
        double precision :: initial_coefficients(Number_of_molecule/2+2)
        double precision :: exponential_Dt_real(Number_of_molecule/2+2, Number_of_molecule/2+2)
        double precision :: exponential_Dt_imaginary(Number_of_molecule/2+2, Number_of_molecule/2+2)
        double precision :: sum_coefficient
        double complex :: coefficient_in_k(0:Number_of_molecule/2)
        double complex :: coefficient_in_x(-Number_of_molecule/2: Number_of_molecule/2)
        double precision :: ka
        double precision :: coefficients_squared_x(-Number_of_molecule/2: Number_of_molecule/2) 

        exponential_Dt_real = 0.D0
        exponential_Dt_imaginary = 0.D0 
        sum_coefficient = 0.D0

        eig_vector_matrix_t = TRANSPOSE(eig_vector_matrix) 
        
        do i=1, Number_of_molecule/2 + 2
          exponential_Dt_real(i,i) = COS(eig_value(i)*time)
          exponential_Dt_imaginary(i,i) = SIN(-eig_value(i)*time)
        end do

        coefficients_real = MATMUL(  eig_vector_matrix, MATMUL( exponential_Dt_real, 
     &                          MATMUL(eig_vector_matrix_t, initial_coefficients) )  )

        coefficients_imaginary = MATMUL(  eig_vector_matrix, MATMUL( exponential_Dt_imaginary, 
     &                          MATMUL(eig_vector_matrix_t, initial_coefficients) )  )

        coefficients_squared = coefficients_real**2 + coefficients_imaginary**2
        sum_coefficient = SUM(coefficients_squared) 
        coefficients_squared = coefficients_squared/sum_coefficient 

        ! wavefunction in k space
        coefficient_in_k(0:Number_of_molecule/2) = coefficients_real(2:Number_of_molecule/2+2) 
     &                      + coefficients_imaginary(2:Number_of_molecule/2+2)*(0.D0,1.D0)  

        ! convert from k space to x space
        do i = -Number_of_molecule/2, Number_of_molecule/2 ! i ---- molecule index
          do j = 1, Number_of_molecule/2 ! ka>0
            ka = j*pi/(Number_of_molecule/2)
            coefficient_in_x(i) = coefficient_in_x(i) + 2*coefficient_in_k(j)*exp(ka*i*(0.D0,1.D0))
            ! 2 accouts for ka<0
          end do
          coefficient_in_x(i) = coefficient_in_x(i) + coefficient_in_k(0)
        end do
        coefficients_squared_x = ABS(coefficient_in_x)**2
        sum_coefficient = SUM(coefficients_squared_x)
        coefficients_squared_x = coefficients_squared_x/sum_coefficient 
      end subroutine probability_dis








      subroutine biexciton_ham_K(Number_of_molecule, electric_field, theta, ham, pair_number, k1_array, k2_array, K_n)
        use exciton
        implicit none
        double precision, parameter :: Pi=3.141592653589793115997963468d0
        integer :: pair_number
        integer :: K_n
        integer :: k1_array(pair_number), k2_array(pair_number)  
        double precision :: electric_field
        integer :: Number_of_molecule
        double precision :: dyn_interact
        double precision :: ham(pair_number-1,pair_number-1)
        double precision :: twoexciton_energy_array(pair_number)
        integer :: ka1n, ka2n, ka3n, ka4n
        integer :: ka3n_last, ka4n_last
        double precision :: E_min, E_max
        double precision :: theta
        integer :: i, j
        integer :: ndim
        double precision :: delta_E
        common /energy_difference/ delta_E
                
        ndim = pair_number -1
        ham = 0.D0
        ka3n_last = k1_array(pair_number)
        ka4n_last = k2_array(pair_number)
!----------------------------------------------------------------------------------------------
! consider the dynamical terms
        do i=1, ndim
           ka1n = k1_array(i)
           ka2n = k2_array(i)

          do j = 1, ndim
             ka3n = k1_array(j)
             ka4n = k2_array(j)
             ham(i,j) = dyn_interact(Number_of_molecule,electric_field,theta, ka1n,ka2n,ka3n,ka4n)
     &                -dyn_interact(Number_of_molecule,electric_field,theta, ka1n,ka2n,ka3n_last,ka4n_last)  
          end do
        end do
!------------------------------------------------------------------------------------------------


 ! include the kinematic interaction
        do i=1, ndim
           ka1n = k1_array(i)
           ka2n = k2_array(i)
!
          do j = 1, ndim
             ka3n = k1_array(j)
             ka4n = k2_array(j)

            ham(i,j) = ham(i,j)- 2*( E_10(electric_field,number_of_molecule,theta,ka3n)
     &                 + E_10(electric_field,number_of_molecule,theta,ka4n) )/number_of_molecule 
     &                         + 2*( E_10(electric_field,number_of_molecule,theta,ka3n_last)
     &                 + E_10(electric_field,number_of_molecule,theta,ka4n_last) )/number_of_molecule
          end do
        end do


! consider the noninteracting two-exciton terms
        call  Twoexciton_energy_angle(electric_field, Number_of_molecule,theta, 
     &                                   pair_number,k1_array, k2_array, twoexciton_energy_array) 
       
        do i =1, ndim
          ham(i,i) = ham(i,i) + twoexciton_energy_array(i)
        end do


! record the mimum and maxium of energy for two free N=1 excitons
        E_min = MINVAL(twoexciton_energy_array)
        E_max = MAXVAL(twoexciton_energy_array)
        open(unit=10,file='free_excitons_Emin.dat',position='append')
        write(10,*)  K_n*pi/(Number_of_molecule/2), E_min + delta_E
        close(10)

        open(unit=20,file='free_excitons_Emax.dat',position='append')
        write(20,*)  K_n*pi/(Number_of_molecule/2), E_max + delta_E
        close(20)

      end subroutine biexciton_ham_K











! calculate the matrix element corresponding to the dynamical interaction 
! < ka1, ka2 || H_dyn || ka3, ka4 >
! |ka1, ka2>  ka1 >= ka2
! |ka3, ka4>  ka3 >= ka4
! note the electric field can only be specified once in the main program
! otherwise it will lead to errors
! electric field in atomic unit
      Function dyn_interact(Number_of_molecule,electric_field, theta, ka1n, ka2n, ka3n, ka4n)
        use unit_conversion
        use sums
        use dipole_dipole
        implicit none
        double precision, parameter :: Pi=3.141592653589793115997963468d0
        double precision :: dyn_interact
        double precision :: electric_field 
        integer :: ka1n, ka2n, ka3n, ka4n
        double precision :: ka1, ka2, ka3, ka4
!        double precision,external :: DDInteract_in_Field_angle2  
        integer :: Number_of_molecule
      	DOUBLE PRECISION :: RotationalConstant
        double precision :: summation, prefactor
      	DOUBLE PRECISION :: Dipole, Lattice_Constant
        double precision :: theta
        double precision, save :: D_1111
        integer, save :: counter
!        double precision, external :: half_sum4
        data counter / 0 /


        if (counter /= 1) then
          counter =1
       	  RotationalConstant=11.7998D9/2.D0  !unit=Hz 1GHz=10^9 Hz
      	  Dipole=5.529D0 ! Debye
          Lattice_Constant=4.D-7 ! 400nm         
        
          CALL Dipole_in_atomic_unit(Dipole)
          CALL Length_in_Bohr(Lattice_Constant)
          CALL Hz_to_atomic_unit(RotationalConstant)
          ! calculate the dynamical term <ff|V|ff> + <gg|V|gg> -2*<fg|V|fg>
      	  D_1111 = DDInteract_in_Field_angle2(Dipole,Electric_Field,RotationalConstant,Lattice_Constant,
     &                             1,0,  1,0,  1,0,  1,0, theta)
     &  + DDInteract_in_Field_angle2(Dipole,Electric_Field,RotationalConstant,Lattice_Constant,0,0, 0,0, 0,0, 0,0, theta)
     &  -2*DDInteract_in_Field_angle2(Dipole,Electric_Field,RotationalConstant,Lattice_Constant,1,0,  0,0,  1,0,  0,0, theta)
          
          D_1111 = D_1111/2.D0 !divided by 2 because of 1/2 in front of the double summation

          call Energy_in_kHz(D_1111) !convert atomic unit to kHz 
          write(*,*) "D_1111 =" , D_1111
          open(unit=40,file='D_1111.dat')
          write(40,*) "According to the notation in Vektaris' paper on biexciton"
          write(40,*) "mu =", 2*D_1111, "kHz"
          close(40)
        endif

        ka1 = ka1n*pi/(Number_of_molecule/2) 
        ka2 = ka2n*pi/(Number_of_molecule/2)
        ka3 = ka3n*pi/(Number_of_molecule/2)
        ka4 = ka4n*pi/(Number_of_molecule/2)  
     


! using another method
! nearest-neighbor-approximation
!         summation = 2*half_sum1(2, ka1-ka3) + 2*half_sum1(2, ka1-ka4)

         summation = 2*half_sum1(Number_of_molecule, ka1-ka3) + 2*half_sum1(Number_of_molecule, ka1-ka4)
         prefactor = 2.D0



         
!        summation = 2*half_sum4(Number_of_molecule, ka1, ka2, ka4) ! times 2 because it is half of the sum
! using nearest-neighbor-approximation, N =2
!         summation = 2*half_sum4(2, ka1, ka2, ka4)
!
!        if ( (ka1n+ka2n==ka3n+ka4n).or.(ka1n+ka2n==ka3n+ka4n+2*(Number_of_molecule/2))  
!     &    .or.(ka1n+ka2n==ka3n+ka4n-2*(Number_of_molecule/2)) ) then
!
!          prefactor = 2.D0
!!---------------------------------------------------------------------
!          if ( (ka1n==ka2n).or.(ka1n==ka2n+2*(Number_of_molecule/2))
!     &         .or.(ka1n==ka2n-2*(Number_of_molecule/2)) ) then
!            prefactor = prefactor/sqrt(2.D0)
!          else
!            prefactor = prefactor 
!          endif
!
!          if ( (ka3n==ka4n).or.(ka3n==ka4n+2*(Number_of_molecule/2))
!     &         .or.(ka3n==ka4n-2*(Number_of_molecule/2)) ) then
!            prefactor = prefactor/sqrt(2.D0)
!          else
!            prefactor = prefactor
!          endif
!!---------------------------------------------------------------------
!        else
!          prefactor = 0.D0
!        endif
        
        dyn_interact = prefactor*D_1111*summation


        dyn_interact = dyn_interact/Number_of_molecule
        return
      end function dyn_interact

! normalize an eigvector matrix
      subroutine normalize_eigvector_matrix(eig_vector_matrix, N)
        implicit none
        integer :: N
        double precision :: eig_vector_matrix(N,N)
        double precision :: column_squared(N)
        double precision :: column_sum
        double precision :: normalization_factor
        integer :: i, j
        
        do j = 1, N
          column_squared(1:N) = eig_vector_matrix(1:N,j)
          column_squared = column_squared**2
          column_sum = SUM(column_squared) 
          normalization_factor = sqrt(column_sum)
          eig_vector_matrix(1:N,j) = eig_vector_matrix(1:N,j)/normalization_factor
        end do
      end subroutine

 
      subroutine kick_one_out(old, N, element_index, new)
        implicit none
        double precision :: old(N)
        double precision :: new(N-1)
        integer :: N 
        integer :: element_index ! the index for the element you want to kick out
        
        if (element_index==1) then
          new(1:N-1) = old(2:N)
        elseif (element_index==N) then
          new(1:N-1) = old(1:N-1)
        else
          new(1:element_index-1) = old(1:element_index-1)
          new(element_index: N-1) = old(element_index+1:N)
        endif
        
      end subroutine
      
      
   
   
   
   
   
   
   
   
   
   
      
      
