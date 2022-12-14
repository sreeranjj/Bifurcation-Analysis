function F= NewInterpol(alpha,coeff,xval)

%alpha=ALPHA_BREAK ->Vector
%coeff=Aerodynamic Coefficient which needs to be interpolated-> Vector
%xval=Alpha value in degree

global ALPHA_BREAK clift clift_de cd0 cd_de cm cm_de m g S rho Tm V gamma 
% This function interpolates the angle of attack data points
% This function returns the interpolated value of the aerodynamic
% coefficients corresponding to a particular angle of attack



%SUBROUTINE SPLNE(FI,XI,Xval,val)

% F1 can be any of the aerodynamic coefficient
% X1 is alpha or ALPHA_BREAK
% Xval is the current value of alpha
% val is the interpolated value of the coefficient

F1=coeff;
X1=alpha;
Xval=xval;

F=val;
% An example of creating cubic-spline approximation of
% a discrete function fi=f(xi).

  
  %INTEGER, PARAMETER :: N=26, M=100
  %INTEGER :: I, K,ST,SIZ,J,loca
  %DOUBLE PRECISION :: X, F, DX, H, ALPHA, BETA, GAMM, ETA,Xval,val,mini,M1,M2
  %DOUBLE PRECISION, DIMENSION (N+1) :: XI, FI, P2
  %DOUBLE PRECISION, DIMENSION (M+1) :: Xnew,Fnew

  
 N=26; % Number of data points minus one
 M=100; % Number of data points using which the coefficients will be interpolated
 
 Xnew=zeros(M+1);
 Fnew=zeros(M+1);
 P2=zeros(N+1);
  %CALL CUBIC_SPLINE(N, XI, FI, P2)
  
  %SUBROUTINE CUBIC_SPLINE (N, XI, FI, P2)

% Function to carry out the cubic-spline approximation
% with the second-order derivatives returned.

  %INTEGER :: I
  %INTEGER, INTENT (IN) :: N
  %DOUBLE PRECISION, INTENT (IN), DIMENSION (N+1):: XI, FI
  %DOUBLE PRECISION, INTENT (OUT), DIMENSION (N+1):: P2
  %DOUBLE PRECISION, DIMENSION (N):: G, H
  %DOUBLE PRECISION, DIMENSION (N-1):: D, B, C

% Assign the intervals and function differences

G=zeros(N);
H=zeros(N);
D=zeros(N-1);
B=zeros(N-1);
C=zeros(N-1);

for I=1:N
    
  %DO I = 1, N
    H(I) = XI(I+1) - XI(I);
    G(I) = FI(I+1) - FI(I);
  
    %END DO
end

% Evaluate the coefficient matrix elements
 for I=1:N-1
     
%DO I = 1, N-1
    D(I) = 2*(H(I+1)+H(I));
    B(I) = 6*(G(I+1)/H(I+1)-G(I)/H(I));
    C(I) = H(I+1);
  %END DO
 end

% Obtain the second-order derivatives
%Cm_de=interp1(ALPHA_BREAK,cm_de,x(2));

  %CALL TRIDIAGONAL_LINEAR_EQ (N-1, D, C, C, B, G)
  
  %SUBROUTINE TRIDIAGONAL_LINEAR_EQ (L, D, E, C, B, Z)

% Functione to solve the tridiagonal linear equation set.

  %INTEGER, INTENT (IN) :: L
  %INTEGER :: I
  %%DOUBLE PRECISION, INTENT (IN), DIMENSION (L):: D, E, C, B
  %DOUBLE PRECISION, INTENT (OUT), DIMENSION (L):: Z
  %DOUBLE PRECISION, DIMENSION (L):: Y, W
  %DOUBLE PRECISION, DIMENSION (L-1):: V, T

 %Evaluate the elements in the LU decomposition
 W=zeros(L);
 T=zeros(L-1)
  W(1) = D(1);
  V(1)  = C(1);
  T(1)  = E(1)/W(1);
  for I=2:L-1
  %DO I = 2, L - 1;
    W(I) = D(I)-V(I-1)*T(I-1);
    V(I) = C(I);
    T(I) = E(I)/W(I);
  end
  
    %END DO
  W(L) = D(L)-V(L-1)*T(L-1);

% Forward substitution to obtain y

  Y(1) = B(1)/W(1)
  for I=2:L
  %DO I = 2, L
    Y(I) = (B(I)-V(I-1)*Y(I-1))/W(I);
  %END DO
  end
% Backward substitution to obtain z
  Z(L) = Y(L);
  DO I = L-1, 1, -1
  for I=L-1:1:-1
    Z(I) = Y(I) - T(I)*Z(I+1)
  %END DO
  end
END SUBROUTINE TRIDIAGONAL_LINEAR_EQ
  
  
  
  
  
  
  
  
  
  
  
  
  
  P2(1) = 0
  P2(N+1) = 0
  DO I = 2, N 
    P2(I) = G(I-1)
  END DO
END SUBROUTINE CUBIC_SPLINE
!
  
  
  
  
  

% Find the approximation of the function
  ST=1;
  H = (XI(N+1)-XI(1))/M;
  X = XI(1);
  
 for I= 1:M-1
     X=X+H;
 
  %DO I = 1, M-1
  %  X = X + H

% Find the interval that x resides
    K = 1;
    DX = X-XI(1);
    while (DX>=0)
        K=K+1;
        DX=X-X1(K);
    end
    K=K-1;
    
        
    %DO WHILE (DX .GE. 0) %.GE.          >=   (greater than or equal to)
      %K = K + 1
      %DX = X-XI(K)
    %END DO
    %K = K - 1

% Find the value of function f(x)

    DX = XI(K+1) - XI(K);
    ALPHA = P2(K+1)/(6*DX);
    BETA = -P2(K)/(6*DX);
    GAMM = FI(K+1)/DX - DX*P2(K+1)/6;
    ETA = DX*P2(K)/6 - FI(K)/DX;
    F = ALPHA*(X-XI(K))*(X-XI(K))*(X-XI(K)) +BETA*(X-XI(K+1))*(X-XI(K+1))*(X-XI(K+1)) +GAMM*(X-XI(K))+ETA*(X-XI(K+1));
    Xnew(ST)=X;
    Fnew(ST)=F;
    
    ST=ST+1;
 end 

SIZ=length(Xnew);

mini=0.5;
loca=0;

for J=1:SIZ-2
%DO J=1,SIZ-2
   if ((Xnew(J))>Xval)
    loca=J;
    break;
    
    %EXIT
   end
end


M1=Xnew(J)-Xval;
M2=Xval-Xnew(J-1);
val=(M1*Fnew(J-1)+M2*Fnew(J))/(M1+M2);


  




%END SUBROUTINE SPLNE




SUBROUTINE TRIDIAGONAL_LINEAR_EQ (L, D, E, C, B, Z)
!
! Functione to solve the tridiagonal linear equation set.
!
  INTEGER, INTENT (IN) :: L
  INTEGER :: I
  DOUBLE PRECISION, INTENT (IN), DIMENSION (L):: D, E, C, B
  DOUBLE PRECISION, INTENT (OUT), DIMENSION (L):: Z
  DOUBLE PRECISION, DIMENSION (L):: Y, W
  DOUBLE PRECISION, DIMENSION (L-1):: V, T
!
! Evaluate the elements in the LU decomposition
!
  W(1) = D(1)
  V(1)  = C(1)
  T(1)  = E(1)/W(1)
  DO I = 2, L - 1
    W(I) = D(I)-V(I-1)*T(I-1)
    V(I) = C(I)
    T(I) = E(I)/W(I)
  END DO
  W(L) = D(L)-V(L-1)*T(L-1)
!
! Forward substitution to obtain y
!
  Y(1) = B(1)/W(1)
  DO I = 2, L
    Y(I) = (B(I)-V(I-1)*Y(I-1))/W(I)
  END DO
!
! Backward substitution to obtain z
  Z(L) = Y(L)
  DO I = L-1, 1, -1
    Z(I) = Y(I) - T(I)*Z(I+1)
  END DO
END SUBROUTINE TRIDIAGONAL_LINEAR_EQ