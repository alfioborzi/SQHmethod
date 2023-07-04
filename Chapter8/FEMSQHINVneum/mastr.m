function [INTPT,W,PHI,PHIX,PHIY, ...
          INTPTE1,WE1,PHIE1, ...
          INTPTE2,WE2,PHIE2, ...
          INTPTE3,WE3,PHIE3] = mastr(n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% purpose:
%
%    define the master element and master element basis function info.
%
% notation:
%
%    T  = number of nodes in a single element
%    Q  = number of quadrature points in the interior of the master element
%    QE = number of quadrature points along the master element edges
%
%    INTPT   = zeros(Q,2)  = master element volume integration points
%    W       = zeros(Q,1)  = master element volume integration weights
%    PHI     = zeros(T,Q)  = values of local basis funcs at integ pts
%    PHIX    = zeros(T,Q)  = values of x-derivs of local basis funcs integ pts
%    PHIY    = zeros(T,Q)  = values of y-derivs of local basis funcs integ pts
%
%    INTPTE1 = zeros(QE,2) = type (2,3) EDGE integration points
%    WE1     = zeros(QE,1) = type (2,3) EDGE integration weights
%    PHIE1   = zeros(T,QE) = type (2,3) basis function values at EDGE integ pts
%
%    INTPTE2 = zeros(QE,2) = type (3,1) EDGE integration points
%    WE2     = zeros(QE,1) = type (3,1) EDGE integration weights
%    PHIE2   = zeros(T,QE) = type (3,1) basis function values at EDGE integ pts
%
%    INTPTE3 = zeros(QE,2) = type (1,2) EDGE integration points
%    WE3     = zeros(QE,1) = type (1,2) EDGE integration weights
%    PHIE3   = zeros(T,QE) = type (1,2) basis function values at EDGE integ pts
%
% author:   michael holst 101095
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% initialize

   T  = 3;   %%% triangles
   Q  = 1;   %%% one volume integration point
   QE = 2;   %%% two edge integration points

   INTPT   = zeros(Q,2);
   W       = zeros(Q,1);
   PHI     = zeros(T,Q);
   PHIX    = zeros(T,Q);
   PHIY    = zeros(T,Q);

   INTPTE1 = zeros(QE,2);
   WE1     = zeros(QE,1);
   PHIE1   = zeros(T,QE);

   INTPTE2 = zeros(QE,2);
   WE2     = zeros(QE,1);
   PHIE2   = zeros(T,QE);

   INTPTE3 = zeros(QE,2);
   WE3     = zeros(QE,1);
   PHIE3   = zeros(T,QE);

%%% define volume integration pts, weights, and master element basis functions

   INTPT(1:Q,1:2) = [1/3 1/3];
   W(1:Q,1:1)     = 1/2;
   PHI(1:T,1:Q)   = [1/3 ; 1/3 ; 1/3];
   PHIX(1:T,1:Q)  = [-1 ; 1 ; 0];
   PHIY(1:T,1:Q)  = [-1 ; 0 ; 1];

%%% define edge integration pts, weights, and master element basis functions

   c1 = 1/2;
   c2 = 1/(2*sqrt(3));

   INTPTE1(1:QE,1:2) = [ (1-c1+c2) (c1-c2) ; (1-c1-c2) (c1+c2) ];
   WE1(1:QE,1:1)     = [ c1 ; c1 ];
   PHIE1(1:T,1:QE)   = [ 0 0 ; (c1+c2) (c1-c2) ; (c1-c2) (c1+c2) ];

   INTPTE2(1:QE,1:2) = [ 0 (c1+c2) ; 0 (c1-c2) ];
   WE2(1:QE,1:1)     = [ c1 ; c1 ];
   PHIE2(1:T,1:QE)   = [ (c1-c2) (c1+c2) ; 0 0 ; (c1+c2) (c1-c2) ];

   INTPTE3(1:QE,1:2) = [ (c1-c2) 0 ; (c1+c2) 0 ];
   WE3(1:QE,1:1)     = [ c1 ; c1 ];
   PHIE3(1:T,1:QE)   = [ (c1+c2) (c1-c2) ; (c1-c2) (c1+c2) ; 0 0 ];
end

