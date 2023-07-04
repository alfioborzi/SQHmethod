function [Jt]=jfunctional(ELNODE,NODECO,BONODE,EDGE, ...
                     INTPT,W,PHI,atilde,Ut,ydtilde)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% notation:
%
%    T  = number of nodes in a single element
%    Q  = number of quadrature points in the interior of the master element
%    QE = number of quadrature points along the master element edges
%
%    L  = number of elements in mesh
%    M  = number of nodes in mesh
%    M0 = number of nodes on boundary with essential BC
%    N  = number of interior points + number of neumann BC nodes
%       = M - M0 = total number of unknowns
%    L1 = total number of element edges on boundary (BOTH diri. and neum.)
%
%    ELNODE  = zeros(L,3)  = elements; node numbers given counter-clockwise
%    NODECO  = zeros(M,2)  = nodes; (x,y)-coordinates of all nodes
%    BONODE  = zeros(M0,1) = boundary nodes; essential boundary conditions
%    EDGE    = zeros(L1,3) = list of all edges along boundary;
%              (:,1) = triangle which has this external edge
%              (:,2) = local node number opposite the edge
%              (:,3) = type: 1=dirichlet, 2=neumann
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
% Output:
%
%    J = Tikhonov functional
%
% author:   michael holst 101095
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


   %%% First recover various problem dimensions.
   [L,T]      = size(ELNODE);
   [M,two]    = size(NODECO);
   [M0,one]   = size(BONODE);
   [L1,three] = size(EDGE);
   [Q,two]    = size(INTPT);
   
   N = M - M0;

   %%% Initialize global vector.
   
   F = zeros(M,1);

   %%% Get coordinates of vertices of each element.
   Xi0(1:L,1) = NODECO(ELNODE(1:L,1),1);
   Xi1(1:L,1) = NODECO(ELNODE(1:L,2),1);
   Xi2(1:L,1) = NODECO(ELNODE(1:L,3),1);
   Yi0(1:L,1) = NODECO(ELNODE(1:L,1),2);
   Yi1(1:L,1) = NODECO(ELNODE(1:L,2),2);
   Yi2(1:L,1) = NODECO(ELNODE(1:L,3),2);

   %%% Build affine transformations (from master element to arbitrary).
   f00(1:L,1) = Xi1(1:L,1) - Xi0(1:L,1);
   f01(1:L,1) = Xi2(1:L,1) - Xi0(1:L,1);
   f10(1:L,1) = Yi1(1:L,1) - Yi0(1:L,1);
   f11(1:L,1) = Yi2(1:L,1) - Yi0(1:L,1);
   b0(1:L,1)  = Xi0(1:L,1);
   b1(1:L,1)  = Yi0(1:L,1);

   %%% Build the jacobians.
   D  (1:L,1) = abs( f00 .* f11 - f01 .* f10 );
   Di (1:L,1) = 1.0 ./ D(1:L,1);

   %%% Build inverse transformations (from arbitrary element to master).
   g00(1:L,1) =  Di(1:L,1) .* f11(1:L,1);
   g01(1:L,1) = -Di(1:L,1) .* f01(1:L,1);
   g10(1:L,1) = -Di(1:L,1) .* f10(1:L,1);
   g11(1:L,1) =  Di(1:L,1) .* f00(1:L,1);

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Build the load vector at the same time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 %%% the computed state on the quad pts / elements
   ytilde(1:L,1) = zeros(L,1);
   for r=1:T 
      i(1:L,1) = ELNODE(1:L,r);
      for l=1:L         
          ytilde(l,1) = ytilde(l,1) + Ut(i(l))/3.;       
      end 
   end
   
   a0diff = ones(L,1);
   
   %%% the integrand 
   ftilde(1:L,1) = zeros(L,1);
   ftilde(1:L,1) = 0.5*(ytilde(1:L,1)-ydtilde(1:L,1)).^2 ... 
       + gcost(atilde(1:L,1)-a0diff(1:L,1));

   %%% At each r loop iter, compute el load vec entries
   for r=1:T

      %%% initialize element load entry fr
      fr(1:L,1) = zeros(L,1);

         %%% Go thru quad pts m=1:Q, computing contrib. to element matrix.
         for m=1:Q

            %%% Get quad pts by mapping master Gauss pts to our element
            xm(1:L,1)=f00(1:L,1)*INTPT(m,1)+f01(1:L,1)*INTPT(m,2)+b0(1:L,1);
            ym(1:L,1)=f10(1:L,1)*INTPT(m,1)+f11(1:L,1)*INTPT(m,2)+b1(1:L,1);
                      
               %%% We must evaluate function at current quadrature point m
               %% ftilde(1:L,1) = ff(xm,ym);
               
               thetam(1:L,1) = ftilde(1:L,1) * PHI(r,m) .* D(1:L,1);

               %%% Do quad using weights W, evaluating function at quad pts
               fr(1:L,1) = fr(1:L,1) + W(m) * thetam(1:L,1);
            
         end

      %%% the weithed integrand on nodes
      i(1:L,1) = ELNODE(1:L,r);
      for l=1:L
         F(i(l)) = F(i(l)) + fr(l);
      end

   end
   
  
   % the integral
   Jt = sum(F);


end

