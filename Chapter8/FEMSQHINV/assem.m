function [A,F]=assem(ELNODE,NODECO,BONODE,EDGE, ...
                     INTPT,W,PHI,PHIX,PHIY, ...
                     INTPTE1,WE1,PHIE1, ...
                     INTPTE2,WE2,PHIE2, ...
                     INTPTE3,WE3,PHIE3,atilde,ftilde)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% purpose:
%
%    Finite element matrix assembly; our notation and approach is similar to
%    Axelsson and Barker, Chapter 5.  The main noticable difference is our
%    use of an edge array for ALL elements along the boundary (not just
%    Neumann edges); this is for our refinement and domain decomposition
%    routines.  The edge array is given an additional (third) column, which
%    specifies the edge type.
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
%    A = zeros(M,M) = stiffness matrix
%    F = zeros(M,1) = load vector
%
% author:   michael holst 101095
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   t7 = clock;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


   %%% First recover various problem dimensions.
   [L,T]      = size(ELNODE);
   [M,two]    = size(NODECO);
   [M0,one]   = size(BONODE);
   [L1,three] = size(EDGE);
   [Q,two]    = size(INTPT);
   [QE,two]   = size(INTPTE1);
   N = M - M0;

   %%% Initialize global stiffness matrix and global load vector.
   A = spalloc(M,M,0);  
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
%%% Build the stiffness matrix and load vector at the same time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %  fprintf('   building stiffness matrix...')
   

   %%% (Go through all elements AT ONCE and do element integrals)

   %%% At each r loop iter, compute el stiffness matrix and el load vec entries
   for r=1:T

      %%% initialize element load entry fr
      fr(1:L,1) = zeros(L,1);

      %%% At each s loop iter, compute el stiffness matrix entries
      for s=1:r

         %%% initialize element matrix entry ars 
         ars(1:L,1) = zeros(L,1);

         %%% Go thru quad pts m=1:Q, computing contrib. to element matrix.
         for m=1:Q

            %%% Get quad pts by mapping master Gauss pts to our element
            xm(1:L,1)=f00(1:L,1)*INTPT(m,1)+f01(1:L,1)*INTPT(m,2)+b0(1:L,1);
            ym(1:L,1)=f10(1:L,1)*INTPT(m,1)+f11(1:L,1)*INTPT(m,2)+b1(1:L,1);
            
            %%% do the load vector while we are here...
            if (s==1)
               %%% We must evaluate function at current quadrature point m
               %% ftilde(1:L,1) = ff(xm,ym);
               thetam(1:L,1) = ftilde(1:L,1) * PHI(r,m) .* D(1:L,1);

               %%% Do quad using weights W, evaluating function at quad pts
               fr(1:L,1) = fr(1:L,1) + W(m) * thetam(1:L,1);
            end

            %%% Evaluate PDE & control coefficients at the quadrature points.
            %  atilde(1:L,1) = aa(xm,ym);
            btilde(1:L,1) = bb(xm,ym);

            %%% We must evaluate function at the current quadrature point m
            prx(1:L,1) = g00(1:L,1) * PHIX(r,m) + g10(1:L,1) * PHIY(r,m);
            pry(1:L,1) = g01(1:L,1) * PHIX(r,m) + g11(1:L,1) * PHIY(r,m);
            psx(1:L,1) = g00(1:L,1) * PHIX(s,m) + g10(1:L,1) * PHIY(s,m);
            psy(1:L,1) = g01(1:L,1) * PHIX(s,m) + g11(1:L,1) * PHIY(s,m);
            thetam(1:L,1) = D(1:L,1) .* ( ...
                  atilde(1:L,1) .* ( ...
                      prx(1:L,1).*psx(1:L,1) + pry(1:L,1).*psy(1:L,1) ) ...
                + btilde(1:L,1) * PHI(r,m) * PHI(s,m) ...
            );

            %%% Perform quadrature using the weights in W.
            ars(1:L,1) = ars(1:L,1) + W(m) * thetam(1:L,1);
         end

         %%% Add contribution to global stiffness matrix.
         i(1:L,1) = ELNODE(1:L,r);
         j(1:L,1) = ELNODE(1:L,s);
         for l=1:L
            if ( i(l) <= j(l) )
               A(i(l),j(l)) = A(i(l),j(l)) + ars(l);
            else
               A(j(l),i(l)) = A(j(l),i(l)) + ars(l);
            end
         end

      end

      %%% Add contribution to global load vector
      i(1:L,1) = ELNODE(1:L,r);
      for l=1:L
         F(i(l)) = F(i(l)) + fr(l);
      end

   end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Handle the natural (Neumann) boundary conditions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %  fprintf('   imposing natural conditions...');
   

   %%% Cycle through (sub)list of natural edges and do edge integrals
   EDGE_N = EDGE(find(EDGE(:,3)==2),:);
   [LL1,three] = size(EDGE_N); 
   for edge=1:LL1

      %%% Recover element num and local node numbers in element forming edge 
      el      = EDGE_N(edge,1); % Element containing the edge
      theEdge = EDGE_N(edge,2); % Local node in el opposite edge
      if (theEdge == 1)     %% Edge opposite vertex 1 ("Edge 0" in the notes)
         n1_loc = 2;
         n2_loc = 3;
      elseif (theEdge == 2) %% Edge opposite vertex 2 ("Edge 1" in the notes)
         n1_loc = 3;
         n2_loc = 1;
      elseif (theEdge == 3) %% Edge opposite vertex 3 ("Edge 2" in the notes)
         n1_loc = 1;
         n2_loc = 2;
      else 
         fprintf('problem in assembly...\n');
      end

      %%% global node numbers and corresponding coordinates
      n1 = ELNODE(el,n1_loc);
      n2 = ELNODE(el,n2_loc);
      X1 = NODECO(n1,1); Y1 = NODECO(n1,2);
      X2 = NODECO(n2,1); Y2 = NODECO(n2,2);

      %%% edge length as jacobian in transformation to master element
      ddd = sqrt( (X2 - X1)^2 + (Y2 - Y1)^2 );

      %%% For each r, compute single element load vector entry g_r.
      for rr=1:2
         ffrr = 0.0;
         if (rr==1)
            r=n1_loc;
         else
            r=n2_loc;
         end

         %%% Cycle thru the quadrature points (3 diff edge integ possible)
         for m=1:QE
            if (theEdge == 1)
               xm = f00(el) * INTPTE1(m,1) + f01(el) * INTPTE1(m,2) + b0(el);
               ym = f10(el) * INTPTE1(m,1) + f11(el) * INTPTE1(m,2) + b1(el);
               ffrr = ffrr + WE1(m) * hh(xm,ym) * PHIE1(r,m) * ddd;
            elseif (theEdge == 2)
               xm = f00(el) * INTPTE2(m,1) + f01(el) * INTPTE2(m,2) + b0(el);
               ym = f10(el) * INTPTE2(m,1) + f11(el) * INTPTE2(m,2) + b1(el);
               ffrr = ffrr + WE2(m) * hh(xm,ym) * PHIE2(r,m) * ddd;
            elseif (theEdge == 3)
               xm = f00(el) * INTPTE3(m,1) + f01(el) * INTPTE3(m,2) + b0(el);
               ym = f10(el) * INTPTE3(m,1) + f11(el) * INTPTE3(m,2) + b1(el);
               ffrr = ffrr + WE3(m) * hh(xm,ym) * PHIE3(r,m) * ddd;
            else
               fprintf('problem in assembly...\n');
            end
         end

         %%% Add contribution to global load vector
         i = ELNODE(el,r);
         F(i) = F(i) + ffrr;

      end

   end

   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Handle the essential (Dirichlet) boundary conditions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %  fprintf('   imposing essential conditions...')
  

   %%% First, create a list of non-boundary nodes.
   BONODE = sort(BONODE);
   tmp = [1:M];
   tmp(BONODE) = zeros(length(BONODE),1);
   NOTBONODE = find(tmp);

   %%% Now take single pass through boundary node array.
   for jj=1:M0

      %%% The current boundary node (global node number) we are handling.
      j = BONODE(jj);

      %%% The coordinates of this boundary node. 
      Xj(1:2) = NODECO(j,1:2);

      %%% Find all non-boundary nodes which interact with the boundary node.
      IR = find(A(j,NOTBONODE));
      [mm,nn]=size(IR); 
      JIR=NOTBONODE(IR);

      %%% If we found some, take some actions.
      if ( (mm > 0) & (nn > 0) )

         %%% Modify the load vector appropriately.
         F(JIR) = F(JIR) - gg(Xj(1),Xj(2)) * A(j,JIR)';

         %%% Decouple equation for this "known" value from other equations.
         if ( (mm == 1) & (nn == 1) )
            A(j,JIR) = 0;
         else
            A(j,JIR) = zeros(size(A(j,JIR)));
         end
      end

      %%% Find all non-boundary nodes which interact with the boundary node.
      IC = find(A(NOTBONODE,j));
      [mm,nn]=size(IC);
      JIC = NOTBONODE(IC);

      %%% If we found some, take some actions.
      if ( (mm > 0) & (nn > 0) )

         %%% Modify the load vector appropriately.
         F(JIC) = F(JIC) - gg(Xj(1),Xj(2)) * A(JIC,j);

         %%% Decouple equation for this "known" value from other equations.
         if ( (mm == 1) & (nn == 1) )
            A(JIC,j) = 0;
         else
            A(JIC,j) = zeros(size(A(JIC,j)));
         end
      end

      %%% Find all boundary nodes which interact with the boundary node.
      IR = find(A(j,BONODE(1:(jj-1))));
      [mm,nn]=size(IR);

      %%% Simply decouple the equations if cases were found.
      if ( (mm > 0) & (nn > 0) )
         JIR = BONODE(IR);
         if ( (mm == 1) & (nn == 1) )
            A(j,JIR) = 0;
         else
            A(j,JIR) = zeros(size(A(j,JIR)));
         end
      end

      %%% Create an identity equation for the known boundary node.
      F(j) = gg(Xj(1),Xj(2));
      A(j,j) = 1;

   end

  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Now create full system from symmetric construction above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %  fprintf('   forming full system from symmetric construction...')
   
   A_s = A;
   A_upper = triu(A,1);
   A = A_s + A_upper';
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute total time and return.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   time = etime(clock,t7);
%   fprintf('   total assembly time = %g\n',time);

end

