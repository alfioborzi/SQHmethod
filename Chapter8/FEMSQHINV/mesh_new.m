function [ELNODE,NODECO,BONODE,EDGE] = mesh(n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% purpose:
%
%    define a simplical finite element mesh.
%
% domain:
%
%    fixed size uniform mesh on the unit square.
%
% notation:
%
%    L  = number of elements in mesh
%    M  = number of nodes in mesh
%    M0 = number of nodes on boundary with essential BC
%    N  = number of interior points + number of neumann BC nodes
%       = M - M0 = total number of unknowns
%    L1 = total number of element edges on boundary
%
%    ELNODE = zeros(L,3)  = elements; node numbers given counter-clockwise
%    NODECO = zeros(M,2)  = nodes; (x,y)-coordinates of all nodes
%    BONODE = zeros(M0,1) = boundary nodes; essential boundary conditions
%    EDGE   = zeros(L1,3) = list of all edges along boundary;
%            (:,1) = triangle which has this external edge
%            (:,2) = local node number opposite the edge
%            (:,3) = type: 1=dirichlet, 2=neumann
%
% all counterclockwise
%
% author:   michael holst 101095
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% initialization

%% a square domain & Dirichlet boundary conditions

% number of nodes on one side 
   S = 10; 
% number of nodes
   M  = S * S;
% number of elements
   L  = (S-1)*(S-1)*2;
% number of bc Dirichlet nodes
   M0 = 4*(S-1);  
% number of edges on boundary   
   L1 = 4*(S-1);
   
   %   M0 = 10;  % mix

%%% more initialization 

   ELNODE = zeros(L,3);
   NODECO = zeros(M,2);
   BONODE = zeros(M0,1);
   EDGE   = zeros(L1,3);

%%% node array
   IC=0; 
   for J=1:S
   for I=1:S
       IC = IC+1;
       NODECO(IC,1:2)  = [I-1 J-1] / (S-1) ;
   end
   end

 %  NODECO = NODECO / (S-1);  %%% so we are on the unit square...

%%% element array
  
   IC=0;
   JC=0;
%  corresp. J  make a strip of lower & upper triangles
   for J=1:S-1
       for I=1:S-1
       IC=IC+1;
       ELNODE(IC,1:3)    = [ I+JC  I+1+JC  I+S+JC];
       IC=IC+1;
       ELNODE(IC,1:3)  = [ I+1+JC  I+1+S+JC  I+S+JC];
       end
       JC=JC+S;
   end

   %%% edge array 
   % south 
   IC=0;
   for I=1:S-1
       IC=IC+1;
       IE=2*I-1;
   EDGE(IC,1:3)  = [ IE 3 1];
   end
   % east
  
   for J=1:S-1
       IC=IC+1;
       IE=2*(S-1)*J;
   EDGE(IC,1:3)  = [ IE 3 1];
   end
   % north
   
   for I=1:S-1
       IC=IC+1;
       IE=2*(S-1)*(S-1)-(I-1)*2;
   EDGE(IC,1:3)  = [ IE 1 1];
   end
   % west 
    for J=1:S-1
       IC=IC+1;
       IE=2*(S-1)*(S-1)-2*(S-1)+1-(J-1)*(S-1)*2;
   EDGE(IC,1:3)  = [ IE 2 1];
   end
  
    
   %%% dirichlet boundary array
   IC=0;
   % south 
   for I=1:S
       IC=IC+1;
   BONODE(IC)  =  I;
   end
   % west and east
   for J=2:S-1
       % west 
   IF=(J-1)*S+1;
   IC=IC+1;
   BONODE(IC)  =  IF;  
   %east
   IF=IF+S-1;
   IC=IC+1;
   BONODE(IC)  =  IF;
   end
   % north  
   for I=1:S
   IF=(S-1)*S+I;
   IC=IC+1;
   BONODE(IC)  =  IF;
   end
   

end
