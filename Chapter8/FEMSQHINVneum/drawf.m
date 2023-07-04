function  drawf(ELNODE,NODECO,BONODE,USOL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% purpose:
%
%    draw the finite element function
%
% author:   michael holst 101095
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% recover various problem dimensions

   [L,T]    = size(ELNODE);
   [M,two]  = size(NODECO);
   [M0,one] = size(BONODE);
   
%    ELNODE  = zeros(L,3)  = elements; node numbers given counter-clockwise
%    NODECO  = zeros(M,2)  = nodes; (x,y)-coordinates of all nodes
%    BONODE  = zeros(M0,1) = boundary nodes; essential boundary conditions

%%% setup for plot

   clf;
   hold off;
   
   
   X = NODECO(1:M,1);
   Y = NODECO(1:M,2);
   Z = USOL(1:M,1);

grid on
xv = linspace(min(X), max(X));
yv = linspace(min(Y), max(Y));
[XX,YY] = meshgrid(xv, yv);
ZZ = griddata(X,Y,Z,XX,YY,'linear');
ZZ = reshape(ZZ,size(XX));

%contour(X,Y,Z)

surface(XX,YY,ZZ)
axis on;

%%% plot the essential boundary nodes with circles

% plot3(NODECO(BONODE(1:M0),1),NODECO(BONODE(1:M0),2),USOL(BONODE(1:M0)),'or');

end

