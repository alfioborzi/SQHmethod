function  drawa(ELNODE,NODECO,BONODE,INTPT,adiff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% purpose:
%
%    draw the gradient of the finite element function
%
% author:   michael holst 101095
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% recover various problem dimensions

   [L,T]    = size(ELNODE);
   [M,two]  = size(NODECO);
   [M0,one] = size(BONODE);
   [Q,two]  = size(INTPT);
   
   
%    ELNODE  = zeros(L,3)  = elements; node numbers given counter-clockwise
%    NODECO  = zeros(M,2)  = nodes; (x,y)-coordinates of all nodes
%    BONODE  = zeros(M0,1) = boundary nodes; essential boundary conditions
%    INTPT   = zeros(Q,2)  = master element volume integration points

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
   
%%% Get quad pts & position where gradient sol is evaluated
    for m=1:Q
    xm(1:L,1)=f00(1:L,1)*INTPT(m,1)+f01(1:L,1)*INTPT(m,2)+b0(1:L,1);
    ym(1:L,1)=f10(1:L,1)*INTPT(m,1)+f11(1:L,1)*INTPT(m,2)+b1(1:L,1);
    end
 
%%% setup for plot

   clf;
   hold off;
   
   %%% plot 
       
grid on
xv = linspace(min(xm), max(xm));
yv = linspace(min(ym), max(ym));
[X,Y] = meshgrid(xv, yv);
Z = griddata(xm,ym,adiff,X,Y,'cubic');
Z = reshape(Z,size(X));

%contour(X,Y,Z)

surface(X,Y,Z)

xlim([0 1]);
ylim([0 1]);
zlim([0.5 1.5]);

colorbar


end

