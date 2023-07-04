function  [xm,ym] = draw(ELNODE,NODECO,BONODE,INTPT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% purpose:
%
%    draw the finite element mesh
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

%%% plot all nodes with dots

   plot(NODECO(1:M,1),NODECO(1:M,2), '.b')
   axis off;
   hold on;

%%% cycle through the elements and draw the edges

   for ell = 1:L
      Xil(1:T,1:2) = NODECO(ELNODE(ell,1:T),1:2);
      plot([Xil(1:T,1)' Xil(1,1)], [Xil(1:T,2)' Xil(1,2)]);
   end

%%% plot the essential boundary nodes with circles

   plot(NODECO(BONODE(1:M0),1), NODECO(BONODE(1:M0),2), 'or');

end

