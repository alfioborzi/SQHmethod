function  [Zm,Zx,Zy] = drawg(ELNODE,NODECO,BONODE,INTPT,USOL)
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

%%% cycle through the elements and draw the polygons

   C = [0.5 0.5 0.5];
   X = zeros(T,L);
   Y = zeros(T,L);
   Z = zeros(T,L);
   Zm = zeros(L,1);
   Zx = zeros(L,1);
   Zy = zeros(L,1);
   for element = 1 : L
      X(1:T,element) = NODECO(ELNODE(element,1:T),1);
      Y(1:T,element) = NODECO(ELNODE(element,1:T),2);
      Z(1:T,element) = USOL(ELNODE(element,1:T));
  %%% 3 vertices     
  a = [ X(1,element) Y(1,element) Z(1,element)];
  b = [ X(2,element) Y(2,element) Z(2,element)];
  c = [ X(3,element) Y(3,element) Z(3,element)]; 
  %%% area of triangle - element      
  % area = .5*abs(a(1)*(b(2)-c(2)) + b(1)*(c(2)-a(2)) + c(1)*(a(2)-b(2))); 
  % normal = cross(a-b, a-c);
   
   %%% use z = A x + B y + C to get Zx = A and Zy = B
   
   M = [ a(1) a(2) 1 ; b(1) b(2) 1; c(1) c(2) 1];
   F = [ a(3) b(3) c(3)]' ; 
   V = M \ F ;
   
   %%% compute sol at quad pts
   
   Zm(element,1) = V(1)*xm(element,1) + V(2)*ym(element,1) + V(3);
   
   %%% Check by comparison 
   % [Zm(element,1) (a(3)+b(3)+c(3))/3.]
  
   %%% The gradient components at (xm(element,1),ym(element,1))
   
   Zx(element,1) = V(1);
   Zy(element,1) = V(2);
 
   end
   
   %%% plot the gradient field (gx,gy)
   
   [xx,yy] = meshgrid(xm(1:L,1),ym(1:L,1));
   [gx,gy] = meshgrid(Zx(1:L,1),Zy(1:L,1));
   
   quiver(xx,yy,gx,gy)
   
   axis off;
   hold on;

end

