function p = adjoint_function(Fp,nu,y,yd,pT,v,x)
%
%   p = adjoint_function(Fp,nu,y,yd,pT,v,x)
%
%  Calculates the adjoint function with backward integration
%
%  p = calculated adjoint function, as column vector
%
%  Fp = function of the rhs of the adjoint equation
%  nu = matrix nu(Nx,Nv) of the PDF
%  y  = vector of the evalued state function
%  yd = vecotr of the evalued desired function
%  pT = terminal condition function evaluated in y(end)
%  v  = vector of the control values
%  x  = grid of the independent variable, as row vector


% Fp is the rhs of the equation p' = Fp(x,y,yd,p,v) 
% x,y,yd,v are evaluated. p is symbolic scalar evaluated in the ode1 function
%

% FP is a function of the averaged rhs of the adjoint equation
% w.r. to v
%
FP = @(xx,p) [ interp1(x,Avg_FP(Fp.n1,nu(:,:,1),x,y,yd,v(1,:),p),xx) ...
               interp1(x,Avg_FP(Fp.n2,nu(:,:,2),x,y,yd,v(2,:),p),xx) ]';

p = ode1(FP,flip(x),pT(y(end,:)));  % flip because the EDO is solved backward
% [~,p] = ode23(FP,flip(x),pT(y(end)));

p = flip(p);

% p is two column ?

end



% --- function for averaging ---
function E_Fp = Avg_FP(Fp,nu,x,y,yd,v,p)
%
%    E_Fp = Avg_FP(Fp,nu,x,y,yd,v,p)
%
%  Calculates the average of the rhs of the adjoint equation
%
%  Fp = function of the rhs of the adjoint equation
%  nu = matrix nu(Nx,Nv) of the PDF
%  x = state domain 
%  y = values of the state variable
%  yd= values of the desired function
%  v = values of the controls (row vector)
%  p = scalar of the unknown adjoint function

dv = v(2)-v(1);

E_Fp = x;  % allocates the output

% make p a row vector (of 2 elements)
% y = y(:)';

p = p(:)';

% average point-by-point as a scalar product with respect to v
% v must be set as column vector
%
for i=1:length(x)
    % !!!! adesso Fp e' una matrice. 
%     val_Fp = Fp(x(i),y(i),yd(i),p,v');
%     for n=1:2
%         E_Fp(i,n) = nu(i,:,n)*val_Fp(n)*dv(n);
%     end 

   val_Fp = Fp(x(i),y(i,:),yd(i,:),p,v)';
    
  E_Fp(i) =  nu(i,:)*val_Fp*dv;    
    
    
end

end
