function J = functional_J(L,gb,nu,y,yd,x,v)
%
%    J = functional_J(L,gb,nu,y,yd,x,v)
%
%   Evaluates the cost function 
%
%   J = value of the functional
%
%   L  = cost model function
%   gb = is the terminal cost function evaluated in y(end)
%   nu = matrix nu(Nx,Nv) of the PDF
%   y  = vector of the state function
%   x  = grid vector of the independent variable
%   v  = grid (row) vector of the control values 
%
dx = x(2)-x(1);
dv = v(2)-v(1);

% Evaluate the functional according to the model
% averages in v with nu, then quadrature in v and x
% v must be transposed in order to set as column vector
% 
J=0;

for i=1:length(x)
    J = J + nu(i,:)*L(x(i),y(i),yd(i),v')*dv;
end

J = J*dx + gb(y(end));   % add the terminal cost


