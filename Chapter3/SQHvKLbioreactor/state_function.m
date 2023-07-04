function [y,u] = state_function(F,nu,ya,x,v)
%
%    [y,u] = state_function(F,nu,ya,x,v)
%
% Calculates the state function by solving the 1st-order EDO with 
% averaged control with respect to the measure nu on the space v
%
%   y = value of the state function (two column vectors)
%   u = averaged control  (2 column vectors)
%
%   nu = 
%   doppia matrice !!!! matrix (Nx,Nv,2) of the PDF. Nv is the lenght of v, Nx the length
%        of the grid vector of independent variable x
%   ya = initial conditions of the state variable [ y1a  y1b]
%   x  = (row?) column vector of the grid of independent variable
%   v  = row vector of the grid value of the state space

dv = v(:,2)-v(:,1);

u(:,1) = dv(1)*nu(:,:,1)*v(1,:)';   % averaged control
u(:,2) = dv(2)*nu(:,:,2)*v(2,:)';

% Defines the function for the ODE integrator. y is treated as a scalar
% Makes the average w.r. to v, point-by-point in x

% FARE 2 INTERP1 IN UN UNICO HANDLE
FY = @(xx,y) [ interp1(x,Avg_F(F.n1,nu(:,:,1),x,y,v(1,:)),xx) ...
               interp1(x,Avg_F(F.n2,nu(:,:,2),x,y,v(2,:)),xx) ]';

% calculates the state function

% FY sono due vettori colonna uno per ogni stato ****
y = ode1(FY,x,ya);

% [~,y] = ode23(FY,x,ya);

% y is a two column vector 

end


% --- function for averaging ---
function E_Fy = Avg_F(F,nu,x,y,v)
%
%    E_Fp = Avg_F(F,nu,x,y,v)
%
%  Calculates the average of the rhs of the state equation
%
%  E_Fy = two column vector of the resulting averages
%
%  F = handle to the rhs of the state equation
%  nu= matrix nu(Nx,Nv,2) of the PDF
%  x = state domain 
%  y = scalar of the unknown state function
%  v = values of the controls (row vector)
%

dv = v(2)-v(1);

E_Fy = x;

% average point-by-point as a scalar product with respect to v.
% v must transposed in order to be set as column vector

% make y a row vector (of 2 elements)
y = y(:)';

for i=1:length(x) 
    
%     val_F = F(x(i),y,v)';
  E_Fy(i) =  nu(i,:)*F(x(i),y,v)'*dv;
  
end

end
