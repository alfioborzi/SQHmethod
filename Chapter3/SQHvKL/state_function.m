function [y,u] = state_function(F,nu,ya,x,v)
%
%    [y,u] = state_function(F,nu,ya,x,v)
%
% Calculates the state function by solving the 1st-order EDO with 
% averaged control with respect to the measure nu on the space v
%
%   y = value of the state function (column vector)
%   u = averaged control  (column vector)
%
%   nu = matrix (Nx,Nv) of the PDF. Nv is the lenght of v, Nx the length
%        of the grid vector of independent variable x
%   ya = initial condition of the state variable
%   x  = row vector of the grid of independent variable
%   v  = row vector of the grid value of the state space

dv = v(2)-v(1);

u = dv*nu*v';   % averaged control

% Defines the function for the ODE integrator. y is treated as a scalar
% Makes the average w.r. to v, point-by-point in x
FY = @(xx,y) interp1(x,Avg_F(F,nu,x,y,v),xx);

% calculates the state function
y = ode1(FY,x,ya);

% y is a column vector

end


% --- function for averaging ---
function E_Fy = Avg_F(F,nu,x,y,v)
%
%    E_Fp = Avg_F(F,nu,x,y,v)
%
%  Calculates the average of the rhs of the state equation
%
%  E_Fy = column vector of the resulting average
%
%  F = handle to the rhs of the state equation
%  nu= matrix nu(Nx,Nv) of the PDF
%  x = state domain 
%  y = scalar of the unknown state function
%  v = values of the controls (row vector)
%

dv = v(2)-v(1);

E_Fy = x;  % allocates the output

% average point-by-point as a scalar product with respect to v.
% v must transposed in order to be set as column vector
%
for i=1:length(x)
    E_Fy(i) =  nu(i,:)*F(x(i),y,v')*dv;
end

end
