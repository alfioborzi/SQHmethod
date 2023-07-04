%
%The function solves the adjoint equation of the elliptic problem in the corresponding paper,
% see Subsection 4.3 for details

function [p] = backward(y,yd,A,OCP )

n=OCP.N+1;          %number of intervals
N=(n+1)*(n+1)-4*n;  %Number of grid points where y is unknown

%Adjoint equation as a linear system of equation
p=zeros(n+1,n+1);
A=A+spdiags([reshape((3*(y(2:n,2:n).^2))',[(n-1)*(n-1),1])],[0],N,N); %Include the part that comes from the nonlinearity y^3
y_temp=y-yd;
y_temp=y_temp(2:n,2:n);
y_temp=reshape(y_temp',[(n-1)*(n-1),1]);
p_temp=A\y_temp;
p_temp=reshape(p_temp,[n-1,n-1])';
p(2:n,2:n)=p_temp;
end

