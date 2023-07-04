%
%The function solves the state equation of an elliptic linear equation with zero boundary conditions, 
%see corresponding paper, Subsection 4.1 for details

function [y] = forward_y(u,A,OCP)

n=OCP.N+1;          %number of intervals
N=(n+1)*(n+1)-4*n;  %Number of grid points where y is unknown

%State equation as a linear system of equation
y=zeros(n+1,n+1);
u=reshape(u',[N,1]);
y_temp=A\u;
y_temp=reshape(y_temp,[n-1,n-1])';
y(2:n,2:n)=y_temp;
end

