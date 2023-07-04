%
%The function solves the state equation of an elliptic bilinear equation with zero boundary conditions, 
%see corresponding paper, Subsection 4.2 for details

function [y] = forward_y(u,f,A,OCP)

n=OCP.N+1;          %number of intervals
N=(n+1)*(n+1)-4*n;  %Number of grid points where y is unknown

%State equation as a linear system of equation
y=zeros(n+1,n+1);
f=reshape(f',[N,1]);
u=reshape(u',[N,1]);
A=A+spdiags(u,0,N,N); %Include the bilinear term into the discretised equation
y_temp=A\f;
y_temp=reshape(y_temp,[n-1,n-1])';
y(2:n,2:n)=y_temp;
end

