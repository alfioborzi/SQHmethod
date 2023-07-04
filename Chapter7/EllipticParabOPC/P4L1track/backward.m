%
%The function solves the adjoint equation of the elliptic problem in the corresponding paper,
% see Subsection 4.4 for details

function [p ] = backward(y,yd,A,OCP)

%Assemble Laplacian
h=(OCP.b-OCP.a)/(OCP.N+1);  
n=OCP.N+1;
N=(n+1)*(n+1)-4*n;
% v=ones(N,1);
% v_block=[ones(n-2,1);0];
% v_block=kron(ones(N/(n-1),1),v_block); 
% A_hilf=spdiags([(-1/h^2)*v_block (-1/h^2)*[v_block(N); v_block(1:N-1)]],[ -1  1 ],N,N);
% A=spdiags([(-1/h^2)*v  (4/h^2)*v  (-1/h^2)*v],[-(n-1)  0  n-1],N,N);
% A=A+A_hilf;

p=zeros(OCP.N+2,OCP.N+2);
h1=zeros(OCP.N,OCP.N);
h2=zeros(OCP.N,OCP.N);

for i=2:OCP.N+1
    %Calculate the "derivative" of the absolute value and the max function, 
    %see corresponding paper Subsection 4.4 for details
    for j=2:OCP.N+1
        if(y(i,j)>=yd(i,j))
            h1(i-1,j-1)=1;
        else
            h1(i-1,j-1)=-1;
        end
        
        if(y(i,j)>=0)
            h2(i-1,j-1)=1;
        end
    end
end

%Include the terms from the absolute value (L^1 trackig term) of the cost functional and max(y,0) 
%into the linear system of equations
A=A+spdiags([reshape(h2',[(n-1)*(n-1),1])],[0],N,N); 
h1_temp=reshape(h1',[(n-1)*(n-1),1]);
p_temp=A\h1_temp;
p_temp=reshape(p_temp,OCP.N,OCP.N)';

p(2:OCP.N+1,2:OCP.N+1)=p_temp;

end
