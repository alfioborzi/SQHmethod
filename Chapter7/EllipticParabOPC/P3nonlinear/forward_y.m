%
%The function solves the state equation of an elliptic linear equation with zero boundary conditions 
%and a nonlineartiy y^3

function [ yit ] = forward_y(y,f,A,OCP)
% f is u 
%Assemble the Laplacian
h=(OCP.b-OCP.a)/(OCP.N+1); 
n=OCP.N+1;
N=(n+1)*(n+1)-4*n;
% v=ones(N,1);
% v_block=[ones(n-2,1);0];
% v_block=kron(ones(N/(n-1),1),v_block);
% A_hilf=spdiags([(-1/h^2)*v_block (-1/h^2)*[v_block(N); v_block(1:N-1)]],[ -1  1 ],N,N);
% A=spdiags([(-1/h^2)*v  (4/h^2)*v  (-1/h^2)*v],[-(n-1)  0  n-1],N,N);
% A=A+A_hilf;

yit=zeros(OCP.N+2,OCP.N+2);
yit = y; 

y_temp=y(2:n,2:n);

 for k=1:10 %Inner loop of Picard iterations
       
    %State equation as a linear system of equation   - Picard
    u_temp = f - y_temp.^3; 
    ut = reshape(u_temp',[(n-1)*(n-1),1]);
    y_temp=A\ut;
    y_temp=reshape(y_temp,[n-1,n-1])';
    
    yit(2:n,2:n) = y_temp; 

    %Determine residum of the left hand side and the right hand side f
    res=norm(A*reshape(yit(2:OCP.N+1,2:OCP.N+1)',[OCP.N*OCP.N,1])+reshape(((yit(2:OCP.N+1,2:OCP.N+1)).^3)',[OCP.N*OCP.N,1])...
        -reshape(f',[OCP.N*OCP.N,1]))*h; 
    
    if(res<10^-6)   %If residuum is less than 10^-6, then stop
        break;
    end
 end
end

