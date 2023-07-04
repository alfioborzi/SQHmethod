%
%The function solves the state equation of an elliptic linear equation with zero boundary conditions 
%and with max(y,0) in the state equation


function [ yit ] = forward_y(y,f,A,OCP)

%
 h=(OCP.b-OCP.a)/(OCP.N+1); 
 n=OCP.N+1;
 N=(n+1)*(n+1)-4*n;
% 

 yit=zeros(OCP.N+2,OCP.N+2);
 yit = y; 
 
 y_temp=y(2:n,2:n);
 
 for k=1:10
 
    %State equation as a linear system of equation   - Picard
    u_temp = f - max(0,y_temp); 
    ut = reshape(u_temp',[(n-1)*(n-1),1]);
    y_temp=A\ut;
    y_temp=reshape(y_temp,[n-1,n-1])';
    
    yit(2:n,2:n) = y_temp; 
 
 
 
 %Determine the residuum of the state equation
 res=norm(A*reshape(yit(2:OCP.N+1,2:OCP.N+1)',[OCP.N*OCP.N,1])+max(0,reshape(yit(2:OCP.N+1,2:OCP.N+1)',[OCP.N*OCP.N,1]))...
     -reshape(f',[OCP.N*OCP.N,1]))*h; 
     if(res<10^-6)  
         %Stop if the residuum of the state equation is less than 10^-8
         break;
     end
 end
 
end

