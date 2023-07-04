function [y] = forward(A,B,y0,u,OCP)

T =OCP.T;
Nt=OCP.Nt;
dt=OCP.T/OCP.Nt;

N=max(size(y0));
y=zeros(N,Nt+1);
y(1:N,1)=y0;
I=eye(N,N);




% General Midpoint implementation
% initial conditions 
q = y(1:N,1); 
q1 = q ; 
q2 = q1 ;

for i = 1:Nt 
u_mid = ((u(i)+u(i+1)))/2; 
% implicit solve by fix point iteration 
% first guess by Euler explicit
q2 = q1 + dt * stateDyn(A,B,q1,u_mid) ; 
for k = 1:100 
q_mid = (q2+q1)/2;
q2 = q1 + dt * stateDyn(A,B,q_mid,u_mid) ; 
end

q1 = q2; 
y(1:N,i+1) = q1;
end

end
