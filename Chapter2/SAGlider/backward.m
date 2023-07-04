function [p] = backward(y,u,T,OCP,DYN)

global A lambda

Nt=OCP.Nt;
dt=T/Nt;

NNt = Nt;

N=DYN.ND;
p=zeros(N,OCP.Nt+1);


% Terminal condition
p(1,NNt+1)= 1;
p(2,NNt+1)= - A*y(2,Nt+1) - lambda ; 
%p(2,NNt+1)= - cot(y(4,NNt+1)) ; 
p(3,NNt+1)=0;
p(4,NNt+1)=0;


% General Midpoint implementation
% initial conditions 
q = p(1:N,NNt+1); 
q1 = q ; 


for i = NNt:-1:1
u_mid = ((u(:,i)+u(:,i+1)))/2; 
y_mid = (y(1:N,i)+y(1:N,i+1))/2;
% implicit solve by fix point iteration 
% first guess by Euler explicit
q2 = q1 + dt * adjointDyn(DYN,q1,y_mid,u_mid);
for k = 1:100 
q_mid = (q2+q1)/2;
q2 = q1 + dt * adjointDyn(DYN,q_mid,y_mid,u_mid) ; 
end

q1 = q2; 
p(1:N,i) = q1;
end

end

