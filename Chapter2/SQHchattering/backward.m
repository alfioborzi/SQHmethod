function [p] = backward(A,B,y,yd,u,OCP)

T =OCP.T;
Nt=OCP.Nt;
dt=OCP.T/OCP.Nt;

N=1;
p(1:N,Nt+1)=0;
I=eye(N,N);


% General Midpoint implementation
% initial conditions 
q = p(1:N,Nt+1); 
q1 = q ; 
q2 = q1 ;

for i = Nt:-1:1
u_mid = ((u(i)+u(i+1)))/2; 
y_mid = (y(1:N,i)+y(1:N,i+1))/2;
ydm   = (yd(1:N,i)+yd(1:N,i+1))/2;
% implicit solve by fix point iteration 
% first guess by Euler explicit
q2 = q1 + dt * adjointDyn(A,B,q1,y_mid,ydm,u_mid);
for k = 1:100 
q_mid = (q2+q1)/2;
q2 = q1 + dt * adjointDyn(A,B,q_mid,y_mid,ydm,u_mid) ; 
end

q1 = q2; 
p(1:N,i) = q1;
end

end

