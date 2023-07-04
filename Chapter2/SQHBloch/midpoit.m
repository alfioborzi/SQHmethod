function [y] = midpoit(A,B,B2,y0,u,OCP)
dt=OCP.dt;
Nt=OCP.Nt;

N=max(size(y0));
y=zeros(N,Nt+1);
y(1:N,1)=y0;


% initial conditions 
q = y(1:N,1); 
q1 = q ; 
q2 = q1 ;

for i = 1:Nt 
u_mid = ((u(i)+u(i+1)))/2; 
% implicit solve by fix point iteration 
for k = 1:100 
q_mid = (q2+q1)/2;
q2 = q1 + dt * stateDyn(A,B,B2,q_mid,u_mid) ; 
end

q1 = q2; 
y(1:N,i+1) = q1;
end

end