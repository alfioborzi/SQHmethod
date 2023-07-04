function y = forward(y0,u,T,OCP,DYN)

Nt=OCP.Nt;
dt=T/Nt;

TT=T;
NNt=Nt;

N=DYN.ND;
y=zeros(N,Nt+1);
y(1:N,1)=y0;

% General Midpoint implementation
% initial conditions 
q = y(1:N,1); 
q1 = q ; 

for i = 1:Nt 
    tt = (i-1)*dt; 
u_mid = ((u(:,i)+u(:,i+1)))/2; 

% implicit solve by fix point iteration 
% first guess by Euler explicit
q2 = q1 + dt * stateDyn(DYN,q1,u_mid) ; 
for k = 1:100 
q_mid = (q2+q1)/2;
q2 = q1 + dt * stateDyn(DYN,q_mid,u_mid) ; 
end

q1 = q2; 
y(1:N,i+1) = q1;

    
end

end
