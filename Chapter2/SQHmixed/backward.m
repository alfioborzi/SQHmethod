function [p] = backward(A,B1,B2,y,yd,u,ct,OCP)

T =OCP.T;
Nt=OCP.Nt;
dt=OCP.T/OCP.Nt;
alpha=OCP.alpha;

N=max(size(yd));
p(1:N,Nt+1)=-(y(:,end)-yd);
I=eye(N,N);

% Midpoint in the bilinear case 
for n=Nt:-1:1
    C=transpose( 2*A + (u(1,n)+u(1,n+1))*B1 + (u(2,n)+u(2,n+1))*B2);
    CT = 0.5*dt*alpha ... 
        *(max(0,gc(y(:,n),u(:,n))+ct(1,n)/alpha) * dgdy(y,u)' ... 
          + max(0,gc(y(:,n+1),u(:,n+1))+ct(1,n+1)/alpha) * dgdy(y,u)' );
    p(1:N,n)=(I-0.25*dt*C)\((I+0.25*dt*C)*p(1:N,n+1)-CT);
end


end

