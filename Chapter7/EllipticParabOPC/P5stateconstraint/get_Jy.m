%
%Calculates the cost functional value for given y,u and yd

function [ J ] = get_Jy(u,y,yd,OCP)

y=reshape(y,[(OCP.N+2)*(OCP.N+2),1]);
yd=reshape(yd,[(OCP.N+2)*(OCP.N+2),1]);
c=OCP.c;
gamma=OCP.gamma;
u=reshape(u,[OCP.N*OCP.N,1]);
h=(OCP.b-OCP.a)/(OCP.N+1);

J=(0.5*((y-yd)')*(y-yd)+0.5*OCP.alpha*(u')*u+OCP.beta*sum(abs(u).*(abs(u)>OCP.s))...
    +gamma*sum((max(0,y-c)).^3))*h*h;
end


