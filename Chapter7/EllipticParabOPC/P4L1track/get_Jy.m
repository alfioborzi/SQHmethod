%
%Calculates the cost functional value for given y,u and yd

function [ J ] = get_Jy(u,y,yd,OCP)

y=reshape(y,[(OCP.N+2)*(OCP.N+2),1]);
yd=reshape(yd,[(OCP.N+2)*(OCP.N+2),1]);
u=reshape(u,[OCP.N*OCP.N,1]);
h=(OCP.b-OCP.a)/(OCP.N+1);

J=(sum(abs(y-yd))+0.5*OCP.alpha*(u')*u+OCP.beta*sum(log(1+abs(u))))*h*h;
end


