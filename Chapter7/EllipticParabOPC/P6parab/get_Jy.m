%
%Calculates the cost functional value for given y,u and yd

function [ J ] = get_Jy(y,yd,u,OCP )

N=OCP.N;
Nt=OCP.Nt;
d=(OCP.b-OCP.a)*OCP.T/(N*Nt);

J=(0.5*d*transpose((y-yd))*(y-yd))+(0.5*OCP.alpha*d*(transpose(u)*u))+OCP.beta*d*sum(abs(u).*(abs(u)>OCP.s));
end
