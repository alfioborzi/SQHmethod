function [ J ] = get_J(y,yd,u,OCP )
J=0;
N=OCP.N;
Nt=OCP.Nt;
d=(OCP.b-OCP.a)*OCP.T/(N*Nt);

J=(0.5*d*transpose((y-yd))*(y-yd))+(0.5*OCP.alpha*d*(transpose(u)*u))+OCP.beta*d*sum(abs(u).*(abs(u)>OCP.s));
end


