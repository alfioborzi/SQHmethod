function [ v ] = argminH(i,j,dKU,y,yd,p,u,epsilon,f,OCP )

posi=(j-1)*(OCP.N+1)+1+i;


% Disc_KU=[-10,-6,-2,0,2,6,10];
Du = (OCP.u_up - OCP.u_lo)/dKU;
Disc_KU=(OCP.u_lo:Du:OCP.u_up)';

[~,pos]=min((OCP.alpha/2)*Disc_KU.^2+OCP.beta*abs(Disc_KU).*(abs(Disc_KU)>OCP.s)+(Disc_KU+f(posi)).*p(posi)+epsilon*((Disc_KU-u).^2));
    
v=Disc_KU(pos);

end

