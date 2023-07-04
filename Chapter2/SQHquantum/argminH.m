function [ v ] = argminH(y,u,p,A,B,epsilon,OCP )
beta=OCP.beta;
s=OCP.s;

u(1)=min(max((-p'*B*y-beta)/OCP.nu,s),OCP.umax);
u(2)=min(max((-p'*B*y+beta)/OCP.nu,OCP.umin),-s);
u(3)=min(max((-p'*B*y)/OCP.nu,-s),s);
H(1)=OCP.nu*0.5*u(1)^2 + OCP.beta*(abs(u(1))>s).*abs(u(1))+p'*B*y*u(1);
H(2)=OCP.nu*0.5*u(2)^2 + OCP.beta*(abs(u(2))>s).*abs(u(2))+p'*B*y*u(2);
H(3)=OCP.nu*0.5*u(3)^2 +p'*B*y*u(3);
[~,pos]=min(H);
v=u(pos);




% 
%     num_Int=100;
%     Disc_KU=linspace(OCP.umin,OCP.umax,num_Int);
%     Disc_KU=sort([0,Disc_KU]);
%     dist=inf;
% 
%     while(dist>10^-14)
%     [~,pos]=min((OCP.nu/2)*Disc_KU.^2+OCP.beta*(abs(Disc_KU)>0) + (p')*A*y+(p')*B*y.*Disc_KU+epsilon*((Disc_KU-u).^2));
%     v=Disc_KU(pos);
% 
%     if (pos==1)
%     Disc_KU=linspace(Disc_KU(pos),Disc_KU(pos+1),num_Int);
%     dist=(Disc_KU(pos+1)-Disc_KU(pos));
%     elseif  (pos==max(size(Disc_KU)))
%     Disc_KU=linspace(Disc_KU(pos-1),Disc_KU(pos),num_Int);
%     dist=(Disc_KU(pos)-Disc_KU(pos-1));
%     else
%     Disc_KU=linspace(Disc_KU(pos-1),Disc_KU(pos+1),num_Int);
%     dist=Disc_KU(pos+1)-Disc_KU(pos-1);
%     end
%     end



end


