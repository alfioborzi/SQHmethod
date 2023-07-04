function [] = numeric_optimality(u,y0,yd,A,B,OCP)
N=max(size(y0));
Nt=round((OCP.T/OCP.dt))+1;
y=zeros(N,Nt);
p=zeros(N,Nt);
s=OCP.s;
y=forward(A,B,y0,u,OCP);
p=backward(A,B,y(:,end),yd,u,OCP);



for i=1:Nt
       Hopt(i)=(OCP.nu/2)*u(1,i)^2 + OCP.beta*(abs(u(1,i))>s).*abs(u(1,i)) + (p(:,i)')*(A+u(1,i)*B)*y(:,i);         
end

u_old=ones(1,Nt+1);
for i=1:Nt
    u(1,i)=argminH(y(:,i),u(i),p(:,i),A,B,0,OCP);
    Hnew(i)=(OCP.nu/2)*u(1,i)^2 + OCP.beta*(abs(u(1,i))>s).*abs(u(1,i)) + (p(:,i)')*(A+u(1,i)*B)*y(:,i);
    
%         epsi=0;
%         uz(1)=min(max(0,(2*epsi*u_old(1,i)-(p(:,i)')*B*y(:,i)-OCP.beta)/(2*epsi+OCP.nu)),OCP.umax);
%         uz(2)=min(max(OCP.umin,(2*epsi*u_old(1,i)-(p(:,i)')*B*y(:,i)+OCP.beta)/(2*epsi+OCP.nu)),0);
%         H(1)=(OCP.nu/2)*uz(1)^2 + OCP.beta*abs(uz(1)) + (p(:,i)')*(A+uz(1)*B)*y(:,i) + epsi*(uz(1)-u_old(1,i))^2;
%         H(2)=(OCP.nu/2)*uz(2)^2 + OCP.beta*abs(uz(2)) + (p(:,i)')*(A+uz(2)*B)*y(:,i) + epsi*(uz(2)-u_old(1,i))^2;
%         [Hnew(i),pos]=min(H);
%         u(1,i)=uz(pos);
      
end

%Determine the maximum difference 
fprintf('Max DH: %d\n',max(max(Hopt-Hnew)));
fprintf('Fraction of grid points where 0<=DH<=10^-l is fulfilled:\n')
for l=[2,4,6,8,10]
count=0;
for i=1:Nt
    %Count the numer of grid points where the solution u is PMP optimal up to the tolerance 10^-l
    if (Hopt(i)-Hnew(i)<=10^(-l))      
        count=count+1;
    end
end

ratio=count/Nt;
fprintf('l=%i: %f\t',l,ratio);
end
fprintf('\n');
end

