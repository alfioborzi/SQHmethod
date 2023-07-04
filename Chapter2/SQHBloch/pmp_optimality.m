function [] = pmp_optimality(u,y0,yd,A,B1,B2,gv,OCP)
N=max(size(y0));

T =OCP.T;
Nt=OCP.Nt;
dt=OCP.T/OCP.Nt;

y=zeros(N,Nt);
p=zeros(N,Nt);

y=forward(A,B1,B2,gv,y0,u,OCP);
p=backward(A,B1,B2,y,yd,u,OCP);


% The HP function with the SQH solution u - using H_epsilon - 
for i=1:Nt
    uz(1)=u(1,i);
    uz(2)=u(2,i);
    Hopt(i)=HPfunction(A,B1,B2,gv,y(:,i),p(:,i),uz',OCP);        
end

% The HP function with u=argmax - epsilon =0 - 
for i=1:Nt
    uz=argmaxH(y(:,i),p(:,i),A,B1,B2,gv,OCP);
    Hnew(i)= HPfunction(A,B1,B2,gv,y(:,i),p(:,i),uz',OCP);   
end

%Determine the maximum difference 
fprintf('Max DH: %d\n',max(max(Hnew-Hopt)));
fprintf('Fraction of grid points where 0<=DH<=10^-l is fulfilled:\n')
for l=[2,3,4,5,8,15]
count=0;
for i=1:Nt
    %Count the numer of grid points where the solution u is PMP optimal up to the tolerance 10^-l
    if (Hnew(i)-Hopt(i)<=10^(-l))      
        count=count+1;
    end
end

ratio=count/Nt;
fprintf('l=%i: %f\t',l,ratio);
end
fprintf('\n');
end

