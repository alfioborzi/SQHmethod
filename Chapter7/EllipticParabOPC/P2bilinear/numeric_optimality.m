%
% Check PMP optimality of the SQH solution 
%

function [] = numeric_optimality(u,f,yd,OCP)

%Discretised Laplacian
h=(OCP.b-OCP.a)/(OCP.N+1);
n=OCP.N+1;
N=(n+1)*(n+1)-4*n;
v=ones(N,1);
v_block=[ones(n-2,1);0];
v_block=kron(ones(N/(n-1),1),v_block);
A_hilf=spdiags([(-1/h^2)*v_block (-1/h^2)*[v_block(N); v_block(1:N-1)]],[ -1  1 ],N,N);
A=spdiags([(-1/h^2)*v  (4/h^2)*v  (-1/h^2)*v],[-(n-1)  0  n-1],N,N);
A=A+A_hilf;

%Calculate the state variable with the calculated control u
y=forward_y(u,f,A,OCP);
%Calculate the adjoint variable with the control u
p=backward(u,y,yd,A,OCP);

n=OCP.N; 
for i=1:n
    for j=1:n
        %Calculate pointwise the value of the Hamiltonian where y,u,p returned by the SQH method are utilised
        Hopt(i,j)=(OCP.alpha/2)*u(i,j).^2+OCP.beta*abs(u(i,j)).*(abs(u(i,j))>OCP.s)-u(i,j).*y(i+1,j+1)*p(i+1,j+1); 
        %Determine for y and p and the returned u again the value for the
        %control where the Hamiltonian takes its minimum
        unew(i,j)=argMinH(u(i,j),y(i+1,j+1),p(i+1,j+1),0,OCP); 
    end
end

for i=1:n
    for j=1:n
        %Calculate pointwise the value of the Hamiltonian at unew
        Hnew(i,j)=(OCP.alpha/2)*unew(i,j).^2+OCP.beta*abs(unew(i,j)).*(abs(unew(i,j))>OCP.s)-unew(i,j).*y(i+1,j+1)*p(i+1,j+1); %Calculate pointwise the value of the Hamiltonian at the minimum at unew
    end
end

%Determine the maximum difference 
fprintf('Max DH: %d\n',max(max(Hopt-Hnew)));
fprintf('Fraction of grid points where 0<=DH<=10^-l is fulfilled:\n')
for l=2:2:12
count=0;
for i=1:n
    for j=1:n
         %Count the numer of grid points where the solution u is PMP optimal up to the tolerance 10^-l
        if (Hopt(i,j)-Hnew(i,j)<=10^(-l))      
            count=count+1;
        end
    end
end

ratio=count/(n*n);
fprintf('l=%i: %f\t',l,ratio);
end
fprintf('\n');
end

