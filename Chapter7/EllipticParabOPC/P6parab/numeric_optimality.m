%
% Check PMP optimality of the SQH solution 
%

function [maxDiff,ratio] = numeric_optimality( u,y0,f,yd,OCP )

%Calculate the state variable with the calculated control u
y=forward_y(u,y0,f,OCP);  
%Calculate the adjoint variable with the control u 
p=backward(y,yd,u,OCP); 

for j=1:OCP.Nt
         for i=1:OCP.N-1 
             uzw=u((j-1)*(OCP.N-1)+i,1);
             %Calculate pointwise the value of the Hamiltonian where y,u,p returned by the SQH method are utilised
             Hopt((j-1)*(OCP.N-1)+i)=(OCP.alpha/2)*uzw.^2+OCP.beta*(abs(uzw)*(abs(uzw)>OCP.s))-uzw.*p((j-1)*(OCP.N+1)+1+i)*y((j-1)*(OCP.N+1)+1+i); 
         end
end

for j=1:OCP.Nt
    for i=1:OCP.N-1
        %Determine for y and p the value u for the
        %control where the Hamiltonian takes its minimum
        s=OCP.s;
        pint=p((j-1)*(OCP.N+1)+1+i);
        yint=y((j-1)*(OCP.N+1)+1+i);
        u1=min(max(0,pint*yint/(2*(0.5*OCP.alpha))),s);
        u2=min(max(s,-(-yint*pint+OCP.beta)/(2*(0.5*OCP.alpha))),OCP.u_up);
        H1=OCP.alpha*u1- pint*yint*u1;
        H2=OCP.alpha*u2- pint*yint*u2+OCP.beta*u2;
        [~,pos]=min([H1,H2]);
        switch pos
            case 1
                u_new((j-1)*(OCP.N-1)+i,1)=u1;
            case 2
                u_new((j-1)*(OCP.N-1)+i,1)=u2;
        end
    end
end

for j=1:OCP.Nt
         for i=1:OCP.N-1 
             uzw=u_new((j-1)*(OCP.N-1)+i,1);
             %Calculate pointwise the value of the Hamiltonian at uzw
             Hnew((j-1)*(OCP.N-1)+i)=(OCP.alpha/2)*uzw.^2+OCP.beta*(abs(uzw)*(abs(uzw)>OCP.s))...
                                    -uzw.*p((j-1)*(OCP.N+1)+1+i)*y((j-1)*(OCP.N+1)+1+i); 
         end
end

%Determine the maximum difference 
fprintf('Max DH: %d\n',max(Hopt-Hnew));
fprintf('Fraction of grid points where 0<=DH<=10^-l is fulfilled:\n')

for l=2:2:12 
count=0;
for j=1:OCP.Nt
    for i=1:OCP.N-1
        %Count the numer of grid points where the solution u is PMP optimal up to the tolerance 10^-l
        if(Hopt((j-1)*(OCP.N-1)+i)-Hnew((j-1)*(OCP.N-1)+i)<10^(-l))
            count=count+1;
        end
    end
end

ratio=count/(OCP.Nt*(OCP.N-1));
fprintf('l=%i: %f\t',l,ratio);
end
fprintf('\n');
end



