%
%Function to set values for the desired state yd

function [ y ] = desiredState(OCP)

N=OCP.N;    %Number of intervals in space
Nt=OCP.Nt;  %Number of intervals in time
a=OCP.a;
b=OCP.b;

y=0*ones((Nt+1)*(N+1),1);
dx=(b-a)/(N+1);
for n=2:Nt+1
    x_quer=0.5+0.2*sin(2*pi*(n-1)/Nt);
    
    
    counter=0;
    for x=a+dx:dx:b-dx
        counter=counter+1;
        if (x_quer-0.1*(b-a)<=x && x<=x_quer+0.1*(b-a))
            y((n-1)*(N+1)+counter)=0.5;
        end
    end
end
end

