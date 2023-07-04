function [ y ] = desiredState( OCP )
N=OCP.N;
Nt=OCP.Nt;
a=OCP.a;
b=OCP.b;

y=zeros((Nt+1)*(N+1),1);
dx=(b-a)/N;
for n=2:Nt
    x_quer=(b+a)/2;
    counter=0;
    for x=a+dx:dx:b-dx
        counter=counter+1;
        if (x_quer-0.1<=x && x<=x_quer+0.1)
            y((n-1)*(N+1)+counter)=5*sin(2*pi*n/Nt);
        end
    end
end


end

