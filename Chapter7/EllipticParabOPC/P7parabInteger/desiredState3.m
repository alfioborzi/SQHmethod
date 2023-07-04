function [ y ] = desiredState3( OCP )
N=OCP.N;
Nt=OCP.Nt;

for j=1:Nt+1
    for i=1:N+1
        y((j-1)*(N+1)+i,1)=2*sin(2*pi*(j-1)/Nt)*exp(-0.5*(i-N/2)^2/N);%sin(2*pi*i/N)*sin(2*pi*j/Nt);
    end
end

end


