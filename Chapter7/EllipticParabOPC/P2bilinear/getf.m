%
%Sets the value of f to a consant value

function [ f ] = getf( OCP )

N=OCP.N;
for i=1:N
    for j=1:N
        f(i,j)=1;
    end
end


end

