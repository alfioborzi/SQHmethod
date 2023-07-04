% Function defining the augmented Hamiltonian function 
function Out = Hamiltonian(X,P,u,U,eps)

N = size(X,1);   
u = reshape(u,[N,N+1]);


RegPar1 = 1.0e-4;
RegPar2 = 0.0;
H = -sum( dot( P , tanh( u(:,1:end-1) * X + u(:,end) ) ) ) ... 
    + RegPar1 .* 0.5 .* norm(u,'fro').^2 ... 
    + RegPar2 .* norm(u,1);

Aug = eps .* norm(u - U,'fro').^2;

Out = H + Aug;
end