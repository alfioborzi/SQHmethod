% Function defining the augmented Hamiltonian function 
function Out = Hamiltonian(X,P,u,U,eps)

N = size(X,1);   
u = reshape(u,[N,N+1]);


RegPar = 0.0001;
H = -sum( dot( P , tanh( u(:,1:end-1) * X + u(:,end) ) ) ) + RegPar .* 0.5 .* norm(u,'fro').^2; 
Aug = eps .* norm(u - U,'fro').^2;

Out = H + Aug;
end