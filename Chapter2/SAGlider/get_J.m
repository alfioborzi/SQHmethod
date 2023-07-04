function [J] = get_J(y,u,OCP,NE)

global A lambda

J= -y(1,NE) + lambda*y(2,NE) + 0.5*A*y(2,NE).^2 ;  
end

