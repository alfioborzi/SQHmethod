% Function defining the Hmax method for one layer
function Out = H_max(X,P,U,eps)

% Initiazizing the Hamiltonian as a function of the Weight U
H_eps = @(u) Hamiltonian(X,P,u,U,eps);
N = size(U,1);
u_pm = zeros(size(U));
RegPar = 1.0e-4; % Regularization parameter for the regularization term, 
                 % its value must fit the value in LossL2 and Hamiltonian

beta = 1./(RegPar + (2 .* eps)); 
alpha = (2 .* eps) .* beta;

O = ones(1,size(X,2));
Xbar = vertcat(X,O);

df_dy = 1-tanh(U * Xbar).^2;

    for i = 1:N
        
        u_pm(i,:) = alpha .* U(i,:) + beta .* sum(P(i,:) .* df_dy(i,:) .* Xbar, 2).';
       
    end

    % If the determined weight does not reduce the augmented Hamiltonian
    % with respect to the prior weight, the new weight is searched among
    % the convex combinations of u_pm and the old weight U
    if H_eps(u_pm) <= H_eps(U)
   
        Out = u_pm;
        
    
        
    else 
        u_temp0 = U;
        for lambda = 0.9:-0.1:0
        
            u_temp1 = lambda.^2 .* u_pm + (1 - lambda.^2) .* U; 
       
            
            if H_eps(u_temp1) <= H_eps(u_temp0)
          
               u_temp0 = u_temp1;
          
            end       
        end
        Out = u_temp0;
    end
end       
