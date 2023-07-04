% Main Function for the SQH training algorithm, nessecary input is
% restricted to the Initial weights for each Layer as a 3D Tensor, as well
% as the preprocessed training Data and Training Labels.
function [Out,Jh] = SQH(U,Xtrain,Ytrain,iHB)
% Call all nessecary functions and assign the Parameters for running the
% algorithm
p = ParameterSQH();
NN = NNStructure();
S = size(Xtrain(1,:),2);

[eps,stop,increps,decreps,suffdecr,kmax] = deal(p(1),p(2),p(3),p(4),p(5),p(6));
[N,L,delta] = deal(NN(1),NN(2),NN(3));

% Extract the number of samples from the Input data of the algorithm

% Prepare two 3D tensors for saving the states and co-states
X = zeros(N,S,L);
Xtemp = zeros(N,S,L);
X(:,:,1) = Xtrain; 
Xtemp(:,:,1) = Xtrain;

P = zeros(N,S,L);

% For-Loop over all Layers of the neural network (Forwardpropagation)
for l = 1:(L-1)
    X(:,:,l+1) = ResNet(X(:,:,l), U(:,:,l), delta);
end    

% Calculate the Inital value for the backpropagation 
P(:,:,L) = -dPhi_dx(X(:,:,L),Ytrain);

% For-Loop over all layers to generate the adjoint states (Backpropagation)
for l = 1:(L-1)
    P(:,:,L-l) = Adjoint(P(:,:,(L-l)+1),X(:,:,L-l), U(:,:,L-l), delta);
end    
iteration = 1;
% While-Loop over all iterations including a stopping condition 
while iteration < kmax
    
    % Initialize temporary tensor for weights    
    Utemp = zeros(N,N+1,L-1);
    
    if iHB == 1 
    % Maximization step with for-loop over all layers using H-Max update
    for l = 1:L-1     
        Utemp(:,:,l) = H_max(X(:,:,l),P(:,:,l+1),U(:,:,l),eps);     
    end
    elseif iHB==2 
    % Maximization step with for-loop over all layers using BFGS update
    
    options = setoptimoptions(...
    'algorithm', 'fminlbfgs',...       
    'goalsexactachieve', 0,...
    'hessupdate', 'bfgs',...
    'storeN', 10,...
        'display', 'iter',...        
        'TolCon', 1e-3,...        
        'TolX' , 1e-4,...
        'TolFun' ,1e-4,...
        'Gradobj', 'off',...
             'Display', 'off',...
        'GradConstr', 'off',...
        'MaxIter', 10,...
        'MaxFunEvals', 5e2,...
        'popsize', 10);
      
     

       for l = 1:L-1         
       Utemp(:,:,l) = fminlbfgs(@(u)Hamiltonian(X(:,:,l),P(:,:,l+1),u,U(:,:,l),eps), ... 
           U(:,:,l),options); 
       end
     end
    
    % For-Loop over all Layers of the neural network (Forwardpropagation)
    for l = 1:(L-1)
        Xtemp(:,:,l+1) = ResNet(Xtemp(:,:,l), Utemp(:,:,l), delta);
    end    
    
    % Norm difference between the newly generated control and the previous
    % one 
    DeltaU = delta * norm(Utemp(:) - U(:),'fro').^2;
    
    % Check for the suffiecient decrease condition 
    % If the sufficient decrease condition is fullfilled decrease epsilon
    % and proceed with the backwardpropagation. Else increase epsilon
    % proceed from the beginning of the while loop, without changing the
    % iteration count
    if (LossL2(Xtemp(:,:,L),Ytrain,Utemp,delta) - LossL2(X(:,:,L),Ytrain,U,delta)) > -(suffdecr * DeltaU)
       
        eps = increps .* eps;
        
    else
        Jh(iteration)=LossL2(X(:,:,L),Ytrain,U,delta);
        eps = decreps .* eps;
        U = Utemp;
        X = Xtemp;   
    
        % Calculate the Inital value for the backpropagation 
        P(:,:,L) = -dPhi_dx(X(:,:,L),Ytrain);

        % For-Loop over all layers to generate the adjoint states (Backpropagation)
        for l = 1:(L-1)
            P(:,:,L-l) = Adjoint(P(:,:,(L-l)+1),X(:,:,L-l), U(:,:,L-l), delta);
        end    
        
        % Stopping condition
        if DeltaU < stop
            Out = U;
            fprintf('--------------- Iteration: %d --------------\n', iteration);
            fprintf('Value of reg. loss functional: %f \n', LossL2(X(:,:,L),Ytrain,U,delta));
            fprintf('Value of augmentation: %f \n', eps./decreps);
            fprintf('Norm of (u - w): %.10f \n', DeltaU);
            break;
        else
            fprintf('--------------- Iteration: %d --------------\n', iteration);
            fprintf('Value of reg. loss functional: %f \n', LossL2(X(:,:,L),Ytrain,U,delta));
            fprintf('Value of augmentation: %f \n', eps./decreps);
            fprintf('Norm of (u - w): %.10f \n', DeltaU);
            iteration = iteration + 1;
        end
    end
    
end
Out = U;
end



