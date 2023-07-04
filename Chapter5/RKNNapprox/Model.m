function Out = Model(X0,U)

NN = NNStructure();
[N,L,delta] = deal(NN(1),NN(2),NN(3));

S = size(X0(1,:),2);
X = zeros(N,S,L);
X(:,:,1) = X0; 

% For-Loop over all Layers of the neural network (Forwardpropagation)
for l = 1:(L-1)
    X(:,:,l+1) = ResNet(X(:,:,l), U(:,:,l), delta);
end    

Out = sum(X(:,:,L),1);

end