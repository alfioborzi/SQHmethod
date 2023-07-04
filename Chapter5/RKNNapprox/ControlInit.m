% Function for the initialization of the weights for the neural network. 
function Out = ControlInit()

NN = NNStructure();
N = NN(1);
L = NN(2)-1;
rng(0);
W0 = randn(N,N,L);
b0 = zeros(N,1,L);
Out = horzcat(W0,b0);

end