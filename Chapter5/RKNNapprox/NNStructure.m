% Function defining the Structure of the neural network including the
% number of noded, Layers and stepsize delta.
function Out = NNStructure()
Nodes = 5;
FinalTime = 5;
Layer = 20;
delta = FinalTime./Layer;
Out = [Nodes, Layer, delta];
end