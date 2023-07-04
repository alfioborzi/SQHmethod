% Function defining the Structure of the neural network including the
% number of noded, Layers and stepsize delta.
function Out = NNStructure()
Nodes = 8;
FinalTime = 10;
Layer = 40;
delta = FinalTime./Layer;
Out = [Nodes, Layer, delta];
end