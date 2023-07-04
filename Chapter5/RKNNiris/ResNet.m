% Function defining the forward propagation process through one layer of
% ResNet.
function Out = ResNet(Xl,Ul,delta)
Out = Xl + delta .* tanh( Ul(:,1:end-1) * Xl + Ul(:,end) );
end