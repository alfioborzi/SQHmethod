% Function defining the Initial value for the Backpropagation, that is the
% the derivative of the costfunctional with respect to the state variable,
% evaluated at the last layer.
function Out = dPhi_dx(XL,Y)
h = sum(XL,1);
Out = ((2.*(h - Y))/size(XL,2) ).* ones(size(XL));
end