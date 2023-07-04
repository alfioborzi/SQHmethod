% Function defining the Lossfunction with L2 refulaization.
% Regularizationparameter (RegPar) must coincide with the one defined in the
% function Hamiltonian.
function Out = LossL2(XL,Y,U,delta)

RegPar = 0.0001;
h = sum(XL,1);

Phi = mean((h - Y).^2);
R = RegPar .* 0.5 .* norm(U(:),'fro').^2;

Out = Phi + delta .* R;
end