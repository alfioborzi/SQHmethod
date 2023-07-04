% Function defining the Lossfunction with L2 refulaization.
% Regularizationparameter (RegPar) must coincide with the one defined in the
% function Hamiltonian.
function Out = LossL2(XL,Y,U,delta)

RegPar1 = 1.0e-4;
RegPar2 = 0.0;

h = sum(XL,1);

Phi = mean((h - Y).^2);
R1 = RegPar1 .* 0.5 .* norm(U(:),'fro').^2;
R2 = RegPar2 .* norm(U(:),1);
R = R1 + R2; 

Out = Phi + delta .* R;
end