% Function defining the Adjoint equation, resempbling the backpropagation
% process through one layer of a ResNet architecture
function Out = Adjoint(Pl,Xl,Ul,delta)

W = Ul(:,1:end-1);
df_dy = 1-tanh(W * Xl + Ul(:,end)).^2;
Out = Pl + delta .* (W * (df_dy .* Pl)); 

end
