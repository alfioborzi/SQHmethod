% Function for generating training and test Data, for a 
% univariate approximation problem.
function [Out1, Out2, Out3, Out4, Out5, Out6] = Generate_Data()
%        [Xtrain,xtrain,ytrain,Xtest,xtest,ytest] = Generate_Data();
% Predefine the size of the uniform grid for the approximation task
a = 0;
b = 1;
step = 0.0002;
StepTest = 0.00123;



% Generate datapoints for training
x = (a:step:b);
% H = heaviside(x) evaluates the Heaviside step function at x. 
% The Heaviside function is a discontinuous function that 
% returns 0 for x < 0, 1/2 for x = 0, and 1 for x > 0.
y = 3*heaviside(x - 0.313) + heaviside(x - 0.747) + 2*cos(4*pi*x) ;

xtest = (a:StepTest:b);
ytest = 3*heaviside(xtest - 0.313) + heaviside(xtest - 0.747) + 2*cos(4*pi*xtest) ;

% fit the data to the shape of the neural network
NN = NNStructure();
N = NN(1);

Out1 = repmat(x,N,1);
Out2 = x;
Out3 = y;
Out4 = repmat(xtest,N,1);   
Out5 = xtest;
Out6 = ytest;
end