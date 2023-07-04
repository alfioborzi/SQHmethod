% Function for generating training and test Data
% Iris data-set
function [Out1, Out2, Out3, Out4, Out5, Out6] = Generate_Data()
%        [Xtrain,xtrain,ytrain,Xtest,xtest,ytest] = Generate_Data();

% Generate datapoints for training

load('datatrain2.mat'); % training data stored in arrays X, y - 120

% normalizing the data ?
% X = (X - min(X)) ./ ( max(X) - min(X) );


Xtrain = X'; 
ytrain = y'; 
nstep  = size(y,1)-1;
step = 1./nstep; 
xtrain = (0:step:1);

% test part

load('datatest2.mat'); % test X y - 30

% normalizing the data ?
% X = (X - min(X)) ./ ( max(X) - min(X) );

Xtest = X';
ytest = y'; 
nstep  = size(y,1)-1;
step = 1./nstep; 
xtest = (0:step:1);

% fit the data to the shape of the neural network
NN = NNStructure();
N = NN(1);


% size(X)
Out1 = repmat(Xtrain,2,1);  % Number of nodes is 2 times number of known attributes
%Out1 = Xtrain;
Out2 = xtrain;
Out3 = ytrain;

Out4 = repmat(Xtest,2,1);   
%Out4 = Xtest;
Out5 = xtest;
Out6 = ytest;
end

