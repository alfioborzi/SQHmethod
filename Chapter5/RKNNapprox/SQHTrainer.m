%
disp('A. Borzì');
disp('The Sequential Quadratic Hamiltonian Method');
disp('1st Edition, CRC Press, 2023.');
disp(' ');
%
% SQH Method for Runge-Kutta Neural Networks 
%
close all;
clear all; 

ltype = {'r-*','b-','m-x','b-*','r:','m-.'};

% First the necessary data for the desired task needs to be generated
[Xtrain,xtrain,ytrain,Xtest,xtest,ytest] = Generate_Data();

% Next the initial weights need to be determmined 
U0 = ControlInit();

% The generated training data and Control initialization can now be used to
% train ResNet with the SQH algorithm 
%
% Set iHB=1 to choose H-max or iHB=2 to choose BFGS 

tStart = tic;  
iHB = 1; 
[U,Jh] = SQH(U0,Xtrain,ytrain,iHB);
tEnd = toc(tStart); 

% The trained model can be validated on the test set using the Mean-Squared
% Error as performance measure.
% Here yhat is the output produced by the neural network after training
% with the SQH method

yhat = Model(Xtest, U);
MSE = mean((yhat-ytest).^2);


fprintf('\n CPU time %d and MSE %d \n',tEnd,MSE);


figure(1)
plot(xtest,yhat,ltype{1},'Linewidth',2); hold on; 
plot(xtrain,ytrain,ltype{2},'Linewidth',2);
legend({'$NN$','$Train \, data set$'},'Interpreter','Latex','FontSize',12, ...
    'Location','northwest')

print('-depsc2', 'approxRKNN01.eps','-b0'); 
print('-dpdf', 'approxRKNN01.pdf','-b0');

figure(2)
plot(Jh,ltype{2},'Linewidth',2);

print('-depsc2', 'approxRKNN02.eps','-b0'); 
print('-dpdf', 'approxRKNN02.pdf','-b0');

