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

ltype = {'b-','r--','m-x','b-*','r:','m-.'};

% First the necessary data for the desired task needs to be generated
[Xtrain,xtrain,ytrain,Xtest,xtest,ytest] = Generate_Data();

% Initial weights 
U0 = ControlInit();


% The generated training data and Control initialization can now be used to
% train ResNet with the SQH algorithm with H-max 
tStart = tic;  
[U,Jh] = SQH(U0,Xtrain,ytrain);
tEnd = toc(tStart); 

% The trained model can be validated on the test set using the Mean-Squared
% Error as performance measure.
% Here yhat is the output produced by the neural network after training
% with the SQH method

yhat = Model(Xtest, U);
MSE = mean((yhat-ytest).^2);

fprintf('\n CPU time %d and MSE %d \n',tEnd,MSE);


figure(1)
stairs(xtest,yhat,ltype{1},'Linewidth',2); hold on;
plot(xtest,ytest,ltype{3},'Linewidth',2); hold on; 
% plot(xtrain,ytrain,ltype{2},'Linewidth',2);
legend({'$NN$','$Train \, data \, set$'},'Interpreter','Latex','FontSize',12, ...
    'Location','northwest')


print('-depsc2', 'irisRKNN01.eps','-b0'); 
print('-dpdf', 'irisRKNN01.pdf','-b0');

% figure(2)
% plot(xtest,yint,ltype{4},'Linewidth',2); hold on;
% plot(xtest,ytest,ltype{3},'Linewidth',2); hold on; 

%  figure(2)
%  plot(Jh,ltype{2},'Linewidth',2);



