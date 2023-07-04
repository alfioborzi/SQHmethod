% Function defining the Parameters for the SQH learning algorithm 
function  Out = ParameterSQH()
eps = 1;
stop = 10e-10;
increps = 1.1;
decreps = 0.9;
suffdecr = 10e-5;
kmax = 5000;
Out = [eps,stop,increps,decreps,suffdecr,kmax];
end