function [a,b,m,h,x,y,Nt,t,dt,alpha,beta,nu] = parameters(kt,OCPA)

% Domain definition O = (a, b) x (a, b)
a = OCPA.a; 
b = OCPA.b;

dt = OCPA.T/OCPA.Nt;
dx = ((OCPA.b-OCPA.a)/OCPA.Nx);
h = dx; 

T  = OCPA.T;
Nt = OCPA.Nt;
N  = OCPA.Nx;

% Number of subintervals
m = N;

% Generating mesh
x = linspace(a,b,m+1);
y = linspace(a,b,m+1);



% Time grid of a time window
t = linspace(0,T,Nt+1);




% Values of weights in J 
alpha = OCPA.alpha ;
beta = OCPA.beta ;
nu = 0;


