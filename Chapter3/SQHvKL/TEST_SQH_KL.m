% TEST SCRIPT FOR THE CALCULATION MEASURE WITH THE SQH_KL METHOD
close all;
clear all; 

ltype = {'b-','r--','m-.'};

% ----  Define the model functions of the optimization problem  ----
% - Use the dot operators
%
% % Tracking - regular
% Kad = [-2, 2 ];
% Iab = [ 0, 1 ];
% y_b = 0.0;
% model.ya = 0.5;                                      % Initial condition for y(x)
% model.pT = @(yb) 0;                                  % Terminal condition for p(x)
% model.F  = @(x,y,v)     y+v;                         % The rhs state equation
% model.Fp = @(x,y,yd,p,v)  2*(y-yd)-p+0*v;            % The rhs adjoint equation
% model.L  = @(x,y,yd,v)  (y-yd).^2 + 0.01*v.^2;       % The cost function
% model.gb = @(yb) 0;                                  % Terminal cost
% model.yD = @(x)         0.1*sin(2*pi*x) + 0.5;       % The desired function
%
% % Example - oscillations 1 & 2
% Kad = [-2, 2 ];
% % num exp Kad = [-1, 1 ];
% Iab = [ 0, 1 ];
% model.ya = 0.5; 
% % num exp model.ya = 0.0;                                      % Initial condition for y(x)
% model.pT = @(yb) 0;                                  % Terminal condition for p(x)
% model.F  = @(x,y,v)     v+0*x;                       % The rhs state equation
% model.Fp = @(x,y,yd,p,v)  2*(y-yd)+0*v;              % The rhs adjoint equation
% model.L  = @(x,y,yd,v)  (y-yd).^2 + (v.^2-1).^2;     % The cost function
% model.gb = @(yb) 0;                                  % Terminal cost
% model.yD = @(x)         0.0*sin(10*x);               % The desired function

% Example - concentrations
% Kad = [ 0, 20 ];
% Iab = [ 0, 1 ];
% y_b = 1.0;
% model.ya = 0.0;                                         % Initial condition for y(x)
% model.pT = @(yb) -2*(yb-y_b);                           % Terminal condition for p(x)
% model.F  = @(x,y,v)     0*x+v;                          % The rhs state equation
% model.Fp = @(x,y,yd,p,v)  0*x+0*v;                      % The rhs adjoint equation
% model.L  = @(x,y,yd,v)  v.*(x-0.5).^2;                  % The cost function
% model.gb = @(yb)  (yb-y_b)^2;                           % Terminal cost
% model.yD = @(x)         0.0*sin(10*x);                  % The desired function


% % Tracking 
Kad = [-2, 2 ];
Iab = [ 0, 1 ];
model.ya = 0.5;                                      % Initial condition for y(x)
model.pT = @(yb) 0;                                  % Terminal condition for p(x)
model.F  = @(x,y,v)     y+v;                         % The rhs state equation
model.Fp = @(x,y,yd,p,v)  2*(y-yd)-p+0*v;            % The rhs adjoint equation
model.L  = @(x,y,yd,v)  (y-yd).^2 + (v.^2-1).^2;     % The cost function
model.gb = @(yb) 0;                                  % Terminal cost
model.yD = @(x)         0.1*sin(2*pi*x) + 0.5;       % The desired function

% % Bang-bang (at x=1)
% Kad = [0, 1 ];
% Iab = [0, 2 ];
% model.ya = 0.5;                                      % Initial condition for y(x)
% model.pT = @(yb) 0;                                  % Terminal condition for p(x)
% model.F  = @(x,y,v)       y.*v;                      % The rhs state equation
% model.Fp = @(x,y,yd,p,v)  v-1-p.*v;                  % The rhs adjoint equation
% model.L  = @(x,y,yd,v)    -y.*(1-v);                 % The cost function
% model.gb = @(yb) 0;                                  % Terminal cost
% model.yD = @(x)           0.0*sin(2*pi*x);           % The desired function

% % Singular
% Kad = [-1 , 1]; 
% Iab = [ 0, 1 ];
% model.ya = 1.0;                                      % Initial condition for y(x)
% model.pT = @(yb) 0;                                  % Terminal condition for p(x)
% model.F  = @(x,y,v)       v+0*x;                     % The rhs state equation
% model.Fp = @(x,y,yd,p,v)  v+0*x;                     % The rhs adjoint equation
% model.L  = @(x,y,yd,v)    y.*v;                      % The cost function
% model.gb = @(yb) 0;                                  % Terminal cost
% model.yD = @(x)           0.0*sin(2*pi*x);           % The desired function

%
% SQH parameters into a structure
%
sqh_p.kmax   = 10000;
sqh_p.epsi   = 100.0;
sqh_p.sigma  = 1.1;
sqh_p.zeta   = 0.9;
sqh_p.eta    = 1e-3; 
sqh_p.rkappa = 1e-4;


% --- define the domain interval ---

% Iab = [0 1]; 
Nx  = 201;     % grid points on Iab incl. init and final point
x   = linspace(Iab(1),Iab(2),Nx);
dx  = x(2)-x(1);


% --- define K_ad the adimissible set value ---

% Kad = [-1 1] ;

% - mesh size and grid on Kad 
Nv = 401;                       % number of grid points on Kad incl. boundary
v  = linspace(Kad(1),Kad(2),Nv); 
dv = v(2)-v(1);

% --- prepare the initial PDF, uniform distribution for all x
mu = ones(size(v));
% mu = exp(-(v-0.5).^2/0.2);
mu = mu/( dv * sum(mu));
nu0 = repmat(mu,Nx,1);     % use mu to initialize PDF for all x  



% *****  Call the function of the SQH with KL distance ****
%
%
[nu, y, p, u, yd, err, info] = SQH_KL(nu0,v,x,model,sqh_p);



if err>0
    fprintf('\n WARINING: ACCURACY NOT SATISFIED. Error code = %i \n',err);
end

disp(info.msg);

fprintf('\n eps = %g \n',info.eps);

% plot the YM at all x 

figure(1)
mesh(v,x,nu)
view(-15,30)
xlabel('v'); ylabel('x');
print('-depsc2', 'fig_meshNU.eps'); 

figure(2)
plot(x,y,ltype{1},'Linewidth',2)
xlabel('x'); ylabel('y');
print('-depsc2', 'fig_y.eps'); 


 
  figure(3)
  plot(x,p,ltype{1},'Linewidth',2)
 xlabel('x'); ylabel('p');
print('-depsc2', 'fig_p.eps'); 


  
  figure(4)
  plot(x,u,ltype{1},'Linewidth',2)
  xlabel('x'); ylabel('u mean');
  print('-depsc2', 'fig_urel.eps'); 

figure(5)
% num exp semilogx(info.hist_J,ltype{1},'Linewidth',2); 
plot(info.hist_J,ltype{1},'Linewidth',2); 
xlabel('iter'); ylabel('J');
print('-depsc2', 'fig_J.eps'); 




