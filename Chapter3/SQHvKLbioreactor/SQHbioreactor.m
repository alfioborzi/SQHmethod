%
disp('A. Borzì');
disp('The Sequential Quadratic Hamiltonian Method');
disp('1st Edition, CRC Press, 2023.');
disp(' ');
%
% SQH Method for optimal relaxed controls of a bioreactor
%
% **** EXTENDED VERSION FOR THE VECTOR MODEL ****
close all;
clear all; 



% ----  Define the model functions of the optimization problem  ----
% - Use the dot operators
%
% Example - bioreactor
G=1; D=1; K=3; L=1; alfa=0.1; beta=0.1; gamma=0.01; 

Kad = [0, 2 ];                         
                                     
Iab = [ 0, 3 ];
model.ya = [ 0.5   0.5];                                % Initial condition for y1(x) and y2(x)
%model.pT = @(yb) [ 0  0];                               % Terminal condition for p1(x) ans p2(x)

model.pT = @(yb) [ 0 -gamma*(yb(2))];   

% y and p are a double coloumn vector
% v is a double row vector, v1 and v2 are single row vectors

model.F.n1 = @(x,y,v1)   G*y(:,1).*v1-D*y(:,1).^2  ;     % State equation n. 1                
model.F.n2 = @(x,y,v2)  -K*y(:,1).*y(:,2)+L*v2;         % State equation n. 2

model.Fp.n1 = @(x,y,yd,p,v1)  -G*p(:,1).*v1 + 2*D*y(:,1).*p(:,1) + ... 
                                K*p(:,2).*y(:,2)+(y(:,1)-yd);    % Adjoint equation n. 1
                            
model.Fp.n2 = @(x,y,yd,p,v2)    K*y(:,1).*p(:,2)+0*v2 ;          % Adjoint equation n. 2               

%model.L.n1  = @(x,y1,yd,v1)     0.5*(y1-yd).^2 + alfa*(v1.^2-1).^2;
model.L.n1  = @(x,y1,yd,v1)     0.5*(y1-yd).^2 + alfa*abs(v1.*(1-v1));

%model.L.n2  = @(x,y1,yd,v2)     0*y1 + beta*abs(v2).*(1-abs(v2));  % The cost function
model.L.n2  = @(x,y1,yd,v2)     0*y1 + beta*abs(v2.*(1-v2));  % The cost function
                
model.gb = @(yb)                0.5*gamma*yb(2)^2;                       % Terminal cost
model.yD = @(x)                 0.0*sin(10*x)+0.6;                 % The desired function


% SQH parameters into a structure
%
sqh_p.kmax   = 5000;
sqh_p.epsi   = 1.0; %10
sqh_p.sigma  = 1.1;
sqh_p.zeta   = 0.9;
sqh_p.eta    = 1e-2; 
sqh_p.rkappa = 1e-6;


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
mu = mu/( dv * sum(mu));
nu0 = repmat(mu,Nx,1,2);     % use mu to initialize PDF for all x  
                            % **** there are two densities nu0(:,:,1), nu0(:,:,1)

v = [v ; v];    % ***** raddoppia v

% *****  Call the function of the SQH with KL distance ****
%
[nu, y, p, u, yd, err, info] = SQH_KL_v(nu0,v,x,model,sqh_p);



if err>0
    fprintf('\n WARINING: ACCURACY NOT SATISFIED. Error code = %i \n',err);
end

disp(info.msg);

fprintf('\n eps = %g \n',info.eps);

%% Plot section


% --- plot results

figure(1)

subplot(1,2,1)
mesh(v(1,:),x,nu(:,:,1))
view(15,15)
xlabel('v'); ylabel('t');
title('$\nu_1$','Interpreter','latex');

subplot(1,2,2)
mesh(v(2,:),x,nu(:,:,2))
view(15,15)
xlabel('v'); ylabel('t');
title('$\nu_2$','Interpreter','latex');

print('-depsc2', 'meshNU12reactor.eps'); 

% ----- plot the state variables ----------

figure(2)

subplot(1,2,1);
 plot(x,y(:,1),'LineWidth',2)
xlabel('t'); ylabel('y1');

subplot(1,2,2);
 plot(x,y(:,2),'LineWidth',2)
xlabel('t'); ylabel('y2');

print('-depsc2', 'y12reactor.eps'); 

% ---- plot the adjont variables ----
 
  figure(3)
  
  subplot(1,2,1);
  plot(x,p(:,1),'LineWidth',2)
  xlabel('t'); ylabel('p1');

 subplot(1,2,2);
  plot(x,p(:,2),'LineWidth',2)
  xlabel('t'); ylabel('p2');

print('-depsc2', 'p12reactor.eps'); 


% ---------- plot the u mean ----------

figure(4)

subplot(1,2,1);
plot(x,u(:,1),'LineWidth',2)
xlabel('t'); ylabel('u1 mean');

subplot(1,2,2);
plot(x,u(:,2),'LineWidth',2)
xlabel('t'); ylabel('u2 mean');

print('-depsc2', 'u12meanReactor.eps'); 

% -------- plot the functional ----

figure(5)
plot(info.hist_J,'LineWidth',2); 
xlabel('SQH iterations'); ylabel('J');
print('-depsc2', 'Jreactor.eps'); 


