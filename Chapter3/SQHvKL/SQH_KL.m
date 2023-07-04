function [nu, y, p, u, yd, err, info] = SQH_KL(nu,v,x,Model,sqh_params)
%
%   [nu, y, p, u, yd, err, info] = SQH_KL(nu,v,x,Model,sqh_params)
%
%   Calculates the relaxed measure nu(x,v) related to the optimization 
%   Model by using a modifed version of the Sequential Quadratic 
%   Hamiltonian with the Kullback-Leibler divergence.
%   This algorithm is based on the paper
%   M. Annunziato, A. Borzi', "On a ...."
%   (Release beta 0.83)
%    
%   OUTPUT:
%
%   nu  = calculated measure
%   y   = state function
%   p   = adjoint function
%   u   = averaged control
%   yd  = desired trajectory
%   err = output error:  0: no error 
%                        1: ended with kmax steps, i.e. k=kmax
%                        2: epsi > MAX_EPS 
%                        3: epsi < MIN_EPS
%
%  info: structure with some useful information. Field list:
%
%  info.tau         = estimated error
%  info.tot_iter    = total iterations
%  info.fail_iter   = number of failed iterations
%  info.hist_J      = history of the cost function 
%  info.eps         = last value of the penalty epsi
%  info.msg         = string with an error message
%
%  INPUT:
%
%   nu = starting pdf distribution
%   v = grid points of the space value control
%   x = grid points of the independent variable
%   Model: structure with the functions of the optimization problem.
%   The assignements are as an example. Do not change the arguments of 
%   the functions and the name of the variables. Use the point operator '.'
%   for the algebraic operations.
%
%       Model.ya = 0.5;                                      % Initial condition for y(x)
% 	    Model.pT = @(yb) -2*(yb-1);                          % Terminal condition for p(x)
%       Model.F  = @(x,y,v)     y+v;                         % The rhs state equation
%       Model.Fp = @(x,y,yd,p,v)  2*(y-yd)-p;                % The rhs adjoint equation
%       Model.L  = @(x,y,yd,v)  (y-yd).^2 + (v.^2-1).^2;     % The cost function
%       Model.gb = @(yb)  (yb-1)^2;                          % Terminal cost (yb is scalar)
%       Model.yD = @(x)         0.1*sin(10*x);               % The desired function
%
%       Note: to define a zero desired function,  a vector of 0 function
%       must be defined like this: Model.yD=@(x) 0*x.
%       The variables x,y,yd,p lay on the numerical grid of x, i.e. Iab.
%       The variable v lay on the grid space of v, i.e. Kad.
%       In the defintion of the functions both grid space should be
%       defined. If there is no dependence on one of these, that should
%       appears like a vector of zeros. 
%       Example 1: if Fp(y) = y, then the definition is
%       Model.Fp = @(x,y,yd,p,v)  y+0*v;
%       so that both grid spaces are defined
%       Example 2: if F(v) = v, the the definition is
%       Model.F = @(x,y,v) v+0*x;
%       For the terminal cost and the terminal condition, simply set 0 in 
%       case of vanishing conditions: Model.pT = @(yb) 0; Model.gb = @(yb) 0;
%
%   sqh_params: structure with parameters for the algorithm. It contains:
%    
%       kmax    = maximum of iterations
%       epsi    = initial penalisation parameter        
%       sigma   = amplification factor of the penalisation
%       zeta    = contraction factor of the penalisation
%       eta     = tolerance of the discrimination test
%       rkappa  = tolerance of the stop criteria
%         

% defines some constants 
MAX_EPS  = 1E+10;    % max value for the penalty 
MIN_EPS  = 1E-10;    % min value for the penalty
last_tau = Inf;

% error variable
% this error is cleared in case of regular exit
err = 1;    % output with k=kmax, i.e. has reached the maximun number of loops 

% defines the Hamiltonian function
H = @(x,y,yd,p,v) p.*Model.F(x,y,v) - Model.L(x,y,yd,v) ;

% local parameters name
kmax    = sqh_params.kmax;
epsi    = sqh_params.epsi;
sigma   = sqh_params.sigma;
zeta    = sqh_params.zeta;
eta     = sqh_params.eta;
rkappa  = sqh_params.rkappa;


% ---- check the format of the variables ----
% set v as row vector ---
%
if iscolumn(v), v=v'; end
 
% set x as column vector 
%
if isrow(x),  x=x'; end
%---------

% --- Some numerical setting ---

ya = Model.ya;  % initial condition for the state equation
pT = Model.pT;           % terminal condition for this model

yd = Model.yD(x);    % Evaluates the desired function 

dx=x(2)-x(1);    % step size of the independent variable grid
dv=v(2)-v(1);    % step size of the control value grid





%  --- SQH loop ---


k_fail = 0;         % counter of the failed steps
kount = 0;         % counter of the successful steps
   
JJ = zeros(kmax,1); % history of cost values
JJ(:) = -Inf;

% for the loop number 0
[y, u] = state_function(Model.F,nu,ya,x,v); % calculates the state function
Jold   = functional_J(Model.L,Model.gb,nu,y,yd,x,v); % first value of J
p      = adjoint_function(Model.Fp,nu,y,yd,pT,v,x); % calculates the adjoint function
    
% -- the loop ----     

for k = 1:kmax

% --- SQH Step 2. ---  attempt to update nu 

% new pdf from the gradient of H_eps with the Kullback-Leibler distance
% x,y,yd,p are column vectors, v is row vector.
nu_temp = nu.*exp(H(x,y,yd,p,v)/epsi);

% normalize along the variable v
nu_temp = nu_temp./(sum(nu_temp,2)*dv);

% tolerance 
tau = dx*dv*sum(sum(abs(nu-nu_temp)));
tau = tau^2;

% calculate the state function for the new attempt of the PDF measure
[y_temp,u_temp] = state_function(Model.F,nu_temp,ya,x,v);

% cost function
J = functional_J(Model.L,Model.gb,nu_temp,y_temp,yd,x,v);


% --  SQH adaptivity --
chk = (J - Jold) + eta*tau ;

 %  -- check for the failed step --
 if  chk > 0    
       
     if epsi > MAX_EPS  %  Irregular exit on failed step
         err=2;
         break
     end
     epsi = sigma * epsi;  % update eps with a bigger value
       
     k_fail = k_fail+1;    % increase the counter of failed steps
%      fprintf(' FAILED \n');
     continue              % if failed, jump to the beginning of the loop !
 end     
 
 %  -------- %
 
 kount = kount + 1; 
 
 fprintf(' \n k=%3i ko=%3i   tau=%d   J=%d    epsi=%6.5g  ',k,kount,tau,J,epsi);
 
% [k kount]
 
 % **** if the step is accepted updates the new values !!  ******
 y  = y_temp;
 u  = u_temp;
 nu = nu_temp;
 
 Jold  = J;      % save J for the next loop
 JJ(kount) = J;

 % when accepted updates p. After a failed loop, no need to update p
 p = adjoint_function(Model.Fp,nu,y,yd,pT,v,x);

 
 last_tau = tau;   % save the last value of accepted tolerance
 
 % --  check for the exit --
 if tau < rkappa        % tolerance reached. Regular exit
     err = 0;           % clear the error flag
     break;
 end
 
 
 if epsi < MIN_EPS       % irregular exit on the accepted step
     err=3;
     break;
 end
 
 % when the step is accepted, shrink the constant penalty eps
 epsi = zeta * epsi;    

end   
% ----- END OF THE SQH LOOP ----


% erases the failed steps from the history of J, i.e. J=-Inf
JJ = JJ(JJ > -Inf);

% save output informations into a structure
info.tau        = last_tau;
info.tot_iter   = k;
info.fail_iter  = k_fail;
info.hist_J     = JJ;
info.eps        = epsi;

% set error messages
switch err
    
    case 0
        info.msg = [newline ...
            'The execution is terminated without errors.'];
    case 1
        info.msg = [newline 'The execution is terminated with the maximum '...
        'number of loops,' newline ...
        ' but the required precision has not been reached.' newline ...
        ' Try increasing kmax'];
    case 2
        info.msg = [newline 'The execution is terminated with eps > ' num2str(MAX_EPS,1) ...
            ', ' newline 'but the required precision has not been reached.'];
    case 3
        info.msg = [newline 'The execution is terminated with eps z ' num2str(MIN_EPS,1) ...
            ', ' newline 'but the required precision has not been reached.'];

end

end
