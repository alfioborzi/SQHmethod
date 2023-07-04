
% Riccati equation (open loop)- two players 
% dX1/d tau = A' X1 + X1 A + X1 E1 X1 + X1 E2 X2  - Q1; X1(0)=XT1
% dX2/d tau = A' X2 + X2 A + X2 E2 X2 + X2 E1 X1  - Q2; X2(0)=XT2 

close all; clear all;

u_min = [-200;-200];%[-5;-5];
u_max = [200;200]; %[5;5];

Id = eye(2,2);
 
A = [1, 0; 0,2]; %1; 

B1 = [1,0;0,-1] ; %1; 
Q1 = 0.1* Id; % alpha
R1 = 0.01*Id ; % nu
D1 = 0.01*Id; % gamma
 
B2 = [2,-1;0,2]; %2; 
Q2 = 1.0* Id; % alpha
R2 = 0.01*Id; % nu
D2 = 0.01*Id; % gamma


% Assemble XT for Nash
XT = [D1(:); D2(:)];  % 2n^2 x 1 - vector


% Final time, t0=0 initial time
T  = 0.250; t0=0; 
% Let S0 be the initial state 
S0 = [2;1];  %y0


% Time grid from t0 to T
dt = 1e-04;
Time = t0:dt:T;
iTime=size(Time,2);

%----------- determine the Nash point --------

% Nash: For efficiency precompute E = B R^-1 B'
R1i = inv(R1);
E1 = B1*R1i*(B1');

R2i = inv(R2);
E2 = B2*R2i*(B2');

% ode5 for Nash
[X] = ode5(@(t,X)mRiccatiNash(t, X, A, E1, E2, Q1, Q2), Time, -XT);



%---------------------------------------
% Nash
% assemble A + E * X 

nX  =size(X,2);
nXm =nX/2; 


% ACN=zeros(iTime,1); 
for i=1:1:iTime
X1 = reshape(X(iTime-i+1,1:nXm),size(A));
X2 = reshape(X(iTime-i+1,nXm+1:nX),size(A));
AMN = A + E1*X1 + E2*X2;
ACN(i,:) = AMN(:);
end

% State with initial condition y0 at t0 and Nash controls
y0=S0;
[y] = ode5(@(t,y)statedyn(t, y, ACN, A, t0, dt), Time, y0);
y=y';



% Nash controls corresponding to y 
u1 = zeros(2,iTime);
u2 = zeros(2,iTime);
for i=1:1:iTime
    
XX1 = reshape(X(iTime-i+1,1:nXm),size(A));
XX2 = reshape(X(iTime-i+1,nXm+1:nX),size(A));
YS  = y(:,i);

% u1(:,i) = R1i*B1'*XX1*YS;
% u2(:,i) = R2i*B2'*XX2*YS;   %unconstrained


%with control constraints:
u1(:,i) = max(min(R1i*B1'*XX1*YS, u_max), u_min);
u2(:,i) = max(min(R2i*B2'*XX2*YS, u_max), u_min);

end


% Nash: compute the cost functionals
% for initial state S0 at t0 
J1N = 0; 
for i=1:1:iTime-1
   
J1N = J1N + 0.5*(y(:,i)'*Q1*y(:,i) + u1(:,i)'*R1*u1(:,i));
end
J1N = dt * J1N + 0.5*y(:,iTime)'*D1*y(:,iTime);
%  J1 = 0.5* 0*sum((y(:,end)).^2) + 0.5*dt*0.1*sum(sum(u1.^2)) + 0.5*1*dt*sum(sum(y.^2))


J2N = 0; 
for i=1:1:iTime-1
J2N = J2N + 0.5*(y(:,i)'*Q2*y(:,i) + u2(:,i)'*R2*u2(:,i));
end
J2N = dt * J2N + 0.5*y(:,iTime)'*D2*y(:,iTime);
%  J2 = 0.5* 0*sum((y(:,end)).^2) + 0.5*dt*1*sum(sum(u2.^2)) + 0.5*1*dt*sum(sum(y.^2))

% (J1N,J2N) is the Nash point
NE = [J1N,J2N];
%------------------------------------ 

save('uRiccati.mat', 'u1','u2');
save('yRiccati.mat', 'y');







%--------- Set up ODEs -----------
function dydt = statedyn(t, y, AC, AF, t0, dt)
% State equation
% dy/dt = A(t) * y ; 
% time step
jf = ((t-t0)/dt);
j = uint32(jf)+1;
% A(t)
AT = AC(j,:);
AD = reshape(AT, size(AF)); %Convert from "n^2"-by-1 to "n"-by-"n"
dydt = AD*y;                %Determine derivative

end

%--------- 


function dXdt = mRiccatiNash(t, X, A, E1, E2, Q1, Q2)
% Our transformed coupled Riccati equations
% dX1/d tau = A' X1 + X1 A + X1 E1 X1 + X1 E2 X2 - Q1 ; 
% dX2/d tau = A' X2 + X2 A + X2 E2 X2 + X2 E1 X1 - Q2 ; 

X1 = reshape(X(1:4,:), size(A));
X2 = reshape(X(5:8,:), size(A));

dX1dt = A'*X1 + X1*A + X1*E1*X1 + X1*E2*X2 - Q1; %Determine derivative
dX2dt = A'*X2 + X2*A + X2*E2*X2 + X2*E1*X1 - Q2; %Determine derivative

dXdt = [dX1dt(:); dX2dt(:)] ; %Convert from 2 "n"-by-"n" to "2 n^2"-by-1
end


