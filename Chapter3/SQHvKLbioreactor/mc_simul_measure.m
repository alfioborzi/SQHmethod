function  y = mc_simul_measure(F,ya,nu,t,v)
%
%    y = mc_simul_measure(F,ya,nu,t,v)
%
%

% sz = size(nu);

nT=length(t); 
%  dt=t(2)-t(1);
% nV=length(v);   
dv=v(2)-v(1);

nu = cumsum(nu,2)*dv;

w = rand(size(t));

% 
% y=zeros(size(t));
% y(1)=ya;


val=zeros(size(t));


for n=1:nT
    
    [~,idv]=max((nu(n,:)-w(n))>0);
    
    val(n)=v(idv); 
           
end

FY = @(xx,y) interp1(t,F(t,y,val),xx);

% calculates the state function
y = ode1(FY,t,ya);

end