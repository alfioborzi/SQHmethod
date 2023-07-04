function ff = initpT(x,t)
global L T v
%ff = 0.*exp(-(x(:)-x(mpx)).^2/omega);

ff = 0*(pi*v/L)*sin(pi*x/L).*sin(pi*v*t/L);
end