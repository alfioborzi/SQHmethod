function ff = inityt(x)
global L T v
%ff = 0.*exp(-(x(:)-x(mpx)).^2/omega);

ff =  0*sin(pi*x(:)/L);
end