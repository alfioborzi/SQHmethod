function ff = frhs(x,t)
global L T v
ff = 0*x(:).*(1-x(:)/L)/L;
end