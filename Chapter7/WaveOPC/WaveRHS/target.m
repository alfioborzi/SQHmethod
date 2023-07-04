function ff = target(x,t)
global L T v
% 
ff = (1-t)*sin(pi*x(:)/L) + t*(sin(2*pi*x(:)/L));
end