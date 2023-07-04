function ff = target(x,t)
global L T v

 ff = 0.1*cos(3*pi*x(:)/L) ;
% ff = 0.1*x(:).*(1-x(:)/L) ;

end