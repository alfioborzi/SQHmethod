function ff = initp(x,t)
global L T v

ff = 0*sin(pi*x/L).*cos(pi*v*t/L);
end