function ff = exact(x,t)
global L T v

ff = sin(pi*x/L).*cos(pi*v*t/L);
end