function [hxy] = hh(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% purpose:
%
%    natural boundary condition function
%
% author:   michael holst 101095
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   [L,one] = size(x);
   hxy = zeros(L,1);

   u_x = pi * cos(pi*x) .* sin(pi*y);

   axy = aa(x,y);
   hxy = axy .* u_x ;

end

