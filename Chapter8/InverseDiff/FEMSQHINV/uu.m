function [uxy] = uu(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% purpose:
%
%    true solution (for testing)
%
% author:   michael holst 101095
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   [L,one] = size(x);
   uxy = zeros(L,1);

   uxy = sin(pi*x) .* sin(pi*y);

end

