function [gxy] = gcost(axy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% purpose:
%
%    essential boundary condition function
%
% author:   michael holst 101095
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     global alpha ulo uup

%    [L,one] = size(x);
%    gxy = zeros(L,1);

   
   gxy = 0.5 * alpha * axy.^2 ;

end

