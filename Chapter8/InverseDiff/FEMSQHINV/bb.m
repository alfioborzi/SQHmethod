function [bxy] = bb(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% purpose:
%
%    the zero-order coefficient function
%
% author:   michael holst 101095
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   [L,one] = size(x);
   
   bxy = ones(L,1);
   bxy = x.*y;
end

