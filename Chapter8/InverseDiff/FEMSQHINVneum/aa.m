function [axy] = aa(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% purpose:
%
%    the first order coefficient function
%
% author:   michael holst 101095
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   [L,one] = size(x);
   axy = ones(L,1);
   
  for l=1:L
   z=0;
   if 0.4 < x(l) 
       if x(l) < 0.6 
       if 0.4 < y(l)
           if y(l) < 0.6 
           z=0.3;       
           end
       end
       end
   end
       axy(l) = 1 + z ;
  end
  

end

