function [ v ] = argmaxH(y,p,A,B1,B2,gv,ct,OCP)

    nInt=200;
    
    % admissible control values [umin,umax]
    umin =OCP.umin;
    umax =OCP.umax;
    du   = (umax-umin)/nInt ; 
    ux = [umin:du:umax] ; % defines the discrete set in each component
    Nu = length(ux) ;
    
         HPMax = -1.0e+10; 
        % search for max in the discrete set 
        for ii=1:Nu
        for jj=1:Nu
                
        uus = [ux(ii);ux(jj)]' ; 
        HP=HPfunction(A,B1,B2,gv,y,p,uus,ct,OCP);
      
        if ( HP > HPMax ) 
           HPMax = HP;
           iim=ii; jjm=jj; 
        end
      
        end
        end

        v(1)=ux(iim);
        v(2)=ux(jjm);
end


