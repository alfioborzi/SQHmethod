function [v] =get_xd(T,dt)

Nt=T/dt;




for i=1:Nt
    x_quer(1,i)= 2*((i-1)/Nt) - 1 ;
    y_quer(1,i)=sin(pi*(i-1)/Nt);
    
end

v=[x_quer;y_quer];

end
