function [v] = argMinH(u,y,p,epsilon,OCP)
%Calculates the minimum of the augmented Hamiltonian in a certain grid
%point with a secant method

if (OCP.beta==0)     %If there is a smooth cost functional, calculate the minimum directly   
    v=min(max((2*epsilon*u-p)/(OCP.alpha+2*epsilon),OCP.u_lo),OCP.u_up);
else                %Perform a secant method to find a minimum of the augmented Hamiltonian
    num_Int=100;    %Number of intervals K_U is intersected in
    Disc_KU=linspace(OCP.u_lo,OCP.u_up,num_Int);    %Discretise K_U into num_Int intervals
    dist=inf;
    while(dist>10^-12)                              %If the interval within the minimum is is smaller than 10^-12 stop
        [~,pos]=min((OCP.alpha/2)*Disc_KU.^2+OCP.beta*abs(Disc_KU).*(abs(Disc_KU)>OCP.s)-Disc_KU.*y*p+epsilon*((Disc_KU-u).^2)); %Find the position of a minimum on the disretised Disc_KU
        v=Disc_KU(pos);         %Set the argument that minimised the augmented Hamiltonian as the current value for the control in the corresponding grid point
        if (pos==1)                                 %If the interval is on the lower boundary of Disc_KU
            Disc_KU=linspace(Disc_KU(pos),Disc_KU(pos+1),num_Int); %Discretise the interval around the minimum 
            dist=(Disc_KU(pos+1)-Disc_KU(pos));                    %New distance of an interval 
        elseif  (pos==max(size(Disc_KU)))                          %Case if minimum is on the upper bound
            Disc_KU=linspace(Disc_KU(pos-1),Disc_KU(pos),num_Int);
            dist=(Disc_KU(pos)-Disc_KU(pos-1));
        else
            Disc_KU=linspace(Disc_KU(pos-1),Disc_KU(pos+1),num_Int);    %Case if the minimum is in the inner of Disc_KU, discretise around the minimum
            dist=Disc_KU(pos+1)-Disc_KU(pos-1);
        end
    end

end 
end

