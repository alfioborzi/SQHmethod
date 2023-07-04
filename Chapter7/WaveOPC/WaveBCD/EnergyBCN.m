function J = EnergyBCN(y,j,Nx)
global L T v

dx = L/Nx; 

J1 = 0;
for i=1:Nx
J1 = J1 + (y(j,i+1)-y(j,i)).^2;
end
J1 = 0.5*J1/dx;

J2 = 0;
for i=1:Nx
J2 = J2 + (y(j,i)-y(j-1,i)).^2;
end
J2 = 0.5*J2/dx;



J = J1 + J2 ;

end