%Calculating the solution to the adjoint equation with a Chang-Cooper space
%discretization and implicit Euler time discretization

function f1 = adjfpsolve(u_1,u_2,Vp,xyd,OCPA)

dt = OCPA.T/OCPA.Nt;
dx = ((OCPA.b-OCPA.a)/OCPA.Nx);
h = dx; 

T  = OCPA.T;
Nt = OCPA.Nt;
N  = OCPA.Nx;

x = linspace(OCPA.a,OCPA.b,N+1);
y = linspace(OCPA.a,OCPA.b,N+1);
t = linspace(0,T,Nt+1);



%Right hand-side of the adjoint equation

F=zeros(N+1,N+1,Nt+1);

for i=1:Nt+1
    [xt,yt] = xyd(t(i));
   
   for k=1:N+1
    for j=1:N+1
            dkji=0.5*OCPA.alpha*(u_1(k,j,i).^2+u_2(k,j,i).^2);
    F(k,j,i)= Vp(x(k),y(j),xt,yt) + dkji;
    end
   end
end





% Terminal condition for the adjoint
[xT,yT] = xyd(T);
for i=1:N+1
    for j=1:N+1
                f1(i,j,Nt+1) = Vp(x(i),y(j),xT,yT);
    end
end


    
    m=N;


% Using Euler for the first time step
k = Nt;

% Interpolating u along x-axis
for j = 1:m+1
    for i = 1:m
        u_11(i,j) = inter(u_1(:,:,k),i,j,1);
        u_21(i,j) = inter(u_2(:,:,k),i,j,1);
    end
end

% Interpolating u along y-axis
for i = 1:m+1
    for j = 1:m
        u_12(i,j) = inter(u_1(:,:,k),i,j,2);
        u_22(i,j) = inter(u_2(:,:,k),i,j,2);
    end
end 

% Mapping the 2D solution matrix to 1D solution vector.
fm1 = mat_to_vec(f1(2:m,2:m,Nt+1));

% Number of grid points along each side is m+1.
% Number of boundary points along each side is 2.
% Hence total number of degrees of freedom along each side is m-1.
d = (m-1)^2;

A = sparse(d,d);

t1  = t(k);


for i = 2:m
    for j = 2:m
        K_1 = (1-delta_1(t1,x(i)+h/2,y(j),h,u_11(i,j),u_21(i,j),OCPA))*B_1(t1,x(i)+h/2,y(j),h,u_11(i,j),u_21(i,j)) + (1/h)*C_1(t1,x(i)+h/2,y(j),h,OCPA);
        K_2 = (1/h)*C_1(t1,x(i)+h/2,y(j),h,OCPA) - delta_1(t1,x(i)+h/2,y(j),h,u_11(i,j),u_21(i,j),OCPA)*B_1(t1,x(i)+h/2,y(j),h,u_11(i,j),u_21(i,j));
        K_3 = (1-delta_1(t1,x(i)-h/2,y(j),h,u_11(i-1,j),u_21(i-1,j),OCPA))*B_1(t1,x(i)-h/2,y(j),h,u_11(i-1,j),u_21(i-1,j)) + (1/h)*C_1(t1,x(i)-h/2,y(j),h,OCPA);
        K_4 = (1/h)*C_1(t1,x(i)-h/2,y(j),h,OCPA) - delta_1(t1,x(i)-h/2,y(j),h,u_11(i-1,j),u_21(i-1,j),OCPA)*B_1(t1,x(i)-h/2,y(j),h,u_11(i-1,j),u_21(i-1,j));
        
        L_1 = (1-delta_2(t1,x(i),y(j)+h/2,h,u_12(i,j),u_22(i,j),OCPA))*B_2(t1,x(i),y(j)+h/2,h,u_12(i,j),u_22(i,j)) + (1/h)*C_2(t1,x(i),y(j)+h/2,h,OCPA);
        L_2 = (1/h)*C_2(t1,x(i),y(j)+h/2,h,OCPA) - delta_2(t1,x(i),y(j)+h/2,h,u_12(i,j),u_22(i,j),OCPA)*B_2(t1,x(i),y(j)+h/2,h,u_12(i,j),u_22(i,j));
        L_3 = (1-delta_2(t1,x(i),y(j)-h/2,h,u_12(i,j-1),u_22(i,j-1),OCPA))*B_2(t1,x(i),y(j)-h/2,h,u_12(i,j-1),u_22(i,j-1)) + (1/h)*C_2(t1,x(i),y(j)-h/2,h,OCPA);
        L_4 = (1/h)*C_2(t1,x(i),y(j)-h/2,h,OCPA) - delta_2(t1,x(i),y(j)-h/2,h,u_12(i,j-1),u_22(i,j-1),OCPA)*B_2(t1,x(i),y(j)-h/2,h,u_12(i,j-1),u_22(i,j-1));
        
        l = index_map(i,j,m);
        
        % Building the matrix
        A(l,l) = 1 + (dt/h)*(K_3+K_2+L_2+L_3);
        
        if(i == 2)
            A(l,l+m-1) = -(dt/h)*K_2;
            if(j==2)
                A(l,l+1) = -(dt/h)*L_2;
            elseif(j == m)
                A(l,l-1) = -(dt/h)*L_3;
            else
                A(l,l-1) = -(dt/h)*L_3;
                A(l,l+1) = -(dt/h)*L_2;
            end
        elseif(i==m)
            A(l,l-(m-1)) = -(dt/h)*K_3;
            if(j==2)
                A(l,l+1) = -(dt/h)*L_2;
            elseif(j == m)
                A(l,l-1) = -(dt/h)*L_3;
            else
                A(l,l-1) = -(dt/h)*L_3;
                A(l,l+1) = -(dt/h)*L_2;
            end
        else
            A(l,l+m-1) = -(dt/h)*K_2;
            A(l,l-(m-1)) = -(dt/h)*K_3;
            A(l,l-1) = -(dt/h)*L_3;
            A(l,l+1) = -(dt/h)*L_2;
        end
        
        % Building the right hand side
        
        b1(l) = fm1(l) + dt * F(i,j,k); 
    end
    
end

fmn = A\b1';

f1(:,:,k) = vec_to_mat(fmn);

%Boundary condition.
for i=1:m+1
    f1(i,1,k) = 0;
    f1(i,m+1,k) = 0;
    f1(1,i,k) = 0;
    f1(m+1,i,k) = 0;
end



% Using BDF2 for the time discretization.
for k = Nt-1:-1:1
    
    t1 = t(k);
    
    
    % Interpolating u along x-axis
    for j = 1:m+1
        for i = 1:m
            u_11(i,j) = inter(u_1(:,:,k),i,j,1);
            u_21(i,j) = inter(u_2(:,:,k),i,j,1);
        end
    end
    
    % Interpolating u along y-axis
    for i = 1:m+1
        for j = 1:m
            u_12(i,j) = inter(u_1(:,:,k),i,j,2);
            u_22(i,j) = inter(u_2(:,:,k),i,j,2);
        end
    end
    
    fm1 = mat_to_vec(f1(2:m,2:m,k+1));
    fm2 = mat_to_vec(f1(2:m,2:m,k+2));
    d = (m-1)^2;
    
    A = sparse(d,d);
    %b = zeros(d,1);
    for i = 2:m
        for j = 2:m
             K_1 = (1-delta_1(t1,x(i)+h/2,y(j),h,u_11(i,j),u_21(i,j),OCPA))*B_1(t1,x(i)+h/2,y(j),h,u_11(i,j),u_21(i,j)) + (1/h)*C_1(t1,x(i)+h/2,y(j),h,OCPA);
             K_2 = (1/h)*C_1(t1,x(i)+h/2,y(j),h,OCPA) - delta_1(t1,x(i)+h/2,y(j),h,u_11(i,j),u_21(i,j),OCPA)*B_1(t1,x(i)+h/2,y(j),h,u_11(i,j),u_21(i,j));
             K_3 = (1-delta_1(t1,x(i)-h/2,y(j),h,u_11(i-1,j),u_21(i-1,j),OCPA))*B_1(t1,x(i)-h/2,y(j),h,u_11(i-1,j),u_21(i-1,j)) + (1/h)*C_1(t1,x(i)-h/2,y(j),h,OCPA);
             K_4 = (1/h)*C_1(t1,x(i)-h/2,y(j),h,OCPA) - delta_1(t1,x(i)-h/2,y(j),h,u_11(i-1,j),u_21(i-1,j),OCPA)*B_1(t1,x(i)-h/2,y(j),h,u_11(i-1,j),u_21(i-1,j));
             
             L_1 = (1-delta_2(t1,x(i),y(j)+h/2,h,u_12(i,j),u_22(i,j),OCPA))*B_2(t1,x(i),y(j)+h/2,h,u_12(i,j),u_22(i,j)) + (1/h)*C_2(t1,x(i),y(j)+h/2,h,OCPA);
             L_2 = (1/h)*C_2(t1,x(i),y(j)+h/2,h,OCPA) - delta_2(t1,x(i),y(j)+h/2,h,u_12(i,j),u_22(i,j),OCPA)*B_2(t1,x(i),y(j)+h/2,h,u_12(i,j),u_22(i,j));
             L_3 = (1-delta_2(t1,x(i),y(j)-h/2,h,u_12(i,j-1),u_22(i,j-1),OCPA))*B_2(t1,x(i),y(j)-h/2,h,u_12(i,j-1),u_22(i,j-1)) + (1/h)*C_2(t1,x(i),y(j)-h/2,h,OCPA);
             L_4 = (1/h)*C_2(t1,x(i),y(j)-h/2,h,OCPA) - delta_2(t1,x(i),y(j)-h/2,h,u_12(i,j-1),u_22(i,j-1),OCPA)*B_2(t1,x(i),y(j)-h/2,h,u_12(i,j-1),u_22(i,j-1));
             
            l = index_map(i,j,m);
            
            % Building the matrix 
            A(l,l) = 3 + 2*(dt/h)*(K_3+K_2+L_2+L_3);
            
            if(i == 2)
                A(l,l+m-1) = -2*(dt/h)*K_2;
                if(j==2)
                    A(l,l+1) = -2*(dt/h)*L_2;
                elseif(j == m)
                    A(l,l-1) = -2*(dt/h)*L_3;
                else
                    A(l,l-1) = -2*(dt/h)*L_3;
                    A(l,l+1) = -2*(dt/h)*L_2;
                end   
            elseif(i==m)
                A(l,l-(m-1)) = -2*(dt/h)*K_3;
                if(j==2)
                    A(l,l+1) = -2*(dt/h)*L_2;
                elseif(j == m)
                    A(l,l-1) = -2*(dt/h)*L_3;
                else
                    A(l,l-1) = -2*(dt/h)*L_3;
                    A(l,l+1) = -2*(dt/h)*L_2;
                end
            else
                 A(l,l+m-1) = -2*(dt/h)*K_2;
                 A(l,l-(m-1)) = -2*(dt/h)*K_3;
                 A(l,l-1) = -2*(dt/h)*L_3;
                 A(l,l+1) = -2*(dt/h)*L_2;                    
            end
            
            % Building the right hand side
            b1(l) = 4*fm1(l) - fm2(l) + 2*dt*F(i,j,k); 
        end
        
    end
    size(b1);    
    fmn = A\b1';
    f1(:,:,k) = vec_to_mat(fmn);
     %Boundary condition.
    for i=1:m+1
        f1(i,1,k) = 0;
        f1(i,m+1,k) = 0;
        f1(1,i,k) = 0;
        f1(m+1,i,k) = 0;
    end
    
end
