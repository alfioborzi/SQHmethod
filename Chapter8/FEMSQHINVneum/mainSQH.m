%
disp('A. Borzì');
disp('The Sequential Quadratic Hamiltonian Method');
disp('1st Edition, CRC Press, 2023.');
disp(' ');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    linear finite elements and sequential quadratic hamiltonian method 
%    for estimation of diffusion coefficient 'a'
%
%    the governing model:
%
%       - \nabla \cdot (a \nabla u) + b u = f in Omega in R^2
%
%    boundary conditions:
%                               u = g   on \Gamma_D 
%              a \nabla u \cdot n = h   on \Gamma_N
% 
%    where \partial \Omega = \Gamma = \Gamma_D \cup \Gamma_N.
%
%
%    author : Alfio Borzi
%    author of the FEM structure : Michael Holst
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all;

ltype = {'b-','r--','m-x','b-*','r:','m-.'};

global alpha ulo uup

      alpha = 0.0;
      ulo =  0.1;
      uup =  10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% set problem and generate measured data


   %%% controlling parameters
   level = 1 ;    %%% parameter for mesh definition etc.
         %%% create the mesh
      %   fprintf('defining mesh...\n');
         [INTPT,W,PHI,PHIX,PHIY, ...
          INTPTE1,WE1,PHIE1, ...
          INTPTE2,WE2,PHIE2, ...
          INTPTE3,WE3,PHIE3] = mastr(level);
         [ELNODE,NODECO,BONODE,EDGE] = mesh(level);

 %%% plot the mesh and get quad points
    
    %   fprintf('plotting mesh...\n');
         figure(1)
         %%% Get quad pts & position where control is defined
         [xm,ym] = draw(ELNODE,NODECO,BONODE,INTPT);
         
         print('-depsc2', 'meshInv01.eps','-b0'); 
       % print('-dpdf', 'meshInv01.pdf','-b0');


  %%% get sizes, create an array noting the non-boundary points
      [L_elements,three] = size(ELNODE);
      [M_nodes,two]      = size(NODECO);
      [M0_nodes,one]     = size(BONODE);
      N_nodes = M_nodes - M0_nodes;
      TMP = ones(M_nodes,1); 
      TMP(BONODE) = zeros(size(BONODE));
      NOTBONODE = find(TMP);
      
  %%% data and initialization
      [L,T] = size(ELNODE);
      
      frhs(1:L,1) = ff(xm,ym);     % rhs f 
 
      atrue(1:L,1) = aa(xm,ym);    % diff coeff to be reconstructed
    
      %%% assemble model problem 
%      fprintf('assembling linear system for state equation ...\n');
      [A,F]=assem(ELNODE,NODECO,BONODE,EDGE, ...
                     INTPT,W,PHI,PHIX,PHIY, ...
                     INTPTE1,WE1,PHIE1, ...
                     INTPTE2,WE2,PHIE2, ...
                     INTPTE3,WE3,PHIE3,atrue,frhs);

      %%% solve discrete system 
%      fprintf('solving linear system for state ...\n');
      U_F = A \ F;
      
       [nodes,two] = size(NODECO);
      
%       %%% discretization error test if y is known (uu)
%       U_t = zeros(nodes,1);
%       U_t(1:nodes,1) = uu(NODECO(1:nodes,1), NODECO(1:nodes,2));   
%       error = norm(U_t - U_F) / norm(U_t);
%       fprintf('discretization error is = %g\n',error); 

   %%% add noise 
   % Spread of the Gaussians be R% of the U_F values
   sigmas = 0.05 * U_F; 
   % Create the noise values to add 
 %  randNoise = rand(length(U_F),1) .* sigmas;
   
      randNoise = sin(4*pi*(1:nodes)/nodes) .* sigmas;

   % Add noise 
   U_F = U_F + randNoise; 

      
      ymeas(1:L,1) = zeros(L,1);
      for r=1:T 
         i(1:L,1) = ELNODE(1:L,r);
         for l=1:L         
             ymeas(l,1) = ymeas(l,1) + U_F(i(l))/3.;       
         end 
      end            
%%%      ymeas(1:L,1) = yd(xm,ym);    % measured state

      
      
   fprintf('----------------------------------------------------\n'); 
   fprintf('STARTING SQH\n');
   fprintf('----------------------------------------------------\n'); 
   fprintf('\n'); 
      
      %%% initialization
      epsilon = 1.0;
      %Tolerance for convergence
      kappa=10^-7;                
      % SQH Algorithm paramters
      eta  = 10e-8;
      zeta = 0.9;
      sigma= 1.1;
      kmax = 10000;
      
      %%% reference value for diff coefficient 
      a0diff = ones(L,1);
      
      acontrol(1:L,1) = a0diff(1:L,1);      % diff coeff initialization 
      anew(1:L,1) = zeros(L,1);        
      count_updates=1;  
      
      %%% initial calculations 
      
      %%% assemble model problem 
%      fprintf('assembling linear system for state equation ...\n');
      [A,F]=assem(ELNODE,NODECO,BONODE,EDGE, ...
                     INTPT,W,PHI,PHIX,PHIY, ...
                     INTPTE1,WE1,PHIE1, ...
                     INTPTE2,WE2,PHIE2, ...
                     INTPTE3,WE3,PHIE3,acontrol,frhs);

      %%% solve discrete system
%      fprintf('solving linear system for state ...\n');
      U_F = A \ F;
            
      
      %%% Tikhonv functional 
     [Jt] = jfunctional(ELNODE,NODECO,BONODE,EDGE, ...
                     INTPT,W,PHI,acontrol,U_F,ymeas);
     Jk(1)=Jt;
     fprintf(' Initial J = %g\n',Jt);

     %%% gradient of U_F 
     [Zx,Zy] = gradient(ELNODE,NODECO,BONODE,INTPT,U_F);
     
     %%% SQH loop
      for iter = 1:kmax
          
      %%% assemble adjoint problem
      
      %%% construct rhs for adjoint problem with U_F 
      ytilde(1:L,1) = zeros(L,1);
      for r=1:T 
         i(1:L,1) = ELNODE(1:L,r);
         for l=1:L         
             ytilde(l,1) = ytilde(l,1) + U_F(i(l))/3.;       
         end 
      end       
      frhsa(1:L,1) = ytilde(1:L,1) - ymeas(1:L,1);     % rhs (y-y_d) 
      
          
        
 %     fprintf('assembling linear system for adjoint equation ...\n');
      [A,F]=assem(ELNODE,NODECO,BONODE,EDGE, ...
                     INTPT,W,PHI,PHIX,PHIY, ...
                     INTPTE1,WE1,PHIE1, ...
                     INTPTE2,WE2,PHIE2, ...
                     INTPTE3,WE3,PHIE3,acontrol,frhsa);

      %%% solve discrete system
 %     fprintf('solving linear system for adjoint ...\n');
      U_A = A \ F;          
    
      %%% gradient of U_A
      [Qx,Qy] = gradient(ELNODE,NODECO,BONODE,INTPT,U_A);
      
      %%% update the diffusion coefficient 
      
     
             
%      fprintf('update the diffusion coefficient ...\n');
      [anew]=minHP(ELNODE,NODECO,BONODE,EDGE, ...
             INTPT,acontrol,U_F,U_A,Zx,Zy,Qx,Qy,epsilon);
              
      %%% tacitaly is the size of the volume element in epsilon
      tau = norm((acontrol - anew),2)^2 ;
 %        fprintf(' tau = %g\n',tau);
            
      %%% assemble model problem with anew
 %     fprintf('assembling linear system for state equation ...\n');
      [A,F]=assem(ELNODE,NODECO,BONODE,EDGE, ...
                     INTPT,W,PHI,PHIX,PHIY, ...
                     INTPTE1,WE1,PHIE1, ...
                     INTPTE2,WE2,PHIE2, ...
                     INTPTE3,WE3,PHIE3,anew,frhs);

      %%% solve discrete system
%      fprintf('solving linear system for state ...\n');
      U_Fnew = A \ F;
                 
     %%% Tikhonv functional 
     [Jt] = jfunctional(ELNODE,NODECO,BONODE,EDGE, ...
                     INTPT,W,PHI,anew,U_Fnew,ymeas);
 %    fprintf(' J = %g\n',Jt);
                 
    if(Jt-Jk(count_updates) > -eta * tau )  
        % The improvemnt to the cost functional value is not sufficiently large, 
        % increase epsilon
        epsilon=epsilon*sigma;   
     %   fprintf('Update step fails\n');
    else                 
        % The improvment to the cost functional value is sufficently large, 
        % decrease epsilon
        count_updates=count_updates+1;
        Jk(count_updates)=Jt;        
        acontrol = anew; 
        U_F = U_Fnew;          
        %%% gradient of U_F 
        [Zx,Zy] = gradient(ELNODE,NODECO,BONODE,INTPT,U_F);
        epsilon=epsilon*zeta;  
        
        % Print current values      
        fprintf('k %i | k up %i | J %e | eps %e | tau %e\n',iter,... 
            count_updates-1,Jt,epsilon,tau) 
        
    end 
    
        if(tau < kappa && iter > 100 )      
        fprintf('Convergence achieved | k %i | tau %e\n',iter,tau);
        break;
        end
             
      end


      %%% plot the solution
      fprintf('plotting solution ...\n');           
      figure(2)
      drawf(ELNODE,NODECO,BONODE,U_F);
      view(30,30);
      axis on; 
      
        print('-depsc2', 'ysolInv01.eps','-b0'); 
       % print('-dpdf', 'ysolInv01.pdf','-b0');
      

      
%       %%% plot the adjoint 
%       fprintf('plotting adjoint ...\n');           
%       figure(3)
%       drawf(ELNODE,NODECO,BONODE,U_A);
%       view(30,30);
%     
%%% plot the gradient of the solution at quadrature points 
%       fprintf('plotting gradient of solution ...\n');          
%       figure(4)
%       %%% sol Zm and the gradient of sol (Zx,Zy) at quad pts
%       Zm = zeros(L,1);
%       Zx = zeros(L,1);
%       Zy = zeros(L,1);
%       [Zm,Zx,Zy] = drawg(ELNODE,NODECO,BONODE,INTPT,U_F);
%       view(30,30);
%     
%%% plot the gradient of the adjoint at quadrature points 
%       fprintf('plotting gradient of adjoint ...\n');          
%       figure(5)
%       %%% sol Zm and the gradient of sol (Zx,Zy) at quad pts
%       Qm = zeros(L,1);
%       Qx = zeros(L,1);
%       Qy = zeros(L,1);
%       [Qm,Qx,Qy] = drawg(ELNODE,NODECO,BONODE,INTPT,U_A);
%       view(30,30);      
      
      figure(6)
      drawa(ELNODE,NODECO,BONODE,INTPT,anew);
      view(30,30);
      axis on; 

      print('-depsc2', 'acompInv01.eps','-b0'); 
    % print('-dpdf', 'acompInv01.pdf','-b0');

      
      figure(7)
      drawa(ELNODE,NODECO,BONODE,INTPT,atrue);
      view(30,30);
      
      print('-depsc2', 'atrueInv01.eps','-b0'); 
    % print('-dpdf', 'atrueInv01.pdf','-b0');

      
      %Plot convergence history
      figure(8)
      loglog(Jk,ltype{3},'Linewidth',2)
      % title('Convergence history of J');
       print('-depsc2', 'invJhist01.eps','-b0'); 
       % print('-dpdf', 'invJhist01.pdf','-b0');


