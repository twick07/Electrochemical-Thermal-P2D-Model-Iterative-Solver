%%%% Main Simulation Loop %%%%

for i = 1 : time_max
    %% Set timestep for current iterations
    %If externally driven
    if COMSOL_Driven == 1

        time_step = time_step_vector(i);
    end

    if time_variable == 1
        
        time_step = time_step_vector(i);

    end
   
    %% Get C_e for current timestep
   %Get Reduced C_e for solver based on C_e calculation approach
   if mass_con == 0
     C_e_reduced = C_e(i,:);
   end

   if mass_con == 1
    C_e_reduced = [C_e(i,1:N_neg) C_e(i,N_neg+2:N_neg+N_sep+1) C_e(i,N_neg+N_sep+3:N_tot_active+2)];
   end

   if mass_con == 2
     C_e_reduced = C_e(i,:);
   end

   if mass_con == 3
    C_e_reduced = [C_e(i,1:N_neg) C_e(i,N_neg+2:N_neg+N_sep+1) C_e(i,N_neg+N_sep+3:N_tot_active+2)];
   end

   if mass_con == 4
     C_e_reduced = C_e(i,:);
   end
 
   if mass_con == 5
     C_e_reduced = C_e(i,:);
   end
   
   %% Get Temperature Dependent Parameters For Current time step
   %Evaluate Temperature Dependent Material Parameters 
   %kappa, kd, D_s, D_e,Sigma,k_0
   Temperature_Material_Parameters;
    
   %% Get Open Circuit potential and temperature for current time step
   %Calculate current stoichiometry in the model 
   stoic(i,:) = [C_se(i,1:N_neg)./C_max(1:N_neg) zeros(1,N_sep) C_se(i,N_neg+N_sep+1:N_tot_active)./C_max(N_neg+N_sep+1:N_tot_active) ];
    
   x = stoic(i,1:N_neg);
   y =  stoic(i,N_neg+N_sep+1:N_tot_active);
    
   %Calculate OCPs
   Electrode_OCPs_Temperature;
        
   %Calculate exchange current
   ex_current(i,: ) = [ (K_0(1:N_neg).*sqrt(C_e_reduced(1:N_neg).*(C_max_neg*ones(1,N_neg) - C_se(i,1:N_neg)).*C_se(i,1:N_neg)))  zeros(1,N_sep) ...
          (K_0(N_neg+N_sep+1:N_tot_active).*sqrt(C_e_reduced(N_neg+N_sep+1:N_tot_active).*(C_max_pos*ones(1,N_pos) - C_se(i,N_neg+N_sep+1:N_tot_active)).* C_se(i,N_neg+N_sep+1:N_tot_active)))] ;

   %% Iterative Solver
   % Iterative Solver Algorithm For All Ageing Forms
   
     Iterative_Solver;
   
   
   %% Check if variable time_step and adjust
   
    if variation == 2
           
       variation_neg = (abs(max(j_flux(i,1:N_neg)) - min(j_flux(i,1:N_neg)))/abs(j_ave_neg))*100;
       variation_pos = (abs(max(j_flux(i,N_neg+N_sep+1:N_tot_active)) - min(j_flux(i,N_neg+N_sep+1:N_tot_active)))/abs(j_ave_pos))*100;
       
       if or((variation_neg > thresh_neg) , (variation_pos > thresh_pos)) == 1
           
           time_step = time_step_fine ;
           
       else
           
           time_step = time_step_norm ;
           
       end
           
    end
   
    %% Calculate psi_e, psi_s, eta, i_e and i_s and differentials of each
   %%% Calculate i_e and i_s
   
   A = tril(ones(N_tot_active));
   B = [zeros(1,N_tot);A];
   i_e(i,:) = (zeros(N_tot_active+1,1) + B*((F*As.*del.*j_tot(i,:))'))' ;
   
   %i_s(i,:) = I_app(i) - i_e(i,:); % Should be changed to use j_flux

   i_s(i,:) = (I_app(i)*ones(N_tot_active+1,1) - B*((F*As.*del.*j_tot(i,:))'))' ; 
   %%% Calculate psi_e

   psi_e(i,:) = psi_e(i,1) + (A*(diff_psi_e_m(i,1:end-1).*del)')';

   %%% Calculate eta

   eta_tot(i,:) = (R_const/(alpha*F))*T(i,x_neg_start:x_pos_end).*asinh((1/2).*(j_tot(i,:)./ex_current(i,:))) ;

   eta_tot(i,N_neg+1:N_neg+N_sep) = zeros(1,N_sep);

   %%% Calculate psi_s
   psi_s(i,1:N_neg) = eta(i,1:N_neg) + psi_e(i,1:N_neg) + OCP(i,1:N_neg)+ F*R_SEI_m(i,1:N_neg).*j_tot(i,1:N_neg);
    
   psi_s(i,N_neg+N_sep+1:N_tot_active) = eta(i,N_neg+N_sep+1:N_tot_active) + psi_e(i,N_neg+N_sep+1:N_tot_active) + OCP(i,N_neg+N_sep+1:N_tot_active);


   %%% Calculate eta_SEI

   if SEI > 0
       eta_SEI(i,1:N_neg) = psi_s(i,1:N_neg) - psi_e(i,1:N_neg) - OCP_SEI - F*R_SEI_m(i,1:N_neg).*j_tot(i,1:N_neg);
    
       j_SEI(i,1:N_neg)= -(I_app(i)<0)*(j_tot(i,1:N_neg)<0).*(i_0_SEI_v(i,1:N_neg)/F).*exp( -(alpha_SEI*F./(R_const.*T(i,1:N_neg))).*(eta_SEI(i,1:N_neg)) ) ;
   end
    
   %%% Calculate j_flux
   
   j_flux(i,:) = j_tot(i,:) - j_SEI(i,:);

   eta(i,:) = (R_const*T(i,:)./(alpha*F)).*asinh((1/2)*(j_flux(i,:)./ex_current(i,:))) ;

   eta(i,N_neg+1:N_neg+N_sep) = zeros(1,N_sep);


   %% Calculate Time Derivatives %%%
   %% Calculate C_e %%%

   if mass_con == 0

   %%% Use Standard FVM %%%
   % Calculate C_e using harmonic mean approach
   diff_C_e_1 = [ 0 diff(C_e(i,:),1) 0 ];
   
    
   del_ave_neg_sep = (del(N_neg) +  del(N_neg+1))/2 ;
   del_ave_sep_pos = (del(N_neg + N_sep) +  del(N_neg+N_sep+1))/2 ;

   del_extended = [del(1:N_neg) del_ave_neg_sep del(N_neg+2:N_neg+N_sep) del_ave_sep_pos del(N_neg+N_sep+1:N_tot_active)];
   
   % Creat Harmonic Mean at anode/seperator interface and seperator/cathode
   % interface
   beta_neg_sep = del(N_neg)/(del(N_neg) + del(N_neg+1));
   beta_sep_pos = del(N_neg+N_sep)/(del(N_neg+N_sep) + del(N_neg+N_sep+1));
   
   D_HM_neg_sep = (D_e(N_neg)*D_e(N_neg+1))/( beta_neg_sep*D_e(N_neg+1) + (1 - beta_neg_sep)*D_e(N_neg) );
   
   D_HM_sep_pos = (D_e(N_neg+N_sep)*D_e(N_neg+N_sep+1))/( beta_sep_pos*D_e(N_neg+N_sep+1) + (1 - beta_sep_pos)*D_e(N_neg+N_sep) );
   
   D_e_eff_ext = [D_e(1:N_neg) D_HM_neg_sep D_e(N_neg+2:N_neg+N_sep) D_HM_sep_pos D_e(N_neg+N_sep+1:N_tot_active)];
   
   %caclulate C_e time derivative
   diff_C_e_2 = (1./del).*diff((D_e_eff_ext.*diff_C_e_1)./del_extended,1); 
   difft_C_e =  (1./(eps_elec)).*diff_C_e_2 + (1-t_li)*(As./(eps_elec)).*j_tot(i,:);

 
   %Store time derivatives 
   difft_C_e_m(i,:) = difft_C_e;

   % Calculate C_e for next time_step
   C_e(i+1,:) = C_e(i,:) + time_step*difft_C_e ; 

    
   
   end 


   if mass_con == 1
   
   %%% Use FVM set up at boundary %%%

   eps_del_extended = [eps_neg_elec*del(1:N_neg) (eps_sep_elec*del(N_neg+1) + eps_neg_elec*del(N_neg))/2 eps_sep_elec*del(N_neg+1:N_neg+N_sep) (eps_sep_elec*del(N_neg+N_sep) + eps_pos_elec*del(N_neg+N_sep+1))/2 eps_pos_elec*del(N_neg+N_sep+1:N_tot)];
   D_e_extended = [D_e(1:N_neg) mean([D_e(N_neg) D_e(N_neg+1)]) D_e(N_neg+1:N_neg+N_sep) mean([D_e(N_neg+N_sep) D_e(N_neg+N_sep+1)]) D_e(N_neg+N_sep+1:N_tot_active)];
   A = eye([N_tot_active+2+1,N_tot_active+2]);
   A = 0.5*A + 0.5*circshift(A,[1 0]);
   A(1,1) = 1;
   A(end,end) = 1;
   D_e_extended = A*D_e_extended';
   D_e_extended = D_e_extended';
   del_extended = [del(1:N_neg) mean([del(N_neg) del(N_neg+1)]) del(N_neg+1:N_neg+N_sep) mean([del(N_neg+N_sep) del(N_neg+N_sep+1)]) del(N_neg+N_sep+1:N_tot_active)];
   del_extended = A*del_extended';
   del_extended = del_extended';
   
   diff_C_e_1 = [ 0 diff(C_e(i,:),1,2) 0 ]./(del_extended);

   C_e_flux_temp = D_e_extended.*diff_C_e_1;
   C_e_flux(i,:) = (diff(C_e_flux_temp,1,2));
   
  
   
   C_e_add_extended(i,:) = [(1-t_li)*del(1:N_neg).*As(1:N_neg).*j_tot(i,1:N_neg) 0.5*del(N_neg)*(1-t_li)*As_neg.*j_tot(i,N_neg) (1-t_li)*del(N_neg+1:N_neg+N_sep).*As(N_neg+1:N_neg+N_sep).*j_tot(i,N_neg+1:N_neg+N_sep) ...
       0.5*del(N_neg+N_sep+1)*(1-t_li)*As_pos.*j_tot(i,N_neg+N_sep+1) (1-t_li)*del(N_neg+N_sep+1:N_tot_active).*As(N_neg+N_sep+1:N_tot_active).*j_tot(i,N_neg+N_sep+1:N_tot_active)];
   
   difft_C_e_m(i,:) = (1./eps_del_extended).*( C_e_flux(i,:)  + C_e_add_extended(i,:) );

   
    
   C_e(i+1,:) = C_e(i,:) + time_step*difft_C_e_m(i,:) ; 

   end

   if mass_con == 2
   
   %%% Use Standard FDM implementation %%%

   %Create Differentiation Matrix
   A = eye(N_tot_active);
   B = -2*A + circshift(A,[1,0]) + circshift(A,[-1,0]);
   B(end,1) = 0;
   B(1,end) = 0;
   B(1,1) = -1;
   B(end,end) = -1;
    
   %Find time change vector

   difft_C_e_m(i,:) = (1./eps_elec).*( (D_e./(del.^2)).*( (B*(C_e(i,:)'))')  + (1-t_li)*As.*j_tot(i,:) )  ;
   

   C_e(i+1,:) = C_e(i,:) + time_step*difft_C_e_m(i,:) ; 


   end

   

   if mass_con == 3

   %%% Use FDM with FVM at at boundarys %%%
       

       %Create Differentiation Matrix
       A = eye(N_tot_active+2);
       B = -2*A + circshift(A,[1,0]) + circshift(A,[-1,0]);
       B(end,1) = 0;
       B(1,end) = 0;
       B(1,1) = -1;
       B(end,end) = -1;
        
       %Find time change vector
    
       difft_C_e_m(i,1:N_neg) = (1./eps_elec(1:N_neg)).*( (D_e(1:N_neg)./(del(1:N_neg).^2)).*( (B(1:N_neg,1:N_neg+1)*(C_e(i,1:N_neg+1)'))')  + (1-t_li)*As(1:N_neg).*j_tot(i,1:N_neg) )  ;

       difft_C_e_m(i,N_neg+1) = (2/(eps_elec(N_neg)*del(N_neg) + eps_elec(N_neg+1)*del(N_neg+1)))*( -(D_e(N_neg)./(del(N_neg)))*diff(C_e(i,N_neg:N_neg+1),1,2) ...
            + (D_e(N_neg+1)./(del(N_neg+1)))*diff(C_e(i,N_neg+1:N_neg+2),1,2) + 0.5*del_neg*(1-t_li)*As_neg.*j_tot(i,N_neg));

       difft_C_e_m(i,N_neg+2:N_neg + N_sep+1) = (1./eps_elec(N_neg+1:N_neg+N_sep)).*( (D_e(N_neg+1:N_neg+N_sep)./(del(N_neg+1:N_neg+N_sep).^2)).*( (B(N_neg+2:N_neg+N_sep+1,N_neg+1:N_neg+N_sep+2)*(C_e(i,N_neg+1:N_neg+N_sep+2)'))')  + (1-t_li)*As(N_neg+1:N_neg+N_sep).*j_tot(i,N_neg+1:N_neg+N_sep) )  ;

       difft_C_e_m(i,N_neg+N_sep+2) = (2/(eps_elec(N_neg+N_sep)*del(N_neg+N_sep) + eps_elec(N_neg+N_sep+1)*del(N_neg+N_sep+1)))*( -(D_e(N_neg+N_sep)./(del(N_neg+N_sep)))*diff(C_e(i,N_neg+N_sep+1:N_neg+N_sep+2),1,2) ...
            + (D_e(N_neg+N_sep+1)./(del(N_neg+N_sep+1)))*diff(C_e(i,N_neg+N_sep+2:N_neg+N_sep+3),1,2) + 0.5*(1-t_li)*del_pos*As_pos.*j_tot(i,N_neg+N_sep+1));

       difft_C_e_m(i,N_neg+N_sep+3:N_tot_active+2) = (1./eps_elec(N_neg+N_sep+1:N_tot_active)).*( (D_e(N_neg+N_sep+1:N_tot_active)./(del(N_neg+N_sep+1:N_tot_active).^2)).*( (B(N_neg+N_sep+3:N_tot_active+2,N_neg+N_sep+2:N_tot_active+2)*(C_e(i,N_neg+N_sep+2:N_tot_active+2)'))')  + ...
           (1-t_li)*As(N_neg+N_sep+1:N_tot_active).*j_tot(i,N_neg+N_sep+1:N_tot_active) )  ;

       C_e(i+1,:) = C_e(i,:) + time_step*difft_C_e_m(i,:) ; 


   end

   
   if mass_con == 4

   %%% Use standard FDM with ode15s %%%
   

       %Create Differentiation Matrix
       A = eye(N_tot_active);
       B = -2*A + circshift(A,[1,0]) + circshift(A,[-1,0]);
       B(end,1) = 0;
       B(1,end) = 0;
       B(1,1) = -1;
       B(end,end) = -1;
    
       op = odeset('reltol',tol);
       [t,C_e_temp] = ode15s(@(t,C_e_temp) ( (1./eps_elec').*( ((D_e./(del.^2)).*( (B*(C_e_temp))'))' + ((1-t_li)*(As).*j_tot(i,:))' ) ), [0 time_step], C_e(i,:)', op);
       C_e(i+1,:) = C_e_temp(end,:); % update
    
       difft_C_e_m(i,:) = (1/time_step)*(C_e(i+1,:)-C_e(i,:));
   
   end 

   if mass_con == 5
  
   %%% Use standard FVM with ode15s %%%
  
  
   % Make Extended Del Matrix
   del_ave_neg_sep = (del(N_neg) +  del(N_neg+1))/2 ;
   del_ave_sep_pos = (del(N_neg + N_sep) +  del(N_neg+N_sep+1))/2 ;

   del_extended = [del(1:N_neg) del_ave_neg_sep del(N_neg+2:N_neg+N_sep) del_ave_sep_pos del(N_neg+N_sep+1:N_tot_active)];
   
   % Creat Harmonic Mean at anode/seperator interface and seperator/cathode
   % interface

   beta_neg_sep = del(N_neg)/(del(N_neg) + del(N_neg+1));
   beta_sep_pos = del(N_neg+N_sep)/(del(N_neg+N_sep) + del(N_neg+N_sep+1));
   
   D_HM_neg_sep = (D_e(N_neg)*D_e(N_neg+1))/( beta_neg_sep*D_e(N_neg+1) + (1 - beta_neg_sep)*D_e(N_neg) );
   
   D_HM_sep_pos = (D_e(N_neg+N_sep)*D_e(N_neg+N_sep+1))/( beta_sep_pos*D_e(N_neg+N_sep+1) + (1 - beta_sep_pos)*D_e(N_neg+N_sep) );
   
   D_e_eff_ext = [D_e(1:N_neg) D_HM_neg_sep D_e(N_neg+2:N_neg+N_sep) D_HM_sep_pos D_e(N_neg+N_sep+1:N_tot_active)];
   
   %Caclulate using FVM
   op = odeset('reltol',tol);
   [t,C_e_temp] = ode15s(@(t,C_e_temp) ( (1./eps_elec').*( ((1./del).*diff((D_e_eff_ext.*[ 0 ; diff(C_e_temp,1) ; 0 ]')./del_extended,1))' + ((1-t_li)*(As).*j_tot(i,:))' ) ), [0 time_step], C_e(i,:)', op);
   
   
   C_e(i+1,:) = C_e_temp(end,:); % update

   
   difft_C_e_m(i,:) = (1/time_step)*(C_e(i+1,:)-C_e(i,:));

   

   end

   %% Calculate solid phase (C_s) diffusion
    if order == 1
      
        
        difft_C_s = [-3*(1./Rs(1:N_neg)) zeros(1,N_sep) -3*(1./Rs(N_neg+N_sep+1:N_tot_active))].*j_flux(i,:) ; 

        C_s(i+1,:) = C_s(i,:) + time_step*difft_C_s;

        C_s_ave(i+1,:) = C_s(i+1,:);
    
    
    
        % Calculate C_se
    
        %C_se(i+1,:) = C_s(i+1,:) -  (1/5).*[Rs(1:N_neg)./D_s_eff(1:N_neg) zeros(1,N_sep) Rs(N_neg+N_sep+1:N_tot_active)./D_s_eff(N_neg+N_sep+1:N_tot_active)].*j_flux(i,:);
    
        C_se(i+1,:) = C_s(i+1,:) -  (1/5).*[Rs(1:N_neg)./D_s(1:N_neg) zeros(1,N_sep) Rs(N_neg+N_sep+1:N_tot_active)./D_s(N_neg+N_sep+1:N_tot_active)].*j_flux(i,:);

    end
    
    if order == 2
        
        %Calculate C_s
        difft_C_s = [-3*(1./Rs(1:N_neg)) zeros(1,N_sep) -3*(1./Rs(N_neg+N_sep+1:N_tot_active))].*j_flux(i,:) ; 

        C_s(i+1,:) = C_s(i,:) + time_step*difft_C_s;
        
        %Calculate q_x
        difft_q_x = -30*[D_s_eff(1:N_neg)./(Rs(1:N_neg).^2) zeros(1,N_sep) D_s_eff(N_neg+N_sep+1:N_tot_active)./(Rs(N_neg+N_sep+1:N_tot_active).^2)].*q_x(i,:) - (45/2)*[(1./(Rs(1:N_neg).^2)) zeros(1,N_sep) (1./(Rs(N_neg+N_sep+1:N_tot_active).^2))].*j_flux(i,:);

        q_x(i+1,:) = q_x(i,:) + time_step*difft_q_x ;

        %Calculate C_se
        C_se(i+1,:) = C_s(i+1,:) -  (1/35).*[Rs(1:N_neg)./D_s_eff(1:N_neg) zeros(1,N_sep) Rs(N_neg+N_sep+1:N_tot_active)./D_s_eff(N_neg+N_sep+1:N_tot_active)].*j_flux(i,:) + 8*Rs.*q_x(i,:);

        
    end
    
    if order == 3
        for k = 1 : N_tot_active
            C_s_shell = (C_s_full(i,:,k));
             
            
            if k <N_neg+1
                N_flux_neg(i,:,k) = -D_s(k)*diff(C_s_shell)/dR_neg; % flux density at surfaces between "bins" [mols m^{-2} s^{-1}]
                M_neg = N_flux_neg(i,:,k).*Sa_neg(1:end-1); % total moles crossing surface between bins [mols s^{-1}]
                C_s_shell = C_s_shell + ([0 M_neg] - [M_neg 0])*time_step./dV_neg; % conc. change via diffusion
                C_s_shell(end) = C_s_shell(end) - j_flux(i,k)*Sa_neg(end)*time_step/dV_neg(end); % at boundary
                
                C_se(i+1,k) = C_s_shell(end); %Set Surface Concentration
                C_s_full(i+1,:,k) = C_s_shell;

                C_s_ave(i,k) = mean(C_s_shell);


                %op = odeset('reltol',tol);

                %[t,C_s_shell_temp] = ode15s(@(t,C_s_shell_temp) [ -time_step*(D_s(k)*diff(C_s_shell_temp(1:2))/dR_neg).*(Sa_neg(1))./dV_neg(1) ;diff(-(D_s(k)*diff(C_s_shell_temp(2:end-1))/dR_neg).*(Sa_neg(2:end-1)'))*time_step./(dV_neg(2:end-1)') ;time_step*(D_s(k)*diff(C_s_shell_temp(end-1:end))/dR_neg).*(Sa_neg(end-1))./dV_neg(end) - j_flux(i,k)*Sa_neg(end)*time_step/dV_neg(end)]  ,[0 time_step], C_s_shell',op);
   
                %[t,C_s_shell_temp] = ode15s(@(t,C_s_shell_temp) (D_s(k)./dV_neg').*[Sa_neg(1)*diff(C_s_shell_temp(1:2))./dR_neg ;diff(diff(Sa_neg'.*C_s_shell_temp./dR_neg)) ; -Sa_neg(end)*j_flux(i,k)./D_s(k) - Sa_neg(end-1)*diff(C_s_shell_temp(end-1:end))./dR_neg ],[0 time_step], C_s_shell',op);

            end

            if k>N_neg+N_sep
                N_flux_pos(i,:,k) = -D_s(k)*diff(C_s_shell)/dR_pos; % flux density at surfaces between "bins" [mols m^{-2} s^{-1}]
                M_pos = N_flux_pos(i,:,k).*Sa_pos(1:end-1); % total moles crossing surface between bins [mols s^{-1}]
                C_s_shell = C_s_shell + ([0 M_pos] - [M_pos 0])*time_step./dV_pos; % conc. change via diffusion
                C_s_shell(end) = C_s_shell(end) - j_flux(i,k)*Sa_pos(end)*time_step/dV_pos(end); % at boundary
                
                C_se(i+1,k) = C_s_shell(end); %Set Surface Concentration
                C_s_full(i+1,:,k) = C_s_shell;

                C_s_ave(i,k) = mean(C_s_shell);

            end
            

        end
    end

    if order == 4
       
       %Calculate Solid Concentration Change In Anode
       N_flux_neg(i,:,:) = diff(C_s_full(i,:,1:N_neg),1,2)/dR_neg;
       

       N_flux_neg_temp = -squeeze(N_flux_neg(i,:,:))*diag(D_s(1:N_neg));
       M_neg =   diag(Sa_neg(1:end-1),0)*N_flux_neg_temp ;
       
       difft_C_s_neg = diag(1./dV_neg,0)*([zeros(1,N_neg) ; M_neg] - [M_neg ; zeros(1,N_neg)]); 
       
       %Calculate Solid Concentration Change In Cathode
       N_flux_pos(i,:,:) = diff(C_s_full(i,:,N_neg+N_sep+1:N_tot_active),1,2)/dR_pos;

       N_flux_pos_temp = -squeeze(N_flux_pos(i,:,:))*diag(D_s(N_neg+N_sep+1:N_tot_active));
       
       M_pos =   diag(Sa_pos(1:end-1),0)*N_flux_pos_temp ;

       difft_C_s_pos = diag(1./dV_pos,0)*([zeros(1,N_pos) ; M_pos] - [M_pos ; zeros(1,N_pos)]); 
        
       %Calculate Solid Concentration Change in Entire Cell
       difft_C_s = [difft_C_s_neg zeros(N_shell,N_sep) difft_C_s_pos] ;

       
       
       C_s_shells_next = difft_C_s*time_step;
       
       % Add j_flux effect
       C_s_shells_next(end,:) = C_s_shells_next(end,:) - time_step*j_flux(i,:).*[(Sa_neg(end)/dV_neg(end))*ones(1,N_neg) zeros(1,N_sep) (Sa_pos(end)/dV_pos(end))*ones(1,N_pos)];
        
       C_s_full(i+1,:,:) = C_s_full(i,:,:) +  permute(C_s_shells_next,[3 1 2]);
       
       C_se(i+1,:) = squeeze(C_s_full(i+1,N_shell,:))'; %Set Surface Concentration
      
       C_s_ave(i+1,:) = squeeze(mean(C_s_full(i+1,:,:)))' ;

       % If heat of mixing is to be calculated
       if heat_of_mixing == 1
            
           Q_mix_calculation;

       end
                    

    end

    if order == 5
        
        op = odeset('reltol',tol);
   
        [t,C_s_ave_temp] = ode15s(@(t,C_s_ave_temp) ([-3*(1./Rs(1:N_neg)) zeros(1,N_sep) -3*(1./Rs(N_neg+N_sep+1:N_tot_active))].*j_flux(i,:))' , [0 time_step], C_s_ave(i,:)',op);

        C_s_ave(i+1,:) = C_s_ave_temp(end,:); % update

        difft_C_s = (1/time_step)*(C_s_ave(i+1,:) - C_s_ave(i,:));

        C_se(i+1,:) = C_s_ave(i+1,:) -  (1/5).*[Rs(1:N_neg)./D_s(1:N_neg) zeros(1,N_sep) Rs(N_neg+N_sep+1:N_tot_active)./D_s(N_neg+N_sep+1:N_tot_active)].*j_flux(i,:);
   


    end

    if order == 6

        for k = 1 : N_tot_active
            C_s_shell = (C_s_full(i,:,k));

            if k <N_neg+1

                op = odeset('reltol',tol);
                [t,C_s_shell_temp] = ode15s(@(t,C_s_shell_temp) (D_s(k)./dV_neg').*[Sa_neg(1)*diff(C_s_shell_temp(1:2))./dR_neg ;diff(Sa_neg(1:end-1)'.*diff(C_s_shell_temp)./dR_neg) ; -Sa_neg(end)*j_flux(i,k)./D_s(k) - Sa_neg(end-1)*diff(C_s_shell_temp(end-1:end))./dR_neg ],[0 time_step], C_s_shell',op);
                
                C_s_shell_next = C_s_shell_temp(end,:);

                C_se(i+1,k) = C_s_shell_next(end); %Set Surface Concentration
                C_s_ave(i,k) = mean(C_s_shell_next);
                C_s_full(i+1,:,k) = C_s_shell_next;
            end


            if k>N_neg+N_sep
                [t,C_s_shell_temp] = ode15s(@(t,C_s_shell_temp) (D_s(k)./dV_pos').*[Sa_pos(1)*diff(C_s_shell_temp(1:2))./dR_pos ;diff(Sa_pos(1:end-1)'.*diff(C_s_shell_temp)./dR_pos) ; -Sa_pos(end)*j_flux(i,k)./D_s(k) - Sa_pos(end-1)*diff(C_s_shell_temp(end-1:end))./dR_pos ],[0 time_step], C_s_shell',op);
                
                C_s_shell_next = C_s_shell_temp(end,:);

                C_se(i+1,k) = C_s_shell_next(end); %Set Surface Concentration
                C_s_ave(i,k) = mean(C_s_shell_next);
                C_s_full(i+1,:,k) = C_s_shell_next;

            end

        end
    end

    %Special Solid Diffusion Case For Long Cycles 

    if order == 7
        %Calculate Solid Concentration Change In Anode
       N_flux_neg = diff(C_s_shells(:,1:N_neg),1,1)/dR_neg;
       


       N_flux_neg_temp = -N_flux_neg*diag(D_s(1:N_neg));
       M_neg =   diag(Sa_neg(1:end-1),0)*N_flux_neg_temp ;
       
       difft_C_s_neg = diag(1./dV_neg,0)*([zeros(1,N_neg) ; M_neg] - [M_neg ; zeros(1,N_neg)]); 
       
       %Calculate Solid Concentration Change In Cathode
       N_flux_pos = diff(C_s_shells(:,N_neg+N_sep+1:N_tot_active),1,1)/dR_pos;
       

       N_flux_pos_temp = -N_flux_pos*diag(D_s(N_neg+N_sep+1:N_tot_active));
       
       M_pos =   diag(Sa_pos(1:end-1),0)*N_flux_pos_temp ;

       difft_C_s_pos = diag(1./dV_pos,0)*([zeros(1,N_pos) ; M_pos] - [M_pos ; zeros(1,N_pos)]); 
        
       %Calculate Solid Concentration Change in Entire Cell
       difft_C_s = [difft_C_s_neg zeros(N_shell,N_sep) difft_C_s_pos] ;

       
       
       C_s_shells_next = difft_C_s*time_step;
       
       % Add j_flux effect
       C_s_shells_next(end,:) = C_s_shells_next(end,:) - time_step*j_flux(i,:).*[(Sa_neg(end)/dV_neg(end))*ones(1,N_neg) zeros(1,N_sep) (Sa_pos(end)/dV_pos(end))*ones(1,N_pos)];
        

       C_s_shells = C_s_shells + C_s_shells_next;
       
       C_se(i+1,:) = C_s_shells(N_shell,:); %Set Surface Concentration
        
       C_s_ave(i+1,:) = mean(C_s_shells);
       


    end
    
    %% SEI Degredation Modeling
    
    %%% Kinetic Limited SEI Growth
    if SEI ==1
        %%% Add SEI Growth
        %j_SEI(i,:) = j_tot(i,:) - j_flux(i,:);
        L_SEI_m(i+1,:) = L_SEI_m(i,:) - time_step*(M_SEI_v./row_SEI_v).*j_SEI(i,:);
        R_SEI_m(i+1,:) = R_SEI_m(i,:) - time_step*(M_SEI_v./(row_SEI_v.*kappa_SEI_v)).*j_SEI(i,:)  ;
        %c_SEI_m(i+1,:) = c_SEI_m(i,:) - j_SEI(i,:);
        c_SEI_m(i+1,:) = c_SEI_m(i,:) - time_step*As.*j_SEI(i,:);

    end


    %%% Diffusion Limited SEI Growth

    if SEI ==2
        
        %Determine  R_SEI, L_SEI and c_SEI for next time step
        L_SEI_m(i+1,:) = L_SEI_m(i,:) - time_step*(M_SEI_v./row_SEI_v).*j_SEI(i,:);
        R_SEI_m(i+1,:) = R_SEI_m(i,:) - time_step*(M_SEI_v./(row_SEI_v.*kappa_SEI_v)).*j_SEI(i,:)  ;
        c_SEI_m(i+1,:) = c_SEI_m(i,:) - time_step*As.*j_SEI(i,:);
        %Determine C_sol for next time step
        %Condition for diffusion during charge
        %c_sol_m(i+1,1:N_neg) = ((L_SEI_m(i,1:N_neg).*j_SEI(i,1:N_neg))./(D_SEI_v(1:N_neg))) + epss_SEI*C_sol_0;

        
         %c_sol_charge       = ((L_SEI_m(i,1:N_neg).*j_SEI(i,1:N_neg))./(D_SEI_v(1:N_neg))) + epss_SEI*C_sol_0;

        %c_sol_m(i+1,1:N_neg) = max([zeros(1,N_neg) ; c_sol_m(i+1,1:N_neg)],[],1 );
        %Condition for diffusion during charge
        %diff_c_sol_surf = (D_SEI_v(1:N_neg)./L_SEI_m(i,1:N_neg)).*( ((epss_SEI*C_sol_0 - c_sol_m(i,1:N_neg))./L_SEI_m(i,1:N_neg))) ;

        %c_sol_m(i+1,1:N_neg) = c_sol_m(i,1:N_neg) + time_step*diff_c_sol_surf;

        %c_sol_discharge= c_sol_m(i,1:N_neg) + time_step*diff_c_sol_surf;

        if I_app(i) < 0
            %c_sol_m(i+1,1:N_neg) = max([zeros(1,N_neg) ; c_sol_m(i+1,1:N_neg)],[],1 );

            c_sol_temp = ((L_SEI_m(i,1:N_neg).*j_SEI(i,1:N_neg))./(D_SEI_v(1:N_neg))) + epss_SEI*C_sol_0;

            %0.5*epss_SEI*C_sol_0;

            c_sol_temp = (c_sol_m(i,1:N_neg) > 100).*max([zeros(1,N_neg) ; c_sol_temp],[],1 );

            c_sol_m(i+1,1:N_neg) = c_sol_temp ;

        end
        if I_app(i) > 0

            c_sol_m(i+1,1:N_neg) = epss_SEI*C_sol_0*ones(1,N_neg);

        end

        %c_sol_m(i+1,1:N_neg) = (I_app(i)<0)*max([zeros(1,N_neg) ; c_sol_m(i+1,1:N_neg)],[],1 );

        %c_sol_m(i+1,1:N_neg) = (I_app(i)>0)*epss_SEI*C_sol_0;

        %diff_c_sol_surf = (D_SEI_v(1:N_neg)./L_SEI_m(i,1:N_neg)).*( ((epss_SEI*C_sol_0 - c_sol_m(i,1:N_neg))./L_SEI_m(i,1:N_neg))  + j_SEI(i,1:N_neg)./D_SEI_v(1:N_neg) ) ;

        %c_sol_m(i+1,1:N_neg) = c_sol_m(i,1:N_neg) + time_step*diff_c_sol_surf;

        %c_sol_m(i+1,1:N_neg) = max([zeros(1,N_neg) ; c_sol_m(i+1,1:N_neg)],[],1 );
        %Determine i_0_SEI for next time step
        i_0_SEI_v(i+1,1:N_neg) = F*k_SEI_v(i,1:N_neg).*c_sol_m(i+1,1:N_neg);

    end

    %%% Lithium plating calculation

    if li_pl == 1
        c_li_m(i+1,:) = c_li_m(i,:) - time_step*As.*j_li_pl(i,:) - time_step*As.*j_li_st(i,:);
        L_SEI_m(i+1,:) =  L_SEI_m(i,:) - time_step*(M_li_v./row_Li_v).*j_li_pl(i,:) - time_step*(M_li_v./row_Li_v).*j_li_st(i,:);
        R_SEI_m(i+1,:) = L_SEI_m(i+1,:)./kappa_SEI_v ;
        i_o_pl_st_v(i+1,:) = F*k_pl_v.*(C_e(i,:).^(alpha_pl_neg)) ;
        

    end
    %% Calculate Temperature Change (T)

    %% Isothermal case %%
    if thermal == 0
        T(i+1,:) = T(i,:); 
    end

    %% Temperature Models for case of no current collector %%
    if (thermal >0) && (current_collector == 0) 
        % Create del_extended 
        del_ave_neg_sep = (del_neg + del_sep)/2 ;
        del_ave_sep_pos = (del_sep + del_pos)/2 ;
   
        del_extended = [del(1:N_neg) del_ave_neg_sep del(N_neg+2:N_neg+N_sep) del_ave_sep_pos del(N_neg+N_sep+1:N_tot_active) ];


        %Calculate Harmonic Mean for Lambda at Boundaries
        
        beta_neg_sep = del_neg/(del_neg + del_sep);
        beta_sep_pos = del_sep/(del_sep + del_pos);
        
   
        
        lambda_HM_neg_sep = (lambda(x_neg_end)*lambda(x_neg_end+1))/( beta_neg_sep*lambda(x_neg_end+1) + (1 - beta_neg_sep)*lambda(x_neg_end) );
   
        lambda_HM_sep_pos = (lambda(x_sep_end)*lambda(x_sep_end+1))/( beta_sep_pos*lambda(x_sep_end+1) + (1 - beta_sep_pos)*lambda(x_sep_end) );
   
        lambda_mean_neg_sep = mean([lambda(x_neg_end) lambda(x_neg_end+1)]);

        lambda_mean_sep_pos = mean([lambda(x_sep_end) lambda(x_sep_end+1)]);

        %lambda_HM_ncc_neg = mean([lambda_ncc lambda_neg]);
        %lambda_HM_neg_sep = mean([lambda_neg lambda_sep]);
        %lambda_HM_sep_pos = mean([lambda_sep lambda_pos]);
        %lambda_HM_pos_pcc = mean([lambda_pos lambda_acc]);
        %D_e_eff_ext = [D_e(1:N_neg) D_HM_neg_sep D_e(N_neg+2:N_neg+N_sep) D_HM_sep_pos D_e(N_neg+N_sep+1:N_tot_active)];
        
        %lambda_ext = [lambda(x_neg_start:x_neg_end) lambda_HM_neg_sep lambda(x_sep_start+1:x_sep_end) lambda_HM_sep_pos lambda(x_pos_start:x_pos_end)];

        %above works

        lambda_ext = [lambda(x_neg_start:x_neg_end) lambda_mean_neg_sep lambda(x_sep_start+1:x_sep_end) lambda_mean_sep_pos lambda(x_pos_start:x_pos_end)];

        %Calculate temperature change in time in active domains
        Q_rxn(i,:) = F*As.*j_flux(i,:).*eta(i,:);

        Q_rev(i,:) = F*As.*j_flux(i,:).*T(i,x_neg_start:x_pos_end).*diff_OCP_T(1,:);
    
        % Get diff_psi_s_m , diff_psi_e_m and diffx_C_e_log in centre node
        % points
        A = eye([N_tot_active,N_tot_active+1]);
        A = A + circshift(A,[-1 0]);
        A(end,end) = 1;
        A(end,1) = 0;
        A = 0.5*A;

        diff_psi_s_mid = A*diff_psi_s_m(i,:)';
        diff_psi_e_mid = A*diff_psi_e_m(i,:)';
        diffx_C_e_log_mid = A*diffx_C_e_log';
    
        %Changed
        diff_psi_s_mid = diff_psi_s_mid';
        diff_psi_e_mid = diff_psi_e_mid';
        diffx_C_e_log_mid = diffx_C_e_log_mid';


        %Q_ohm(i,:) = sigma(x_neg_start:x_pos_end).*(diff_psi_s_mid.^2) + kappa(i,1:N_tot_active).*(diff_psi_e_mid.^2) +...
            %kd(i,1:N_tot_active).*T(i,x_neg_start:x_pos_end).*diffx_C_e_log_mid.*diff_psi_e_mid;

        Q_ohm(i,:) = -A*(i_s(i,:).*diff_psi_s_m(i,:) + i_e(i,:).*diff_psi_e_m(i,:))';

        if heat_of_mixing == 0
    
        Q_gen(i,:) = [Q_ohm(i,:)+Q_rxn(i,:)+Q_rev(i,:)]; %
   
        end

        if heat_of_mixing == 1
    
        Q_gen(i,:) = [Q_ohm(i,:)+Q_rxn(i,:)+Q_rev(i,:)+Q_mix(i,:)]; %
   
        end

        if thermal == 5
            del_extended_mid = [del_neg*ones(1,N_neg) del_sep*ones(1,N_sep) del_pos*ones(1,N_pos) ];

            A = 0.5*eye(N_tot+1,N_tot);
        
            B = A + circshift(A,[1 0]);

            del_extended = B*del_extended_mid';

            del_extended = del_extended';

            %lambda_ext = B*lambda';

            %lambda_ext = lambda_ext';
   
            lambda_diff_T_x = [ -h_ncc*(T_amb-T(i,1))  lambda_ext(2:end-1).*diff(T(i,:),1,2)  -h_pcc*(T(i,end) - T_amb) ];

            %above works

            %lambda_diff_T_x = [ -h_ncc*(T_amb-T(i,1))  lambda_ext(2:end-1).*diff(T(i,:),1,2)./del_extended(2:end-1)  -h_pcc*(T(i,end) - T_amb) ]; 

            %lambda_diff_T_x = [ -h_ncc*(T_amb-T(i,1))  lambda_ext(2:end-1).*diff(T(i,:),1,2)  -h_pcc*(T(i,end) - T_amb) ];
        
            %./del_extended(2:end-1)
            diff_T_x_2_comp(i,:) = diff(lambda_diff_T_x,1,2) ; 
            %diff_T_x_2_comp(i,:) = diff(lambda_diff_T_x,1,2)./del_extended_mid ; 
            %diff_T_x_2_comp(i,:) = diff(lambda_diff_T_x,1,2) ;
            %./del_extended_mid
            %diff_T_t(i,:) = (1./(row.*C_p)).*( diff_T_x_2_comp(i,:)+ Q_gen(i,:) );
        
            diff_T_t(i,:) = (1./(row.*C_p.*(del_extended_mid))).*(diff_T_x_2_comp(i,:)+ Q_gen(i,:).*(del_extended_mid) ); 

            %above is correct

            %diff_T_t(i,:) = (1./(row.*C_p.*(del_extended_mid.^2))).*(diff_T_x_2_comp(i,:).*(del_extended_mid) + Q_gen(i,:).*(del_extended_mid.^2) ); 
            %diff_T_t(i,:) = (1./(row.*C_p.*(del_extended_mid))).*( diff_T_x_2_comp(i,:)+ Q_gen(i,:).*(del_extended_mid.^2) );
            T(i+1,:) = T(i,:) + time_step*diff_T_t(i,:);
    
            T_ave(i) = mean(T(i,:));
            T_ave_neg(i) = mean(T(i,x_neg_start:x_neg_end));
            T_ave_sep(i) = mean(T(i,x_sep_start:x_sep_end));
            T_ave_pos(i) = mean(T(i,x_pos_start:x_pos_end));

            %T_ave_2(i,:) = (1/(L_tot))*sum( T(i,:).*del );
            


        end

        % Use ode15s to calculate tempearature variation

        if thermal == 9
            
            A = 0.5*eye(N_tot+1,N_tot);
        
            B = A + circshift(A,[1 0]);

            del_extended = B*del';

            del_extended = del_extended';
            
            op = odeset('reltol',tol);
            %[t,T_temp] = ode15s(@(t,T_temp) ( (1./(row.*C_p.*(del))').*( (diff([-h_ncc*(T_amb-T_temp(1)) ; (lambda_ext(2:end-1)').*diff(T_temp,1) ; -h_pcc*(T_temp(end) - T_amb)]',1))' + (del.*Q_gen(i,:))' ) ) , [0 time_step], T(i,:)' , op);
            %[t,T_temp] = ode15s(@(t,T_temp) ( (1./(row.*C_p.*(del))').*( (diff( (([-h_ncc*(T_amb-T_temp(1)) ; (lambda_ext(2:end-1)').*diff(T_temp,1) ; -h_pcc*(T_temp(end) - T_amb)])./del_extended')',1))' + (del.*Q_gen(i,:))' ) ) , [0 time_step], T(i,:)' , op);

            %[t,T_temp] = ode15s(@(t,T_temp) ( (1./(row.*C_p)').*( (diff( (([-h_ncc*(T_amb-T_temp(1)) ; (lambda_ext(2:end-1)').*diff(T_temp,1) ; -h_pcc*(T_temp(end) - T_amb)])./del_extended')',1)./(del))' + (Q_gen(i,:))' ) ) , [0 time_step], T(i,:)' , op);

            [t,T_temp] = ode15s(@(t,T_temp) ( (1./(row.*C_p)').*( (diff( (([-h_ncc*(T_amb-T_temp(1)) ; (lambda_ext(2:end-1)').*diff(T_temp,1)./del_extended(2:end-1)' ; -h_pcc*(T_temp(end) - T_amb)]))',1)./(del))' + (Q_gen(i,:))' ) ) , [0 time_step], T(i,:)' , op);


            %[t,T_temp] = ode15s(@(t,T_temp) ( ((1./(row.*C_p.*(del)))').*( (diff( (([-h_ncc*(T_amb-T_temp(1)) ; (lambda_ext(2:end-1)').*diff(T_temp,1) ; -h_pcc*(T_temp(end) - T_amb)]))',1)./(del))' + (Q_gen(i,:))' ) ) , [0 time_step], T(i,:)' , op); %works

            %[t,T_temp] = ode15s(@(t,T_temp) ( (1./(row.*C_p)').*( (diff( (([-h_ncc*(T_amb-T_temp(1)) ; (lambda_ext(2:end-1)').*diff(T_temp,1) ; -h_pcc*(T_temp(end) - T_amb)])./del_extended')',1)./(del))' + (Q_gen(i,:)./(del))' ) ) , [0 time_step], T(i,:)' , op);

            %[t,T_temp] = ode15s(@(t,T_temp) ( (1./(row.*C_p)').*( (diff( (([-h_ncc*(T_amb-T(i,1)) ; (lambda_ext(2:end-1)').*diff(T_temp,1) ; -h_pcc*(T(i,end) - T_amb)])./del_extended')',1)./(del))' + (Q_gen(i,:))' ) ) , [0 time_step], T(i,:)' , op);

            T(i+1,:) = T_temp(end,:); % update

            diff_T_t(i,:) = (1/time_step)*(T(i+1,:)-T(i,:));

            T_ave(i) = mean(T(i,:));
            T_ave_neg(i) = mean(T(i,x_neg_start:x_neg_end));
            T_ave_sep(i) = mean(T(i,x_sep_start:x_sep_end));
            T_ave_pos(i) = mean(T(i,x_pos_start:x_pos_end));


        end

    end
    
    %% Temperature Models for case of current collector %%
    if (thermal >0) && (current_collector == 1) 

        % Create del_extended 
        del_ave_ncc_neg = (del_cc+del_neg)/2;
        del_ave_neg_sep = (del_neg + del_sep)/2 ;
        del_ave_sep_pos = (del_sep + del_pos)/2 ;
        del_ave_pos_pcc = (del_pos+del_cc)/2;
   
        del_extended = [del_cc*ones(1,N_cc) del_ave_ncc_neg del(1:N_neg-1) del_ave_neg_sep del(N_neg+2:N_neg+N_sep) del_ave_sep_pos del(N_neg+N_sep+2:N_tot_active) del_ave_pos_pcc del_cc*ones(1,N_cc)];

    
        diff_T_x(i,:) = [-(h_ncc/(lambda(1)))*(T_amb-T(i,1)) diff(T(i,:),1,2)./del_extended(2:N_tot)  -(h_pcc/(lambda(end)))*(T(i,end) - T_amb)]; %changed
    
        %Apply Contiunity at interfaces
        %diff_T_x(i,x_ncc_end+1) = 0;%mean([diff_T_x(i,x_ncc_end) diff_T_x(i,x_ncc_end+2)]);
        %diff_T_x(i,x_neg_end+1) = 0;%mean([diff_T_x(i,x_neg_end) diff_T_x(i,x_neg_end+2)]);
        %diff_T_x(i,x_sep_end+1) = 0;%mean([diff_T_x(i,x_sep_end) diff_T_x(i,x_sep_end+2)]);
        %diff_T_x(i,x_pos_end+1) = 0;%mean([diff_T_x(i,x_pos_end) diff_T_x(i,x_pos_end+2)]);

        %Calculate Harmonic Mean for Lambda at Boundaries
        beta_ncc_neg = del_cc/(del_cc + del_neg);
        beta_neg_sep = del_neg/(del_neg + del_sep);
        beta_sep_pos = del_sep/(del_sep + del_pos);
        beta_pos_pcc = del_pos/(del_pos + del_cc);
   
        lambda_HM_ncc_neg = (lambda(N_cc)*lambda(N_cc+1))/( beta_ncc_neg*lambda(N_cc+1) + (1 - beta_ncc_neg)*lambda(N_cc) );

        lambda_HM_neg_sep = (lambda(x_neg_end)*lambda(x_neg_end+1))/( beta_neg_sep*lambda(x_neg_end+1) + (1 - beta_neg_sep)*lambda(x_neg_end) );
   
        lambda_HM_sep_pos = (lambda(x_sep_end)*lambda(x_sep_end+1))/( beta_sep_pos*lambda(x_sep_end+1) + (1 - beta_sep_pos)*lambda(x_sep_end) );
   
        lambda_HM_pos_pcc = (lambda(x_pos_end)*lambda(x_pos_end+1))/( beta_ncc_neg*lambda(x_pos_end+1) + (1 - beta_ncc_neg)*lambda(x_pos_end) );

        %lambda_HM_ncc_neg = mean([lambda_ncc lambda_neg]);
        %lambda_HM_neg_sep = mean([lambda_neg lambda_sep]);
        %lambda_HM_sep_pos = mean([lambda_sep lambda_pos]);
        %lambda_HM_pos_pcc = mean([lambda_pos lambda_acc]);

        lambda_ext = [lambda(1:N_cc) lambda_HM_ncc_neg lambda(x_neg_start+1:x_neg_end) lambda_HM_neg_sep lambda(x_sep_start+1:x_sep_end) lambda_HM_sep_pos lambda(x_pos_start+1:x_pos_end) lambda_HM_pos_pcc lambda(x_pcc_start:N_tot)];
   
   
    

        %Calculate temperature change in time in each domain
        %Calculate temperature change in time in current collectors
        diff_T_t_ncc = (1./(row(1:N_cc).*C_p(1:N_cc))).*( (1./(del_extended(1:N_cc))).*((diff(lambda_ext(1:N_cc+1).*diff_T_x(i,1:N_cc+1),1,2))) + (I_app(i)^2)./(sigma(1:N_cc)) ) ;

        diff_T_t_pcc = (1./(row(x_pcc_start:x_end).*C_p(x_pcc_start:x_end))).*( (1./(del_extended(x_pcc_start:x_end))).*(diff(lambda_ext(x_pcc_start:x_end+1).*diff_T_x(i,x_pcc_start:x_end+1),1,2)) + (I_app(i)^2)./(sigma(x_pcc_start:x_end))) ;
    

        %Calculate temperature change in time in active domains
        Q_rxn(i,:) = F*As.*j_flux(i,:).*eta(i,:);

        Q_rev(i,:) = F*As.*j_flux(i,:).*T(i,x_neg_start:x_pos_end).*diff_OCP_T(1,:);
    
        % Get diff_psi_s_m , diff_psi_e_m and diffx_C_e_log in centre node
        % points
        A = eye([N_tot_active,N_tot_active+1]);
        A = A + circshift(A,[-1 0]);
        A(end,end) = 1;
        A(end,1) = 0;
        A = 0.5*A;

        diff_psi_s_mid = A*diff_psi_s_m(i,:)';
        diff_psi_e_mid = A*diff_psi_e_m(i,:)';
        diffx_C_e_log_mid = A*diffx_C_e_log';
    
        %Changed
        diff_psi_s_mid = diff_psi_s_mid';
        diff_psi_e_mid = diff_psi_e_mid';
        diffx_C_e_log_mid = diffx_C_e_log_mid';


        Q_ohm(i,:) = sigma(x_neg_start:x_pos_end).*(diff_psi_s_mid.^2) + kappa(i,1:N_tot_active).*(diff_psi_e_mid.^2) +...
            kd(i,1:N_tot_active).*T(i,x_neg_start:x_pos_end).*diffx_C_e_log_mid.*diff_psi_e_mid;

        Q_ncc = ((I_app(i)^2)./(sigma(1:N_cc))).*ones(1,N_cc);
        Q_pcc = ((I_app(i)^2)./(sigma(x_pcc_start:x_end))).*ones(1,N_cc);
    
        Q_gen(i,:) = [Q_ncc  Q_ohm(i,:)+Q_rxn(i,:)+Q_rev(i,:)  Q_pcc]; %
    
        %Caculate Active T change
        Q_diff_T_x = (1./(del_extended(1:x_end))).*(diff(lambda_ext(1:x_end+1).*diff_T_x(i,1:x_end+1),1,2));

        diff_T_t_active = (1./(row(x_neg_start:x_pos_end).*C_p(x_neg_start:x_pos_end))).*( (1./(del_extended(x_neg_start:x_pos_end))).*(diff(lambda_ext(x_neg_start:x_pos_end+1).*diff_T_x(i,x_neg_start:x_pos_end+1),1,2)) + ...
            + Q_rxn(i,:) + Q_rev(i,:) + Q_ohm(i,:));

        del_extended_mid = [del_cc*ones(1,N_cc) del_neg*ones(1,N_neg) del_sep*ones(1,N_sep) del_pos*ones(1,N_pos) del_cc*ones(1,N_cc)];
    
        diff_T_x_2_comp(i,:) = (1./(del_extended_mid)).*(diff(lambda_ext(1:x_end+1).*diff_T_x(i,1:x_end+1),1,2));

        %(1./(del_extended(1:x_end))).*
        %Modifications
        %diff_T_x_2_comp(i,x_neg_start) = 0;
        %diff_T_x_2_comp(i,x_ncc_end) = diff_T_x_2_comp(i,x_neg_start);
        %diff_T_x_2_comp(i,x_sep_start) = 0;
        %diff_T_x_2_comp(i,x_neg_end) = diff_T_x_2_comp(i,x_sep_start);
        %diff_T_x_2_comp(i,x_pos_start) = 0;
        %diff_T_x_2_comp(i,x_sep_end) = diff_T_x_2_comp(i,x_pos_start);
        %diff_T_x_2_comp(i,x_pos_end) = 0;
        %diff_T_x_2_comp(i,x_pcc_start) = diff_T_x_2_comp(i,x_pos_end);

        %Calculate Temperature change
        %diff_T_t(i,:) = [diff_T_t_ncc diff_T_t_active diff_T_t_pcc];
        %diff_T_t(i,:) = (1./(row(1:x_end).*C_p(1:x_end))).*(diff_T_x_2_comp(i,:) + Q_gen(i,:));
    
   
        %Nonisothermal case
        if thermal == 1
            diff_T_t(i,:) = [diff_T_t_ncc diff_T_t_active diff_T_t_pcc];
            T(i+1,:) = T(i,:) + time_step*diff_T_t(i,:);
        %Add BC for CCs
            %   T(i+1,x_ncc_end+1) = mean([T(i+1,x_ncc_end) T(i+1,x_ncc_end+1)]);
            %T(i+1,x_ncc_end) = T(i+1,x_ncc_end+1);
            %T(i+1,x_neg_end) = mean([T(i+1,x_neg_end) T(i+1,x_neg_end+1)]);
            %T(i+1,x_sep_start) = T(i+1,x_neg_end);
            %T(i+1,x_pos_start) = mean([T(i+1,x_sep_end) T(i+1,x_pos_start)]);
            %T(i+1,x_sep_end) = T(i+1,x_pos_start);
            %T(i+1,x_pcc_start) = mean([T(i+1,x_pos_end) T(i+1,x_pcc_start)]);
            %T(i+1,x_pos_end) = T(i+1,x_pcc_start);
            %T(i+1,x_neg_end) = T(i+1,x_neg_end+1);
            %T(i+1,x_sep_end) = T(i+1,x_sep_end+1);
            %T(i+1,x_pos_end) = T(i+1,x_pos_end+1);
            %T(i+1,x_pcc_start) = T(i+1,x_pcc_start-1);
        end
    
        if thermal == 2
            %diff_T_x_2_comp(i,x_neg_start) = 0;
            %diff_T_x_2_comp(i,x_neg_start-1) = diff_T_x_2_comp(i,x_neg_start);
            %diff_T_x_2_comp(i,x_sep_start-1) = diff_T_x_2_comp(i,x_sep_start);
            %diff_T_x_2_comp(i,x_sep_start) = 0;
            %diff_T_x_2_comp(i,x_pos_start-1) = diff_T_x_2_comp(i,x_pos_start);
            %diff_T_x_2_comp(i,x_pos_start) = 0;
            %diff_T_x_2_comp(i,x_pcc_start) = diff_T_x_2_comp(i,x_pcc_start-1);
            %diff_T_x_2_comp(i,x_pcc_start-1) = 0;

            diff_T_t(i,:) = (1./(row.*C_p)).*(diff_T_x_2_comp(i,:)  + Q_gen(i,:));
            T(i+1,:) = T(i,:) + time_step*diff_T_t(i,:);
            T_ave(i) = mean(T(i,:));
            %diff_T_x_2_comp(i,:)

        end
    
        if thermal == 3
        A = eye([N_tot_active,N_tot_active+2]);
            B = A + -2*circshift(A,[-1 0])+ circshift(A,[-1 1]);
            B(end,1) = 0;
            B(end,2) = 0;
            B(end,end) = 1;
            B(end,end-1) = -2;
            T_x = B*T(i,x_neg_start-1:x_pos_end+1)';
            T_x = T_x';
            diff_T_t_active = (1./(row(x_neg_start:x_pos_end).*C_p(x_neg_start:x_pos_end)).*(del)).*( ...
                lambda(x_neg_start:x_pos_end).*T_x./(del) + Q_gen(i,x_neg_start:x_pos_end).*del_extended_mid(x_neg_start:x_pos_end) );
    
            %diff_T_t(i,:) = [diff_T_t_ncc diff_T_t_active diff_T_t_pcc];

            diff_T_t(i,:) = [del_extended(1:N_cc).*diff_T_t_ncc diff_T_t_active del_extended(1:N_cc).*diff_T_t_pcc];

            T(i+1,:) = T(i,:) + time_step*diff_T_t(i,:);
    
            T_ave(i) = mean(T(i,:));
            %Apply BCs
            %T(i+1,x_ncc_end) = T(i+1,x_ncc_end+1);
            %T(i+1,x_neg_end+1) = T(i+1,x_neg_end);
            %T(i+1,x_sep_end) = T(i+1,x_sep_end+1);
            %T(i+1,x_pos_end+1) = T(i+1,x_pos_end);
    
        end
        if thermal == 4

            Q_gen(i,:) = [Q_ncc  Q_ohm(i,:)+Q_rxn(i,:)+Q_rev(i,:)  Q_pcc]; %

            A = 0.5*eye(N_tot+1,N_tot);
        
            B = A + circshift(A,[1 0]);

            del_extended = B*del_extended_mid';

            del_extended = del_extended';

            %lambda_ext = B*lambda';

            %lambda_ext = lambda_ext';
   

            lambda_diff_T_x = [ h_ncc*(T_amb-T(i,1))  lambda_ext(2:end-1).*diff(T(i,:),1,2)./del_extended(2:end-1)  -h_pcc*(T(i,end) - T_amb) ];

            diff_T_x_2_comp(i,:) = diff(lambda_diff_T_x,1,2) ;

            diff_T_t_ncc = (1./(row(1:N_cc).*C_p(1:N_cc))).*(diff_T_x_2_comp(i,1:N_cc) + Q_gen(i,1:N_cc));

            diff_T_t_pcc = (1./(row(x_pcc_start:x_end).*C_p(x_pcc_start:x_end))).*(diff_T_x_2_comp(i,x_pcc_start:x_end) + Q_gen(i,x_pcc_start:x_end));

            A = eye([N_tot_active,N_tot_active+2]);
            B = A + -2*circshift(A,[-1 0])+ circshift(A,[-1 1]);
            B(end,1) = 0;
            B(end,2) = 0;
            B(end,end) = 1;
            B(end,end-1) = -2;
            T_x = B*T(i,x_neg_start-1:x_pos_end+1)';
            T_x = T_x';
    
   

            diff_T_t_active = (1./(row(x_neg_start:x_pos_end).*C_p(x_neg_start:x_pos_end).*(del))).*( ...
                lambda(x_neg_start:x_pos_end).*T_x + Q_gen(i,x_neg_start:x_pos_end).*(del));
            %.*del_extended_mid(x_neg_start:x_pos_end)
            %.*del_extended_mid
            diff_T_t(i,:) = [diff_T_t_ncc diff_T_t_active diff_T_t_pcc];

            T(i+1,:) = T(i,:) + time_step*diff_T_t(i,:);
    
            T_ave(i) = mean(T(i,:));




        end


        if thermal == 5

            %Calculate temperature change in time in active domains
            Q_rxn(i,:) = F*As.*j_flux(i,:).*eta(i,:);
    
            Q_rev(i,:) = F*As.*j_flux(i,:).*T(i,x_neg_start:x_pos_end).*diff_OCP_T(1,:);
        
            % Get diff_psi_s_m , diff_psi_e_m and diffx_C_e_log in centre node
            % points
            A = eye([N_tot_active,N_tot_active+1]);
            A = A + circshift(A,[-1 0]);
            A(end,end) = 1;
            A(end,1) = 0;
            A = 0.5*A;
    
            diff_psi_s_mid = A*diff_psi_s_m(i,:)';
            diff_psi_e_mid = A*diff_psi_e_m(i,:)';
            diffx_C_e_log_mid = A*diffx_C_e_log';
        
            %Changed
            diff_psi_s_mid = diff_psi_s_mid';
            diff_psi_e_mid = diff_psi_e_mid';
            diffx_C_e_log_mid = diffx_C_e_log_mid';
    
    
            Q_ohm(i,:) = sigma(x_neg_start:x_pos_end).*(diff_psi_s_mid.^2) + kappa(i,1:N_tot_active).*(diff_psi_e_mid.^2) +...
                kd(i,1:N_tot_active).*T(i,x_neg_start:x_pos_end).*diffx_C_e_log_mid.*diff_psi_e_mid;
    
            %Q_ohm(i,:) = -A*(i_s(i,:).*diff_psi_s_m(i,:) + i_e(i,:).*diff_psi_e_m(i,:));
            Q_ncc = ((I_app(i)^2)./(sigma(1:N_cc))).*ones(1,N_cc);
            Q_pcc = ((I_app(i)^2)./(sigma(x_pcc_start:x_end))).*ones(1,N_cc);
    
            Q_gen(i,:) = [Q_ncc  Q_ohm(i,:)+Q_rxn(i,:)+Q_rev(i,:)  Q_pcc]; %
           

        

            A = 0.5*eye(N_tot+1,N_tot);
        
            B = A + circshift(A,[1 0]);

            del_extended = B*del_extended_mid';

            del_extended = del_extended';

            %lambda_ext = B*lambda';

            %lambda_ext = lambda_ext';
   

            lambda_diff_T_x = [ -h_ncc*(T_amb-T(i,1))  lambda_ext(2:end-1).*diff(T(i,:),1,2)  -h_pcc*(T(i,end) - T_amb) ];

            %lambda_diff_T_x = [ -h_ncc*(T_amb-T(i,1))  lambda_ext(2:end-1).*diff(T(i,:),1,2)  -h_pcc*(T(i,end) - T_amb) ];
        
            %./del_extended(2:end-1)
            diff_T_x_2_comp(i,:) = diff(lambda_diff_T_x,1,2) ; 
            %diff_T_x_2_comp(i,:) = diff(lambda_diff_T_x,1,2)./del_extended_mid ; 
            %diff_T_x_2_comp(i,:) = diff(lambda_diff_T_x,1,2) ;
            %./del_extended_mid
            %diff_T_t(i,:) = (1./(row.*C_p)).*( diff_T_x_2_comp(i,:)+ Q_gen(i,:) );
        
            diff_T_t(i,:) = (1./(row.*C_p.*del_extended_mid)).*( diff_T_x_2_comp(i,:)+ Q_gen(i,:).*del_extended_mid );

            %diff_T_t(i,:) = (1./(row.*C_p.*(del_extended_mid))).*( diff_T_x_2_comp(i,:)+ Q_gen(i,:).*(del_extended_mid.^2) );
            T(i+1,:) = T(i,:) + time_step*diff_T_t(i,:);
    
            T_ave(i) = mean(T(i,:));
            T_ave_ncc(i) = mean(T(i,1:N_cc));
            T_ave_neg(i) = mean(T(i,x_neg_start:x_neg_end));
            T_ave_sep(i) = mean(T(i,x_sep_start:x_sep_end));
            T_ave_pos(i) = mean(T(i,x_pos_start:x_pos_end));
            T_ave_pcc(i) = mean(T(i,x_pcc_start:x_end));


        end


        if thermal == 6
            A = 0.5*eye(N_tot+1,N_tot);
        
            B = A + circshift(A,[1 0]);

            del_extended = B*del_extended_mid';

            del_extended = del_extended';

            %lambda_ext = B*lambda';

            %lambda_ext = lambda_ext';
   

            lambda_diff_T_x = [ -h_ncc*(T_amb-T(i,1))  (1./del_extended(2:end-1)).*lambda_ext(2:end-1).*diff(T(i,:),1,2) -h_pcc*(T(i,end) - T_amb) ];

            %lambda_diff_T_x = [ -h_ncc*(T_amb-T(i,1))  lambda_ext(2:end-1).*diff(T(i,:),1,2)  -h_pcc*(T(i,end) - T_amb) ];
        
            %./del_extended(2:end-1)
            diff_T_x_2_comp(i,:) = diff(lambda_diff_T_x,1,2) ; 
            %./del_extended_mid
            %diff_T_t(i,:) = (1./(row.*C_p)).*( diff_T_x_2_comp(i,:)+ Q_gen(i,:) );
        
            diff_T_t(i,:) = (1./(row.*C_p.*del_extended_mid)).*( diff_T_x_2_comp(i,:)+ Q_gen(i,:).*del_extended_mid  );
            T(i+1,:) = T(i,:) + time_step*diff_T_t(i,:);
    
            T_ave(i) = mean(T(i,:));


        end
    
    %Lumped Element Case
        if thermal == 7
            T_ave(i) = mean(T(i,:));
            diff_T_t_ave = (1./(mean_row*mean_C_p)*L_tot)*( L_tot_active*mean(Q_gen(i,x_neg_start:x_pos_end)) - mean([ h_ncc h_pcc])*(T_ave(i) - T_amb));
    
            %L_tot_active*sum
            T_ave(i+1)  = T_ave(i) + time_step*diff_T_t_ave;

            T(i+1,:) = T_ave(i+1)*ones(1,N_tot);

        end
    %lumped sum with individual domain temperatures
        if thermal == 8
            T_ncc = mean(T(i,1:N_cc));
            T_neg = mean(T(i,x_neg_start:x_neg_end));
            T_sep = mean(T(i,x_sep_start:x_sep_end));
            T_pos = mean(T(i,x_pos_start:x_pos_end));
            T_pcc = mean(T(i,x_pcc_start:x_end));

            diff_T_t_ncc = (1./(row_ncc*C_p_ncc))*(L_ncc*sum(Q_gen(i,1:x_ncc_end))+h_ncc*(T_amb - 2*T_ncc));
            diff_T_t_neg = (1./(row_neg*C_p_neg))*(L_neg*sum(Q_gen(i,1:x_neg_end))+lambda_neg*(T_ncc - 2*T_neg+T_sep));
            diff_T_t_sep = (1./(row_sep*C_p_sep))*(L_sep*sum(Q_gen(i,1:x_neg_end))+h_ncc*(T_ncc - 2*T_neg+T_sep));

        end

        

    end
    

    

    %% Verlet Time Integration For C_e and Temperature

    %%% Verlet Integration
    if time_int == 2
        if i >2

            %%% For C_e
            
            C_e_accel = (1/time_step)*(C_e(i+1,:) - 2*C_e(i,:) + C_e(i-1,:)) ;

            C_e(i+1,:) = 2*C_e(i,:) - C_e(i-1,:) + (((time_step^2)/2))*C_e_accel;

           

            %%% For Temperature
            T_accel = (1/time_step)*(T(i+1,:) - 2*T(i,:) + T(i-1,:)) ;

            T(i+1,:) = 2*T(i,:) - T(i-1,:) + ((time_step)^2)*T_accel;

            %%% For C_s
            
            if order == 5

                C_se_accel = (1/time_step)*(C_se(i+1,:) - 2*C_se(i,:) + C_se(i-1,:)) ;

                C_se(i+1,:) = 2*C_se(i,:) - C_se(i-1,:) + ((time_step)^2)*C_se_accel;

                C_s_ave_accel = (1/time_step)*(C_s_ave(i+1,:) - 2*C_s_ave(i,:) + C_s_ave(i-1,:)) ;

                C_s_ave(i+1,:) = 2*C_s_ave(i,:) - C_s_ave(i-1,:) + ((time_step)^2)*C_s_ave_accel;

                

            end
            %}

        end

    end
    
    %% Simulation End Conditions
    V_cell(i) = psi_s(i,N_tot_active) - psi_s(i,1);
    
    capacity = 100*mean(stoic(i,N_neg+N_sep+1:N_tot_active));
    time = time + time_step  ;
    time_vector(i+1) = time;

    %Conditions for cycling
    if COMSOL_Driven == 0
        if I_app(i) ~= 0
            I_app(i+1) = I_app(i);
    
        end
    
        %Simulation Stop Condition
         if I_app(i) > 0
            if or(100*max(stoic(i,N_neg+N_sep+1:N_tot_active))> Cut_off_capacity,V_cell(i)< V_cut_off_low) %98.3 for 0.1 timestep
                I_app(i+1) = -I_app(i); % Change applied current
                %break;
            end
    
        end
    
         if I_app(i) < 0
            if or(100*max(stoic(i,N_neg+N_sep+1:N_tot_active))> Cut_off_capacity,V_cell(i)> V_cut_off_high) %98.3 for 0.1 timestep
                I_app(i+1) = -I_app(i); % Change applied current
                %break;
            end
    
         end
    
         %%% Cycle_counter
         if I_app(i)*I_app(i+1) < 0
            cycle_point = cycle_point+1;
            cycle_count(cycle_point) = i;

            if mod(cycle_point,2) ~= 0
                output_cycle_count = [num2str(cycle_point/2 - 0.5), ' number of cycles completed'];
                disp(output_cycle_count);
            end
    
            use_prev_j = 0;
            rep_j = 0;
    
            if cycle_point > 2*No_of_cycles
                break;
    
            end
    
         else
    
             use_prev_j = 1;
    
         end

    end
    
    %Conditions for a drive cycle
    if COMSOL_Driven ~= 0

        if i == time_max

            break;
        end


        use_prev_j = 0;
        rep_j = 0;
        


    end
    %}
    
end

%%% End of Main Simulation

%%%% Post Processing %%%%%
%%% Calculate Cell Capacity During Simulation
Q_tot = (1000/3600)*cumsum(I_app)*time_step; % Capacity change over the simulation [mAh/m^2]