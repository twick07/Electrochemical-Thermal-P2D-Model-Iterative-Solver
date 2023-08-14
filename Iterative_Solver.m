%% Beginning of Iterative Solver
if I_app(i) ~= 0
    %% Iterative Solver For Negative Electrode
    
    if use_prev_j == 1
    if I_app(i) == 0
    
        rep_j = 0;
        j_ave_neg = I_app(i)/(As_neg*F*L_neg)  ;
    
    else
    
        rep_j = rep_j +1 ;
    
        if rep_j == 1
        
            j_ave_neg = I_app(i)/(As_neg*F*L_neg)  ;
    
        end
    
       
    
    end
    end
    
    if use_prev_j == 0
    j_ave_neg = I_app(i)/(As_neg*F*L_neg)  ;
    end
    
    %Set inital guesses
    beta1 = 0.1;
    beta2 = 1.5;
    beta3 = (beta1 + beta2)/2 ;
    
    
    j_sol = zeros(3,N_tot_active);
    i_e_sol = zeros(3,N_tot_active);
    i_s_sol = zeros(3,N_tot_active);
    
    %Set Initial j guesses
    %j1 = beta1*j_ave_neg;
    %j2 = beta2*j_ave_neg;
    %beta3 = (beta1 +beta2)/2;
    %j3 = beta3*j_ave_neg;
    
    %j_attempts = [j1 j2 j3];
    
    err = ones(3,1);
    counter_convergence = 0;
    counter_boundary_finder = 0;
    ITP_count = 0;
    loop = 0;
    while norm(err(3),2) > tol
    
    if rep_j < 2
        j1 = beta1*j_ave_neg;
        j2 = beta2*j_ave_neg;
        %beta3 = (beta1 +beta2)/2;
        j3 = beta3*j_ave_neg;
    
    else
        j1 = beta1*j_flux(i-1,1);
        j2 = beta2*j_flux(i-1,1);
        %beta3 = (beta1 +beta2)/2;
        j3 = beta3*j_flux(i-1,1);
    end
    
    j_attempts = [j1 j2 j3];
    for k = 1:length(j_attempts)
    
       % Use old J
        if i == 1
            j_init = 0;
            j_R_tot = zeros(1,N_tot_active);
            
    
        else
    
          j_init   =  j_tot(i-1,1);
          j_R_tot = j_tot(i-1,1:N_tot_active);
        end
       
    
        %Set first guess
        j_flux(i,1) = j_attempts(k);
        
        %Set Reference Potential
        psi_e(i,1) = 0;
    
        %Caclulate eta for j_flux guess 
        eta(i,1) = (R_const*T(i,x_neg_start)/(alpha*F))*asinh((1/2)*(j_flux(i,1)./ex_current(i,1))) ;
    
        %Calculate psi_s
        psi_s(i,1) = eta(i,1) + psi_e(i,1) + OCP(i,1) + F*R_SEI_m(i,1)*j_init;
    
        % Set Boundary Conditions
        diff_psi_s_m(i,1) = -I_app(i)/(sigma_neg_eff) ;
    
        diff_psi_e_m(i,1) = 0;
        
        %%% Try different
        
            A = 0.5*eye(N_tot_active+1,N_tot_active);
        
            B = A + circshift(A,[1 0]);
            
            del_mid = B*del';
            
            del_mid = del_mid';
            
            diffx_C_e_log = diff(log(C_e_reduced),1,2) ;
    
            diffx_C_e_log = [0 diffx_C_e_log 0]./del_mid;
    
       
    
        %%% end different
        %diffx_C_e_log = diff(log(C_e_reduced),1,2)./del_tot ;
    
        %diffx_C_e_log = [0 diffx_C_e_log 0];
    
        for j = 1:N_neg-1
            
            % Calculate psi_s for j+1 element
           
    
            
            diff_psi_s_m(i,j+1) = diff_psi_s_m(i,j) + (F*As(j)*del(j)/sigma(x_ncc_end+j))*j_flux(i,j);
            
            
            psi_s(i,j+1) = psi_s(i,j) + del(j)*diff_psi_s_m(i,j+1);
            
            %Calculate psi_e for j+1 element
            kappa_k_plus_half = (kappa(i,j+1) + kappa(i,j))/2 ;
            if j == 1
                kappa_k_minus_half = kappa_k_plus_half;
            else
                kappa_k_minus_half = (kappa(i,j) + kappa(i,j-1))/2 ;
            end
    
            
            T_temp = mean([T(i,x_ncc_end+j+1) T(i,x_ncc_end+j)]);
          
            kd_k_plus_half = ((2*kappa_k_plus_half*R_const*T_temp)/(F))*(1-t_li) ;
    
            if j == 1
                kd_k_minus_half = kd_k_plus_half;
            else
                T_temp = mean([T(i,x_ncc_end+j-1) T(i,x_ncc_end+j)]);
                kd_k_minus_half = ((2*kappa_k_minus_half*R_const*T_temp)/(F))*(1-t_li) ;
            end
    
            diff_psi_e_m(i,j+1) = (1/kappa_k_plus_half)*( (kappa_k_minus_half)*diff_psi_e_m(i,j) + ( kd_k_plus_half*diffx_C_e_log(j+1) - kd_k_minus_half*diffx_C_e_log(j) ) - (F*As(j)*del(j)*j_flux(i,j)) );
    
            psi_e(i,j+1) = psi_e(i,j) + del(j)*diff_psi_e_m(i,j+1);
            
            % Determine j_fux for j+1 element
            eta(i,j+1) = psi_s(i,j+1) - psi_e(i,j+1) - OCP(i,j+1) - F*R_SEI_m(i,j+1).*j_R_tot(1,j+1);
    
            j_flux(i,j+1) = 2*ex_current(i,j+1).*sinh((alpha*F/(R_const*T(i,x_neg_start+j)))*eta(i,j+1)) ;
        end
        
        j_sol(k,:) = j_flux(i,:);
    end
    
    %%% err is caluclated wrt to j_ave_neg
    err =  1  - mean(j_sol(:,1:N_neg),2)./j_ave_neg;
    err;
    if sign(err(1)*err(3)) == sign(err(2)*err(3))
    
        %error("Doesnt Converge")
        beta1 = beta1 -0.1;
        beta2 = beta2 + 0.1;
        beta3 = (beta1+beta2)/2 ; 
    
        counter_boundary_finder = counter_boundary_finder +1;
    else
    
    Root_Finding_Algorithm;
    
    end
    
    end
    count_convergence(1,i) = counter_convergence;
    count_boundary_finding(1,i) = counter_boundary_finder;
    %[Min Ind] = min(abs(err));
    j_flux(i,:) = j_sol(3,:);
    
    %% Caclulate psi_e from current collector to first element of cathode
    
    for j = 1:N_neg+N_sep
    %Calculate psi_e for j+1 element
    kappa_k_plus_half = (kappa(i,j+1) + kappa(i,j))/2 ;
    
    if j == 1
       
        kappa_k_minus_half = kappa_k_plus_half;
    else
        kappa_k_minus_half = (kappa(i,j) + kappa(i,j-1))/2 ;
    end
            
    T_temp = mean([T(i,x_ncc_end+j+1) T(i,x_ncc_end+j)]);
          
    kd_k_plus_half = ((2*kappa_k_plus_half*R_const*T_temp)/(F))*(1-t_li) ;
    
    if j == 1
          kd_k_minus_half = kd_k_plus_half;
    else
          T_temp = mean([T(i,x_ncc_end+j-1) T(i,x_ncc_end+j)]);
          kd_k_minus_half = ((2*kappa_k_minus_half*R_const*T_temp)/(F))*(1-t_li) ;
    end
    
         
    
    diff_psi_e_m(i,j+1) = (1/kappa_k_plus_half)*((kappa_k_minus_half)*diff_psi_e_m(i,j) + (kd_k_plus_half*diffx_C_e_log(j+1) - kd_k_minus_half*diffx_C_e_log(j)) - (F*As(j)*del(j)*j_flux(i,j)));
    
    psi_e(i,j+1) = psi_e(i,j) + del(j)*diff_psi_e_m(i,j+1);
    
    end
    
    %% Iterative Solver For Positive Electrode
    j_ave_pos = -I_app(i)/(As_pos*F*L_pos) ;
    beta1 = 0.1;
    beta2 = 1.9;
    beta3 = (beta1 +beta2)/2;
    
    err = ones(3,1);
    counter_convergence = 0;
    counter_boundary_finder = 0;
    ITP_count = 0;
    loop = 0;
    while norm(err(3),2) > tol
    if rep_j < 2
        j1 = beta1*j_ave_pos;
        j2 = beta2*j_ave_pos;
        %beta3 = (beta1 +beta2)/2;
        j3 = beta3*j_ave_pos;
    
    else
        j1 = beta1*j_flux(i-1,N_neg+N_sep+1);
        j2 = beta2*j_flux(i-1,N_neg+N_sep+1);
        %beta3 = (beta1 +beta2)/2;
        j3 = beta3*j_flux(i-1,N_neg+N_sep+1);
    end
    
    
    
    j_attempts = [j1 j2 j3];
    
    for k = 1:length(j_attempts)
        %Set first guess
        j_flux(i,N_neg+N_sep+1) = j_attempts(k);
    
        %Calculate eta
    
        eta(i,N_neg+N_sep+1) = (R_const*T(i,x_pos_start)/(alpha*F))*asinh((1/2)*(j_flux(i,N_neg+N_sep+1)./ex_current(i,N_neg+N_sep+1))) ;
        
        %Calculate psi_s
        psi_s(i,N_neg+N_sep+1) = eta(i,N_neg+N_sep+1) + psi_e(i,N_neg+N_sep+1) + OCP(i,N_neg+N_sep+1);
    
        % Set Boundary Conditions
        diff_psi_s_m(i,N_neg+N_sep+1) = 0 ;
    
        
        
         for j = N_neg+N_sep+1:N_tot_active-1
            
            % Calculate psi_s for j+1 element
            diff_psi_s_m(i,j+1) = diff_psi_s_m(i,j) + (F*As(j)*del(j)/sigma(x_ncc_end+j))*j_flux(i,j);
            
            psi_s(i,j+1) = psi_s(i,j) + del(j)*diff_psi_s_m(i,j+1);
            
            %Calculate psi_e for j+1 element
            kappa_k_plus_half = (kappa(i,j+1) + kappa(i,j))/2 ;
            if j == 1
                kappa_k_minus_half = kappa_k_plus_half;
            else
                kappa_k_minus_half = (kappa(i,j) + kappa(i,j-1))/2 ;
            end
            
            kd_k_plus_half = ((2*kappa_k_plus_half*R_const*mean([T(i,x_ncc_end+j+1) T(i,x_ncc_end+j)]))/(F))*(1-t_li) ;
    
            kd_k_minus_half = ((2*kappa_k_minus_half*R_const*mean([T(i,x_ncc_end+j-1) T(i,x_ncc_end+j)]))/(F))*(1-t_li) ;
    
         
    
            diff_psi_e_m(i,j+1) = (1/kappa_k_plus_half)*((kappa_k_minus_half)*diff_psi_e_m(i,j) + (kd_k_plus_half*diffx_C_e_log(j+1) - kd_k_minus_half*diffx_C_e_log(j)) - (F*As(j)*del(j)*j_flux(i,j)));
    
            psi_e(i,j+1) = psi_e(i,j) + del(j)*diff_psi_e_m(i,j+1); % changed
            
            % Determine j_fux for j+1 element
            eta(i,j+1) = psi_s(i,j+1) - psi_e(i,j+1) - OCP(i,j+1);
    
            j_flux(i,j+1) = 2*ex_current(i,j+1).*sinh((alpha*F/(R_const*T(i,x_ncc_end+j+1)))*eta(i,j+1)) ;
        end
        
        j_sol(k,:) = j_flux(i,:);
        
    end
    err =  1  - mean(j_sol(:,N_neg+N_sep+1:N_tot_active),2)./j_ave_pos;
     
    if sign(err(1)*err(3)) == sign(err(2)*err(3))
    
        %error("Doesnt Converge")
        %beta1 = beta1/beta1;
        beta1 = beta1 - 3;%Changed from 0.1
        beta2 = beta2 + 3;
        beta3 = (beta1+beta2)/2 ; 
    
        counter_boundary_finder = counter_boundary_finder +1;
     else
    
     Root_Finding_Algorithm;
    
     end
    
    end
    count_convergence(2,i) = counter_convergence;
    count_boundary_finding(2,i) = counter_boundary_finder;
    %[Min Ind] = min(abs(err));
    j_flux(i,:) = j_sol(3,:);
    j_tot(i,:) = j_flux(i,:);
    
    %Get Last diff_psi_e and diff_psi_s
    j = N_tot_active;
    diff_psi_s_m(i,j+1) = diff_psi_s_m(i,j) + (F*As(j)*del(j)/sigma(x_ncc_end+j))*j_flux(i,j);
    
    kd_k_minus_half = ((2*kappa_k_minus_half*R_const*mean([T(i,x_ncc_end+j-1) T(i,x_ncc_end+j)]))/(F))*(1-t_li) ;
    
    kd_k_plus_half = kd_k_minus_half ;
    
         
    diff_psi_e_m(i,j+1) = (1/kappa_k_plus_half)*((kappa_k_minus_half)*diff_psi_e_m(i,j) + (kd_k_plus_half*diffx_C_e_log(j+1) - kd_k_minus_half*diffx_C_e_log(j)) - (F*As(j)*del(j)*j_flux(i,j)));
    
    %%% End of iterative solver

end

if I_app(i) == 0
    
    j_flux(i,:) = zeros(1,N_tot_active);

    j_tot(i,:) = j_flux(i,:);

    i_e(i,:) = zeros(1,N_tot_active+1);

    i_s(i,:) = zeros(1,N_tot_active+1);

    A = eye([N_tot_active+1,N_tot_active]);
    A = A + circshift(A,[1 0]);
    B = A*0.5;
    B(1,1) = 1;
    B(end,end) = 1;

     
            
    del_mid = B*del';
            
    del_mid = del_mid';
            
    diffx_C_e_log = diff(log(C_e_reduced),1,2) ;
    
    diffx_C_e_log = [0 diffx_C_e_log 0]./del_mid;

    


    diff_psi_e_m(i,:) = -i_e(i,:) + (2*R_const/F)*(((B)*(T(i,:).*kappa(i,:))')').*diffx_C_e_log  ;

    diff_psi_s_m(i,:) = zeros(1,N_tot_active+1);
    

end