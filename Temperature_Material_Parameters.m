%%% Script to choose conductivity model equation
%% Taken from Liu et al
if param == 1



    kappa(i,:) = 1.0793e-2 + (6.7461e-4).*(C_e_reduced) - (5.2245e-7).*(C_e_reduced.^2) + (1.3605e-10).*(C_e_reduced.^3) - (1.1724e-14).*(C_e_reduced.^4);

    if kappa_correction == 1
        

        kappa(i,:) = (eps_elec.^brugg).*kappa(i,:);

    end

    kd(i,:) = (2*(1-t_li)*R_const/F).*kappa(i, :);

    %K_0 = K_0_init;

    K_0 = K_0_init.*exp( -(E_a_K/R_const).*( (1./T(i,x_neg_start:x_pos_end)) - (1/T_amb) ) ); %Changed


    %D_s = D_s_eff;

    D_s = D_s_eff.*exp( -(E_a_D/R_const).*( (1./T(i,x_neg_start:x_pos_end)) - (1/T_amb) ) ); %Changed


    %D_e = D_e_eff;
    
    D_e =  D_e_eff.*exp( -(E_a_D_e/R_const).*( (1./T(i,x_neg_start:x_pos_end)) - (1/T_amb) ) ); %Changed


    if SEI == 2
        k_SEI_v(i,:) = k_SEI.*exp( -(E_a_k_SEI/R_const).*( (1./T(i,x_neg_start:x_pos_end)) - (1/T_amb) ) );
    end




end





%% All Equations are from Hosseinzadeh et al, Liebing et al (NMC111/C6 Cell)
if param == 2

%Electrolyte Conductivity (kappa) & kd

%kappa(i, :) = (15.8e-2)*(eps_elec.^brugg).*C_e_reduced.*exp(-13472.*C_e_reduced.^(1.4)).*exp( -(E_a_kappa/R_const).*( (1./T(i,x_neg_start:x_pos_end)) - (1/T_amb) ) );

kappa(i, :) = (1.0793e-2 + (6.7461e-4).*(C_e_reduced) - (5.2245e-7).*(C_e_reduced.^2) + (1.3605e-10).*(C_e_reduced.^3) - (1.1724e-14).*(C_e_reduced.^4)).*exp( -(E_a_kappa/R_const).*( (1./T(i,x_neg_start:x_pos_end)) - (1/T_amb) ) );

kd(i,:) = (2*(1-t_li)*R_const/F).*kappa(i, :);

%Reaction Rate Constant (K_0)
K_0 = K_0_init.*exp( -(E_a_K/R_const).*( (1./T(i,x_neg_start:x_pos_end)) - (1/T_amb) ) );

% Diffusion Coeffisients (D_s & D_e)
D_s = D_s_eff.*exp( -(E_a_D/R_const).*( (1./T(i,x_neg_start:x_pos_end)) - (1/T_amb) ) );

D_e =  D_e_eff.*exp( -(E_a_D_e/R_const).*( (1./T(i,x_neg_start:x_pos_end)) - (1/T_amb) ) );





end


%% All Equations are from Torchio et al and Northrop et al
if param ==3

%Electrolyte Conductivity (kappa) & kd
kappa(i, :) = (1e-4)*(eps_elec.^brugg).*C_e_reduced.*((-10.5 + (0.668e-3)*C_e_reduced + (0.494e-6)*(C_e_reduced.^2) + ...
                 (0.074-(1.78e-5)*C_e_reduced - (8.86e-10)*(C_e_reduced.^2)).*T(i,x_neg_start:x_pos_end) + ...
                  (-6.96e-5 + (2.8e-8)*C_e_reduced).*(T(i,x_neg_start:x_pos_end).^2)).^2);

kd(i,:) = (2*(1-t_li)*R_const/F).*kappa(i, :);

%Reaction Rate Constant (K_0)
K_0 = K_0_init.*exp( -(E_a_K/R_const).*( (1./T(i,x_neg_start:x_pos_end)) - (1/T_amb) ) );

% Diffusion Coeffisients (D_s & D_e)
D_s = D_s_eff.*exp( -(E_a_D/R_const).*( (1./T(i,x_neg_start:x_pos_end)) - (1/T_amb) ) );

D_e =  (eps_elec.^brugg)*(1e-4).*(10.^(-4.43 -( (54)./( T(i,x_neg_start:x_pos_end) -229 -(5e-3).*C_e_reduced ) ) ...
          -(0.22e-3).*C_e_reduced ));



end

