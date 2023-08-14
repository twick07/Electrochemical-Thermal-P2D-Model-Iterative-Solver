%%% Full Order P2D Model Simulation

clc;
clear;
%}
close all;




COMSOL_Driven = 0;







%%% Get Model Parameters
%COMSOL_Driven = 0;
%%% param = 1 for temperature model for (LMO/LiC6) from liu et al 
%%% param = 2 for temperature model for NMC111/C6 Cai et al
%%% param = 3 for temperature model with LiCoO2/C6 from tochio et al 

param = 1;

%%% Load Battery Material Parameters


%Load Model Battery Parameters
Battery_Parameters; 



%%% Set Simulation Settings
%%% Time Integration 
%%% Time_int = 1 for eurler only
%%% Time_int = 2 for verlet integration
time_int = 2;


%%% Set number of elements
N_tot = 30;
del_tot = L_tot/N_tot ;


sim_time = (116)*(60); %Total sim time in seconds %59.25 for solver 1, order 3, %8.5 mins for solver 2 ,order 3, 86.4

%set up variable time_vector
time_variable = 0 ;% Parameter required for simulation code
time_stop = 3600; % Stop time

%time_step = 1.0*time_step_initial;
if param == 1
    time_step = 1.3; % 1.3 works for thermal = 9  mass_con = 5, 2.6 works for thermal = 0
end


if param == 2
    time_step = 1.3;% 1.3 works for thermal = 9  mass_con = 5, 2.6 works for thermal = 0
end

if param ==3
    time_step = 0.45; %Works for thermal ==0 and thermal == 9
    
end

time_max = round(time_stop/time_step);
%time_max = 1;
No_of_cycles = 1;

time_max = 5*2*No_of_cycles*time_max;
time_max = round(time_max);



C_rate = 1.0;
I_const =  C_rate*I_1c; 



I_app =I_const*ones(1,time_max);


%%% Set Simulation Limits
Cut_off_capacity = 99.9; % Maximum Capacity in electrode
if param == 2
    V_cut_off_low = 3.1; % Lowest voltage in simulation 2.9 for param 1 and 3, 3.1 for param 2
else
    V_cut_off_low = 2.9; % Lowest voltage in simulation 2.9 for param 1 and 3, 3.1 for param 2
end
V_cut_off_high = 4.1; % Highest voltage in simulation

%%% Thermal Charchterisation
%%% thermal = 0 for isothermal simulation
%%% thermal = 3 for non-isothermal simulation
%%% thermal = 7 for lumped sum thermal estimation approach
%%% thermal = 9 for ode15s implementation of temperature
thermal = 9;

%%% Error Condition
%%% err_con = 1 for using j_ave
%%% err_con = 2 for using currents

err_con = 1;


%%% Heat of Mixing
%%% heat_of_mixing = 0 for no heat of mixing
%%% heat_of_mixing = 1 for including heat of mixing

heat_of_mixing = 0;

%%% Ageing Model Choice

%%% SEI growth model
%%% SEI = 0, no SEI growth
%%% SEI = 1, For Kinetics limited SEI growth model
%%% SEI = 2, For Diffusion limited SEI growth model
%%% n_0 is the number of points to find the root for j_flux)i,j)
SEI = 1;



%%% Lithium Plating
%%% li_pl = 0, no plating
%%% li_pl = 1, Aurora et al Plating Model
li_pl = 0;


%%% Discretize model

%%% set element size in x
%Mesh = 1 for proportional Discretisation
%Mesh = 2 for COMSOL Mesh Import
%Mesh = 3 for seperator mesh being fixed
Mesh = 1;

N_shell = 10;

del_Rs_neg = Rs_neg/N_shell;
del_Rs_pos = Rs_pos/N_shell;
%%% set element size in r
if Mesh == 1
    if current_collector == 0
        
        %%% Set Model Discretisation
        %%% Linear Discretisation
        
    
        del_tot = L_tot/N_tot;
    
        N_neg = round(L_neg/del_tot);                                                                                                                                                                                                                              0;
        del_neg = L_neg/N_neg ;
    
        N_sep = round(L_sep/del_tot);
        del_sep = L_sep/N_sep;
    
        N_pos = round(L_pos/del_tot);
        del_pos = L_pos/N_pos ;
    
    
        N_tot_active = N_neg + N_sep +N_pos;
        N_tot = N_tot_active;
    
        
    
       
    
        % Define Specific Positions in Battery Model
        x_ncc_end = 0;
        x_0 = 1; %First Element in Negative Electrode
        x_neg_start = x_0; %First Element in Negative Electrode
        x_neg_end = N_neg; % Last element of Negative Electrode
        x_sep_start = x_neg_end + 1; %First Element in Seperator
        x_sep_end = N_neg + N_sep; % Last Element of Seperator
        x_pos_start = x_sep_end + 1 ; % First Element of Positive Electrode  
        x_pos_end = N_neg + N_sep + N_pos; % Last Element of Positive Electrode
        x_end = x_pos_end; % Last element of Battery Model
    
    
    end
    
    if current_collector == 1
    
        
    
        %%% Set Model Discretisation
        %%% Linear Discretisation
    
        del_tot = L_tot/N_tot;
    
        N_neg = round(L_neg/del_tot);                                                                                                                                                                                                                              0;
        del_neg = L_neg/N_neg ;
    
        N_sep = round(L_sep/del_tot);
        del_sep = L_sep/N_sep;
    
        N_pos = round(L_pos/del_tot);
        del_pos = L_pos/N_pos ;
    
        N_cc = round(L_ncc/del_tot);
        del_cc = L_ncc/N_cc ;
    
        N_tot_active = N_neg + N_sep +N_pos;
        N_tot = N_cc + N_tot_active + N_cc;
    
        
    
        % Define Specific Positions in Battery Model
    
        x_0 = 1; %First Element in Neg cc
        x_ncc_end = N_cc; % Last element of Neg cc
        x_neg_start = x_ncc_end + 1; %First Element in Negative Electrode
        x_neg_end = N_cc +N_neg; % Last element of Negative Electrode
        x_sep_start = x_neg_end + 1; %First Element in Seperator
        x_sep_end = N_cc + N_neg + N_sep; % Last Element of Seperator
        x_pos_start = N_cc + N_neg + N_sep + 1 ; % First Element of Positive Electrode  
        x_pos_end = N_cc + N_neg + N_sep + N_pos; % Last Element of Positive Electrode
        x_pcc_start = x_pos_end+1;%First Element in Pos cc
        x_end = N_tot; % Last element of Pos cc




end

end

%Keep Mesh In Seperator Fixed
if Mesh ==3


    N_sep = round(N_tot/10);


    del_sep = L_sep/N_sep ; 
    if current_collector == 0

        N_tot_left = N_tot - N_sep ; 

        N_neg = round((L_neg/(L_neg+L_pos))*N_tot_left);

        del_neg = L_neg/N_neg;

        N_pos = N_tot - N_neg - N_sep;

        del_pos = L_pos/N_pos;

        N_tot_active = N_tot;
        
        % Define Specific Positions in Battery Model
        x_ncc_end = 0;
        x_0 = 1; %First Element in Negative Electrode
        x_neg_start = x_0; %First Element in Negative Electrode
        x_neg_end = N_neg; % Last element of Negative Electrode
        x_sep_start = x_neg_end + 1; %First Element in Seperator
        x_sep_end = N_neg + N_sep; % Last Element of Seperator
        x_pos_start = x_sep_end + 1 ; % First Element of Positive Electrode  
        x_pos_end = N_neg + N_sep + N_pos; % Last Element of Positive Electrode
        x_end = x_pos_end; % Last element of Battery Model


    end


end

%%% Solver
%%% Solver 1 for Bisection Method
%%% Solver 2 for False Point Method with bisection method at start
%%% Solver 3 for Illinois Method
%%% Solver 4 for ITP Method
%%% Solver 5 for eSPM Model
%%% Solver 6 for SPM Model
%%% Solver 7 for Anderson–Björck Method 
solver = 2;

tol = 1e-6; % within 0.001 accuracy
tol_change = 10;
% Tolerance for deciding starting point of secant method

%%% Order of Diffusion
%%% Order 1 for reduced order
%%% Order 3 for full order diffusion based solver
%%% Order 4 for linear algebra implementation of full order diffusion
%%% Order 6 for ode15s implementation
order = 4 ;



%%% Use Previous J
%%% use_prev_j = 0 means previous J is not used
%%% use_prev_j = 1 means previous J is used
use_prev_j = 0;
rep_j = 0;


%%% Parameter Corrections according to different papers 
%%% correction = 3 for LMO Corrections
correction = param;

if correction == 1

    D_sep_e_eff = D_elec*(eps_sep_elec^brugg_sep); % Effective Diffusion Coefficient in electrolyte
    D_neg_e_eff = D_elec*(eps_neg_elec^brugg_neg); %Effective Diffusivity in electrode's electrolyte
    D_pos_e_eff = D_elec*(eps_pos_elec^brugg_pos); %Effective Diffusivity in electrode's electrolyte
    sigma_neg_eff = sigma_neg*(eps_neg^brugg_neg); % Effective Electrical Conductivity in negative electrode 
    sigma_pos_eff = sigma_pos*(eps_pos^brugg_pos); % Effective Electrical Conductivity in positive electrode  
    kappa_correction = 1; % Add Correction for Electrolyte Conductivity

end



if correction == 2

    sigma_neg_eff = sigma_neg*(eps_neg); % Effective Electrical Conductivity in negative electrode 
    sigma_pos_eff = sigma_pos*(eps_pos); % Effective Electrical Conductivity in positive electrode  
    D_sep_e_eff = D_elec*(eps_sep_elec^brugg_sep); % Effective Diffusion Coefficient in electrolyte
    D_neg_e_eff = D_elec*(eps_neg_elec^brugg_neg); %Effective Diffusivity in electrode's electrolyte
    D_pos_e_eff = D_elec*(eps_pos_elec^brugg_pos); %Effective Diffusivity in electrode's electrolyte

end


if correction == 3

    
    sigma_neg_eff = sigma_neg*(eps_neg); % Effective Electrical Conductivity in negative electrode 
    sigma_pos_eff = sigma_pos*(eps_pos); % Effective Electrical Conductivity in positive electrode  
    kappa_correction = 1; % Add Correction for Electrolyte Conductivity

end

%%% C_e mass conservation at boundary

%%% mass_con = 0 for standard FVM 
%%% mass_con = 1 for FVM set up at boundary
%%% mass_con = 2 for FDM implementation
%%% mass_con = 3 for FDM with FVM at boundary
%%% mass_con = 4 for standard FDM with ode15s 
%%% mass_con = 5 for standard FVM with ode15s

mass_con = 5;

if mass_con == 0
    C_e = C_e_initial*ones(time_max+1,N_tot_active);
    difft_C_e_m = zeros(time_max,N_tot_active);
end

if mass_con == 1
    C_e = C_e_initial*ones(time_max+1,N_tot_active+2);
    difft_C_e_m = zeros(time_max,N_tot_active+2);
end

if mass_con == 2
    C_e = C_e_initial*ones(time_max+1,N_tot_active);
    difft_C_e_m = zeros(time_max,N_tot_active);
    C_e_ode15s = C_e_initial*ones(time_max+1,N_tot_active);
end

if mass_con == 3
    C_e = C_e_initial*ones(time_max+1,N_tot_active+2);
    difft_C_e_m = zeros(time_max,N_tot_active+2);
end

if mass_con == 4
    C_e = C_e_initial*ones(time_max+1,N_tot_active);
    difft_C_e_m = zeros(time_max,N_tot_active);
    C_e_ode15s = C_e_initial*ones(time_max+1,N_tot_active);
end

if mass_con == 5
    C_e = C_e_initial*ones(time_max+1,N_tot_active);
    difft_C_e_m = zeros(time_max,N_tot_active);
end
%Create Data Containers
%time_vector = zeros(1,time_max);
edge1 = zeros(1,time_max);
edge2 = zeros(1,time_max);
j_flux = zeros(time_max,N_tot_active);
j_tot = zeros(time_max,N_tot_active);
q_x = zeros(time_max+1,N_tot_active);
i_e = zeros(time_max,N_tot_active);
i_s = zeros(time_max,N_tot_active);
C_se = [C_neg_initial*ones(time_max+1,N_neg) zeros(time_max+1,N_sep) C_pos_initial*ones(time_max+1,N_pos)]; 
C_s = [C_neg_initial*ones(time_max+1,N_neg) zeros(time_max+1,N_sep) C_pos_initial*ones(time_max+1,N_pos)]; 
C_max = [C_max_neg*ones(1,N_neg) zeros(1,N_sep) C_max_pos*ones(1,N_pos)];
psi_e = zeros(time_max,N_tot_active);
psi_s = zeros(time_max,N_tot_active);
eta = zeros(time_max,N_tot_active);
eta_tot = zeros(time_max,N_tot_active);
ex_current = zeros(time_max,N_tot_active);
OCP = zeros(time_max,N_tot_active);
stoic = zeros(time_max,N_tot_active);
kappa = zeros(time_max,N_tot_active);
kd = zeros(time_max,N_tot_active);
D_e_eff = [D_neg_e_eff*ones(1,N_neg) D_sep_e_eff*ones(1,N_sep) D_pos_e_eff*ones(1,N_pos)];
D_s_eff = [D_neg_eff_s*ones(1,N_neg) zeros(1,N_sep) D_pos_eff_s*ones(1,N_pos)];
eps_elec = [eps_neg_elec*ones(1,N_neg) eps_sep_elec*ones(1,N_sep) eps_pos_elec*ones(1,N_pos)];
eps_sol = [eps_neg*ones(1,N_neg) zeros(1,N_sep) eps_pos*ones(1,N_pos)];
As = [As_neg*ones(1,N_neg) zeros(1,N_sep) As_pos*ones(1,N_pos)];
Rs = [Rs_neg*ones(1,N_neg) zeros(1,N_sep) Rs_pos*ones(1,N_pos)];
R_SEI_m = R_SEI*ones(time_max+1,N_tot_active);
j_SEI = zeros(time_max,N_tot_active);
eta_SEI = zeros(time_max,N_tot_active);
cycle_count = zeros(1,15);
cycle_point = 1;
cycle_count(cycle_point) = cycle_point;

if current_collector == 0
    sigma = [sigma_neg_eff*ones(1,N_neg) zeros(1,N_sep) sigma_pos_eff*ones(1,N_pos)];

end

if current_collector == 1
    sigma = [sigma_ncc*ones(1,N_cc) sigma_neg_eff*ones(1,N_neg) zeros(1,N_sep) sigma_pos_eff*ones(1,N_pos) sigma_pcc*ones(1,N_cc)];
end
diff_psi_s_m = zeros(time_max,N_tot_active+1);
diff_psi_e_m = zeros(time_max,N_tot_active+1);
del = [del_neg*ones(1,N_neg) del_sep*ones(1,N_sep) del_pos*ones(1,N_pos)];
count_convergence = zeros(2,time_max);
count_boundary_finding = zeros(2,time_max);
brugg = [brugg_neg*ones(1,N_neg) brugg_sep*ones(1,N_sep) brugg_pos*ones(1,N_pos)];
V_cell = zeros(1,time_max);
K_0_init = [K_0_neg*ones(1,N_neg) zeros(1,N_sep) K_0_pos*ones(1,N_pos)];


%Radial Diffusion

C_s_full = zeros(time_max,N_shell,N_tot_active);
C_s_full(1,:,1:N_neg) = C_neg_initial*ones(1,N_shell,N_neg);
C_s_full(1,:,N_neg+N_sep+1:N_tot_active) = C_pos_initial*ones(1,N_shell,N_pos);
C_s_4 = [C_neg_initial*ones(N_shell,N_neg) zeros(N_shell,N_sep) C_pos_initial*ones(N_shell,N_pos)];
C_s_ave = [C_neg_initial*ones(time_max+1,N_neg) zeros(time_max+1,N_sep) C_pos_initial*ones(time_max+1,N_pos)];
dR_pos = Rs_pos/N_shell; % width of each "shell"
dR_neg = Rs_neg/N_shell; % width of each "shell"
Sa_pos = 4*pi*((Rs_pos*(1:N_shell)/N_shell).^2); % outer surface area of each shell
Sa_neg = 4*pi*((Rs_neg*(1:N_shell)/N_shell).^2); % outer surface area of each shell
dV_pos = (4/3)*pi*( ((Rs_pos*(1:N_shell)/N_shell).^3) - ((Rs_pos*(0:N_shell-1)/N_shell).^3) ); % vol. of ea. shell
dV_neg = (4/3)*pi*( ((Rs_neg*(1:N_shell)/N_shell).^3) - ((Rs_neg*(0:N_shell-1)/N_shell).^3) ); % vol. of ea. shell

N_flux_neg = zeros(time_max,N_shell-1,N_neg);
N_flux_pos = zeros(time_max,N_shell-1,N_pos);

%Radial Diffusion for SPM
C_s_SPM_neg = zeros(time_max,N_shell);
C_s_SPM_neg(1,:) = C_neg_initial*ones(1,N_shell);
C_s_SPM_pos = zeros(time_max,N_shell);
C_s_SPM_pos(1,:) = C_pos_initial*ones(1,N_shell);
C_se_SPM_neg = zeros(1,time_max);
C_se_SPM_neg(1,1) = C_neg_initial;
C_se_SPM_pos = zeros(1,time_max);
C_se_SPM_pos(1,1) = C_pos_initial;

% Temperature Specific Modelling
h_ncc =   h_ext; % Heat Exchange Parameter [W/(m^2 K)]
h_pcc = h_ext; % Heat Exchange Parameter [W/(m^2 K)]
T = T_amb*ones(time_max+1,N_tot);
diff_T_x = zeros(time_max,N_tot+1);
diff_T_t = zeros(time_max,N_tot);
if current_collector == 0
    C_p = [C_p_neg*ones(1,N_neg) C_p_sep*ones(1,N_sep) C_p_pos*ones(1,N_pos)];
    row = [row_neg*ones(1,N_neg) row_sep*ones(1,N_sep) row_pos*ones(1,N_pos)];
    lambda = [lambda_neg*ones(1,N_neg) lambda_sep*ones(1,N_sep) lambda_pos*ones(1,N_pos)];
end


if current_collector == 1
    C_p = [C_p_ncc*ones(1,N_cc) C_p_neg*ones(1,N_neg) C_p_sep*ones(1,N_sep) C_p_pos*ones(1,N_pos) C_p_pcc*ones(1,N_cc)];
    row = [row_ncc*ones(1,N_cc) row_neg*ones(1,N_neg) row_sep*ones(1,N_sep) row_pos*ones(1,N_pos) row_pcc*ones(1,N_cc)];
    lambda = [lambda_ncc*ones(1,N_cc) lambda_neg*ones(1,N_neg) lambda_sep*ones(1,N_sep) lambda_pos*ones(1,N_pos) lambda_pcc*ones(1,N_cc)];
end

diff_OCP_T = zeros(1,N_tot_active);
Q_rxn = zeros(time_max,N_tot_active);
Q_rev = zeros(time_max,N_tot_active);
Q_ohm = zeros(time_max,N_tot_active);
Q_gen = zeros(time_max,N_tot);

if heat_of_mixing == 1
    Q_mix = zeros(time_max,N_tot_active);

end
diff_T_x_2_comp = zeros(time_max,N_tot);
E_a_K = [E_a_K_0_neg*ones(1,N_neg) zeros(1,N_sep) E_a_K_0_pos*ones(1,N_pos)];
E_a_D = [E_a_D_s_neg*ones(1,N_neg) zeros(1,N_sep) E_a_D_s_pos*ones(1,N_pos)];
mean_row = ((N_neg*row_neg) + (N_sep*row_sep) + (N_pos*row_pos))/N_tot_active;
mean_C_p = ((N_neg*C_p_neg) + (N_sep*C_p_sep) + (N_pos*C_p_pos))/N_tot_active;
T_ave = zeros(1,time_max);
T_ave_2 = zeros(1,time_max);
T_ave_ncc = zeros(1,time_max);
T_ave_neg = zeros(1,time_max);
T_ave_sep = zeros(1,time_max);
T_ave_pos = zeros(1,time_max);
T_ave_pcc = zeros(1,time_max);




%%% SEI Modelling
if SEI >0 
    
    L_SEI_m = L_SEI*ones(time_max+1,N_tot_active);
    kappa_SEI_v = kappa_SEI*ones(1,N_tot_active);
    M_SEI_v = M_SEI*ones(1,N_tot_active);
    row_SEI_v = rho_SEI*ones(1,N_tot_active);
    j_flux_sol = zeros(3,N_tot_active);
    c_SEI_m = c_SEI_0*ones(time_max+1,N_tot_active);
    c_sol_m = epss_SEI*C_sol_0*ones(time_max+1,N_tot_active);
    i_0_SEI_v = i_0_SEI*ones(time_max+1,N_tot_active);
    if SEI == 2
        i_0_SEI = F*k_SEI*C_sol_0;
        i_0_SEI_v = i_0_SEI*ones(time_max+1,N_tot_active);
        k_SEI_v = k_SEI*ones(time_max+1,N_tot_active);
        D_SEI_v = D_SEI.*ones(1,N_tot_active);
    end

   
end
%%% Plating Modelling
if li_pl >0
    j_li_pl = zeros(time_max,N_tot_active);
    eta_li_pl = zeros(time_max,N_tot_active);
    R_li_m = R_Li*ones(time_max+1,N_tot_active);
    L_li_m = L_li*ones(time_max+1,N_tot_active);
    c_SEI_m = c_Li_0*ones(time_max+1,N_tot_active);
    kappa_Li_v = kappa_li*ones(1,N_tot_active);
    M_SEI_v = M_li*ones(1,N_tot_active);
    row_Li_v = rho_li*ones(1,N_tot_active);
end





%ITP Solver Variables
if solver ==4


k1 = 20.0;
psi = 0.5*(1+sqrt(5));
k2 = 1.3;

eta_zero = 1.0;

end

if solver ~= 5
    i_e = zeros(time_max,N_tot_active+1);
    i_s = zeros(time_max,N_tot_active+1);

end




%%% Set time_step variation
%Variation = 1 - full solver
%Variation = 2 - variable time_step 
variation = 1;

if variation == 2
    
    time_step_fine = time_step;
    time_step_norm = 1.0;
    
    thresh_neg = 35;
    thresh_pos = 35;
end
time = 0;

if (solver < 5 && SEI == 0) 
    tic;
    Main_loop;
    toc;
end

if (solver < 5 && SEI > 0) 
    tic;
    Main_loop;
    toc;
end

if (solver == 7 ) 
    tic;
    Simulation_repeat_j;
    toc;
end
if solver == 5
    tic;
    SPMe_Simulation;
    toc;

end

if solver == 6
    tic;
    SPM_Simulation;
    toc;

end

time_max = i;

if SEI == 0
    Debug_Plotter;
end

if SEI >0
    Debug_Plotter_SEI;
end
