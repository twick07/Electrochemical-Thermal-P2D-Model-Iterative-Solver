%% Script To Choose Model Parameters

%% Taken from Liu et al with no current collector
if param == 1
    

    
    current_collector = 0; % model has current collectors
   

    %Constants
    R_const = 8.314; % Ideal gas constant [J/[K*mol]]
    F = 96487; % Faraday Constant [Coulombs / mol]
    T_amb = 298.15; %Temperature in [K]
    
    %Bruggman constants
    %{
    brugg = 1.7;
    brugg_neg = 1.03; %Bruggeman coefficient for tortuosity in negative electrode
    brugg_sep = brugg; %Bruggeman coefficient for tortuosity in separator
    brugg_pos = brugg; %Bruggeman coefficient for tortuosity in positive electrode
    %}
    
    brugg = 1.5;
    brugg_neg = brugg; %Bruggeman coefficient for tortuosity in negative electrode
    brugg_sep = brugg; %Bruggeman coefficient for tortuosity in separator
    brugg_pos = brugg; %Bruggeman coefficient for tortuosity in positive electrode
    %Cell Dimensions
    L_neg = 100e-6; % Length of negative electrode [m] %added
    L_sep = 52e-6; % Length of separator [m] 
    L_pos = 183e-6; % Length of positive electrode [m]
    L_tot = L_neg + L_sep + L_pos;% [m]


   



    %Electrolyte/Seperator Parameters (1M LiPF6 EC:DMC (2:1))
    C_e_initial = 2000; %Initial Electrolyte Lithium Concentration [mol/[m^3]] 
    D_elec = 7.5e-11; % Diffusion Coefficient for Lithium in electrolyte [[m^2]/s] %Defined in model params
    t_li = 0.363; % Lithium Ion Transference Number [unitless]
    eps_sep_elec = 1.0; %Electrolyte Volume Fraction
    C_p_sep = 700; % Specific Heat Capacity [J/(Kg*K)]
    row_sep = 1.2e3; % Density of Seperator [kg/[m^3]]
    lambda_sep = 1; % Thermal Conductivity [W/(m*K)] 
    E_a_D_e = 35000; % Activation Energy for Conductivity  [J/mol]


    %Negative Electrode
    %soc_initial_neg_charge = 0.02;
    %soc_initial_neg = soc_initial_neg_charge;
    soc_initial_neg = 0.58;
    eps_neg_elec = 0.357; % Electrolyte Volume Fraction (porosity)
    eps_neg = 0.471; % Solid Volume Fraction
    C_max_neg = 26390 ; % Maximum Solid Lithium Concentration [mol/[m^3]]
    C_neg_initial = soc_initial_neg*C_max_neg ;% Initial Lithium Concentration [mol/[m^3]]
    Rs_neg = 12.5e-6; % Electrode Particle Size [m] %changed
    D_neg_s = 3.9e-14; % Solid Diffusion Coefficient [m^2/s]
    D_neg_eff_s = D_neg_s;
    sigma_neg = 100 ; % Electrical Conductivity [S/m]
    alpha = 0.5; % Transfer Coefficient
    row_neg = 2.5e3; % Density of Electrode [kg/[m^3]]
    K_0_neg = 2e-10; % Negative Rate Constant %changed
    As_neg = 3*(eps_neg/Rs_neg) ; % Surface area to volume ration [m-1]
    C_p_neg = 700; % Specific Heat Capacity [J/(Kg*K)]
    E_a_D_s_neg = 4e3; % Activation Energy for Temperature Dependent Solid Phase Diffusion [J/mol]
    E_a_K_0_neg = 30e3; % Activation Energy for Temperature Reaction Constant  [J/mol]
    lambda_neg = 5; % Thermal Conductivity [J/(m*K)]
    diff_U_T_neg =  0 ;%"Entropy Change For Electrode Reaction,Negative Electrode [V/K]"

    %Positive Electrode
    %soc_initial_pos_charge = 0.65;
    %soc_initial_pos = soc_initial_pos_charge;
    soc_initial_pos = 0.19;
    eps_pos_elec = 0.444; %% Electrolyte Volume Fraction 
    eps_pos = 0.297; % Solid Volume Fraction
    C_max_pos = 22860 ; % Maximum Solid Lithium Concentration [mol/[m^3]]
    C_pos_initial = soc_initial_pos*C_max_pos ;% Initial Lithium Concentration [mol/[m^3]]
    Rs_pos = 8.5e-6; % Electrode Particle Size [m] %changed
    D_pos_s = 1e-13; % Solid Diffusion Coefficient [m^2/s]
    D_pos_eff_s = D_pos_s;
    sigma_pos = 3.8 ; % Electrical Conductivity [S/m] From COMSOL
    alpha = 0.5; % Transfer Coefficient
    row_pos = 1.5e3; % Density of Electrode [kg/[m^3]]
    K_0_pos = 2e-10; % Negative Rate Constant %changed
    As_pos = 3*(eps_pos/Rs_pos) ; % Surface area to volume ratio [m-1]
    C_p_pos = 700; % Specific Heat Capacity [J/(Kg*K)]
    E_a_D_s_pos = 20e3; % Activation Energy for Temperature Dependent Solid Phase Diffusion [J/mol]
    E_a_K_0_pos = 30e3; % Activation Energy for Temperature Reaction Constant  [J/mol]
    lambda_pos = 5; % Thermal Conductivity [J/(m*K)]
    diff_U_T_pos = 0; %"Entropy Change For Electrode Reaction,Positive Electrode [V/K]"
        
    %SEI Parameters
    %c_SEI_0 = 4541; % "EC soloution concentration [mol/m^3]"
    D_SEI = 3.5e-20; %1e-13 ;% "EC diffusion coefficient [m^2/s]" -18 works
    E_a_k_SEI  = 200e3 ;% "Activation Energy for reaction coefficient, SEI [J/mol]"
    f_elec = 1; % "Electrolyte Activity coefficient"
    L_SEI =  1e-9;%1e-9 ;%"Initial SEI Thickness [m]"
    M_SEI =  0.162 ;%"SEI Molar Mass [kg/mol]"
    %R_SEI =  1e-4 ;%"SEI initial resistance [ohm*m^2]" 
    alpha_SEI = 0.5 ;%"Charge transfer coefficient, SEI"
    epss_SEI = 0.03;%0.05 ;%"Solid phase volume fraction SEI"
    rho_SEI = 1.69e3 ;% "Separator density [kg/m^3]"
    kappa_SEI = 5e-6 ;%"SEI ionic conductivity [S/m]"
    i_0_SEI = 8e-4 ;%"SEI Exchange Current [A/m^2]"
    OCP_SEI = 0 ; %SEI OCP [V]
    R_SEI =  L_SEI/kappa_SEI;%"SEI initial resistance [ohm*m^2]"
    c_SEI_0 = 0;%4541; % "EC soloution concentration [mol/m^3]" should be = epss_SEI*rho_SEI/M_SEI
    C_sol_0 = 4541; % Concentration of Solvent in bulk electrolyte [mol/m^3]
    SEI_correction = 1.1; %SEI growth correction factor
    i_0_SEI = SEI_correction*i_0_SEI ;
    k_SEI = 3e-10; %1.137e-9; %i_0_SEI/(F*C_sol_0);%1.37e-9"Rate constant, SEI [m/s]" kamyab et al has 4.3e-16 , changed from -5, k_SEI =  i_0_SEI/(F*C_sol_0); for comparing %-9 works well
    
    %Li plating Parameters all taken from Ren et al
    OCP_Li = 0; % "Open Circuit Potential of Lithium [V]"
    alpha_pl_neg = 0.3;%0.5; %0.3; % Anodic Plating Transfer Coefficient , Ren et al
    alpha_pl_pos = 0.7;%0.5;%0.7 % Cathodic Plating Transfer Coefficient , Ren et al
    k_pl = 3.0e-6; % Lithium plating rate constant [mols^0.7/((m^-1.1)(s^-1))] Generic Formula for units = [mols^(1-alpha_pl_an)/m^(-2+3*alpha_pl_an)]
    i_pl_st = F*k_pl*(C_e_initial^alpha_pl_neg); % Lithium plating exchange current, Ren et al, [A/m^2]
    kappa_li = 1/(92.8e-9); % Lithium Conductivity [S/m], taken from wikipedia
    M_li = 6.94e-3; % Molar Mass of Lithium from RSC website [Kg/mol]
    rho_li = 0.534e3; % Density of Lithiun from RSC website [Kg/m^3]
    c_Li_0 = 0; %Initial Plated Lithium [mols/m^3]
    L_li = 0 ; %Initial Plated lithium film [m]
    R_Li = L_li/kappa_li ; %Initial Lithium resistance [ohm*m^2]
    

    % heat transfer coefficient
    h_ext = 5 ; %"heat transfer coefficient [W/(K*m^2)] "
    
  
    %I_1c = 16; % 1C current magnitude for model [A/m^2]
    
    if soc_initial_pos > soc_initial_neg
    Q_pos_initial =  (1/3600)*F*L_pos*eps_pos*C_pos_initial*abs(soc_initial_pos - 0); % Units (Ah/m^2)
    Q_neg_initial = (1/3600)*F*L_neg*eps_neg*C_neg_initial*abs(1.0 - soc_initial_neg ); % Units (Ah/m^2)
    end

    if soc_initial_pos < soc_initial_neg
    Q_pos_initial =  (1/3600)*F*L_pos*eps_pos*C_pos_initial*abs(1.0 - soc_initial_pos); % Units (Ah/m^2)
    Q_neg_initial = (1/3600)*F*L_neg*eps_neg*C_neg_initial*abs(soc_initial_neg - 0); % Units (Ah/m^2)

    end
    I_1c = max(Q_neg_initial,Q_pos_initial); % A/m^2 1c discharge rate;
    I_1c = 17.0;
    
    
    % Generate diff_OCP_C_s 

    x_pchip = linspace(0.1,0.9,1000);
    x = x_pchip;
    
   
    y_pchip = linspace(0.1,0.9,1000);
    y = y_pchip;



   %%% Plot for LixC6 from Liu et al
   OCP_neg_Liu = -0.16 +1.32.*exp(-3.*x) + 10.*exp(-2000.*x);

   %%% Plot for LixMnO2 from Liu et al
   OCP_pos_Liu = 4.19829 +0.0565661.*tanh(-14.5546.*y + 8.60942)  -0.0275479.*(1./((0.998432 - y).^(0.492465)) -1.90111) ...
    -0.157123.*exp(0.04738.*(y.^8)) + 0.810239.*exp(-40.*(y - 0.133875));

   diff_OCP_C_s_neg = diff(OCP_neg_Liu,1,2);
   x_pchip = linspace(0,1,999);

   diff_OCP_C_s_pos = diff(OCP_pos_Liu,1,2);
   y_pchip = linspace(0,1,999);






end


%% Taken from Hosseinzadeh et al, Liebing et al (NMC111/C6 Cell)
if param == 2

current_collector = 0; % model has current collectors
I_1c = 30; % 1C current magnitude for model [A/m^2]

%Constants
R_const = 8.314; % Ideal gas constant [J/[K*mol]]
F = 96485; % Faraday Constant [Coulombs / mol]
T_c = 25.15; %Temperature in Celcius
T = 273 + T_c; %Temperature in kelvin
T_amb = T;
%Bruggman constants

brugg = 1.5;
brugg_neg = 1.5;
brugg_sep = 1.5;
brugg_pos = 1.5;


%Taken from Hosseinzadeh et al and liebig et al
%Cell Dimensions taken from 
L_ncc = 10e-6; % [m] 
L_neg = 110e-6;%75e-6; % [m]
L_sep = 24e-6;%17e-6; % [m] 
L_pos = 129e-6;%40e-6; % [m]
L_pcc = 10e-6; % [m]
L_tot_active = L_neg + L_sep + L_pos; % [m]
L_tot = L_ncc + L_tot_active + L_pcc ;% [m]

%Electrolyte/Seperator Parameters (1M LiPF6 EC:DMC (2:1))
C_e_initial = 1.2e3; %Initial Electrolyte Lithium Concentration [mol/[m^3]] %added
D_elec = 3.8e-10; % Diffusion Coefficient for Lithium in electrolyte [[m^2]/s] %Defined in model params %added
t_li = 0.363; % Lithium Ion Transference Number [unitless]
eps_sep_elec = 0.54; %Electrolyte Volume Fraction
D_sep_e_eff = D_elec*(eps_sep_elec^brugg_sep); % Effective Diffusion Coefficient in electrolyte
C_p_sep = 750; % Specific Heat Capacity [J/(Kg*K)]
row_sep = 1626; % Density of Seperator [kg/[m^3]]
lambda_sep = 0.16; % Thermal Conductivity [J/(m*K)] % changed
E_a_kappa = 20000; % Activation Energy for Conductivity  [J/mol]
E_a_D_e = 35000; % Activation Energy for Conductivity  [J/mol]

%Negative Current Collector (copper)
C_p_ncc = 385; % Specific Heat Capacity [J/(Kg*K)]
row_ncc = 8940; % Density of Seperator [kg/[m^3]]
lambda_ncc = 401; % Thermal Conductivity [J/(m*K)] %changed
sigma_ncc = 3.55e7 ; % Electrical Conductivity [S/m] %Not important in lamorgese model

%Negative Electrode
soc_initial_neg = 0.56; %added
eps_neg_elec = 0.4; % Electrolyte Volume Fraction (porosity)
%eps_neg_fill = 0.0566; % Filler Volume Fraction
eps_neg = 0.58; % Solid Volume Fraction
C_max_neg = 29802 ; % Maximum Solid Lithium Concentration [mol/[m^3]] %added
C_neg_initial = soc_initial_neg*C_max_neg ;% Initial Lithium Concentration [mol/[m^3]] %added
Rs_neg = 26.2e-6; % Electrode Particle Size [m] %changed
D_neg_s = 3e-13; % Solid Diffusion Coefficient [m^2/s] %added
D_neg_eff_s = D_neg_s;
sigma_neg = 100.0 ; % Electrical Conductivity [S/m] %added
sigma_neg_eff = sigma_neg*(eps_neg);  
alpha = 0.5; % Transfer Coefficient
row_neg = 5031.67; % Density of Electrode [kg/[m^3]]
K_0_neg = 5e-10; % Negative Rate Constant
As_neg = 3*(eps_neg/Rs_neg) ; % Surface area to volume ration [m-1]
D_neg_e_eff = D_elec*(eps_neg_elec^brugg_neg); %Effective Diffusivity in electrode's electrolyte
C_p_neg = 750; % Specific Heat Capacity [J/(Kg*K)]
E_a_D_s_neg = 58000; % Activation Energy for Temperature Dependent Solid Phase Diffusion [J/mol]
E_a_K_0_neg = 20000; % Activation Energy for Temperature Reaction Constant  [J/mol]
lambda_neg = 1.7; % Thermal Conductivity [J/(m*K)]



%Positive Current Collector
C_p_pcc = 897; % Specific Heat Capacity [J/(Kg*K)]
row_pcc = 2700; % Density of Seperator [kg/[m^3]]
lambda_pcc = 237; % Thermal Conductivity [J/(m*K)] %changed
sigma_pcc = 3.55e7 ; % Electrical Conductivity [S/m] %Not used in Lamorgesse et al

%Positive Electrode
eps_pos_elec = 0.32; %% Electrolyte Volume Fraction 
%eps_pos_fill = 0.43; % Filler Volume Fraction
eps_pos = 0.43; % Solid Volume Fraction
soc_initial_pos = 0.25; %added
C_max_pos = 87593 ; % Maximum Solid Lithium Concentration [mol/[m^3]] %added
C_pos_initial = soc_initial_pos*C_max_pos ;% Initial Lithium Concentration [mol/[m^3]] %added
Rs_pos = 10.7e-6; % Electrode Particle Size [m] %changed
D_pos_s = 7e-14; % Solid Diffusion Coefficient [m^2/s] %added
D_pos_eff_s = D_pos_s;
sigma_pos = 100; % Electrical Conductivity [S/m] %added
sigma_pos_eff = sigma_pos*(eps_pos);  
alpha = 0.5; % Transfer Coefficient
row_pos = 4670; % Density of Electrode [kg/[m^3]]
K_0_pos = 2.5e-10; % Negative Rate Constant
As_pos = 3*(eps_pos/Rs_pos) ; % Surface area to volume ratio [m-1]
D_pos_e_eff = D_elec*(eps_pos_elec^brugg_pos); %Effective Diffusivity in electrode's electrolyte
C_p_pos = 750; % Specific Heat Capacity [J/(Kg*K)]
E_a_D_s_pos = 29000; % Activation Energy for Temperature Dependent Solid Phase Diffusion [J/mol]
E_a_K_0_pos = 58000; % Activation Energy for Temperature Reaction Constant  [J/mol]
lambda_pos = 2.1; % Thermal Conductivity [J/(m*K)]
diff_U_T_pos = -7.255e-5; % Entropy Change  [V/K] changed from -10/F from comsol



%SEI Parameters
c_SEI_0 = 4541; % "EC soloution concentration [mol/m^3]"
D_SEI =  1e-13 ;% "EC diffusion coefficient [m^2/s]"
E_a_k_SEI  = 200e3 ;% "Activation Energy for reaction coefficient, SEI [J/mol]"
k_SEI = 1.37e-5 ;%"Rate constant, SEI [m/s]"
f_elec = 1; % "Electrolyte Activity coefficient"
L_SEI =  1e-9 ;%"Initial SEI Thickness [m]"
M_SEI =  0.162 ;%"SEI Molar Mass [kg/mol]"
R_SEI =  1e-4 ;%"SEI initial resistance [ohm*m^2]" 
alpha_SEI = 0.5 ;%"Charge transfer coefficient, SEI"
epss_SEI = 0.05 ;%"Solid phase volume fraction SEI"
rho_SEI = 1.69e3 ;% "Separator density [kg/m^3]"
kappa_SEI = 5e-6 ;%"SEI ionic conductivity [S/m]"
i_0_SEI = 8e-4 ;%"SEI Exchange Current [A/m^2]"
OCP_SEI = 0 ; %SEI OCP [V]
R_SEI =  L_SEI/kappa_SEI;%"SEI initial resistance [ohm*m^2]"
c_SEI_0 = 4541; % "EC soloution concentration [mol/m^3]"
C_sol_0 = 4541; % Concentration of Solvent in bulk electrolyte [mol/m^3]
SEI_correction = 1.0; %SEI growth correction factor
i_0_SEI = SEI_correction*i_0_SEI ;


% heat transfer coefficient
h_ext = 5 ; %"heat transfer coefficient [W/(K*m^2)] "


%I_1c = 16; % 1C current magnitude for model [A/m^2]

if soc_initial_pos > soc_initial_neg
    Q_pos_initial =  (1/3600)*F*L_pos*eps_pos*C_pos_initial*abs(soc_initial_pos - 0); % Units (Ah/m^2)
    Q_neg_initial = (1/3600)*F*L_neg*eps_neg*C_neg_initial*abs(1.0 - soc_initial_neg ); % Units (Ah/m^2)
end

if soc_initial_pos < soc_initial_neg
    Q_pos_initial =  (1/3600)*F*L_pos*eps_pos*C_pos_initial*abs(1.0 - soc_initial_pos); % Units (Ah/m^2)
    Q_neg_initial = (1/3600)*F*L_neg*eps_neg*C_neg_initial*abs(soc_initial_neg - 0); % Units (Ah/m^2)

end
I_1c = max(Q_neg_initial,Q_pos_initial); % A/m^2 1c discharge rate;
I_1c = 30.0;
    


end






%% Taken from Torchio and Northrop for temperature modelling
if param == 3

current_collector = 0; % model has current collectors
I_1c = 30; % 1C current magnitude for model [A/m^2]

%Constants
R_const = 8.314; % Ideal gas constant [J/[K*mol]]
F = 96487; % Faraday Constant [Coulombs / mol]
T_c = 25.15; %Temperature in Celcius
T = 273 + T_c; %Temperature in kelvin
T_amb = T;
%Bruggman constants

brugg = 4.0;
brugg_neg = brugg;
brugg_sep = brugg;
brugg_pos = brugg;


%Taken from Torchio and Northrop
%Cell Dimensions
L_pcc = 10e-6; % [m]
L_pos = 80e-6; % [m]
L_sep = 25e-6; % [m] 
L_neg = 88e-6; % [m]
L_ncc = 10e-6; % [m] 
L_tot_active = L_neg + L_sep + L_pos; % [m]
L_tot = L_ncc + L_tot_active + L_pcc ;% [m]

%Electrolyte/Seperator Parameters
C_e_initial = 1.0e3; %Initial Electrolyte Lithium Concentration [mol/[m^3]] %changed
D_elec = 7.5e-10; % Diffusion Coefficient for Lithium in electrolyte [[m^2]/s]
t_li = 0.364; % Lithium Ion Transference Number [unitless]
eps_sep_elec = 0.724; %Electrolyte Volume Fraction
D_sep_e_eff = D_elec*(eps_sep_elec^brugg_sep); % Effective Diffusion Coefficient in electrolyte
C_p_sep = 700; % Specific Heat Capacity [J/(Kg*K)]
row_sep = 1100; % Density of Seperator [kg/[m^3]]
lambda_sep = 0.16; % Thermal Conductivity [J/(m*K)] % changed

%Negative Current Collector
C_p_ncc = 385; % Specific Heat Capacity [J/(Kg*K)]
row_ncc = 8940; % Density of Seperator [kg/[m^3]]
lambda_ncc = 401; % Thermal Conductivity [J/(m*K)] %changed
sigma_ncc = 5.96e7 ; % Electrical Conductivity [S/m]

%Negative Electrode
soc_initial_neg = 0.8551;
%soc_initial_neg = 0.8;
eps_neg_elec = 0.485; % Electrolyte Volume Fraction (porosity)
eps_neg_fill = 0.0326; % Filler Volume Fraction
eps_neg = 1 - eps_neg_elec - eps_neg_fill; % Solid Volume Fraction
C_max_neg = 30555 ; % Maximum Solid Lithium Concentration [mol/[m^3]]
C_neg_initial = soc_initial_neg*C_max_neg ;% Initial Lithium Concentration [mol/[m^3]]
Rs_neg = 2e-6; % Electrode Particle Size [m] %changed
D_neg_s = 3.9e-14; % Solid Diffusion Coefficient [m^2/s]
D_neg_eff_s = D_neg_s;
sigma_neg = 100.0 ; % Electrical Conductivity [S/m]
sigma_neg_eff = sigma_neg*(eps_neg);  
alpha = 0.5; % Transfer Coefficient
row_neg = 2500; % Density of Electrode [kg/[m^3]]
K_0_neg = 5.031e-11; % Negative Rate Constant
As_neg = 3*(eps_neg/Rs_neg) ; % Surface area to volume ration [m-1]
D_neg_e_eff = D_elec*(eps_neg_elec^brugg_neg); %Effective Diffusivity in electrode's electrolyte
C_p_neg = 700; % Specific Heat Capacity [J/(Kg*K)]
E_a_D_s_neg = 5000; % Activation Energy for Temperature Dependent Solid Phase Diffusion [J/mol]
E_a_K_0_neg = 5000; % Activation Energy for Temperature Reaction Constant  [J/mol]
lambda_neg = 1.7; % Thermal Conductivity [J/(m*K)]
%soc_initial_neg = C_neg_initial/C_max_neg;


%Positive Current Collector
C_p_pcc = 897; % Specific Heat Capacity [J/(Kg*K)]
row_pcc = 2700; % Density of Seperator [kg/[m^3]]
lambda_pcc = 237; % Thermal Conductivity [J/(m*K)] %changed
sigma_pcc = 3.55e7 ; % Electrical Conductivity [S/m]

%Positive Electrode
soc_initial_pos = 0.4995;
%soc_initial_pos = 0.5;
eps_pos_elec = 0.385; %% Electrolyte Volume Fraction 
eps_pos_fill = 0.025; % Filler Volume Fraction
eps_pos = 1 - eps_pos_elec - eps_pos_fill; % Solid Volume Fraction
C_max_pos = 51554 ; % Maximum Solid Lithium Concentration [mol/[m^3]]
C_pos_initial = soc_initial_pos*C_max_pos ;% Initial Lithium Concentration [mol/[m^3]]
Rs_pos = 2e-6; % Electrode Particle Size [m] %changed
D_pos_s = 1e-14; % Solid Diffusion Coefficient [m^2/s]
D_pos_eff_s = D_pos_s;
sigma_pos = 100.0 ; % Electrical Conductivity [S/m]
sigma_pos_eff = sigma_pos*(eps_pos);  
alpha = 0.5; % Transfer Coefficient
row_pos = 2500; % Density of Electrode [kg/[m^3]]
K_0_pos = 2.334e-11; % Negative Rate Constant
As_pos = 3*(eps_pos/Rs_pos) ; % Surface area to volume ratio [m-1]
D_pos_e_eff = D_elec*(eps_pos_elec^brugg_pos); %Effective Diffusivity in electrode's electrolyte
C_p_pos = 700; % Specific Heat Capacity [J/(Kg*K)]
E_a_D_s_pos = 5000; % Activation Energy for Temperature Dependent Solid Phase Diffusion [J/mol]
E_a_K_0_pos = 5000; % Activation Energy for Temperature Reaction Constant  [J/mol]
lambda_pos = 2.1; % Thermal Conductivity [J/(m*K)]
%soc_initial_pos = C_pos_initial/C_max_pos;

%SEI Parameters
c_SEI_0 = 4541; % "EC soloution concentration [mol/m^3]"
D_SEI =  1e-13 ;% "EC diffusion coefficient [m^2/s]"
E_a_k_SEI  = 200e3 ;% "Activation Energy for reaction coefficient, SEI [J/mol]"
k_SEI = 1.37e-5 ;%"Rate constant, SEI [m/s]"
f_elec = 1; % "Electrolyte Activity coefficient"
L_SEI =  1e-9 ;%"Initial SEI Thickness [m]"
M_SEI =  0.162 ;%"SEI Molar Mass [kg/mol]"
R_SEI =  1e-4 ;%"SEI initial resistance [ohm*m^2]" 
alpha_SEI = 0.5 ;%"Charge transfer coefficient, SEI"
epss_SEI = 0.05 ;%"Solid phase volume fraction SEI"
rho_SEI = 1.69e3 ;% "Separator density [kg/m^3]"
kappa_SEI = 5e-6 ;%"SEI ionic conductivity [S/m]"
i_0_SEI = 8e-4 ;%"SEI Exchange Current [A/m^2]"
OCP_SEI = 0 ; %SEI OCP [V]
R_SEI =  L_SEI/kappa_SEI;%"SEI initial resistance [ohm*m^2]"
c_SEI_0 = 4541; % "EC soloution concentration [mol/m^3]"
C_sol_0 = 4541; % Concentration of Solvent in bulk electrolyte [mol/m^3]
SEI_correction = 1.1; %SEI growth correction factor
i_0_SEI = SEI_correction*i_0_SEI ;

% heat transfer coefficient
h_ext = 5 ; %"heat transfer coefficient [W/(K*m^2)] "


if soc_initial_pos > soc_initial_neg
    Q_pos_initial =  (1/3600)*F*L_pos*eps_pos*C_pos_initial*abs(soc_initial_pos - 0); % Units (Ah/m^2)
    Q_neg_initial = (1/3600)*F*L_neg*eps_neg*C_neg_initial*abs(1.0 - soc_initial_neg ); % Units (Ah/m^2)
    end

    if soc_initial_pos < soc_initial_neg
    Q_pos_initial =  (1/3600)*F*L_pos*eps_pos*C_pos_initial*abs(1.0 - soc_initial_pos); % Units (Ah/m^2)
    Q_neg_initial = (1/3600)*F*L_neg*eps_neg*C_neg_initial*abs(soc_initial_neg - 0); % Units (Ah/m^2)

    end
    I_1c = max(Q_neg_initial,Q_pos_initial); % A/m^2 1c discharge rate;

end
