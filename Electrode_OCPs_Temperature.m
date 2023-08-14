%%% Script to choose Electrode OCP model equation
%% Taken from Liu et al
if param == 1

   diff_OCP_T_neg = diff_U_T_neg*ones(1,N_neg);
    
   diff_OCP_T_pos =  diff_U_T_pos*ones(1,N_pos);
    
   diff_OCP_T(1,:) = [diff_OCP_T_neg zeros(1,N_sep) diff_OCP_T_pos]; 
   
   %%% Plot for LixC6 from Liu et al
   OCP_neg_Liu = -0.16 +1.32.*exp(-3.*x) + 10.*exp(-2000.*x);

   %%% Plot for LixMnO2 from Liu et al
    OCP_pos_Liu = 4.19829 +0.0565661.*tanh(-14.5546.*y + 8.60942)  -0.0275479.*(1./((0.998432 - y).^(0.492465)) -1.90111) ...
    -0.157123.*exp(0.04738.*(y.^8)) + 0.810239.*exp(-40.*(y - 0.133875));

    OCP_ref = [OCP_neg_Liu  zeros(1,N_sep) OCP_pos_Liu];
   
   OCP(i,:) = OCP_ref +(T(i,x_neg_start:x_pos_end)-T_amb).*diff_OCP_T;

end



%% All Equations are from Hosseinzadeh et al, Liebing et al (NMC111/C6 Cell)
if param == 2
    
    %Calculate OCP Entropy Change
    diff_OCP_T_neg = 0.001*(0.005269056 + 3.299265709*x -91.79325798*(x.^2) + 1004.911008*(x.^3) - 5812.278127*(x.^4) + ...
                    +19329.7549*(x.^5) - 37147.8947*(x.^6) + 38379.18127*(x.^7) - 16515.05308*(x.^8))./(...
                    1 - 48.09287227*(x) + 1017.234804*(x.^2) - 10481.80419*(x.^3) + 59431.3*(x.^4)-195881.6488*(x.^5) + ...
                    374577.3152*(x.^6) - 385821.1607*(x.^7) + 165705.8597*(x.^8));
    
    diff_OCP_T_pos =  diff_U_T_pos*ones(1,N_pos);

    diff_OCP_T(1,:) = [diff_OCP_T_neg zeros(1,N_sep) diff_OCP_T_pos];

    %Calculate Reference OCP

    OCP_ref_neg = 0.6379+0.5416.*exp( -305.5309.*x ) + 0.044.*tanh( -( (x-0.1958)./(0.1088) ) ) ...
                    -0.1978.*tanh( (x-0.0117)./(0.0529) ) - 0.0175.*tanh( (x-0.5692)./(0.0875) );
    

    OCP_ref_pos = -10.72.*(y.^4) + 23.88.*(y.^3) - 16.77.*(y.^2) + 2.595.*y + 4.563;
    
    

    OCP_ref = [OCP_ref_neg zeros(1,N_sep) OCP_ref_pos];

    %Calculate OCPs

    OCP(i,:) = OCP_ref +(T(i,x_neg_start:x_pos_end)-T_amb).*diff_OCP_T;


end


%% Taken From Northop and Torchio

if param == 3
    %Calculate OCP Entropy Change
    diff_OCP_T_neg = 0.001*(0.005269056 + 3.299265709*x -91.79325798*(x.^2) + 1004.911008*(x.^3) - 5812.278127*(x.^4) + ...
                    +19329.7549*(x.^5) - 37147.8947*(x.^6) + 38379.18127*(x.^7) - 16515.05308*(x.^8))./(...
                    1 - 48.09287227*(x) + 1017.234804*(x.^2) - 10481.80419*(x.^3) + 59431.3*(x.^4)-195881.6488*(x.^5) + ...
                    374577.3152*(x.^6) - 385821.1607*(x.^7) + 165705.8597*(x.^8));
    
    diff_OCP_T_pos =  -0.001*(0.199521039-0.928373822*y + 1.364550689000003*(y.^2) - 0.6115448939999998*(y.^3))./(...
                        1-5.661479886999997*y + 11.47636191*(y.^2) - 9.82431213599998*(y.^3) +3.048755063*(y.^4));
    
    diff_OCP_T(1,:) = [diff_OCP_T_neg zeros(1,N_sep) diff_OCP_T_pos];

    %Calculate Reference OCP

    OCP_ref_neg = 0.7222 + 0.1387*x + 0.029*(x.^(0.5)) - (0.0172./x) + (0.0019./(x.^(1.5))) + 0.2808*exp(0.9-15.*x) - 0.7984*exp(0.4465.*x-0.4108);

    OCP_ref_pos = (-4.656 + (88.669.*(y.^2)) - 401.119.*(y.^4) + 342.909*(y.^6) - 462.471.*(y.^8) + 433.434.*(y.^10))./(...
        -1+18.933*(y.^2) - 79.532*(y.^4) + 37.311*(y.^6) - 73.083*(y.^8) + 95.96*(y.^10));

    OCP_ref = [OCP_ref_neg zeros(1,N_sep) OCP_ref_pos];

    %Calculate OCPs

    OCP(i,:) = OCP_ref +(T(i,x_neg_start:x_pos_end)-T_amb).*diff_OCP_T;


end
