% Code to plot data upto error
close all;
stop_time = i;
b = i;
figure(1);
grid on;

plot(linspace(0,stop_time,stop_time)*time_step,[psi_s(1:stop_time,N_tot_active) - psi_s(1:stop_time,1)]');

xlabel('Time (s)');
ylabel('Cell Potential [V]');
if mass_con == 0
figure(2);
hold on;
plot(linspace(0,N_tot_active+2,N_tot_active),C_e(1,:),'r');
plot(linspace(0,N_tot_active+2,N_tot_active),C_e(3,:),'b');
plot(linspace(0,N_tot_active+2,N_tot_active),C_e(round(stop_time/2),:),'g');
plot(linspace(0,N_tot_active+2,N_tot_active),C_e(stop_time,:),'k');

xlabel('Elements');
ylabel('Electrolyte Lithium Ion Concentration');

end
if mass_con == 1
figure(2);
hold on;
plot(linspace(0,N_tot_active+2,N_tot_active+2),C_e(1,:),'r');
plot(linspace(0,N_tot_active+2,N_tot_active+2),C_e(3,:),'b');
plot(linspace(0,N_tot_active+2,N_tot_active+2),C_e(round(stop_time/2),:),'g');
plot(linspace(0,N_tot_active+2,N_tot_active+2),C_e(stop_time,:),'k');

xlabel('Elements');
ylabel('Electrolyte Lithium Ion Concentration');

end


figure(3);
hold on;
plot(linspace(0,N_tot_active,N_tot_active),j_flux(3,:),'r');
plot(linspace(0,N_tot_active,N_tot_active),j_flux(round(stop_time/2),:),'b');
plot(linspace(0,N_tot_active,N_tot_active),j_flux(stop_time,:),'g');

xlabel('Elements');
ylabel('Lithium Ion Transfer Flux in Cathode [mols/(m^2*s)]');

figure(4);
hold on;
plot(linspace(0,N_tot_active,N_tot_active),C_se(1,:),'r');
plot(linspace(0,N_tot_active,N_tot_active),C_se(2,:),'b');
plot(linspace(0,N_tot_active,N_tot_active),C_se(stop_time,:),'g');


xlabel('Elements');
ylabel('Surface Lithium Concentration');

figure(5);
hold on;
plot(linspace(1,N_tot_active,N_tot_active),psi_e(1,:),'r');
plot(linspace(1,N_tot_active,N_tot_active),psi_e(2,:),'b');
plot(linspace(1,N_tot_active,N_tot_active),psi_e(stop_time,:),'g');

xlabel('Elements');
ylabel('Electrolyte Potential (V)');

if current_collector == 1
    if and(thermal > 0,thermal < 10)  
    %figure for overall temperature distribution
    figure(6);
    hold on;
    plot(linspace(0,N_tot,N_tot),T(1,:),'r');
    plot(linspace(0,N_tot,N_tot),T(round(stop_time/4),:),'b');
    plot(linspace(0,N_tot,N_tot),T(round(stop_time/3),:),'g');
    plot(linspace(0,N_tot,N_tot),T(round(stop_time/2),:),'k');
    plot(linspace(0,N_tot,N_tot),T(stop_time,:),'m');
    %plot(linspace(1,N_tot,N_tot),T(stop_time,:),'y');
    
    xlabel('Elements');
    ylabel('Entire Cell Temperature (K)');
    
    xline(N_cc+0.5,'--b','Ncc/Neg');
    xline(x_neg_end+0.5,'--b','Neg/Sep');
    xline(x_sep_end+0.5,'--b','Sep/Pos');
    xline(x_pos_end+0.5,'--b','Pos/Pcc');
    
    
    %figure for temperature distribution in active reigon
    figure(7);
    hold on;
    plot(linspace(1,N_tot_active,N_tot_active),T(1,x_neg_start:x_pos_end),'r');
    plot(linspace(1,N_tot_active,N_tot_active),T(2,x_neg_start:x_pos_end),'b');
    plot(linspace(1,N_tot_active,N_tot_active),T(3,x_neg_start:x_pos_end),'g');
    plot(linspace(1,N_tot_active,N_tot_active),T(4,x_neg_start:x_pos_end),'k');
    plot(linspace(1,N_tot_active,N_tot_active),T(stop_time,x_neg_start:x_pos_end),'m');
    %plot(linspace(1,N_tot_active,N_tot_active),T(stop_time,x_neg_start:x_pos_end),'y');
    
    
    xlabel('Elements');
    ylabel('Active Reigon Temperature (K)');
    
    
    %figure for looking second order differential of Temperature
    figure(8);
    hold on;
    plot(linspace(1,N_tot,N_tot),diff_T_x_2_comp(2,:),'r');
    plot(linspace(1,N_tot,N_tot),diff_T_x_2_comp(b,:),'b');
    %plot(linspace(1,N_tot,N_tot),diff_T_x_2_comp(4,:),'g');
    plot(linspace(1,N_tot,N_tot),diff_T_x_2_comp(stop_time-1,:),'m');
    
    
    xlabel('Elements');
    ylabel('Second Order Differential Temperature * Lambda');
    
    
    %figure for looking a heat generation component.
    figure(9);
    hold on;
    plot(linspace(0,N_tot,N_tot),Q_gen(2,:),'r');
    plot(linspace(0,N_tot,N_tot),Q_gen(round(stop_time/4),:),'b');
    plot(linspace(0,N_tot,N_tot),Q_gen(round(stop_time/3),:),'g');
    plot(linspace(0,N_tot,N_tot),Q_gen(round(stop_time/2),:),'k');
    plot(linspace(0,N_tot,N_tot),Q_gen(stop_time,:),'m');
    
    
    xlabel('Elements');
    ylabel('Heat Generation Component');
    
    xline(N_cc+0.5,'--b','Ncc/Neg');
    xline(x_neg_end+0.5,'--b','Neg/Sep');
    xline(x_sep_end+0.5,'--b','Sep/Pos');
    xline(x_pos_end+0.5,'--b','Pos/Pcc');
    
    
    figure(10);
    hold on;
    plot(linspace(0,stop_time,stop_time)*time_step,T_ave(1:stop_time));
    
    xlabel('Time (s)');
    ylabel('Temperature (K)');
    
    
    figure(12);
    hold on;
    plot(linspace(0,N_tot_active,N_tot_active), Q_gen(b,x_neg_start:x_pos_end),'k');
    plot(linspace(0,N_tot_active,N_tot_active), Q_ohm(b,:),'r');
    plot(linspace(0,N_tot_active,N_tot_active), Q_rxn(b,:),'g');
    plot(linspace(0,N_tot_active,N_tot_active), Q_rev(b,:),'b');
    
    xlabel('Elements');
    ylabel('Heat Generation Component');
    
    legend('Total Heat', 'Ohmic Heating','RXN Heating','Rev Heating');
    
    xline(N_neg+0.5,'--b','Neg/Sep');
    xline(N_neg+N_sep+0.5,'--b','Sep/Pos');
    
    %figure for overall temperature distribution
    figure(13);
    hold on;
    plot(linspace(1,N_tot,N_tot),diff_T_t(1,:),'r');
    %plot(linspace(1,N_tot,N_tot),diff_T_t(round(stop_time/4),:),'b');
    plot(linspace(1,N_tot,N_tot),diff_T_t(b,:),'g');
    %plot(linspace(1,N_tot,N_tot),diff_T_t(round(stop_time/2),:),'k');
    %plot(linspace(1,N_tot,N_tot),diff_T_t(stop_time,:),'m');
    %plot(linspace(1,N_tot,N_tot),T(stop_time,:),'y');
    
    xlabel('Elements');
    ylabel('Temperature (K)');
    
    
    figure(14);
    hold on;
    plot(linspace(0,stop_time,stop_time)*time_step,T_ave(1:stop_time));
    
    xlabel('Time (s)');
    ylabel('Overall Temperature (K)');
    
    figure(15);
    hold on;
    plot(linspace(0,stop_time,stop_time)*time_step,T_ave_ncc(1:stop_time));
    
    xlabel('Time (s)');
    ylabel('NCC Temperature (K)');
    
    figure(16);
    hold on;
    plot(linspace(0,stop_time,stop_time)*time_step,T_ave_neg(1:stop_time));
    
    xlabel('Time (s)');
    ylabel('Neg Temperature (K)');
    
    figure(17);
    hold on;
    plot(linspace(0,stop_time,stop_time)*time_step,T_ave_sep(1:stop_time));
    
    xlabel('Time (s)');
    ylabel('Sep Temperature (K)');
    
    figure(18);
    hold on;
    plot(linspace(0,stop_time,stop_time)*time_step,T_ave_pos(1:stop_time));
    
    xlabel('Time (s)');
    ylabel('Pos Temperature (K)');
    
    figure(19);
    hold on;
    plot(linspace(0,stop_time,stop_time)*time_step,T_ave_pcc(1:stop_time));
    
    xlabel('Time (s)');
    ylabel('PCC Temperature (K)');
    
    %%% Plots at different points
    
    ncc_neg = 1;
    neg_sep = N_neg;
    sep_pos = N_neg+N_sep+1;
    pos_pcc = N_tot_active;
    
    %%%C_e at different points
    figure(20);
    hold on;
    plot(linspace(0,stop_time,stop_time)*time_step,C_e(1:stop_time,ncc_neg)','r');
    plot(linspace(0,stop_time,stop_time)*time_step,C_e(1:stop_time,neg_sep)','b');
    plot(linspace(0,stop_time,stop_time)*time_step,C_e(1:stop_time,sep_pos)','g');
    plot(linspace(0,stop_time,stop_time)*time_step,C_e(1:stop_time,pos_pcc)','k');
    
    ylabel('Electrolyte Concenctration [mols/m^3]');
    xlabel('Time (s)');
    grid on;
    
    xlim([0,3500]);
    ylim([0 1700]);
    
    legend('Electrolyte Concentration at Ncc/Neg Boundary', 'Electrolyte Concentration at Neg/Sep Boundary', 'Electrolyte Concentration at Sep/Pos Boundary', 'Electrolyte Concentration at Pos/Pcc Boundary','Location','best');
    
    %%%C_s_surf at different points
    figure(21);
    hold on;
    plot(linspace(0,stop_time,stop_time)*time_step,C_se(1:stop_time,ncc_neg)','r');
    plot(linspace(0,stop_time,stop_time)*time_step,C_se(1:stop_time,neg_sep)','b');
    plot(linspace(0,stop_time,stop_time)*time_step,C_se(1:stop_time,sep_pos)','g');
    plot(linspace(0,stop_time,stop_time)*time_step,C_se(1:stop_time,pos_pcc)','k');
    
    ylabel('Surface Solid Concenctration [mols/m^3]');
    xlabel('Time (s)');
    grid on;
    
    
    legend('Solid Surface Concentration at Ncc/Neg Boundary', 'Solid Surface Concentration at Neg/Sep Boundary', 'Solid Surface Concentration at Sep/Pos Boundary', 'Solid Surface Concentration at Pos/Pcc Boundary','Location','best');
    
    %C_s_ave at different points
    figure(22);
    hold on;
    plot(linspace(0,stop_time,stop_time)*time_step,C_s_ave(1:stop_time,ncc_neg)','r');
    plot(linspace(0,stop_time,stop_time)*time_step,C_s_ave(1:stop_time,neg_sep)','b');
    plot(linspace(0,stop_time,stop_time)*time_step,C_s_ave(1:stop_time,sep_pos)','g');
    plot(linspace(0,stop_time,stop_time)*time_step,C_s_ave(1:stop_time,pos_pcc)','k');
    
    ylabel('Average Solid Concenctration [mols/m^3]');
    xlabel('Time (s)');
    
    
    legend('Average Solid Concentration at Ncc/Neg Boundary', 'Average Solid Concentration at Neg/Sep Boundary', 'Average Solid Concentration at Sep/Pos Boundary', 'Average Solid Concentration at Pos/Pcc Boundary','Location','best');
    
    grid on;
    
    %psi_e at different points
    figure(23);
    hold on;
    plot(linspace(0,stop_time,stop_time)*time_step,psi_e(1:stop_time,ncc_neg)','r');
    plot(linspace(0,stop_time,stop_time)*time_step,psi_e(1:stop_time,neg_sep)','b');
    plot(linspace(0,stop_time,stop_time)*time_step,psi_e(1:stop_time,sep_pos)','g');
    plot(linspace(0,stop_time,stop_time)*time_step,psi_e(1:stop_time,pos_pcc)','k');
    
    ylabel('Electrolyte Potential [V]');
    xlabel('Time (s)');
    
    
    legend('Electrolyte Potential at Ncc/Neg Boundary', 'Electrolyte Potential at Neg/Sep Boundary', 'Electrolyte Potential at Sep/Pos Boundary', 'Electrolyte Potential at Pos/Pcc Boundary','Location','best');
    
    grid on;
    
    %psi_s at different points
    figure(24);
    hold on;
    plot(linspace(0,stop_time,stop_time)*time_step,psi_s(1:stop_time,ncc_neg)','r');
    plot(linspace(0,stop_time,stop_time)*time_step,psi_s(1:stop_time,neg_sep)','b');
    plot(linspace(0,stop_time,stop_time)*time_step,psi_s(1:stop_time,sep_pos)','g');
    plot(linspace(0,stop_time,stop_time)*time_step,psi_s(1:stop_time,pos_pcc)','k');
    
    ylabel('Solid Potential [V]');
    xlabel('Time (s)');
    
    
    legend('Solid Potential at Ncc/Neg Boundary', 'Solid Potential at Neg/Sep Boundary', 'Solid Potential at Sep/Pos Boundary', 'Solid Potential at Pos/Pcc Boundary','Location','best');
    
    grid on;
    
    
    end

end

if current_collector == 0
    if and(thermal > 0,thermal < 10)  
    %figure for overall temperature distribution
    figure(6);
    hold on;
    plot(linspace(0,N_tot,N_tot),T(1,:),'r');
    plot(linspace(0,N_tot,N_tot),T(round(stop_time/4),:),'b');
    plot(linspace(0,N_tot,N_tot),T(round(stop_time/3),:),'g');
    plot(linspace(0,N_tot,N_tot),T(round(stop_time/2),:),'k');
    plot(linspace(0,N_tot,N_tot),T(stop_time,:),'m');
    %plot(linspace(1,N_tot,N_tot),T(stop_time,:),'y');
    
    xlabel('Elements');
    ylabel('Entire Cell Temperature (K)');
    
   
    xline(x_neg_end+0.5,'--b','Neg/Sep');
    xline(x_sep_end+0.5,'--b','Sep/Pos');
    
    
    
    %figure for temperature distribution in active reigon
    figure(7);
    hold on;
    plot(linspace(1,N_tot_active,N_tot_active),T(1,x_neg_start:x_pos_end),'r');
    plot(linspace(1,N_tot_active,N_tot_active),T(2,x_neg_start:x_pos_end),'b');
    plot(linspace(1,N_tot_active,N_tot_active),T(3,x_neg_start:x_pos_end),'g');
    plot(linspace(1,N_tot_active,N_tot_active),T(4,x_neg_start:x_pos_end),'k');
    plot(linspace(1,N_tot_active,N_tot_active),T(stop_time,x_neg_start:x_pos_end),'m');
    %plot(linspace(1,N_tot_active,N_tot_active),T(stop_time,x_neg_start:x_pos_end),'y');
    
    
    xlabel('Elements');
    ylabel('Active Reigon Temperature (K)');
    
    
    %figure for looking second order differential of Temperature
    figure(8);
    hold on;
    plot(linspace(1,N_tot,N_tot),diff_T_x_2_comp(2,:),'r');
    plot(linspace(1,N_tot,N_tot),diff_T_x_2_comp(b,:),'b');
    %plot(linspace(1,N_tot,N_tot),diff_T_x_2_comp(4,:),'g');
    plot(linspace(1,N_tot,N_tot),diff_T_x_2_comp(stop_time-1,:),'m');
    
    
    xlabel('Elements');
    ylabel('Second Order Differential Temperature * Lambda');
    
    
    %figure for looking a heat generation component.
    figure(9);
    hold on;
    plot(linspace(0,N_tot,N_tot),Q_gen(2,:),'r');
    plot(linspace(0,N_tot,N_tot),Q_gen(round(stop_time/4),:),'b');
    plot(linspace(0,N_tot,N_tot),Q_gen(round(stop_time/3),:),'g');
    plot(linspace(0,N_tot,N_tot),Q_gen(round(stop_time/2),:),'k');
    plot(linspace(0,N_tot,N_tot),Q_gen(stop_time,:),'m');
    
    
    xlabel('Elements');
    ylabel('Heat Generation Component');
    
    
    xline(x_neg_end+0.5,'--b','Neg/Sep');
    xline(x_sep_end+0.5,'--b','Sep/Pos');
   
    
    
    figure(10);
    hold on;
    plot(linspace(0,stop_time,stop_time)*time_step,T_ave(1:stop_time));
    
    xlabel('Time (s)');
    ylabel('Temperature (K)');
    
    
    figure(12);
    hold on;
    plot(linspace(0,N_tot_active,N_tot_active), Q_gen(b,x_neg_start:x_pos_end),'k');
    plot(linspace(0,N_tot_active,N_tot_active), Q_ohm(b,:),'r');
    plot(linspace(0,N_tot_active,N_tot_active), Q_rxn(b,:),'g');
    plot(linspace(0,N_tot_active,N_tot_active), Q_rev(b,:),'b');
    
    xlabel('Elements');
    ylabel('Heat Generation Component');
    
    legend('Total Heat', 'Ohmic Heating','RXN Heating','Rev Heating');
    
    xline(N_neg+0.5,'--b','Neg/Sep');
    xline(N_neg+N_sep+0.5,'--b','Sep/Pos');
    
    %figure for overall temperature distribution
    figure(13);
    hold on;
    plot(linspace(1,N_tot,N_tot),diff_T_t(1,:),'r');
    plot(linspace(1,N_tot,N_tot),diff_T_t(round(stop_time/4),:),'b');
    plot(linspace(1,N_tot,N_tot),diff_T_t(b,:),'g');
    plot(linspace(1,N_tot,N_tot),diff_T_t(round(stop_time/2),:),'k');
    plot(linspace(1,N_tot,N_tot),diff_T_t(stop_time,:),'m');
    plot(linspace(1,N_tot,N_tot),T(stop_time,:),'y');
    
    xlabel('Elements');
    ylabel('Temperature (K)');
    
    
    figure(14);
    hold on;
    plot(linspace(0,stop_time-1,stop_time-1)*time_step,T_ave(1:stop_time-1));
    
    xlabel('Time (s)');
    ylabel('Overall Temperature (K)');
    
    
    
    figure(16);
    hold on;
    plot(linspace(0,stop_time-1,stop_time-1)*time_step,T_ave_neg(1:stop_time-1));
    
    xlabel('Time (s)');
    ylabel('Neg Temperature (K)');
    
    figure(17);
    hold on;
    plot(linspace(0,stop_time-1,stop_time-1)*time_step,T_ave_sep(1:stop_time-1));
    
    xlabel('Time (s)');
    ylabel('Sep Temperature (K)');
    
    figure(18);
    hold on;
    plot(linspace(0,stop_time-1,stop_time-1)*time_step,T_ave_pos(1:stop_time-1));
    
    xlabel('Time (s)');
    ylabel('Pos Temperature (K)');
    
   
    
    %%% Plots at different points
    
    ncc_neg = 1;
    neg_sep = N_neg;
    sep_pos = N_neg+N_sep+1;
    pos_pcc = N_tot_active;
    
    %%%C_e at different points
    figure(20);
    hold on;
    plot(linspace(0,stop_time,stop_time)*time_step,C_e(1:stop_time,ncc_neg)','r');
    plot(linspace(0,stop_time,stop_time)*time_step,C_e(1:stop_time,neg_sep)','b');
    plot(linspace(0,stop_time,stop_time)*time_step,C_e(1:stop_time,sep_pos)','g');
    plot(linspace(0,stop_time,stop_time)*time_step,C_e(1:stop_time,pos_pcc)','k');
    
    ylabel('Electrolyte Concenctration [mols/m^3]');
    xlabel('Time (s)');
    grid on;
    
    xlim([0,3500]);
    ylim([0 1700]);
    
    legend('Electrolyte Concentration at Ncc/Neg Boundary', 'Electrolyte Concentration at Neg/Sep Boundary', 'Electrolyte Concentration at Sep/Pos Boundary', 'Electrolyte Concentration at Pos/Pcc Boundary','Location','best');
    
    %%%C_s_surf at different points
    figure(21);
    hold on;
    plot(linspace(0,stop_time,stop_time)*time_step,C_se(1:stop_time,ncc_neg)','r');
    plot(linspace(0,stop_time,stop_time)*time_step,C_se(1:stop_time,neg_sep)','b');
    plot(linspace(0,stop_time,stop_time)*time_step,C_se(1:stop_time,sep_pos)','g');
    plot(linspace(0,stop_time,stop_time)*time_step,C_se(1:stop_time,pos_pcc)','k');
    
    ylabel('Surface Solid Concenctration [mols/m^3]');
    xlabel('Time (s)');
    grid on;
    
    
    legend('Solid Surface Concentration at Ncc/Neg Boundary', 'Solid Surface Concentration at Neg/Sep Boundary', 'Solid Surface Concentration at Sep/Pos Boundary', 'Solid Surface Concentration at Pos/Pcc Boundary','Location','best');
    
    %C_s_ave at different points
    figure(22);
    hold on;
    plot(linspace(0,stop_time,stop_time)*time_step,C_s_ave(1:stop_time,ncc_neg)','r');
    plot(linspace(0,stop_time,stop_time)*time_step,C_s_ave(1:stop_time,neg_sep)','b');
    plot(linspace(0,stop_time,stop_time)*time_step,C_s_ave(1:stop_time,sep_pos)','g');
    plot(linspace(0,stop_time,stop_time)*time_step,C_s_ave(1:stop_time,pos_pcc)','k');
    
    ylabel('Average Solid Concenctration [mols/m^3]');
    xlabel('Time (s)');
    
    
    legend('Average Solid Concentration at Ncc/Neg Boundary', 'Average Solid Concentration at Neg/Sep Boundary', 'Average Solid Concentration at Sep/Pos Boundary', 'Average Solid Concentration at Pos/Pcc Boundary','Location','best');
    
    grid on;
    
    %psi_e at different points
    figure(23);
    hold on;
    plot(linspace(0,stop_time,stop_time)*time_step,psi_e(1:stop_time,ncc_neg)','r');
    plot(linspace(0,stop_time,stop_time)*time_step,psi_e(1:stop_time,neg_sep)','b');
    plot(linspace(0,stop_time,stop_time)*time_step,psi_e(1:stop_time,sep_pos)','g');
    plot(linspace(0,stop_time,stop_time)*time_step,psi_e(1:stop_time,pos_pcc)','k');
    
    ylabel('Electrolyte Potential [V]');
    xlabel('Time (s)');
    
    
    legend('Electrolyte Potential at Ncc/Neg Boundary', 'Electrolyte Potential at Neg/Sep Boundary', 'Electrolyte Potential at Sep/Pos Boundary', 'Electrolyte Potential at Pos/Pcc Boundary','Location','best');
    
    grid on;
    
    %psi_s at different points
    figure(24);
    hold on;
    plot(linspace(0,stop_time,stop_time)*time_step,psi_s(1:stop_time,ncc_neg)','r');
    plot(linspace(0,stop_time,stop_time)*time_step,psi_s(1:stop_time,neg_sep)','b');
    plot(linspace(0,stop_time,stop_time)*time_step,psi_s(1:stop_time,sep_pos)','g');
    plot(linspace(0,stop_time,stop_time)*time_step,psi_s(1:stop_time,pos_pcc)','k');
    
    ylabel('Solid Potential [V]');
    xlabel('Time (s)');
    
    
    legend('Solid Potential at Ncc/Neg Boundary', 'Solid Potential at Neg/Sep Boundary', 'Solid Potential at Sep/Pos Boundary', 'Solid Potential at Pos/Pcc Boundary','Location','best');
    
    grid on;
    
    
    end

end

if thermal ==7
figure(6);
hold on;
plot(linspace(0,stop_time,stop_time)*time_step,T_ave(1:stop_time));

xlabel('Time (s)');
ylabel('Temperature (K)');
end

