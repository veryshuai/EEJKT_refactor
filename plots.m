function plots(simMoms,Data,Model)

% nonparametric plot of degree distribution 
        figure(1)
        scatter(simMoms.log_matches,simMoms.log_compCDF)      

% plot histogram of frequencies for meeting hazards
        %figure(2)
        %agg_time_gaps = agg_time_gaps(2:size(agg_time_gaps,1),:);
        %histogram(agg_time_gaps(:,3))
        
% plot the meeting hazard and mkt. exit hazard as a functions of # clients and success rate    

        grdsize   = 40;
        b_mkt_exitDAT  = [0.8397,-0.6290,0.1205,0.0211,0.5976,-0.1290]; % from last_matchDAT (and RDC release)
%       b_hazDAT   = [-0.3258+log(12),-0.8181,0.3117,-1.1323,2.4514,-0.7082]; % from match_lag_coefsDAT (and RDC release)
        b_hazDAT   = [-0.3258,-0.8181,0.3117,-1.1323,2.4514,-0.7082]; % from match_lag_coefsDAT (and RDC release)
        
        succ_grid = (min(simMoms.ln_succ_rate):((max(simMoms.ln_succ_rate)-min(simMoms.ln_succ_rate))/grdsize):max(simMoms.ln_succ_rate))';
        csucc_grid = (min(simMoms.ln_csucc):((max(simMoms.ln_csucc)-min(simMoms.ln_csucc))/grdsize):max(simMoms.ln_csucc))';
        
        dat_mkt_exit = zeros(grdsize+1,grdsize+1);
        dat_haz      = zeros(grdsize+1,grdsize+1);
        mod_mkt_exit = zeros(grdsize+1,grdsize+1);
        mod_haz      = zeros(grdsize+1,grdsize+1);
        for i = 1:grdsize+1  % cum success 
            for j = 1: grdsize+1  % success rate 
                dat_mkt_exit(i,j) = b_mkt_exitDAT(1) + b_mkt_exitDAT(2)*csucc_grid(i) ...
                   +b_mkt_exitDAT(3)*csucc_grid(i)^2 + b_mkt_exitDAT(4)*succ_grid(j)...
                    +b_mkt_exitDAT(5)*succ_grid(j)^2 + b_mkt_exitDAT(6)*csucc_grid(i)*succ_grid(j)   ;
                
                mod_mkt_exit(i,j) = simMoms.beta_mkt_exit(1) + simMoms.beta_mkt_exit(2)*csucc_grid(i) ...
                   + simMoms.beta_mkt_exit(3)*csucc_grid(i)^2 + simMoms.beta_mkt_exit(4)*succ_grid(j)...
                    +simMoms.beta_mkt_exit(5)*succ_grid(j)^2 + simMoms.beta_mkt_exit(6)*csucc_grid(i)*succ_grid(j)   ;
                
                dat_haz(i,j) = b_hazDAT(1) + b_hazDAT(2)*csucc_grid(i) ...
                   +b_hazDAT(3)*csucc_grid(i)^2 + b_hazDAT(4)*succ_grid(j)...
                   +b_hazDAT(5)*succ_grid(j)^2  + b_hazDAT(6)*csucc_grid(i)*succ_grid(j)   ;
                
                mod_haz(i,j) =  simMoms.b_haz(1) + simMoms.b_haz(2)*csucc_grid(i) ...
                   +simMoms.b_haz(3)*csucc_grid(i)^2 + simMoms.b_haz(4)*succ_grid(j)...
                   +simMoms.b_haz(5)*succ_grid(j)^2  + simMoms.b_haz(6)*csucc_grid(i)*succ_grid(j)   ;
            end
        end
      

%    x_axis =mm.theta2';
%    y_axis = mm.Phi;
%    surf(x_axis,y_axis, xptr_count)
              

figure(2)
hist3([simMoms.ln_succ_rate,simMoms.ln_csucc],[grdsize+1,grdsize+1]); 
        title('frequency of states--simulated data')
        ylabel('log(cum. success)')
        xlabel('log(1+success rate)')
        zlabel('frequency')        


figure(3)
subplot(2,1,1)

nfreq = hist3([simMoms.ln_csucc,simMoms.succ_rate],[grdsize+1,grdsize+1]) > 0; 
% note: need to switch order of variables to make compatible with dat_haz &
% mod_haz
tdat_haz = dat_haz.*nfreq;
        surf(succ_grid,csucc_grid,tdat_haz)
        title('data-based log match hazard')
        ylabel('log(cum. success)')
        xlabel('log(1+success rate)')
        zlabel('log match hazard')
subplot(2,1,2)        
tmod_haz = mod_haz.*nfreq;        
        surf(succ_grid,csucc_grid,tmod_haz)
        title('model-based log match hazard')
        ylabel('log(cum. success)')
        xlabel('log(1+success rate)')
        zlabel('log match hazard')
        
figure(4)
subplot(2,1,1)

nfreq = hist3([simMoms.ln_csucc,simMoms.succ_rate],[grdsize+1,grdsize+1]) > 0; 
% note: need to switch order of variables to make compatible with dat_haz &
% mod_haz
tdat_haz = dat_haz;
        surf(succ_grid,csucc_grid,tdat_haz)
        title('data-based log match hazard')
        ylabel('log(cum. success)')
        xlabel('log(1+success rate)')
        zlabel('log match hazard')
subplot(2,1,2)        
tmod_haz = mod_haz;        
        surf(succ_grid,csucc_grid,tmod_haz)
        title('model-based log match hazard')
        ylabel('log(cum. success)')
        xlabel('log(1+success rate)')
        zlabel('log match hazard')
        
% figure(6)
% subplot(2,1,1)
% 
%        contour(succ_grid,csucc_grid,tdat_haz)
%         title('data-based log match hazard')
%         ylabel('log(cum. success)')
%         xlabel('log(1+success rate)')
%         zlabel('log match hazard')
% subplot(2,1,2)        
% tmod_haz = mod_haz.*nfreq;        
%         contour(succ_grid,csucc_grid,tmod_haz)
%         title('model-based log match hazard')
%         ylabel('log(cum. success)')
%         xlabel('log(1+success rate)')
%         zlabel('log match hazard')
        

%         tdat_mkt_exit = dat_mkt_exit.*nfreq;
%         tmod_mkt_exit = mod_mkt_exit.*nfreq;
%         figure(5)        
%         subplot(2,1,1)
%         surf(succ_grid,csucc_grid,tdat_mkt_exit)
%         title('data-based mkt. exit')
%         ylabel('cum. success')
%         xlabel('success rate')
%         zlabel('exit probability')
%         
%         subplot(2,1,2)
%         surf(succ_grid,csucc_grid,tmod_mkt_exit)
%         title('model-based mkt. exit')
%         ylabel('cum. success')
%         xlabel('success rate')
%         zlabel('exit probability')
        
% Generate moments vs data plots

% % Table 6: Match hazards, success rates, and endurance
% figure();
% scatter([Data(1:5)';Data(21:26)';Data(33:36)'],[Model(1:5);Model(21:26);Model(33:36)],'filled');
% hold on
% plot(-4:3,-4:3),'red';
% hold off
% xlabel('Data')
% ylabel('Model')
% title('Match Hazards, Success Rates, and Endurance')
% saveas(gcf,'results/pics/MatchHazardsSuccessRatesAndEnduranceFIT.png')
% 
% % Table 7: Client distribution and shipment frequency
% figure();
% scatter(Data(12:14)',Model(12:14),'filled');
% hold on
% plot(-2:1,-2:1),'red';
% hold off
% xlabel('Data')
% ylabel('Model')
% title('Client Distribution and Shipment Frequency')
% saveas(gcf,'results/pics/ClientDistributionAndShipmentFrequencyFIT.png')
% 
% % Table 8: Home and Foreign Sales Regressions
% figure();
% scatter([Data(6:10)';Data(16);Data(20);Data(37:38)'],[Model(6:10);Model(16);Model(20);Model(37:38)],'filled');
% hold on
% plot(-1:12,-1:12),'red';
% hold off
% xlabel('Data')
% ylabel('Model')
% title('Home and Foreign Sales Regressions')
% saveas(gcf,'results/pics/HomeAndForeignSalesRegressionsFIT.png')
%          
% % Table 10: Cohort evolutions
% 
% 
% % Table 11: Exporter client degree dist
% dataClientDegree = [0.778 0.116 0.043 0.021 0.011]';
% modelClientDegree = exp(simMoms.log_compCDF(1:5));
% figure();
% scatter(dataClientDegree,modelClientDegree,'filled');
% hold on
% plot(0:1,0:1,'red');
% hold off
% xlabel('Data')
% ylabel('Model')
% title('Exporter Client Degree Distribution')
% saveas(gcf,'results/pics/ExporterClientDegreeDistributionFIT.png')

        
