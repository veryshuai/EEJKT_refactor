function summary_tables(simMoms,mm)

 %% average log #shipments
 % avg_ln_ships = agg_ln_ships/agg_ship_obs;
   fprintf('\r\n Average log shipments: %.4f\n',simMoms.avg_ln_ships); 

%% match maturation, by type 
[bexit,brooks] = match_summary(simMoms,mm);

%  Data for Brooks table 2 include all generated matches. Data from Brooks table 1 
%  only include types with matches in both the current and the previous year. 
%  Thus the data 1 for table miss some of the marginal firm types, though the discrepanc is minor.

MatchAge={'1-yr old','2-yr old','3-yr old','4-yr old','5+ yr old'};
    format shortG
    quantile_1 = bexit(:,1); quantile_2 = bexit(:,2); quantile_3 = bexit(:,3);
    quantile_4 = bexit(:,4);    
   Match_exit_rates = table(quantile_1,quantile_2,quantile_3,quantile_4,'RowNames',MatchAge)             

 CohortAge = {'1-yr old','2-yr old','3-yr old','4-yr old','5-yr old',...
            '6-yr old','7-yr old','8-yr old','9-yr old','10-yr old'};   
    Firm_count = brooks(2:11,1); Total_exports = brooks(2:11,2); Avg_exports = brooks(2:11,3);  
    % row 1 of brooks() corresponds to active matches with zero
    % shipments--I discard it because they don't show up in the data
    format shortG
    Brooks_table = table(Firm_count,Total_exports,Avg_exports,'RowNames',CohortAge)
    
     Firm_count2 = brooks(2:11,1)./brooks(2,1); Total_exports2 = brooks(2:11,2)./brooks(2,2); Avg_exports2 = brooks(2:11,3)./brooks(2,3);    
     Brooks_table2 = table(Firm_count2,Total_exports2,Avg_exports2,'RowNames',CohortAge)

    %Plots against data
    %Separation rates
%     DataSeps = [.829,.632,.573,.550,.497;...
%                 .756,.584,.494,.468,.437;...
%                 .677,.521,.446,.408,.376;...
%                 .521,.445,.403,.392,.367]';
%     figure();
%     scatter(DataSeps,Match_exit_rates{:,:},'filled',"blue");
%     hold on
%     plot(0:1,0:1,'red');
%     hold off
%     xlabel('Data')
%     ylabel('Model')
%     title('Match Exit Rates by Initial Sales Quartile')
%     saveas(gcf,'results/pics/MatchExitRatesByInitialSalesQuartileFIT.png')
% 
%     %Brooks Table
%     dataBrooks =   [1,1,1;...
%                     0.29,1.11,3.77;...
%                     0.18,0.93,5.03;...
%                     0.14,0.67,4.66;...
%                     0.12,0.63,5.18;...  
%                     0.10,0.51,4.99;...
%                     0.08,0.50,5.72;...
%                     0.08,0.45,5.91;...
%                     0.07,0.39,5.58;...
%                     0.06,0.40,6.58];
%     figure();
%     scatter(dataBrooks,Brooks_table2{:,:},'filled',"blue");
%     hold on
%     plot(0:7,0:7,'red');
%     hold off
%     xlabel('Data')
%     ylabel('Model')
%     title('Average Aggregates by Cohort Age')
%     saveas(gcf,'results/pics/AverageAggregatesByCohortAgeFIT.png') 


      
%% create variables for analysis of degree distribution

% Degree distribution excluding duds
        ff_sim_max      = find(cumsum(simMoms.agg_match_hist)'./sum(simMoms.agg_match_hist)<1);
        log_compCDF     = log(1 - cumsum(simMoms.agg_match_hist(ff_sim_max))'./sum(simMoms.agg_match_hist));
        log_matches     = log(1:1:length(ff_sim_max))';
        xmat            = [ones(length(ff_sim_max),1),log_matches,log_matches.^2];
                 
        temp = cumsum(simMoms.agg_match_hist(ff_sim_max)./sum(simMoms.agg_match_hist(ff_sim_max)));
        ptemp = temp(1:length(ff_sim_max)) - [0,temp(1:length(ff_sim_max)-1)];

% Degree distribution including duds        
        ff_sim_maxD      = find(cumsum(simMoms.agg_match_histD)'./sum(simMoms.agg_match_histD)<1);
        log_compCDF_D    = log(1 - cumsum(simMoms.agg_match_histD(ff_sim_maxD)')./sum(simMoms.agg_match_histD));
        log_matchesD     = log(1:1:length(ff_sim_maxD))';
        xmatD            = [ones(length(ff_sim_maxD),1),log_matchesD,log_matchesD.^2];
        
        tempD = cumsum(simMoms.agg_match_histD(ff_sim_maxD)./sum(simMoms.agg_match_histD(ff_sim_maxD)));
        ptempD = tempD(1:length(ff_sim_maxD)) - [0,tempD(1:length(ff_sim_maxD)-1)];
        
% degree distribution, not counting duds       
        format short
        Category = {'1 buyer','2 buyers','3 buyers','4 buyers','5 buyers','6-10 buyers','11+ buyers'};
        model_share = [ptemp(1:5),sum(ptemp(6:10)),sum(ptemp(11:end-1))]';
        data_share = [0.792,0.112,0.031,0.016,0.009,0.022,0.016]';
        Ergodic_Dist = table(data_share,model_share,'RowNames',Category)
        
% degree distribution, counting duds        
        format short
        Category = {'1 buyer','2 buyers','3 buyers','4 buyers','5 buyers','6-10 buyers','11+ buyers'};
        model_shareD = [ptempD(1:5),sum(ptempD(6:10)),sum(ptempD(11:end-1))]';
        data_shareD = [0.792,0.112,0.031,0.016,0.009,0.022,0.016]';
        Ergodic_Dist = table(data_shareD,model_shareD,'RowNames',Category)
        
        
