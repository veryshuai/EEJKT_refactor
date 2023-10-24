function summary_tables(simMoms,mm)

 %% average log #shipments
 % avg_ln_ships = agg_ln_ships/agg_ship_obs;
   fprintf('\r\n Average log shipments: %.4f\n',simMoms.avg_ln_ships); 

%% match maturation, by type 
[DegreeDistCount,exit_by_age,brooks] = match_summary(simMoms,mm);
%[DegreeDistCount,bexit,brooks] = match_summary_v1(simMoms,mm);

%  Data for Brooks table 2 include all generated matches. Data from Brooks table 1 
%  only include types with matches in both the current and the previous year. 
%  Thus the data for table 1 miss some of the marginal firm types, though the discrepancy is minor.

exit_by_age = exit_by_age(2:6,:);
MatchAge={'1-yr old','2-yr old','3-yr old','4-yr old','5+ yr old'};
    format shortG
    quantile_1 = exit_by_age(:,1); quantile_2 = exit_by_age(:,2); quantile_3 = exit_by_age(:,3);
    quantile_4 = exit_by_age(:,4);    
   Match_exit_rates = table(quantile_1,quantile_2,quantile_3,quantile_4,'RowNames',MatchAge)             

 CohortAge = {'1-yr old','2-yr old','3-yr old','4-yr old','5-yr old',...
            '6-yr old','7-yr old','8-yr old','9-yr old','10-yr old'};   
    Firm_count = brooks(:,1); Total_exports = brooks(:,2); Avg_exports = brooks(:,3);  
    % row 1 of brooks() corresponds to active matches with zero
    % shipments--I discard it because they don't show up in the data
    format shortG
    Brooks_table = table(Firm_count,Total_exports,Avg_exports,'RowNames',CohortAge)
    
     Firm_count2 = brooks(:,1)./brooks(1,1); Total_exports2 = brooks(:,2)./brooks(1,2); Avg_exports2 = brooks(:,3)./brooks(1,3);    
     Brooks_table2 = table(Firm_count2,Total_exports2,Avg_exports2,'RowNames',CohortAge)

      
%% create variables for analysis of degree distribution

     Ntrunc           = 100;
     degreeCDF       = cumsum(DegreeDistCount)'./sum(DegreeDistCount);
     ff_sim_max      = find(degreeCDF <1);
     Nmatches        = length(ff_sim_max);
     Nmatches        = min(Ntrunc,Nmatches);
     log_compCDF     = log(1 - degreeCDF(ff_sim_max(1:Nmatches)))';
     log_matches     = log(1:1:Nmatches)';
        
     figure(3);
     scatter(log_matches,log_compCDF,'filled',"blue");
     xlabel('log # matches')
     ylabel('log(1-CDF)')
     title('Buyers per Seller degree distribution')
     
%      FirmFrac = [DegreeDistCount(1:5)./sum(DegreeDistCount);...
%                  sum(DegreeDistCount(6:10))/sum(DegreeDistCount);...
%                  sum(DegreeDistCount(11:end))/sum(DegreeDistCount)];
%      CumFrac = [degreeCDF(1:5)';degreeCDF(10);degreeCDF(end)];
%                  
%      FirmFreq = {'1 match','2 match','3 matches','4 matches','5 matches',...
%                   '6-10 matches','>10 matches'};
%      DegreeDist = table(FirmFrac,CumFrac,'RowNames',FirmFreq) 
%         
