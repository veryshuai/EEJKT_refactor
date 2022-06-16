
% [Data, W, Model] = read_in_and_organize_moments(simMoms);

std_errors = (diag(W)).^0.5 ;

CI = [Data' - 2*std_errors, Model, Data' + 2*std_errors, Data'];
sortCI = sortrows(CI,4);

figure(1)
hold on
scatter(sortCI(:,1),sortCI(:,4),'filled','blue')
scatter(sortCI(:,2),sortCI(:,4),'filled','black')
scatter(sortCI(:,3),sortCI(:,4),'filled','blue')
xlabel('data-based coefs')
ylabel('model-based coefs')
hold off

figure(2)
scatter(sortCI(:,2),sortCI(:,4),'filled','black')
xlabel('data-based coefs')
ylabel('model-based coefs')

SIMcoefs = cell(12,5);
    SIMcoefs{1,1} = [simMoms.match_exit_rate;simMoms.beta_match_exit(2:5)]; % [match exit rate, 1st yr. dummy, lnXf(ijt), ln(match age), ln(exporter age),mse]
    SIMcoefs{2,1} = [simMoms.ybar_match;simMoms.beta_match(2:4);simMoms.mse_match_ar1]; % [mean ln Xf(ijt), ln Xf(ijt-1), R(ijt-1), ln(exporter age)]   
    SIMcoefs{3,1} = [simMoms.b_degree]; % [intercept, slope, quadratic term]
    SIMcoefs{4,1} = [simMoms.avg_ln_ships]; % average ln(# shipments) 
    SIMcoefs{5,1} = [simMoms.ybar_hfsales;simMoms.beta_hfsales(2);simMoms.mse_hf]; % [mean dep var.,coef,MSE]  
    SIMcoefs{6,1} = [simMoms.ybar_fsales_h;simMoms.beta_fsales_h(2);simMoms.mse_h]; % [mean dep var.,coef,MSE] 
    SIMcoefs{7,1} = [simMoms.mean_ln_haz;simMoms.b_haz(2:6)]; % [mean dep. var, ln(1+a), ln(1+a)^2, ln(1+r), ln(1+r)^2, ln(1+a)*ln(1+r)] 
    SIMcoefs{8,1} = [simMoms.mkt_exit_rate;simMoms.beta_mkt_exit(2:6)]; % [mean dep. var, ln(1+a), ln(1+a)^2, ln(1+r), ln(1+r)^2, ln(1+a)*ln(1+r)]            
    SIMcoefs{9,1}  = [simMoms.mean_succ_rate;simMoms.b_succ_rate(2)]; % [mean succ rate, ln(1+meetings)]
    SIMcoefs{10,1} = [simMoms.mean_usq_succ;simMoms.b_usq_succ(2)]; % [mean dep. var, ln(1+meetings)]
    SIMcoefs{11,1} = [simMoms.avg_expt_rate]; % mean share of exports to U.S. in total sales 
    SIMcoefs{12,1} = [simMoms.share_exptr]; % fraction of firms exporting to U.S.  

    SIMcoefs{1,5} = 'match death coefs';
    SIMcoefs{2,5} = 'match AR1 coefs'  ;
    SIMcoefs{3,5} = 'match degree distribution coefs' ;
    SIMcoefs{4,5} = 'mean average shipment' ;
    SIMcoefs{5,5} = 'exports-domestic sales coefs';  
    SIMcoefs{6,5} = 'domestic AR1 coefs' ;
    SIMcoefs{7,5} = 'ln hazard coefs' ;
    SIMcoefs{8,5} = 'last match coefs';            
    SIMcoefs{9,5} = 'success rate coefs';
    SIMcoefs{10,5} = 'std error var coefs';
    SIMcoefs{11,5} = 'foreign sales share'; 
    SIMcoefs{12,5} = 'exporter fraction'; 

    lb=ones(12,1);
    ub = lb;
    for j = 1:12
       nparam = size(SIMcoefs{j,1},1);
       ub(j) = lb(j) + nparam -1; 
         if j < 12
          lb(j+1) = ub(j)+1;
         end
       SIMcoefs{j,3} = CI(lb(j):ub(j),1);  
       SIMcoefs{j,4} = CI(lb(j):ub(j),2);        
       SIMcoefs{j,2} = Data(lb(j):ub(j))';
    end
    
    for j=1:12
        x = SIMcoefs{j,2};
        y = SIMcoefs{j,1};
        name = SIMcoefs{j,5};
        figdat = sortrows([x,y],1);

        figure(3)
        subplot(4,3,j)
        scatter(figdat(:,1),figdat(:,2))
        title(name);
        xlabel('data-based coefs')
        ylabel('model-based coefs')

    end
         