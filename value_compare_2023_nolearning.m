% %This script calculates values of clients
% clear all;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % First order of business is figuring
% % out typical productivities and success
% % probabilities for "high revenue" and
% % "low revenue" firm sales percentiles
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
rng(80085);
% 
load results/val_dat_nolearning

%create data set containing only variables we need (no information about counterparty)
% all_exporters = [agg_type,firmID,annual_sales,prod_type,theta_type]
%  match_recs: [period, type, firm_ID, sales, shipments, boy Z, eoy Z, match age, firm age]       
match_recs = succ_matches;
positive_sales = match_recs(:,4) > 0; %matches which die during the year sometimes have zero sales
all_matches_one_year = match_recs(match_recs(:,1) >= 588 & match_recs(:,1) < 600 & positive_sales,1:4);
all_matches_one_year(:,3) = all_matches_one_year(:,3) * 2; %some firm IDs are X.5, accumarray does not like this
[unique_firms,uniq_ind,~] = unique(all_matches_one_year(:,3));
annual_sales = accumarray(all_matches_one_year(:,3),all_matches_one_year(:,4));
annual_sales = annual_sales(annual_sales(:,1)>0); %delete zero row for ids which are missing in all_matches_one_year(:,3)
all_exporters = [unique_firms,all_matches_one_year(uniq_ind,2),annual_sales];
all_exporters(:,4:5) = mm.pt_type(all_exporters(:,2),:);
all_exporters = sortrows(all_exporters,3); %sort based on sales

median_prod = [10,floor(prctile(all_exporters(:,4),10));50,floor(prctile(all_exporters(:,4),50));90,floor(prctile(all_exporters(:,4),90))]
median_succ = [10,floor(prctile(all_exporters(:,5),10));50,floor(prctile(all_exporters(:,5),50));90,floor(prctile(all_exporters(:,5),90))]
%DJ NOTE: median_succ does not seem right -- should not be 1.
median_prod = [10,15;50,16;90,17];
median_succ = [10,6;50,6;90,6];

%Generate some figures used in the text
%average sales per shipment
avg_shipment_sales = mean(match_recs(match_recs(:,5)~=0,4)./match_recs(match_recs(:,5)~=0,5))
sd_shipment_sales = (var(match_recs(match_recs(:,5)~=0,4)./match_recs(match_recs(:,5)~=0,5)))^0.5

%Cost function calculations in the paper
cost_one_match_per_yr_no_exp = mm.cost_f(1/12,1)
cost_one_match_two_yrs_no_exp = mm.cost_f(1/24,1)
cost_one_match_per_yr_ten_exp = mm.cost_f(1/12,11)
ratio_cost_one_vs_zero = (mm.cost_f(1/12,11)-mm.cost_f(1/12,1))/mm.cost_f(1/12,1)

%theta stuff
theta_exp = mm.af / (mm.af + mm.bf)
theta_var = mm.af * mm.bf /((mm.af + mm.bf)^2 *(mm.af + mm.bf + 1))
theta_sd = sqrt(theta_var)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Now we have both the median productivity and success probabilities in hand
%simulate value for these different types and plot


% X = [-3.87377400411704,-19.6350579745564,0.141131266607483,0.224048944147957,...
%          0.527541658446327,12.1673847502223,0.0462420845261331,5.13239605341601,...
%          2.38418049835969,15.1614775287665]; % fit PC 11.6912  unix: 11.6935

X = [-7.00428657247973 -23.0915400600166 0.129855165074604 0.147073857937515 0.637554188362434 12.5050138543584 0.0607676627106075 4.33149453421485 2.94198673283020 18.0202455396844]; %BEST FIT NO LEARNING

% %USE FOREIGN PARAMETERS FOR HOME
% %This way we can use value functions to learn about the importance of learning vs network
% X(8) = X(10);
%     
mm = setModelParameters(X);
mm.th_ind_temp = 5;
mm.v_tolerance   = 1e-6; %Need tighter tolerance to get nice looking value plots
policy = generatePolicyAndValueFunctions(mm);

% Key to value functions
% value_f (succ, trial, common succ rate (defunct), network size, prod of firm, F macro shock) 
% value_h (common succ rate (defunct), known succ rate, network size, prod of firm, H macro shock)

%marginal value of another successful foreign client
marg_val_succ_f = zeros(20,3);
val_succ_f = zeros(20,3);
marg_val_succ_f_percent = zeros(20,3);

%marginal value of another unsuccessful foreign client
marg_val_fail_f = zeros(20,3);
val_fail_f = zeros(20,3);
marg_val_fail_f_percent = zeros(20,3);

%marginal value of another successful home client
marg_val_succ_h = zeros(20,3);
val_succ_h = zeros(20,3);
marg_val_succ_h_percent = zeros(20,3);

%marginal value of alternating success and failures for foreign client
marg_val_alt_f = zeros(20,3);
val_alt_f = zeros(20,3);
marg_val_alt_f_percent = zeros(20,3);

%median_succ(:,2) = 7; %only matters for home anyway!  Set to maximum

for type = 1:3

    for k = 1:20
        val_succ_f(k,type) = policy.value_f(k,k,1,k,median_prod(type,2),7); 
        marg_val_succ_f(k,type)  = policy.value_f(k+1,k+1,1,k+1,median_prod(type,2),7) - policy.value_f(k,k,1,k,median_prod(type,2),7);
        marg_val_succ_f_percent(k,type) = (policy.value_f(k+1,k+1,1,k+1,median_prod(type,2),7) - policy.value_f(k,k,1,k,median_prod(type,2),7))/policy.value_f(k,k,1,k,median_prod(type,2),7);
    end
    
    for k = 1:20
        val_fail_f(k,type) = policy.value_f(1,k,1,1,median_prod(type,2),7);
        marg_val_fail_f(k,type)  = policy.value_f(1,k+1,1,1,median_prod(type,2),7) - policy.value_f(1,k,1,1,median_prod(type,2),7);
        marg_val_fail_f_percent(k,type) = (policy.value_f(1,k+1,1,1,median_prod(type,2),7) - policy.value_f(1,k,1,1,median_prod(type,2),7))/policy.value_f(1,k,1,1,median_prod(type,2),7);
    end
    
    for k = 1:20
        val_succ_h(k,type) = policy.value_h(1,7,k,median_prod(type,2),7);
        marg_val_succ_h(k,type)  = policy.value_h(1,7,k+1,median_prod(type,2),7) - policy.value_h(1,7,k,median_prod(type,2),7);
        marg_val_succ_h_percent(k,type) = (policy.value_h(1,7,k+1,median_prod(type,2),7) - policy.value_h(1,7,k,median_prod(type,2),7))/policy.value_h(1,7,k,median_prod(type,2),7);
    end

    for k = 1:20
        val_alt_h(k,type) = policy.value_h(1,7,k,median_prod(type,2),7);
        marg_val_alt_h(k,type)  = policy.value_h(1,7,k+1,median_prod(type,2),7) - policy.value_h(1,7,k,median_prod(type,2),7);
        marg_val_alt_h_percent(k,type) = (policy.value_h(1,7,k+1,median_prod(type,2),7) - policy.value_h(1,7,k,median_prod(type,2),7))/policy.value_h(1,7,k,median_prod(type,2),7);
    end

    for k = 1:20
        succs = floor(k/2) + 1;
        val_alt_f(k,type) = policy.value_f(succs,k,1,succs,median_prod(type,2),7); 
        marg_val_alt_f(k,type)  = policy.value_f(succs-mod(k,2)+1,k+1,1,succs-mod(k,2)+1,median_prod(type,2),7) - policy.value_f(succs,k,1,succs,median_prod(type,2),7);
        marg_val_alt_f_percent(k,type) = (policy.value_f(succs-mod(k,2)+1,k+1,1,succs-mod(k,2)+1,median_prod(type,2),7) - policy.value_f(succs,k,1,succs,median_prod(type,2),7))/policy.value_f(succs,k,1,succs,median_prod(type,2),7);
    end

end

%Plots (TO CREATE THIS PLOT, CHANGE setModelParameters.m:  
%mm.v_tolerance   = 1e-6; % convergence tolerance, value function iterations (WAS 1e-3)
%Important because due to numerical error, there is a slight upward drift in the value of failure.  This is (nearly) eliminated with a lower v_tolerance.
plot(log(val_succ_f(:,3)));
hold on
plot(log(val_alt_f(:,3)));
hold on
plot(log(val_fail_f(:,3)));
xlabel('Matches')
ylabel('1992 USD')
title({'Log value in foreign market','Excludes profit flow from current relationships'})
hold off
saveas(gcf,"results/value_plots/val_f_three_types_no_learning.png");

%Plots
x_vals = 0:size(val_succ_f,1)-1;
plot(x_vals,log(val_succ_f(:,3)));
hold on
plot(x_vals,log(val_alt_f(:,3)));
hold on
plot(x_vals,log(val_fail_f(:,3)));
xlabel('Matches')
ylabel('2023 USD (log scale)')
title({'Value in foreign market','Excludes profit flow from current relationships'})
hold off
ax = gca();
ax.YTickLabel = compose('%0.2g', exp(ax.YTick)');
saveas(gcf,"results/value_plots/val_f_no_learning.png");

plot(marg_val_succ_f(:,2));
hold on
plot(marg_val_alt_f(:,2));
hold on
plot(marg_val_fail_f(:,2));
xlabel('Matches')
ylabel('1992 USD')
title({'Marginal client value in foreign market','Excludes additional profit flow'})
hold off
saveas(gcf,"results/value_plots/marg_val_f_no_learning.png");

plot(marg_val_succ_f_percent(:,2));
hold on
plot(marg_val_alt_f_percent(:,2));
hold on
plot(marg_val_fail_f_percent(:,2));
xlabel('Matches')
ylabel('1992 USD')
title({'Marginal client value in foreign market','fraction of total value'})
hold off
saveas(gcf,"results/value_plots/marg_val_f_percent_no_learning.png");

%Plots
plot(val_alt_h(:,2));
xlabel('Matches')
ylabel('1992 USD')
title({'No learning value in foreign market','Excludes profit flow from current relationships'})
hold on
plot(val_alt_f(:,2));
hold off
saveas(gcf,"results/value_plots/val_h_no_learning.png");

plot(marg_val_alt_h(:,2));
xlabel('Matches')
ylabel('1992 USD')
title({'No learning marginal client value in foreign market','Excludes additional profit flow'})
saveas(gcf,"results/value_plots/marg_val_h_no_learning.png");

plot(marg_val_alt_h_percent(:,2));
xlabel('Matches')
ylabel('1992 USD')
title({'No learning marginal client value in foreign market','fraction of total value'})
hold off
saveas(gcf,"results/value_plots/marg_val_h_percent_no_learning.png");

%%% Time to learn

% In this exercise, we consider how long it will take identical firms to
% reach N matches, conditional on the arrival of successes.  Let us assume
% that the firm will get exactly N/2 successes, the question is when the
% successes arrive.  Just so we have a simple ordering, we assume that the
% successes all arrive consecutively, but start in different years.

cum_year_mat = zeros(3,1);
theta_evolution = zeros(6,1,3); % [match_no, [time, value], first_yr]
for first_succ_yr = 1:3

    succ_seq = zeros(6,1);
    succ_seq(first_succ_yr:first_succ_yr+2) = 1;
    
    cum_years = 0;
    cum_years_lag = 0;
    succs = 1; %first index is zero
    trials = 1; %first index is zero
    for match_no = 0:5
        cum_years_lag = cum_years;
        theta_guess = (mm.af + succs - 1) / (mm.af + mm.bf + trials - 1);
        theta_evolution(match_no+1,1,first_succ_yr) = cum_years;
        theta_evolution(match_no+1,2,first_succ_yr) = theta_guess;
        cum_years = cum_years_lag + 1 /(mm.pd_per_yr * policy.lambda_f(succs,trials,5,succs,median_prod(3,2),7));
        trials = trials + 1;
        succs = succs + succ_seq(match_no+1);
    end
    
    cum_year_mat(first_succ_yr) = cum_years_lag;
end

bar(cum_year_mat);
xlabel('Year of first success');
ylabel('Years to five trials');
ylim([0 40]);
title('three consecutive success in five trials');
saveas(gcf,"results/value_plots/success_order_no_learning.png");

plot(theta_evolution(:,1,3),theta_evolution(:,2,3));
hold on
plot(theta_evolution(:,1,1),theta_evolution(:,2,1));
xlabel('Years')
ylabel('Success probability belief')
title('Effect of Early Discouragement');
legend({'Success last','Success first'},'Location','northeast')
hold off
saveas(gcf,"results/value_plots/success_beliefs_no_learning.png");