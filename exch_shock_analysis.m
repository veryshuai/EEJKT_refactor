%clear;
rng(80085);

load results/val_dat
%load lambdas_temp

files = {'results/exch_shock_plots/baseline_no_shk', ...
    'results/exch_shock_plots/baseline_up_shk', ...
    'results/exch_shock_plots/baseline_down_shk'};

for pol_k = 1:3
    match_recs_appended = [];
    for m=1:200
        filename = sprintf('%d.mat', m);
        load([files{pol_k},filename], 'match_recs');
        %  match_recs: [period, type, firm_ID, sales, shipments, boy Z, eoy Z, match age, firm age]       
        match_recs(:,3) = match_recs(:,3) + 10000 * m; %make firm ids unique across files
        match_recs_appended = [match_recs_appended;match_recs];
    end
    match_recs = match_recs_appended;

    match_recs(:,10) = match_recs(:,1)/mm.pd_per_yr; %years rather than months

    for t=10:max(match_recs(:,10)-1)

        match_recs_one_year = match_recs(match_recs(:,10) == t,:);

        new_id = 0.5*(match_recs_one_year(:,2)+match_recs_one_year(:,3)).*(match_recs_one_year(:,2)+match_recs_one_year(:,3)+1)+match_recs_one_year(:,3);
        [unique_firms,~,~] = unique(new_id);
        total_firms(t,pol_k) = size(unique_firms,1);
        only_new_firms = match_recs_one_year(:,9)/mm.pd_per_yr < t-25;
        only_first_yr_firms = match_recs_one_year(:,9)/mm.pd_per_yr < 1;
        [new_firms,~,~] = unique(new_id(only_new_firms));
        [first_yr_firms,~,~] = unique(new_id(only_first_yr_firms));
        [non_first_yr_firms,~,~] = unique(new_id(~only_first_yr_firms));
        total_new_firms(t,pol_k) = size(new_firms,1);
        total_first_yr_firms(t,pol_k) = size(first_yr_firms,1);
        total_non_first_yr_firms(t,pol_k) = size(non_first_yr_firms,1);

        total_matches(t,pol_k) = size(match_recs_one_year,1);
        only_new_matches = match_recs_one_year(:,8)/mm.pd_per_yr < t-25;
        only_first_yr_matches = match_recs_one_year(:,8)/mm.pd_per_yr < 1;
        total_new_matches(t,pol_k) = size(match_recs_one_year(only_new_firms,:),1);
        total_old_firm_new_matches(t,pol_k) = size(match_recs_one_year(~only_new_firms & only_new_matches,:),1);
        total_old_firm_old_matches(t,pol_k) = size(match_recs_one_year(~only_new_firms & ~only_new_matches,:),1);
        total_first_yr_matches(t,pol_k) = size(match_recs_one_year(only_first_yr_matches,:),1);
        total_non_first_yr_matches(t,pol_k) = size(match_recs_one_year(~only_first_yr_matches,:),1);

        %total_matches_old_firms = size(all_matches_one_year(
        haircut = 1;
        if t>25 && pol_k == 2
            haircut = 6/5;
        end
        if t>25 && pol_k == 3
            haircut = 4/5;
        end
        total_sales(t,pol_k) = sum(haircut * match_recs_one_year(:,4),1);
        total_new_sales(t,pol_k) = sum(haircut * match_recs_one_year(only_new_firms,4),1);
        total_old_firm_new_sales(t,pol_k) = sum(haircut * match_recs_one_year(~only_new_firms & only_new_matches,4),1);
        total_old_firm_old_sales(t,pol_k) = sum(haircut * match_recs_one_year(~only_new_firms & ~only_new_matches,4),1);
        total_first_yr_match_sales(t,pol_k) = sum(haircut * match_recs_one_year(only_first_yr_matches,4),1);
        total_non_first_yr_match_sales(t,pol_k) = sum(haircut * match_recs_one_year(~only_first_yr_matches,4),1);
        total_first_yr_firm_sales(t,pol_k) = sum(haircut * match_recs_one_year(only_first_yr_firms,4),1);
        total_non_first_yr_firm_sales(t,pol_k) = sum(haircut * match_recs_one_year(~only_first_yr_firms,4),1);

    end
    sales_per_firm(:,pol_k)   = total_sales(:,pol_k) ./ total_firms(:,pol_k);
    sales_per_match(:,pol_k)  = total_sales(:,pol_k) ./ total_matches(:,pol_k);
    matches_per_firm(:,pol_k) = total_matches(:,pol_k) ./ total_firms(:,pol_k);
end


% %Value increase after shock (now reported in paper).  In order to be
% %consistent with other results, we need to add the discounted expected
% %value of current clients, which we don't have here.  We take the number from the
% %spreadsheet "value_of_experience_v2.xlsx", which is ultimately taken from EEJKT_baseline_mean_firm_values/mean_firm_match_value.m
% %We assume that these mechanically rise 20% after the exchange rate shock.
% mean_value_before_shock = mean(mean(squeeze(value_per_firm(25:50,:,2)),2)) + 1704000
% median_value_before_shock = mean(mean(squeeze(value_per_firm_med(25:50,:,2)),2)) + 289695
% immediate_m ean_value_jump_depreciation = (mean(squeeze(value_per_firm(51,:,2)),2) + 1.2 * 1704000) / (mean(squeeze(value_per_firm(50,:,2)),2) + 1704000)
% value_net_of_discounted_sales_current_clients = mean(squeeze(value_per_firm(25:50,2)),2)
%
% %Value added per firm
% hold off
% plot(mean(squeeze(value_per_firm_va(11:40,1)),2),'LineWidth',2)
% hold on
% plot(mean(squeeze(value_per_firm_va(11:40,2)),2),'LineWidth',2)
% hold on
% plot(mean(squeeze(value_per_firm_va(11:40,3)),2),'LineWidth',2)
% set(gca,'FontSize',18)
% legend({'Baseline','Favorable Exch shock','Unfavorable Exch shock'},'location','NorthWest','FontSize',20)
% title("Mean value of experience across simulated firms");
% xlim([1,50])
% xticks([6 11 16 21 26 31 36 41 46])
% xticklabels({'-20','-15','-10','-5','0','5','10','15','20'})
% hold off
% saveas(gcf,'results/exch_shock_plots/decomp_plots/intangible_value_per_firm.png')
%
% %Value added per firm median
% hold off
% plot(mean(squeeze(value_per_firm_va_med(11:40,1)),2),'LineWidth',2)
% hold on
% plot(mean(squeeze(value_per_firm_va_med(11:40,2)),2),'LineWidth',2)
% hold on
% plot(mean(squeeze(value_per_firm_va_med(11:40,3)),2),'LineWidth',2)
% set(gca,'FontSize',18)
% legend({'Baseline','Favorable Exch shock','Unfavorable Exch shock'},'location','NorthWest','FontSize',20)
% title("Median value of experience across simulated firms");
% xlim([1,50])
% xticks([6 11 16 21 26 31 36 41 46])
% xticklabels({'-20','-15','-10','-5','0','5','10','15','20'})
% hold off
% saveas(gcf,'results/exch_shock_plots/decomp_plots/intangible_value_per_firm_median.png')
%
% %Value added per firm percent total
% hold off
% plot((mean(squeeze(value_per_firm_pct(11:40,1)),2)),'LineWidth',2)
% hold on
% plot((mean(squeeze(value_per_firm_pct(11:40,2)),2)),'LineWidth',2)
% hold on
% plot((mean(squeeze(value_per_firm_pct(11:40,3)),2)),'LineWidth',2)
% ylim([0.0 1.0])
% xlim([1,50])
% xticks([5 10 15 20 25 30 35 40 45 50])
% xticklabels({'-20','-15','-10','-5','0','5','10','15','20','25'})
% set(gca,'FontSize',18)
% legend({'Baseline','Favorable Exch shock','Unfavorable Exch shock'},'location','NorthWest','FontSize',20)
% title("Value of experience as fraction of total value across simulated firms");
% hold off
% saveas(gcf,'results/exch_shock_plots/decomp_plots/intangible_value_per_firm_pct.png')

%Matches
hold off
area([mean(squeeze(total_new_matches(11:40,2)),2),mean(squeeze(total_old_firm_new_matches(11:40,2)),2),mean(squeeze(total_old_firm_old_matches(11:40,2)),2)])
hold on
plot(mean(squeeze(total_new_matches(11:40,1)),2),'LineWidth',2)
plot(mean(squeeze(total_old_firm_new_matches(11:40,1)),2) + mean(squeeze(total_new_matches(11:40,1)),2),'LineWidth',2)
plot(mean(squeeze(total_old_firm_old_matches(11:40,1)),2) + mean(squeeze(total_old_firm_new_matches(11:40,1)),2) + mean(squeeze(total_new_matches(11:40,1)),2),'LineWidth',2)
set(gca,'FontSize',18)
legend({'New exporter','Incumbent exporter, new importer','Incumbent exporter and importer'},'location','SouthWest','FontSize',20)
title("Total Matches");
xlim([1,30])
xticks([5 10 15 20 25 30])
xticklabels({'-10','-5','0','5','10','15'})
hold off
saveas(gcf,'results/exch_shock_plots/decomp_plots/total_matches.png')

%Matches
hold off
area([mean(squeeze(total_new_matches(11:40,3)),2),mean(squeeze(total_old_firm_new_matches(11:40,3)),2),mean(squeeze(total_old_firm_old_matches(11:40,3)),2)])
hold on
plot(mean(squeeze(total_new_matches(11:40,1)),2),'LineWidth',2)
plot(mean(squeeze(total_old_firm_new_matches(11:40,1)),2) + mean(squeeze(total_new_matches(11:40,1)),2),'LineWidth',2)
plot(mean(squeeze(total_old_firm_old_matches(11:40,1)),2) + mean(squeeze(total_old_firm_new_matches(11:40,1)),2) + mean(squeeze(total_new_matches(11:40,1)),2),'LineWidth',2)
%ylim([0 120])
set(gca,'FontSize',18)
legend({'New exporter','Incumbent exporter, new importer','Incumbent exporter and importer'},'location','SouthWest','FontSize',20)
title("Total Matches");
xlim([1,30])
xticks([5 10 15 20 25 30])
xticklabels({'-10','-5','0','5','10','15'})
hold off
saveas(gcf,'results/exch_shock_plots/decomp_plots/total_matches_unf.png')


A = [mean(total_new_matches(11:40,2),2),mean(total_old_firm_new_matches(11:40,2),2),mean(total_old_firm_old_matches(11:40,2),2)];
area(A./repmat(sum(A,2),1,size(A,2)));
hold on
total_sum = mean(squeeze(total_old_firm_old_matches(11:40,1)),2) + mean(squeeze(total_old_firm_new_matches(11:40,1)),2) + mean(squeeze(total_new_matches(11:40,1)),2);
plot(mean(squeeze(total_new_matches(11:40,1)),2)./total_sum,'LineWidth',2)
plot((mean(squeeze(total_old_firm_new_matches(11:40,1)),2) + mean(squeeze(total_new_matches(11:40,1)),2))./total_sum,'LineWidth',2)
set(gca,'FontSize',18)
legend({'New exporter','Incumbent exporter, new importer','Incumbent exporter and importer'},'location','SouthWest','FontSize',20)
title("Total Matches");
ylim([0.0 1.0])
xlim([1,30])
xticks([5 10 15 20 25 30])
xticklabels({'-10','-5','0','5','10','15'})
hold off
saveas(gcf,'results/exch_shock_plots/decomp_plots/total_matches_pct.png')

%gaps
total_sum_shock = mean(squeeze(total_old_firm_old_matches(11:40,2)),2) + mean(squeeze(total_old_firm_new_matches(11:40,2)),2) + mean(squeeze(total_new_matches(11:40,2)),2);
gap_incum_new_old_matches = mean(squeeze(total_new_matches(11:40,2)),2)./total_sum_shock - mean(squeeze(total_new_matches(11:40,1)),2)./total_sum;
gap_incum_entrant_matches = (mean(squeeze(total_old_firm_new_matches(11:40,2)),2) + mean(squeeze(total_new_matches(11:40,2)),2))./total_sum_shock - (mean(squeeze(total_old_firm_new_matches(11:40,1)),2) + mean(squeeze(total_new_matches(11:40,1)),2))./total_sum;

A = [mean(total_new_matches(11:40,3),2),mean(total_old_firm_new_matches(11:40,3),2),mean(total_old_firm_old_matches(11:40,3),2)];
area(A./repmat(sum(A,2),1,size(A,2)));
hold on
total_sum = mean(squeeze(total_old_firm_old_matches(11:40,1)),2) + mean(squeeze(total_old_firm_new_matches(11:40,1)),2) + mean(squeeze(total_new_matches(11:40,1)),2);
plot(mean(squeeze(total_new_matches(11:40,1)),2)./total_sum,'LineWidth',2)
plot((mean(squeeze(total_old_firm_new_matches(11:40,1)),2) + mean(squeeze(total_new_matches(11:40,1)),2))./total_sum,'LineWidth',2)
set(gca,'FontSize',18)
legend({'New exporter','Incumbent exporter, new importer','Incumbent exporter and importer'},'location','SouthWest','FontSize',20)
title("Total Matches");
ylim([0.0 1.0])
xlim([1,30])
xticks([5 10 15 20 25 30])
xticklabels({'-10','-5','0','5','10','15'})
hold off
saveas(gcf,'results/exch_shock_plots/decomp_plots/total_matches_pct_unf.png')

%Sales
hold off
area([mean(squeeze(total_new_sales(11:40,2)),2),mean(squeeze(total_old_firm_new_sales(11:40,2)),2),mean(squeeze(total_old_firm_old_sales(11:40,2)),2)])
hold on
plot(mean(squeeze(total_new_sales(11:40,1)),2),'LineWidth',2)
plot(mean(squeeze(total_old_firm_new_sales(11:40,1)),2) + mean(squeeze(total_new_sales(11:40,1)),2),'LineWidth',2)
plot(mean(squeeze(total_old_firm_old_sales(11:40,1)),2) + mean(squeeze(total_old_firm_new_sales(11:40,1)),2) + mean(squeeze(total_new_sales(11:40,1)),2),'LineWidth',2)
set(gca,'FontSize',18)
%legend({'New exporter','Incumbent exporter, new importer','Incumbent exporter and importer'},'location','SouthWest','FontSize',20)
%title("Total sales");
xlim([1,30])
xticks([5 10 15 20 25 30])
xticklabels({'-10','-5','0','5','10','15'})
hold off
saveas(gcf,'results/exch_shock_plots/decomp_plots/total_sales.png')

hold off
area([mean(squeeze(total_new_sales(11:40,3)),2),mean(squeeze(total_old_firm_new_sales(11:40,3)),2),mean(squeeze(total_old_firm_old_sales(11:40,3)),2)])
hold on
plot(mean(squeeze(total_new_sales(11:40,1)),2),'LineWidth',2)
plot(mean(squeeze(total_old_firm_new_sales(11:40,1)),2) + mean(squeeze(total_new_sales(11:40,1)),2),'LineWidth',2)
plot(mean(squeeze(total_old_firm_old_sales(11:40,1)),2) + mean(squeeze(total_old_firm_new_sales(11:40,1)),2) + mean(squeeze(total_new_sales(11:40,1)),2),'LineWidth',2)
%ylim([0.0 2e7])
set(gca,'FontSize',18)
legend({'New exporter','Incumbent exporter, new importer','Incumbent exporter and importer'},'location','SouthWest','FontSize',20)
title("Total sales");
xlim([1,30])
xticks([5 10 15 20 25 30])
xticklabels({'-10','-5','0','5','10','15'})
hold off
saveas(gcf,'results/exch_shock_plots/decomp_plots/total_sales_unf.png')

A = [mean(total_new_sales(11:40,2),2),mean(total_old_firm_new_sales(11:40,2),2),mean(total_old_firm_old_sales(11:40,2),2)];
area(A./repmat(sum(A,2),1,size(A,2)));
hold on
total_sum = mean(squeeze(total_old_firm_old_sales(11:40,1)),2) + mean(squeeze(total_old_firm_new_sales(11:40,1)),2) + mean(squeeze(total_new_sales(11:40,1)),2);
plot(mean(squeeze(total_new_sales(11:40,1)),2)./total_sum,'LineWidth',2)
plot((mean(squeeze(total_old_firm_new_sales(11:40,1)),2) + mean(squeeze(total_new_sales(11:40,1)),2))./total_sum,'LineWidth',2)
set(gca,'FontSize',18)
legend({'New exporter','Incumbent exporter, new importer','Incumbent exporter and importer'},'location','SouthWest','FontSize',20)
title("Total sales");
ylim([0.0 1.0])
xlim([1,30])
xticks([5 10 15 20 25 30])
xticklabels({'-10','-5','0','5','10','15'})
hold off
saveas(gcf,'results/exch_shock_plots/decomp_plots/total_sales_pct.png')

%gaps
total_sum_shock = mean(squeeze(total_old_firm_old_sales(11:40,2)),2) + mean(squeeze(total_old_firm_new_sales(11:40,2)),2) + mean(squeeze(total_new_sales(11:40,2)),2);
gap_incum_new_old_sales = mean(squeeze(total_new_sales(11:40,2)),2)./total_sum_shock - mean(squeeze(total_new_sales(11:40,1)),2)./total_sum;
gap_incum_entrant_sales = (mean(squeeze(total_old_firm_new_sales(11:40,2)),2) + mean(squeeze(total_new_sales(11:40,2)),2))./total_sum_shock - (mean(squeeze(total_old_firm_new_sales(11:40,1)),2) + mean(squeeze(total_new_sales(11:40,1)),2))./total_sum;

A = [mean(total_new_sales(11:40,3),2),mean(total_old_firm_new_sales(11:40,3),2),mean(total_old_firm_old_sales(11:40,3),2)];
area(A./repmat(sum(A,2),1,size(A,2)));
hold on
total_sum = mean(squeeze(total_old_firm_old_sales(11:40,1)),2) + mean(squeeze(total_old_firm_new_sales(11:40,1)),2) + mean(squeeze(total_new_sales(11:40,1)),2);
plot(mean(squeeze(total_new_sales(11:40,1)),2)./total_sum,'LineWidth',2)
plot((mean(squeeze(total_old_firm_new_sales(11:40,1)),2) + mean(squeeze(total_new_sales(11:40,1)),2))./total_sum,'LineWidth',2)
set(gca,'FontSize',18)
legend({'New exporter','Incumbent exporter, new importer','Incumbent exporter and importer'},'location','SouthWest','FontSize',20)
title("Total sales");
ylim([0.0 1.0])
xlim([1,30])
xticks([5 10 15 20 25 30])
xticklabels({'-10','-5','0','5','10','15'})
hold off
saveas(gcf,'results/exch_shock_plots/decomp_plots/total_sales_pct_unf.png')

%Number of firms
hold off
area([mean(squeeze(total_new_firms(11:40,2)),2),mean(squeeze(total_firms(11:40,2) - total_new_firms(11:40,2)),2)])
hold on
plot(mean(squeeze(total_new_firms(11:40,1)),2),'LineWidth',2)
plot(mean(squeeze(total_firms(11:40,1)),2),'LineWidth',2)
set(gca,'FontSize',18)
%legend({'New exporter','Incumbent exporter'},'location','SouthWest','FontSize',20)
%title("Total exporters");
xlim([1,30])
xticks([5 10 15 20 25 30])
xticklabels({'-10','-5','0','5','10','15'})
hold off
saveas(gcf,'results/exch_shock_plots/decomp_plots/total_firms.png')

hold off
area([mean(squeeze(total_new_firms(11:40,3)),2),mean(squeeze(total_firms(11:40,3) - total_new_firms(11:40,3)),2)])
hold on
plot(mean(squeeze(total_new_firms(11:40,1)),2),'LineWidth',2)
plot(mean(squeeze(total_firms(11:40,1)),2),'LineWidth',2)
%ylim([0 45])
set(gca,'FontSize',18)
legend({'New exporter','Incumbent exporter'},'location','SouthWest','FontSize',20)
title("Total exporters");
xlim([1,30])
xticks([5 10 15 20 25 30])
xticklabels({'-10','-5','0','5','10','15'})
hold off
saveas(gcf,'results/exch_shock_plots/decomp_plots/total_firms_unf.png')

A = [mean(total_new_firms(11:40,2),2),mean(total_firms(11:40,2) - total_new_firms(11:40,2),2)];
area(A./repmat(sum(A,2),1,size(A,2)));
hold on
total_sum = mean(squeeze(total_firms(11:40,1)),2);
plot(mean(squeeze(total_new_firms(11:40,1)),2)./total_sum,'LineWidth',2)
set(gca,'FontSize',18)
legend({'New exporter','Incumbent exporter'},'location','SouthWest','FontSize',20)
title("Total exporters");
ylim([0.0 1.0])
xlim([1,30])
xticks([5 10 15 20 25 30])
xticklabels({'-10','-5','0','5','10','15'})
hold off
saveas(gcf,'results/exch_shock_plots/decomp_plots/total_firms_pct.png')

%gaps
total_sum_shock = mean(squeeze(total_firms(11:40,2)),2);
gap_incum_entrant_firms = mean(squeeze(total_new_firms(11:40,2)),2)./total_sum_shock - mean(squeeze(total_new_firms(11:40,1)),2)./total_sum;

A = [mean(total_new_firms(11:40,3),2),mean(total_firms(11:40,3) - total_new_firms(11:40,3),2)];
area(A./repmat(sum(A,2),1,size(A,2)));
hold on
total_sum = mean(squeeze(total_firms(11:40,1)),2);
plot(mean(squeeze(total_new_firms(11:40,1)),2)./total_sum,'LineWidth',2)
set(gca,'FontSize',18)
legend({'New exporter','Incumbent exporter'},'location','SouthWest','FontSize',20)
title("Total exporters");
ylim([0.0 1.0])
xlim([1,30])
xticks([5 10 15 20 25 30])
xticklabels({'-10','-5','0','5','10','15'})
hold off
saveas(gcf,'results/exch_shock_plots/decomp_plots/total_firms_pct_unf.png')

elas_var = @(new,base,cov_mat,log_arg) [1 / (log(log_arg)); -1 / (log(log_arg))]' * cov_mat * [1 / (log(log_arg)); -1 / (log(log_arg))];

%elasticities (favorable)
display('Favorable elasticities (Std Err)')
display('SALES');
display([(log(total_sales(26,2)) - log(total_sales(24,1)))/(log(1.2) - log(1)),...
    (log(total_sales(30,2)) - log(total_sales(24,1)))/(log(1.2) - log(1)),...
    (log(total_sales(50,2)) - log(total_sales(24,1)))/(log(1.2) - log(1))]);

display('MATCHES');
display([(log(total_matches(26,2)) - log(total_matches(24,1)))/(log(1.2) - log(1)),...
    (log(total_matches(30,2)) - log(total_matches(24,1)))/(log(1.2) - log(1)),...
    (log(total_matches(50,2)) - log(total_matches(24,1)))/(log(1.2) - log(1))]);

display('FIRMS');
display([(log(total_firms(26,2)) - log(total_firms(24,1)))/(log(1.2) - log(1)),...
    (log(total_firms(30,2)) - log(total_firms(24,1)))/(log(1.2) - log(1)),...
    (log(total_firms(50,2)) - log(total_firms(24,1)))/(log(1.2) - log(1))]);

display('Unfavorable elasticities')
display('SALES')
display([(log(total_sales(26,3)) - log(total_sales(24,1)))/(log(1.2) - log(1)),...
    (log(total_sales(30,3)) - log(total_sales(24,1)))/(log(1.2) - log(1)),...
    (log(total_sales(50,3)) - log(total_sales(24,1)))/(log(1.2) - log(1))]);

display('MATCHES');
display([(log(total_matches(26,3)) - log(total_matches(24,1)))/(log(1.2) - log(1)),...
    (log(total_matches(30,3)) - log(total_matches(24,1)))/(log(1.2) - log(1)),...
    (log(total_matches(50,3)) - log(total_matches(24,1)))/(log(1.2) - log(1))]);

display('FIRMS');
display([(log(total_firms(26,3)) - log(total_firms(24,1)))/(log(1.2) - log(1)),...
    (log(total_firms(30,3)) - log(total_firms(24,1)))/(log(1.2) - log(1)),...
    (log(total_firms(50,3)) - log(total_firms(24,1)))/(log(1.2) - log(1))]);

% display('No shock growth (50-75)?')
% display('SALES');
% display([((mean(squeeze(total_sales(50,1)),2) - mean(squeeze(total_sales(50,:,1)),2)) / mean(squeeze(total_sales(50,:,1)),2))])
% display('MATCHES');
% display([((mean(squeeze(total_matches(50,1)),2) - mean(squeeze(total_matches(50,:,1)),2)) / mean(squeeze(total_matches(50,:,1)),2))])
% display('FIRMS');
% display([((mean(squeeze(total_firms(50,1)),2) - mean(squeeze(total_firms(50,:,1)),2)) / mean(squeeze(total_firms(50,:,1)),2))])
%
% display('No shock growth (25-75)?')
% display('SALES');
% display([((mean(squeeze(total_sales(50,1)),2) - mean(squeeze(total_sales(25,:,1)),2)) / mean(squeeze(total_sales(25,:,1)),2))])
% display('MATCHES');
% display([((mean(squeeze(total_matches(50,1)),2) - mean(squeeze(total_matches(25,:,1)),2)) / mean(squeeze(total_matches(25,:,1)),2))])
% display('FIRMS');
% display([((mean(squeeze(total_firms(50,1)),2) - mean(squeeze(total_firms(25,:,1)),2)) / mean(squeeze(total_firms(25,:,1)),2))])
%
% display('No shock growth (50-99)?')
% display('SALES');
% display([((mean(squeeze(total_sales(99,:,1)),2) - mean(squeeze(total_sales(50,:,1)),2)) / mean(squeeze(total_sales(50,:,1)),2))])
% display('MATCHES');
% display([((mean(squeeze(total_matches(99,:,1)),2) - mean(squeeze(total_matches(50,:,1)),2)) / mean(squeeze(total_matches(50,:,1)),2))])
% display('FIRMS');
% display([((mean(squeeze(total_firms(99,:,1)),2) - mean(squeeze(total_firms(50,:,1)),2)) / mean(squeeze(total_firms(50,:,1)),2))])

% area([mean(total_first_yr_matches(11:40),2),mean(total_non_first_yr_matches(11:40),2)])
% legend({'First year matches','Other matches'},'location','NorthWest','FontSize',14)
% title("First Year Matches");
% saveas(gcf,'results/exch_shock_plots/decomp_plots/total_first_yr_matches.png')
% area([mean(total_first_yr_firms(11:40),2),mean(total_non_first_yr_firms(11:40),2)])
% legend({'First year exporter','Other exporters'},'location','NorthWest','FontSize',14)
% title("First Year Firms");
% saveas(gcf,'results/exch_shock_plots/decomp_plots/total_first_yr_firms.png')
% area([mean(total_first_yr_firm_sales(11:40),2),mean(total_non_first_yr_firm_sales(11:40),2)])
% legend({'First year exporter sales','Other exporter sales'},'location','NorthWest','FontSize',14)
% title("First year firm sales");
% saveas(gcf,'results/exch_shock_plots/decomp_plots/total_first_yr_firm_sales.png')
% area([mean(total_first_yr_match_sales(11:40),2),mean(total_non_first_yr_match_sales(11:40),2)])
% legend({'First year match sales','Other match sales'},'location','NorthWest','FontSize',14)
% title("First year match sales");
% saveas(gcf,'results/exch_shock_plots/decomp_plots/total_first_yr_match_sales.png')
%
% %percentage
% A = [mean(total_first_yr_matches(11:40),2),mean(total_non_first_yr_matches(11:40),2)];
% area(A./repmat(sum(A,2),1,size(A,2)));
% legend({'First year matches','Other matches'},'location','NorthWest','FontSize',14)
% title("First year matches");
% saveas(gcf,'results/exch_shock_plots/decomp_plots/total_first_yr_matches_pct.png')
% A = [mean(total_first_yr_firms(11:40),2),mean(total_non_first_yr_firms(11:40),2)];
% area(A./repmat(sum(A,2),1,size(A,2)));
% legend({'First year exporter','Other exporters'},'location','NorthWest','FontSize',14)
% title("First year Exporters");
% saveas(gcf,'results/exch_shock_plots/decomp_plots/total_first_yr_firms_pct.png')
% A = [mean(total_first_yr_firm_sales(11:40),2),mean(total_non_first_yr_firm_sales(11:40),2)];
% area(A./repmat(sum(A,2),1,size(A,2)));
% legend({'First year exporter sales','Other exporter sales'},'location','NorthWest','FontSize',14)
% title("First year firm sales");
% saveas(gcf,'results/exch_shock_plots/decomp_plots/total_first_yr_firm_sales_pct.png')
% A = [mean(total_first_yr_match_sales(11:40),2),mean(total_non_first_yr_match_sales(11:40),2)];
% area(A./repmat(sum(A,2),1,size(A,2)));
% legend({'First year match sales','Other match sales'},'location','NorthWest','FontSize',14)
% title("First year match sales");
% saveas(gcf,'results/exch_shock_plots/decomp_plots/total_first_yr_match_sales_pct.png')



%     %plots
%     plot(exch_rate_series(25:74,1)/mean(exch_rate_series(49:49,1)))
%     hold on
%     plot(exch_rate_series_shk(25:74,1)/mean(exch_rate_series_shk(49:49,1)))
%     hold on
%     xline(25)
%     ylim([0.85 1.4])
%     xlim([1,50])
%     xticks([0,5,10,15,20,25,30,35,40,45,50])
%     xticklabels({'-25','-20','-15','-10','-5','0','5','10','15','20','25',})
%     legend({'Baseline','Devaluation'},'Location','southeast')
%     hold off
%     saveas(gcf,'results/exch_shock_plots/exch_rate_series_compare.png')
%
%     %title('Total Matches and Exch Rate')
%     plot(exch_rate_series(25:74,1)/mean(exch_rate_series(49:49,1)))
%     hold on
%     plot(total_matches(26:75,1)/mean(total_matches(49:49,1)))
%     xline(25)
%     ylim([0.9 1.15])
%     xlim([1,50])
%     xticks([0,5,10,15,20,25,30,35,40,45,50])
%     xticklabels({'-25','-20','-15','-10','-5','0','5','10','15','20','25',})
%     legend({'Pesos to Dollar Exch Rate','Total matches'},'Location','southeast')
%     hold off
%     saveas(gcf,'results/exch_shock_plots/total_matches_exch_rate_series.png')
%
%
%     %plots
%     plot(total_sales(26:75,1)/mean(total_sales(49:49,1)))
%     title('Total Sales and Exch Rate')
%     hold on
%     plot(exch_rate_series(25:74,1)/mean(exch_rate_series(49:49,1)))
%     xline(25)
%     ylim([0.75 1.2])
%     xlim([1,50])
%     xticks([0,5,10,15,20,25,30,35,40,45,50])
%     xticklabels({'-25','-20','-15','-10','-5','0','5','10','15','20','25',})
%     legend({'Baseline','Pesos to Dollar Exch Rate'},'Location','southeast')
%     hold off
%     saveas(gcf,'results/exch_shock_plots/total_sales_exch_rate_series.png')
%
% plot(mean(total_matches(36:65,:),2)/squeeze(mean(total_matches(49:49,:))))
% title('Total Matches')
% %hold on
% %plot(total_matches(36:65,2)/mean(total_matches(49:49,2)))
% %hold on
% %plot(exch_rate_series_shk(35:64,1)/mean(exch_rate_series_shk(49:49,1)))
% %hold on
% xline(15)
% ylim([0.95 1.25])
% xlim([1,30])
% xticks([0,5,10,15,20,25,30])
% xticklabels({'-15','-10','-5','0','5','10','15'})
% %legend({'Baseline','Devaluation'},'Location','southeast')
% hold off
% saveas(gcf,'results/exch_shock_plots/total_matches.png')
%
% plot(mean(total_firms(36:65,:),2)/mean(total_firms(49:49,:)))
% title('Total Firms')
% hold on
% plot(mean(total_new_firms(36:65,:),2)/mean(total_new_firms(49:49,:)))
% title('Total Firms')
% hold on
% %plot(total_firms(36:65,2)/mean(total_firms(49:49,2)))
% %hold on
% %plot(exch_rate_series_shk(35:64,1)/mean(exch_rate_series_shk(49:49,1)))
% %hold on
% xline(15)
% ylim([0.85 1.15])
% xlim([1,30])
% xticks([0,5,10,15,20,25,30])
% xticklabels({'-15','-10','-5','0','5','10','15'})
% %legend({'Baseline','Devaluation'},'Location','southeast')
% hold off
% saveas(gcf,'results/exch_shock_plots/total_firms.png')
%
% plot(mean(total_sales(36:65,:),2)/mean(total_sales(49:49,:)))
% title('Total Sales')
% %hold on
% %plot(total_sales(36:65,2)/mean(total_sales(49:49,2)))
% %hold on
% %plot(exch_rate_series(35:64,1)/mean(exch_rate_series(49:49,1)))
% %hold on
% ylim([0.7 1.65])
% xlim([1,30])
% xline(15)
% xticks([0,5,10,15,20,25,30])
% xticklabels({'-15','-10','-5','0','5','10','15'})
% %legend({'Baseline','Devaluation'},'Location','southeast')
% hold off
% saveas(gcf,'results/exch_shock_plots/total_sales.png')
%
% plot(sales_per_match(36:65,1)/mean(sales_per_match(49:49,1)))
% title('Sales Per Match')
% hold on
% plot(sales_per_match(36:65,2)/mean(sales_per_match(49:49,2)))
% hold on
% %plot(exch_rate_series(35:64,1)/mean(exch_rate_series(49:49,1)))
% ylim([0.8 1.35])
% xlim([1,30])
% xline(15)
% xticks([0,5,10,15,20,25,30])
% xticklabels({'-15','-10','-5','0','5','10','15'})
% legend({'Baseline','Devaluation'},'Location','southeast')
% hold off
% saveas(gcf,'results/exch_shock_plots/sales_per_match.png')
%
% plot(sales_per_firm(36:65,1)/mean(sales_per_firm(49:49,1)))
% title('Sales Per Firm')
% hold on
% plot(sales_per_firm(36:65,2)/mean(sales_per_firm(49:49,2)))
% ylim([0.75 1.5])
% xlim([1,30])
% xline(15)
% xticks([0,5,10,15,20,25,30])
% xticklabels({'-15','-10','-5','0','5','10','15'})
% legend({'Baseline','Devaluation'},'Location','southeast')
% hold off
% saveas(gcf,'results/exch_shock_plots/sales_per_firm.png')
%
% plot(matches_per_firm(36:65,1)/mean(matches_per_firm(49:49,1)))
% title('Matches Per Firm')
% hold on
% plot(matches_per_firm(36:65,2)/mean(matches_per_firm(49:49,2)))
% xlim([1,30])
% xline(15)
% ylim([0.9 1.15])
% xticks([0,5,10,15,20,25,30])
% xticklabels({'-15','-10','-5','0','5','10','15'})
% legend({'Baseline','Devaluation'},'Location','southeast')
% hold off
% saveas(gcf,'results/exch_shock_plots/matches_per_firm.png')

