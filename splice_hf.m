function [x,y,export_rate,nobs,nobs_f,nobs_h,moms_xy,moms_xx,ysum,nobs_hf] =...
    splice_hf(firm_h_yr_sales,firm_f_yr_sales,typemat,theta_h_firm,periods)

% This function takes panel of realizations on domestic and foreign sales
% for a particular type (theta_f and phi), and splices them together to 
% create synthetic firms

% firm_f_yr_sales: [t,type,firm ID, total exports,total # foreign shipments,firm age in export mkt.]
% firm_h_yr_sales: [t,type,firm ID, total dom. sales, total # dom. shipments,firm age in dom. mkt.]
% typemat: [type, macro state index, theta index, productivity index]

% some firms with active export relationships have zero exports. Weed them out:
exporter = firm_f_yr_sales(:,4)>0;
firm_f_yr_sales = firm_f_yr_sales(exporter,:); % stacks all firm-yr combs. with exports, given micro type
% same deal for domestic market shipments
dom_shipper = firm_h_yr_sales(:,4)>0;  % stacks all firm-yr combs. with dom .sales, given micro type
firm_h_yr_sales = firm_h_yr_sales(dom_shipper,:);

nobs_f = sum(exporter>0);
nobs_h = sum(dom_shipper>0);

type_h = firm_h_yr_sales(:,2);
type_f = firm_f_yr_sales(:,2);
prod_h  = typemat(type_h,4);
prod_f  = typemat(type_f,4);

try % confirm domestic and foreign records are for the same firm productivity type
    minobs = min(size(type_h,1),size(type_f,1));
    assert(sum((prod_h(1:minobs)- prod_f(1:minobs)).^2)==0);
catch
     fprintf('Warning: productivity different in different markets');
end

%% Extract observations on home and foreign sales with same firm_ID and date

% create unique identifiers for each period/firm_ID pair
obs_id_h = firm_h_yr_sales(:,3) + (1/(periods+1))*firm_h_yr_sales(:,1); % micro type and period
obs_id_f = firm_f_yr_sales(:,3) + (1/(periods+1))*firm_f_yr_sales(:,1); % micro type and period
% find home market obs. with identifiers that appear in export mkt. obs.
aa = ismember(obs_id_h,obs_id_f);
temp1 = obs_id_h(aa);
% find foreign market obs. with identifiers that appear in home mkt.
bb = ismember(obs_id_f,temp1);
try
  assert(sum(aa)==sum(bb))
catch
  fprintf('\r Warning: home-foreign merge discrepancy in spice_hf');
end
nobs_hf  = sum(aa);
nobs = nobs_h + nobs_f - nobs_hf; % number of firms with sales in at least one market

% dom_only = ones(nobs_h,1) - aa > 0; % firms with dom. sales only
 ex_only  = ones(nobs_f,1) - bb > 0; % firms with exports only

theta_f = typemat(type_f,3); % Each element of this vector is common to all firms
theta_h = theta_h_firm;      % This is a vector of random draws--one per firm
if sum(bb) == 0
  temp3 = zeros(0,15);
else
  temp3 = [theta_f(bb),theta_h(aa),prod_h(aa),firm_h_yr_sales(aa,:),firm_f_yr_sales(bb,:)];
end
% firm_h_yr_sales: [t,type,firm ID,total dom. sales, total # dom. shipments,firm age in dom. mkt.]
% firm_f_yr_sales: [t,type,firm ID,total exports,total # foreign shipments,firm age in export mkt.]

try % comfirm home and foreign records are for same period and firm ID
assert(sum((temp3(:,4) - temp3(:,10)).^2)==0); % periods match?
assert(sum((temp3(:,6) - temp3(:,12)).^2)==0); % firm ID's match?
sales_splice = [temp3(:,1:4),temp3(:,6),temp3(:,7),temp3(:,13),temp3(:,9),temp3(:,15)];

% sales splice contains records for firms with matched foreign and domestic sales: 
%  [1) theta_f, 2) theta_h, 3) prod_h, 4) t, 5) firm ID, 
%   6) total dom. sales, 7) total exports, 
%   8) age in dom. mkt, 9) age in export mkt.]
catch
   warning('time periods or firm IDs do not match for home versus foreign sales');
   sales_splice = zeros(0,9);
end

% exiting firms are replaced by new ones with the same ID, so need to
% exclude cases where the "firm" age in the export market exceeds the firm
% age in the domestic market to avoid inherited foreign client bases

same_firm   = sales_splice(:,8) >= sales_splice(:,9);
sales_hf    = sales_splice(same_firm,1:7);

% export_rate = [sales_hf(:,7)./(sales_hf(:,6)+sales_hf(:,7));ones(sum(ex_only,1),1)];
% the vector of ones in export_rate accounts for firms that exclusively

export_rate = sales_hf(:,7)./(sales_hf(:,6)+sales_hf(:,7));
% this alternative version excludes firms active solely in the export market

% serve foreign markets

%% moments for regression of log foreign sales on log domestic sales

   nobs_hf = size(sales_hf,1);
   if nobs_hf > 0
    y = log(sales_hf(:,7));
    x = [ones(nobs_hf,1),log(sales_hf(:,6))];
    moms_xx = x'*x;
    moms_xy = x'*y;
    ysum = sum(y);    
   else
        x = zeros(0,2);
        y = zeros(0,1);
        moms_xx = zeros(2,2);
        moms_xy = zeros(2,1);
        ysum    = 0;
        nobs_hf = 0;
    end

end

