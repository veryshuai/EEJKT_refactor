function  [x,y,moms_xx,moms_xy,ysum,n_obs] = firm_reg_h_moms(iterH_in,mm)

% sales and sales_lag:[firm ID, total domestic sales,total domestic shipments,age]

    mat_h_cont_2yr = iterH_in.mat_h_cont_2yr;
    firm_ID        = sort(unique(floor(mat_h_cont_2yr(:,1))));

    sales     = zeros(length(firm_ID),1);
    lag_sales = zeros(length(firm_ID),1);
    for j=1:length(firm_ID)       
       sales(j) =  sum(mat_h_cont_2yr(:,9)...
            .*(mat_h_cont_2yr(:,8)==firm_ID(j)));     
        
        % all sales of previous yr firms that continue to current year
        lag_sales(j) = sum(mat_h_cont_2yr(:,2)...
            .*(mat_h_cont_2yr(:,1)==firm_ID(j)));
    end

    n_obs = length(sales);
    if n_obs>0
        
    y = log(sales+1);
   x = [ones(n_obs,1),log(lag_sales+1)];
    moms_xx = x'*x;
    moms_xy = x'*y;
    ysum = sum(y);
    
    else
        kk = 2; % columns of x matrix
        x = zeros(0,2);
        y = zeros(0,1);
        moms_xx = zeros(2,2);
        moms_xy = zeros(2,1);
        ysum    = 0;
        n_obs = 0;
    end
end


