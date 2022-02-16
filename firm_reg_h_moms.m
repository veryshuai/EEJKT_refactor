function  [x,y,moms_xx,moms_xy,ysum,n_obs] = firm_reg_h_moms(sales,sales_lag,N_firms,t,type)

% sales and sales_lag:[firm ID, total domestic sales,total domestic shipments,age]

    nc = size(sales,2);
  
%   Load populated firm records into associated rows of temp1 and temp2
    temp1 = zeros(N_firms,nc);
    temp2 = zeros(N_firms,nc);
    temp1(sales(:,1),:)= sales; 
    temp2(sales_lag(:,1),:)= sales_lag; 
     
%  Get rid of rows with zeros sales in one or both years and rows where 
%  age doesn't increment. The latter reflects firm exit.
    ff_cont1 = temp1(:,2).*temp2(:,2) > 0;   % positive sales in both years
    ff_cont2 = temp1(:,4) - temp2(:,4) > 0; % age increased from lagged to current year
    ff_cont  = ff_cont1.*ff_cont2 > 0;

    if sum(ff_cont)>0
    
    try
      assert(sum(abs(temp1(ff_cont,1)-temp2(ff_cont,1)))==0)
    catch
      warning('continuing firm mismatch')
      temp1(ff_cont,:)
      temp2(ff_cont,:)
    end
    
    dat = [temp1(ff_cont,:),temp2(ff_cont,:)];

    n_obs = sum(ff_cont);
    y = log(dat(:,2));
    ln_age = log(dat(:,4)+1/4);
 %  x = [ones(n_obs,1),dat(:,8)==0,log(dat(:,6)),ln_age];
   x = [ones(n_obs,1),log(dat(:,6))];
    moms_xx = x'*x;
    moms_xy = x'*y;
    ysum = sum(y);
    
%     if size(x,1) > 40
%      'pause here'
%     end
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


