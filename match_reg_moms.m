function  [x,y,moms_xx,moms_xy,ysum,n_obs] = match_reg_moms(mat_cont,ncols)

     % matches: [firm ID, match-specific sales, shipments, boy Z, eoy Z, match age, firm age]  
     % mat_cont: [matches_lagged, matches] spliced, continuing matches only
     
    ff_obs = find(mat_cont(:,2).*mat_cont(:,ncols+2)>0); % obs. with positive current and lagged sales
    dat = mat_cont(ff_obs,:);
    n_obs = size(ff_obs,1);
    y = log(dat(:,ncols+2));         % log year current sales
    ln_age = log(dat(:,ncols)+1/4);  % log firm age last year
    x = [ones(n_obs,1),log(dat(:,2)),dat(:,4)==0,ln_age];  % 1st yr. match dummy, log(lagged sales), firm age
    moms_xx = x'*x;
    moms_xy = x'*y;
    ysum = sum(y);
    
%     ff2_obs = find(matches(:,3) >0); % obs. with positive current year shipments
%     nship_obs = size(ff2_obs,1);
%     ln_ships = sum(log(matches(ff2_obs,3)));
%     
% %     if nship_obs > 40
% %         'pause here'
% %     end
%     
% %% construct match count vector
%     
%     if sum(matches(:,3))>0       
%         % add up matches with positive shipments, firm by firm:
%         num_match = sum(dummyvar(matches(:,1)).*(matches(:,3)>0)); 
%         % count firms with 0 matches, 1 match, etc.
%         temp = sum(dummyvar([num_match'+1;max_match+1]))';
%         % topcode firms with more than max_match matches
%         temp = min(temp,max_match+1); 
%         match_count = temp(2:max_match+1);
%         match_count(max_match) = match_count(max_match) - 1;
%     else
%         match_count = zeros(max_match,1);
%     end
   
%     if size(x,1) > 40
%      'pause here'
%     end
end


