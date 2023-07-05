function [exit_moms_xx,exit_moms_xy,sum_succ_rate,sum_exits,exit_obs] = mkt_exit_moms(mkt_exit)
% mkt_exit: (1) exit dummy (2) cumulative meetings (3) cumulative successes

ff = mkt_exit(:,2)>0; % obs. with positive number of cumulative meetings
% if sum(ff,1) ~= size(mkt_exit,1)
%     mkt_exit
%     ff
% end

mkt_exit = mkt_exit(ff,:);

exit_obs = size(mkt_exit,1);
x1 = log(1+mkt_exit(:,3));  % log cumulative successes
x2 = x1.^2;
x3 = log(1+(mkt_exit(:,3)./mkt_exit(:,2))); % log match success rate
x4 = x3.^2;
x5 = x1.*x3;
x0 = ones(exit_obs,1);
x  = [x0,x1,x2,x3,x4,x5];
y  = mkt_exit(:,1);
exit_moms_xx = x'*x;
exit_moms_xy = x'*y;
sum_succ_rate = sum(mkt_exit(:,3)./mkt_exit(:,2));
sum_exits = sum(mkt_exit(:,1));

end