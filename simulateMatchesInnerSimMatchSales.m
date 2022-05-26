function [mat_tran,ship_cur,age_vec] = simulateMatchesInnerSimMatchSales(mkt,mm,iter_in,age)

if mkt==1   % foreign market
  scale     = mm.scale_f;
  macro_shk = mm.X_f(iter_in.macro_state_f(iter_in.t));
elseif mkt==2 % home market
  scale     = mm.scale_h;
  macro_shk = mm.X_h(iter_in.macro_state_h(iter_in.t));
end

% This function constucts match transition counts, shipment counts and shipment 
% sales for a specific type of exporter and a time period, t

% mat_tran contains col 1: initial state, col 2: exporter id, 
%       col 3: dest. state (Z index, or 0 for exit); col 4: match rev.
% ship_cur contains number of shipments for each match
% age_vec  contains exporter age for each match


%% create 0/1 matrix "new_st" to track transitions of each firm's individual matches.
[tot_cli, new_st] = simulateMatchesInnerSimMatchTrans(iter_in, mm);
%% construct match sales 
[mat_tran, ship_cur, age_vec] = simulateMatchesInnerSimMatchRevAssign(tot_cli, mm, new_st, scale, iter_in, macro_shk, age);
     

end
