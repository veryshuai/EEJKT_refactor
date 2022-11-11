%%               simulateMatchesInnerSimMatchSales

% This function constucts match transition counts, shipment counts and shipment 
% sales for a specific market (mkt) type of firm (pt_ndx) and time period (t)

%%
function [mat_tran,ship_cur,age_vec] = simulateMatchesInnerSimMatchSales(mkt,mm,iterX_in,age)

% ship_cur contains number of shipments for each match
% age_vec  contains exporter age for each match
% mat_tran contains information on match transitions, period t-1 to t:
%    col 1: initial state, 
%    col 2: exporter id, 
%    col 3: dest. state (Z index, or 0 for exit); 
%    col 4: match revenue 

% create 0/1 matrix "new_st" to track transitions of each firm's individual matches.
[tot_cli, new_st] = simulateMatchesInnerSimMatchTrans(iterX_in.trans_count, mm);

% construct match sales for foreign market (mkt=1) or home market (mkt=2)
  
[mat_tran, ship_cur, age_vec] =...
    simulateMatchesInnerSimMatchRevAssign(tot_cli, mm, new_st, mkt, iterX_in, age);


end
