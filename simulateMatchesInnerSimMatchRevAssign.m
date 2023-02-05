function [mat_tran, ship_cur, age_vec] = simulateMatchesInnerSimMatchRevAssign(tot_cli, mm, new_st, mkt, iterX_in, age)

if mkt==1   % foreign market
  scale     = mm.scale_f;
  macro_shk = mm.X_f(iterX_in.macro_state_f(iterX_in.t));
  poisCDF_shipments = mm.poisCDF_shipmentsF;
  max_ships  = mm.max_shipsF;
elseif mkt==2 % home market
  scale     = mm.scale_h;
  macro_shk = mm.X_h(iterX_in.macro_state_h(iterX_in.t));
  poisCDF_shipments = mm.poisCDF_shipmentsH;
  max_ships  = mm.max_shipsH;
end

   mask = ones(tot_cli,1).*(0:1:size(mm.Z,1));
   t_state = sum(new_st(:,3:size(mm.Z,1)+3).*mask,2); % convert Z_state dummies to indices
% mat_tran col 1: initial state, col 2: exporter id, cols 3: dest. state (Z index, or 0 for exit)
   mat_tran = [new_st(:,1:2),t_state];
 
 % pick off new match shocks (Z's) for all matches, including new ones. 
 % 0's are exits
   new_expZ = sum((ones(tot_cli,1).*exp(mm.Z)').*new_st(:,4:size(mm.Z,1)+3),2);

 % draw random shipment counts for new and continuing matches.
   rr = rand(tot_cli,1);
   select  = ones(tot_cli,1).*poisCDF_shipments > rr; 
   ship_cur = (max_ships.*ones(tot_cli,1) - sum(select,2)).*(new_expZ>0); % shipments, t   
   
   %% modeling choice here
   
   sampl_ship = mat_tran(:,1) == 0; % all new matches generate a sample shipment
%  ship_cur   = ship_cur + 0.5*sampl_ship; % impose sample shipment is 1/2 normal shipment 
   ship_cur   = ship_cur + 0.0*sampl_ship; % impose sample shipment is negligible
   
  % NOTE: the ratio of sample shipment to normal shipment could be estimated, 
  % but for now we fix it at 0.
%%
%  Note: revenue function is based expression in makepie (before Z effect):
%  payoff = 1 / de * exp(sf) * exp((de-1)*st(1,:)+st(2,:)); 
%  Here sf is estimated scalar, st(1,:) is productivity, st(2,:) is macro state  

   mat_rev  = exp(scale + (mm.eta-1)*mm.Phi(mm.pt_type(iterX_in.pt_ndx,1)) + macro_shk).*new_expZ.*ship_cur; 
   
%  [scale,phi,macro_shk,mean(mat_rev),mean(new_expZ)]
   
   mat_tran = [mat_tran,mat_rev];
   age_vec = age(new_st(:,2)); % firm age in periods

   
end