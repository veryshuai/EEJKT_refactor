function [mat_tran,ship_cur,age_vec] = simulateMatchesInnerSimMatchSales(mkt,mm,iter_in,age,pt_ndx,macro_state)
iter_in.trans_count,age,iter_in.pt_ndx,iter_in.macro_state_f(iter_in.t)

if     mkt==1 % foreign market
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

% trans(:,:,:) counts destination z states for all clients of each firm. 
% Dimensions: (1) initial state (2) destination state (3) exporting firm.
% Row 1 contains new entrants in current period; col 1 contains exits.
% State 0 (column 1) identifies exit for incumbent matches and the pre-entry state
% for entrants (row 1). States 1 through N_Z correspond to Z indices.


n_cli2      = zeros(size(iter_in.trans_count,1),size(iter_in.trans_count,3)); 
n_cli2(:,:) = sum(iter_in.trans_count(:,:,:),2); % # clients by initial state (rows 1:N_Z+1) and exporting firm (cols 2:N_firms)
n_cli1      = sum(n_cli2,1)';      % initial client counts, incl. entrants, by exporting firm
  
  [init_st2,expr2] = find(n_cli2>0); % populated (initial states,exporter) addresses, incl. entering exporters
  sz_init_st2 = size(init_st2,1);    % # populated initial state/exporter cells
  tot_cont_cli = sum(n_cli1)-sum(sum(sum(iter_in.trans_count(:,1,:)))); % total # end-of-period clients. trans(:,1,:) contains exits
  %assert(sum(sum(sum(trans(:,2:N_Z+1,:))))==tot_cont_cli) 
  tot_cli = sum(sum(n_cli2)); % total # initial (t-1) clients     
  
  %% create 0/1 matrix "new_st" to track transitions of each firm's individual matches.
  new_st = zeros(tot_cli,size(mm.Z,1)+3); % first 3 cols. give match state, exporter, & exit
  ss = 1;

  for jj=1:1:sz_init_st2 % iterate over populated initial states in tran()
    % # clients of exporter expr(jj) beginning from state init_st2(jj): 
    n_match = sum(iter_in.trans_count(init_st2(jj),:,expr2(jj)),2); % sum over dest. states
    % new state(s) for client(s) of exporter expr(jj) beginning from state init_st2(jj):   
    tloc = find(iter_in.trans_count(init_st2(jj),:,expr2(jj))>0);   
    % count destination states for client(s) of exporter expr(jj) beginning 
    % from state init_st2(jj). Will be < n_match if some share same destination:  
    nloc = size(tloc,2);
    % temp will hold dummies indicating destination states for each of 
    % exporter expr2(jj)'s clients that begin from initial state init_st2(jj)
    temp = zeros(n_match,size(mm.Z,1)+1);
    % # clients transiting from init_st2(jj) to new state(s) tloc, exporter expr(jj):
    bloc = iter_in.trans_count(init_st2(jj),tloc,expr2(jj));
    kk = 1;
    % clients that share initial/destination state pairs treated as separate rows of temp:
    temp(1:bloc(1,1),tloc(1,kk)) = ones(bloc(1,1),1); 
    % populate rows of temp for additional destination states, if any
      if nloc>1
      nn=bloc(1,1);
        for kk=2:nloc
         temp(nn+1:nn+bloc(1,kk),tloc(1,kk)) = ones(bloc(1,kk),1);
         nn = nn + bloc(1,kk);
        end
      end
      % for initial state init_st2(jj), a block of n_match matches in new_st 
      % contains col 1: initial state, col 2: exporter id, cols 3:N_Z+3: destination state dummies
      % Destination state 1 is exit, states 2:N_Z+1 are Z realizations 1:N_Z
      new_st(ss:ss+n_match-1,:) =...
          [ones(n_match,1).*(init_st2(jj)-1),ones(n_match,1).*expr2(jj),temp];
      ss = ss + n_match;
  end
  assert(sum(sum(new_st(:,4:size(mm.Z,1)+3)))==tot_cont_cli);
%% construct match sales 

   mask = ones(tot_cli,1).*(0:1:size(mm.Z,1));
   t_state = sum(new_st(:,3:size(mm.Z,1)+3).*mask,2); % convert Z_state dummies to indices
% mat_tran col 1: initial state, col 2: exporter id, cols 3: dest.
% state (Z index, or 0 for exit)
   mat_tran = [new_st(:,1:2),t_state];
 
 % pick off new match shocks (Z's) for all matches, including new ones. 
 % 0's are exits
   new_expZ = sum((ones(tot_cli,1).*exp(mm.Z)').*new_st(:,4:size(mm.Z,1)+3),2);

 % draw random shipment counts for new and continuing matches.
   rr = rand(tot_cli,1);
   select  = ones(tot_cli,1).*mm.poisCDF_shipments > rr; 
   ship_cur = (mm.max_ships.*ones(tot_cli,1) - sum(select,2)).*(new_expZ>0); % shipments, t   
   
   %% modeling choice here
   
   sampl_ship = mat_tran(:,1) == 0; % all new matches generate a sample shipment
%  ship_cur   = ship_cur + 0.5*sampl_ship; % impose sample shipment is 1/2 normal shipment 
   ship_cur   = ship_cur + 0.0*sampl_ship; % impose sample shipment is negligible
   
  % NOTE: the ratio of sample shipment to normal shipment could be estimated, 
  % but for now I am fixing it at 0.
%%
%  Note: revenue function is based expression in makepie (before Z effect):
%  payoff = 1 / de * exp(sf) * exp((de-1)*st(1,:)+st(2,:)); 
%  Here sf is estimated scalar, st(1,:) is productivity, st(2,:) is macro state  

   mat_rev  = exp(scale + (mm.eta-1)*mm.Phi(mm.pt_type(pt_ndx,1)) + macro_shk).*new_expZ.*ship_cur; 
   
%  [scale,phi,macro_shk,mean(mat_rev),mean(new_expZ)]
   
   mat_tran = [mat_tran,mat_rev];
   age_vec = age(new_st(:,2)); % firm age in periods
     

end
