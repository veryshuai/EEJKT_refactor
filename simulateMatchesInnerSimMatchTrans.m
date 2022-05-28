function [tot_cli, new_st] = simulateMatchesInnerSimMatchTrans(iter_in, mm)

% This function deals with the fact that multiple matches may share the
% same initial and ending Z-state for a given exporter at a given t. The
% matches are separated into distinct records, keeping their common exporter ID.

% iter_in.trans_count(:,:,:) counts destination z states for all clients of each firm. 
% Dimensions: (1) initial state (2) destination state (3) exporting firm.
% Row 1 contains new entrants in current period; col 1 contains exits.
% State 0 (column 1) identifies exit for incumbent matches and the pre-entry state
% for entrants (row 1). States 1 through N_Z correspond to Z indices.


n_cli2      = zeros(size(iter_in.trans_count,1),size(iter_in.trans_count,3)); 
n_cli2(:,:) = sum(iter_in.trans_count(:,:,:),2); % # clients by initial state (rows 1:N_Z+1) and exporting firm (cols 2:N_firms)
n_cli1      = sum(n_cli2,1)';      % initial client counts, incl. entrants, by exporting firm
  
  [init_st2,expr2] = find(n_cli2>0); % populated (initial states/exporter) addresses, incl. entering exporters
  sz_init_st2 = size(init_st2,1);    % # populated initial state/exporter cells
  tot_cont_cli = sum(n_cli1)-sum(sum(sum(iter_in.trans_count(:,1,:)))); % total # end-of-period clients. trans(:,1,:) contains exits
  % assert(sum(sum(sum(trans(:,2:N_Z+1,:))))==tot_cont_cli) 
  tot_cli = sum(sum(n_cli2)); % total # initial (t-1) clients     
  

  new_st = zeros(tot_cli,size(mm.Z,1)+3); % first 3 cols. give match state, exporter, & exit
  ss = 1;

  for jj=1:1:sz_init_st2 % iterate over populated initial states in tran()

    % # clients of exporter expr(jj) beginning from state init_st2(jj): 
    n_match = sum(iter_in.trans_count(init_st2(jj),:,expr2(jj)),2); 
    % sum over dest. states the client(s) of exporter expr(jj) beginning 
    % from state init_st2(jj). 
    
    tloc = find(iter_in.trans_count(init_st2(jj),:,expr2(jj))>0);  
    nloc = size(tloc,2);
    % identify and count destination states for client(s) of exporter expr(jj) beginning 
    % from state init_st2(jj). Will be < n_match if some share same destination:  
    
    temp = zeros(n_match,size(mm.Z,1)+1);
    % temp will hold dummies indicating destination states for each of 
    % exporter expr2(jj)'s clients that begin from initial state init_st2(jj)

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

      new_st(ss:ss+n_match-1,:) =...
          [ones(n_match,1).*(init_st2(jj)-1),ones(n_match,1).*expr2(jj),temp];
      % for initial state init_st2(jj), a block of n_match matches in new_st 
      % contains col 1: initial state, col 2: exporter id, cols 3:N_Z+3: destination state dummies
      % Destination state 1 is exit, states 2:N_Z+1 are Z realizations 1:N_Z

      ss = ss + n_match;
  end
  assert(sum(sum(new_st(:,4:size(mm.Z,1)+3)))==tot_cont_cli);
  
end