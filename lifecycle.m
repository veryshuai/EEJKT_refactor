function [agg_mat_lifecycle,agg_orphan_matches] = lifecycle(NumTF,TF_matdat,mm,max_age)

% This function splices together match lifecycles

  orphan_matches_TF = cell(NumTF,1);
  mat_lifecycle_TF = cell(NumTF,1);

 parfor TF_id=1:NumTF  
    % fprintf('\r firm ID = %0.0f\n',TF_id);  
      all_age = TF_matdat{TF_id};
    % all_age: (1) year, (2) type, (3) firm_ID, (4) sales, (5) shipments, 
    %          (6) boy Z,(7) eoy Z,(8) match age,(9) firm age 
      new_match  = ...
         (all_age(:,6)==0) +... % entered in w/in year period 2 or later
         (all_age(:,6)>0).*(all_age(:,8)==12)+...% entered in period 1 and survived to eoy
         (all_age(:,7)==0).*(all_age(:,8)==1) >0; % entered and lasted 1 period
        % MISSING: entered in period 1, lasted more than 1 period,
        % died before e.o.y. at an age between 2 and 11 months.
        % These are indistinguishable from dying incumbent matches at this
        % point in the code. Those that don't match to previous period
        %  will be treated as new matches later in the code.
        
      old_match  = 1 - new_match;
      residMatch = all_age(logical(old_match),:);
      residMatch = sortrows(residMatch,7,'descend');

      N_TF = sum(new_match); % number of new matches for this firm
      
      mat_lifecycle = zeros(N_TF,5*ceil(max_age./mm.pd_per_yr));
      mat_lifecycle(:,1:5) = [all_age(new_match,1),all_age(new_match,8),...
               all_age(new_match,6:7),all_age(new_match,4)]; 
    % load one-year-olds into first 5 cols of mat_lifecycle       
    % mat_lifecycle(:,1:5): [year, match age (<=1 yr), boy Z, eoy Z, sales]  
       mat_lifecycleC = mat_lifecycle(mat_lifecycle(:,4)>0,:); % continuing new matches
       mat_lifecycleS = mat_lifecycle(mat_lifecycle(:,4)==0,:); % new matches that fail
       N_TFC = size(mat_lifecycleC,1);
   
      for aa = 2:ceil(max_age./mm.pd_per_yr)
          
          stayInLoop = 1;      
          while stayInLoop > 0 && size(residMatch,1) > 0
        
           if N_TFC==0 
           fprintf('\r type-firm ID = %0.0f, number of 1st year matches = %0.0f\n',[TF_id,N_TF]);
           stayInLoop = 0;               
           end
      
              ii = 1; 
              while ii <= N_TFC %looping over matches of age aa for firm-type TF_id
              fprintf('\r type-firm ID = %0.0f, match age = %0.0f, match number = %0.0f',[TF_id,aa,ii]);
              lagFmAge = mat_lifecycleC(ii,5*(aa-2)+1);
              lagMaAge = mat_lifecycleC(ii,5*(aa-2)+2);
              lag_eoyZ = mat_lifecycleC(ii,5*(aa-2)+4);

            % Find compatible matches to splice with last year's match ii

           maYrOld = (residMatch(:,8)  <= lagMaAge + mm.pd_per_yr)...
                    .*(residMatch(:,8) >= lagMaAge+1);
                   
              flg = (lag_eoyZ>0).*(residMatch(:,6)==lag_eoyZ).*...
                    (maYrOld==1).*(residMatch(:,1)==lagFmAge + mm.pd_per_yr);
                  
             % If compatible matches are identified, put first one in the 
             % relevant bloc of mat_lifecycle and remove it from aa_cohort
    
                 if sum(flg)>0
                   mat_cont = find(flg==1,1); % first compatible match in aa_cohort
                   lb = 5*(aa-1)+1;
                   ub = 5*(aa-1)+5;
                   mat_lifecycleC(ii,lb:ub) = [residMatch(mat_cont,1),residMatch(mat_cont,8),...
                       residMatch(mat_cont,6:7),residMatch(mat_cont,4)]; 
              
                % remove matched record from residMatch 
                   keepers = ones(size(residMatch,1),1);
                   keepers(mat_cont,1) = 0;
                   stayInLoop = sum(keepers)>0;
                   residMatch = residMatch(logical(keepers),:); 
                 else
                  stayInLoop = 0;
                 end
              ii = ii + 1;
              end
          end     
      end               
    orphan_matches_TF{TF_id} = [TF_id*ones(size(residMatch,1),1),residMatch];
    mat_lifecycle_TF{TF_id}  = [TF_id*ones(size(mat_lifecycle,1),1), [mat_lifecycleC;mat_lifecycleS]];  
  end  
  
  % Stack match histories for firm-types, putting firm-type ID in col. 1
    agg_mat_lifecycle = double.empty(0,5*ceil(max_age./mm.pd_per_yr)+1);
  %  mat_lifecycle(:,1:6): [TF_id, year, match age, boy Z, eoy Z, sales]    
    agg_orphan_matches = double.empty(0,10);
    for TF_id = 1:NumTF
       agg_mat_lifecycle = [agg_mat_lifecycle; mat_lifecycle_TF{TF_id}];
       agg_orphan_matches =  [agg_orphan_matches;orphan_matches_TF{TF_id}]; 
    end
%% Refinements
    
% Drop matches with 0 sales from orphan list
 PosSale            = agg_orphan_matches(:,6)>0; % No sales and not first yr.
 agg_orphan_matches = agg_orphan_matches(logical(PosSale),:);
 
% Move orphans revealed to be dying new matches to agg_mat_lifecycle
 newDud =...
        logical((agg_orphan_matches(:,9)>1).*(agg_orphan_matches(:,9)<mm.pd_per_yr)==1);
 K2   = size(agg_mat_lifecycle,2);
 temp = [agg_orphan_matches(newDud,1),agg_orphan_matches(newDud,2),agg_orphan_matches(newDud,9),...
         agg_orphan_matches(newDud,7:8),agg_orphan_matches(newDud,5)];          
 agg_mat_lifecycle  = [agg_mat_lifecycle;[temp,zeros(sum(newDud),K2-6)]];
 
 % put remaining orphans in agg_orphan_matches. There should be just a few (<0.1%).
 % They could in principle be matched to the newly identified first-yr
 % matches. Doesn't seem worth the extra code.
 
 realOrphan         = logical(1-newDud); 
 agg_orphan_matches = agg_orphan_matches(realOrphan,:);
 %  agg_orphan_matches: (1) TF_id (2) year, (3) type, (4) firm_ID, (5) sales, 
 %           (6) shipments, (7) boy Z,(8) eoy Z,(9) match age,(10) firm age 
    
end  
    
    
  