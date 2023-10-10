function FirmCount = degree_dist(match_list,mm)
% match_list: [period, type, firm_ID, sales, shipments, boy Z, eoy Z, match age, firm age]
   
   maxYear     = max(match_list(:,1));        
   match_TF    = match_list(:,2) + 0.0001*match_list(:,3); % firm type by firmID identifier
   TF_list     = unique(match_TF);
   NumTF       = size(TF_list,1);    
     
 % count matches for each firm in each year
   firmMatchCount = zeros(NumTF,maxYear);
   for TF_id = 1:NumTF 
     ff = match_TF == ones(size(match_list,1),1).*TF_list(TF_id);
     firmRecs = match_list(ff,:); 
       for yr=1:maxYear
        firmMatchCount(TF_id,yr) = sum(firmRecs(:,1)==yr);
       end    
   end

   maxMatch  = min(mm.max_match, max(max(firmMatchCount)));     
 % convert to frequency counts of # firms with each possible match count
   MatFreq = zeros(maxMatch,maxYear);  
   for yr=1:maxYear
     if sum(firmMatchCount(:,yr))>0 % counts excluding duds
     %the next line replicates firmMatchCount(:,yr) across columns up to
     %max match.  The next term makes a matrix where each row contains the
     %integers from 1 to maxMatch.  Finally, the difference is taken, and
     %the logical is true when the difference is zero.  Finally the sum is
     %taken to make a histogram of firm counts with each possible number of
     %matches. This histogram is created seperately for each year.

     MatFreq(:,yr) = sum(firmMatchCount(:,yr)*ones(1,mm.max_match) - ones(size(firmMatchCount(:,yr),1),1)*(1:mm.max_match)==0,1);
     end
   end 
   
% FirmCount = sum(MatFreq(:,1:end),2); % including burn-in years
  FirmCount = sum(MatFreq(:,mm.burn+1:end),2); % excluding burn-in years

end