% Build life cycle for each new match

function [mat_lifecycle_TF,orphan_matches_TF] = lifecycle(NumTF,TF_matdat,max_age)

orphan_matches_TF = cell(NumTF,1);
mat_lifecycle_TF = cell(NumTF,1);
firstYrCount     = zeros(NumTF,1);

% Build life cycle for each new match

parfor TF_id=1:NumTF

    % fprintf('\r firm ID = %0.0f\n',TF_id);
    all_age = TF_matdat{TF_id};
    % all_age: (1) year, (2) type, (3) firm_ID, (4) sales, (5) shipments,
    %          (6) boy Z,(7) eoy Z,(8) match age,(9) firm age
    new_match_TF = all_age(:,8)<=1;
    max_age_TF = max(all_age(:,8));
    N_TF = sum(new_match_TF); % number of new matches for this firm
    mat_lifecycle = zeros(N_TF,5*max_age);
    mat_lifecycle(:,1:5) = [all_age(new_match_TF,1),all_age(new_match_TF,8),...
        all_age(new_match_TF,6:7),all_age(new_match_TF,4)];
    % load one-year-olds into first 5 cols of mat_lifecycle
    % mat_lifecycle(:,1:5): [year, match age (<=1), boy Z, eoy Z, sales]

    orphan_matches = double.empty(0,9);
    for aa = 2:max_age_TF
        % Grab the aa-year-old matches
        aa_yr_old = all_age(:,8) == aa;
        aa_cohort = all_age(aa_yr_old,:);

        stayInLoop = 1;
        while stayInLoop > 0 && size(aa_cohort,1) > 0

            if firstYrCount(TF_id,1)==0
                fprintf('\r type-firm ID = %0.0f, number of 1st year matches = %0.0f\n',[TF_id,firstYrCount(TF_id,1)]);
                stayInLoop = 0;
            end

            ii =1;
            while ii <= N_TF %looping over matches of age aa for firm-type TF_id
                fprintf('\r type-firm ID = %0.0f, match age = %0.0f, match number = %0.0f',[TF_id,aa,ii]);
                lag_yr   = mat_lifecycle(ii,5*(aa-2)+1);
                lag_age  = mat_lifecycle(ii,5*(aa-2)+2);
                lag_eoyZ = mat_lifecycle(ii,5*(aa-2)+4);

                % Find compatible matches to splice with last year's match ii
                flg = (lag_eoyZ>0).*(aa_cohort(:,6)==lag_eoyZ).*...
                    (aa_cohort(:,8)==lag_age + 1).*(aa_cohort(:,1)==lag_yr + 1);

                % If compatible matches are identified, put first one in the
                % relevant bloc of mat_lifecycle and remove it from aa_cohort
                if sum(flg)>0
                    mat_cont = find(flg==1,1); % first compatible match in aa_cohort
                    lb = 5*(aa-1)+1;
                    ub = 5*(aa-1)+5;
                    mat_lifecycle(ii,lb:ub) = [aa_cohort(mat_cont,1),aa_cohort(mat_cont,8),...
                        aa_cohort(mat_cont,6:7),aa_cohort(mat_cont,4)];

                    % remove matched record from aa_cohort
                    keepers = ones(size(aa_cohort,1),1);
                    keepers(mat_cont,1) = 0;
                    stayInLoop = sum(keepers)>0;
                    aa_cohort = aa_cohort(logical(keepers),:);
                else
                    stayInLoop = 0;
                end
                ii = ii + 1;
            end
        end
        orphans = size(aa_cohort,1) ;
        if orphans > 0 % document leftover matches
            orphan_matches = [orphan_matches;aa_cohort];
            fileID5 = fopen('results/EEJKT_orphan_log.txt','a');
            fprintf(fileID5,'\r\n %0.0f unmatched record(s), %0.0f year-old firm, type %3.6f',[orphans,aa,TF_id,]);
            fclose(fileID5);
        end
    end
    orphan_matches_TF{TF_id} = orphan_matches;
    mat_lifecycle_TF{TF_id} = [TF_id*ones(size(mat_lifecycle,1),1), mat_lifecycle];
end