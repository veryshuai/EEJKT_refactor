function plots_v2(simMoms)

% simMoms.agg_time_gaps:
%          (1) firm_ID, (2) period w/in interval, (3) gap 
%          (4) # new meetings,(5) t (6) cum. meetings, (7) cum succeses

%%      plot degree distribution
        figure(1)      
        scatter(simMoms.log_matchesD,simMoms.log_compCDF_D)
        
%%     plot simulated meeting hazards (based directly on simulated data)

        % Create grid
        grdsize   = 20;
        ln_succ_rate = log(1+(simMoms.agg_time_gaps(:,7)./simMoms.agg_time_gaps(:,6)));
%       ln_csucc     = log(1+simMoms.agg_time_gaps(:,7)); 
        ln_cmeet     = log(1+simMoms.agg_time_gaps(:,6)); 
        succ_grid    = (min(ln_succ_rate):((max(ln_succ_rate)-min(ln_succ_rate))/grdsize):max(ln_succ_rate))';
%       csucc_grid   = (min(ln_csucc):((max(ln_csucc)-min(ln_csucc))/grdsize):max(ln_csucc))';
        cmeet_grid   = (min(ln_cmeet):((max(ln_cmeet)-min(ln_cmeet))/grdsize):max(ln_cmeet))';

        % construct simulated hazards
        MeetHaz    = 1./simMoms.agg_time_gaps(:,3);
        ln_MeetHaz = log(MeetHaz);

        % Move each success rate and cumulative number of successes to
        % closest grid point; construct gridded hazard surface.
        gridVec   = ones(grdsize+1,1);
        gridHaz   = zeros(size(ln_succ_rate,1),3);
        gridLnHaz = zeros(size(ln_succ_rate,1),3);
        for i = 1:size(ln_succ_rate,1)
            % success rate on the grid
            Sndx = find(abs(succ_grid - ln_succ_rate(i)*gridVec) == min(abs(succ_grid - ln_succ_rate(i)*gridVec)));

            % cum successes on the grid
            % Cndx = find(abs(csucc_grid - ln_csucc(i)*gridVec) == min(abs(csucc_grid - ln_csucc(i)*gridVec)));

            % cum meetings on the grid           
            Mndx = find(abs(cmeet_grid - ln_cmeet(i)*gridVec) == min(abs(cmeet_grid - ln_cmeet(i)*gridVec)));

%         gridHaz(i,:)   =  [MeetHaz(i),Sndx,Cndx];
%         gridLnHaz(i,:) =  [ln_MeetHaz(i),Sndx,Cndx];

         gridHaz(i,:)   =  [MeetHaz(i),Sndx,Mndx];
         gridLnHaz(i,:) =  [ln_MeetHaz(i),Sndx,Mndx];
         
        end
        
        % Average (etc.) simulated meeting hazards at each populated grid point
        mean_LnHaz = zeros(grdsize+1,grdsize+1);
        mean_haz   = zeros(grdsize+1,grdsize+1);
        median_haz = zeros(grdsize+1,grdsize+1);
        haz_count  = zeros(grdsize+1,grdsize+1);
        popCell    = zeros(grdsize+1,grdsize+1);
        for i = 1:grdsize+1 % successs rate index
            for j=1:grdsize+1 % cumulative meetings index
                cells_ij = (gridHaz(:,2)==i).*(gridHaz(:,3)==j)==1;
                if sum(cells_ij)>0
                popCell(i,j)   = 1;    
                mean_haz(i,j)   = mean(gridHaz(cells_ij,1));
                mean_LnHaz(i,j) = mean(gridLnHaz(cells_ij,1));
                median_haz(i,j) = median(gridHaz(cells_ij,1));
                haz_count(i,j)  = size(gridHaz(cells_ij,1),1);
                end
            end
        end


         figure(2)
         subplot(3,1,1)
         surf(succ_grid,cmeet_grid,mean_haz); 
         title('mean simulated hazards at populated grid points')
         ylabel('log(cum. meetings)')
         xlabel('log(1+success rate)')
         zlabel('avg. match hazard')
         
         subplot(3,1,2)
         surf(succ_grid,cmeet_grid,median_haz); 
         title('median simulated hazards at populated grid points')
         ylabel('log(cum. meetings)')
         xlabel('log(1+success rate)')
         zlabel('median match hazard')
         
         subplot(3,1,3)
         surf(succ_grid,cmeet_grid,mean_LnHaz); 
         title('mean simulated log hazards at populated grid points')
         ylabel('log(cum. meetings)')
         xlabel('log(1+success rate)')
         zlabel('avg. log match hazard')


         
%% plot histogram of frequencies for simulated meeting hazards
         figure(3)
         surf(succ_grid,cmeet_grid,haz_count); 
       % surf(succ_grid,csucc_grid,popCell); 
         title('frequency of simulated positive match hazards, by state')
         ylabel('log(cum. meetings)')
         xlabel('log(1+success rate)')
         zlabel('number of obs')  
   
%% plot the predicted meeting hazard as a functions of # clients and success rate    

%   without dummy for new exporter (copied from target_stats.m, except intercept)
%    b_hazDAT   = [-0.3528,-0.8181,0.3117,-1.1323,2.4514,-0.7082]; % [mean dep. var, ln(1+a), ln(1+a)^2, ln(1+r), ln(1+r)^2, ln(1+a)*ln(1+r)] 
    b_hazDAT   = [-3.051, 0.837];
    dat_haz    = zeros(grdsize+1,grdsize+1);
    mod_haz    = zeros(grdsize+1,grdsize+1);
    
        for i = 1:grdsize+1  % success rate
            for j = 1: grdsize+1  % cum success 
               
                dat_haz(i,j) = b_hazDAT(1) + b_hazDAT(2)*cmeet_grid(i);
 %              dat_haz(i,j) = b_hazDAT(1) + b_hazDAT(2)*csucc_grid(i) ...                  
%                    +b_hazDAT(3)*succ_grid(i)^2 + b_hazDAT(4)*csucc_grid(j)...
%                    +b_hazDAT(5)*csucc_grid(j)^2  + b_hazDAT(6)*succ_grid(i)*csucc_grid(j)   ;
%   
               mod_haz(i,j) =  simMoms.b_haz(1) + simMoms.b_haz(2)*cmeet_grid(i);
%              mod_haz(i,j) =  simMoms.b_haz(1) + simMoms.b_haz(2)*succ_grid(i) ...
%                    +simMoms.b_haz(3)*succ_grid(i)^2 + simMoms.b_haz(4)*csucc_grid(j)...
%                    +simMoms.b_haz(5)*csucc_grid(j)^2  + simMoms.b_haz(6)*succ_grid(i)*csucc_grid(j)   ;
                    
            end
        end
            

%% constructed data-based versus model-based predicted match hazards

% limit the predictions to states that are populated in the simulated data 
tdat_haz = dat_haz.*popCell;
tmod_haz =  mod_haz.*popCell;
test = tmod_haz - tdat_haz;

figure(5)
subplot(3,1,1)
        surf(succ_grid,cmeet_grid,tdat_haz)
        title('data-based predicted log match hazard')
        ylabel('log(cum. meetings)')
        xlabel('log(1+success rate)')
        zlabel('log match hazard')
subplot(3,1,2)        
           surf(succ_grid,cmeet_grid,tmod_haz)
        title('model-based predicted log match hazard')
        ylabel('log(cum. meetings)')
        xlabel('log(1+success rate)')
        zlabel('log match hazard')       
subplot(3,1,3)  
 surf(succ_grid,cmeet_grid,test)
         title('difference in predicted log match hazards')
        ylabel('log(cum. meetings)')
        xlabel('log(1+success rate)')
        zlabel('model-data difference in log match hazard')
 

end       


        
