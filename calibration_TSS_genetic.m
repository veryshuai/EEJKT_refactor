%This function initiates the genetic algorithm to solve for
%the parameters of the model.  Initial parameter guesses can be set here,
%as well as optional settings for the genetic algorithm.



%% relation between parameters and the X vector:
    
% mm.F_h       = exp(X(1));  % cost of maintaining a client- home 
% mm.scale_h   = X(2);       % Domestic profit function scale parameter
% mm.scale_f   = X(2);       % Export profit function scale parameter (same as home)
% mm.ah        = X(4)*X(3);  % Beta function, home (theta1) success parameter
% mm.bh        = X(4)*(1-X(3));% Beta function, home (theta1) failure parameter
% D_z          = X(5)/mm.pd_per_yr; % Jump size, match productivity shock
% mm.L_bF      = X(6)/mm.pd_per_yr; % Shipment order arrival hazard
% mm.gam       = X(7);       % Network effect parameter
% mm.cs_h      = exp(X(8));  % Cost scaling parameter, home market
% mm.sig_p     = X(9);       %standard deviation of productivity distribution
% mm.F_f       = exp(X(1)); % cost of maintaining a client- foreign
% mm.cs_f      = exp(X(10)); % Cost scaling parameter, foreign market
% mm.optimism  = 0; %parameter on prior distribution (positive means optimistic, negative pessamisitic)


%%

fprintf('\r\n STARTING A NEW RUN: %s\n ', datestr(now,'yyyy_mmdd_HHMM'));
fileID1 = fopen('results/ga_running_output_restricted.txt','a');
 fprintf(fileID1,'\r\n STARTING A NEW RUN %s\n', datestr(now,'yyyy_mmdd_HHMM') );
fclose(fileID1);

fileID2 = fopen('results/ga_fitlog_restricted.txt','a');
  fprintf(fileID2,'\r\n STARTING A NEW RUN %s\n', datestr(now,'yyyy_mmdd_HHMM') );
fclose(fileID2);
          

%  theta = [ -2.88300 -18.24037   0.16754   0.14741   0.25041  27.46668...
%   0.09661   6.40443   2.54228  11.83779];
%  % alt fit metric: 13.137
%  
%  theta = [ -2.83343 -20.26304   0.17207   0.15619   0.29423  23.17569...
%             0.09783   7.18610   2.96590  15.62184 ];
 % alt fit = 12.873 on roar collab
        
 
 theta =[ -2.79107 -20.26304   0.18305   0.15619   0.29423  20.31397...
           0.09783   7.18610   2.96590  15.62184 ];
  % alt fit = 12.69404 on roar collab      
       
% X = theta;     

D0 = distance(theta);
  
% population size
K = length(theta);
PS = 2*K;
disp(['Population size ' num2str(PS)])
H = rand(PS,K);

bounds = [theta - 0.30*abs(theta); theta + 0.30*abs(theta)]';

lb = bounds(:,1);
ub = bounds(:,2);
IR = bounds';
lb = lb(:,ones(PS,1))';
ub = ub(:,ones(PS,1))';
rng = abs(ub-lb);
X0 = lb + rng.*H;
clear lb ub rng
% initial population
population = X0;

cf = 0.30;          % crossover fraction
EC = 4;             % elite count
shrink = 1;
SF = 5; 

GN = 10; % number of generations per run
R =  15; % total number of runs

scale0 = 0.15;
t = 1;
flag = true;
while t<=R
    scale = ((R-t)/R)*scale0;
    mut_fn = @(parents,options,GenomeLength,FitnessFcn,state,thisScore,...
        thisPopulation) ...
        TSS_mutationgaussian(parents,options,GenomeLength,FitnessFcn,...
        state,thisScore,thisPopulation,scale,shrink,t,R,GN);
    disp(['optimization run ',num2str(t),' of ',num2str(R)])
    options = gaoptimset('TolFun',1e-20,'StallGenLimit',800,'MutationFcn',mut_fn,...
        'Generations',GN,'PopulationSize',PS,'InitialPopulation',population,'Display','iter',...
        'UseParallel','always','Vectorized', 'off','CrossoverFraction',cf,'EliteCount',EC,...
        'PopInitRange',IR,'FitnessScaling',@fitscalingrank);
    disp('optimization starting')
    tic
    [x,fval_ga,exitflag,output,population] = ga(@(X) distance(X),K,[],[],[],[],...
        [],[],[],options);
    FileName = ['Output/' 'Optimization','_',datestr(now,'yyyy_mmdd_HHMM'),...
        '_Run',num2str(t),'of',num2str(R)];
    save(FileName,'x','fval_ga','population')
    time_elapsed1 = toc;
    disp(['time elapsed in total: '...
    num2str(time_elapsed1/60) ...
    ' minutes'])
    disp(['time elapsed per generation: ' ...
    num2str(time_elapsed1/(GN+1)) ...
    ' seconds'])
   pop_range_t = [min(population',[],2) max(population',[],2)];
   IRt = IR';
   pop_range0 = round(100*mean(abs(pop_range_t(:,2)-pop_range_t(:,1))./abs(IRt(:,2)-IRt(:,1))));
   disp(['Average (max minus min) population range as percentage of initial range: ' num2str(pop_range0)])
   pop_range1 = round(100*(sum(pop_range_t(:,1)<IRt(:,1)) + sum(pop_range_t(:,2)>IRt(:,2)))/(2*K));
   disp(['Percentage of initial bounds that are violated by at least one parameter vector in the current population: ' ...
       num2str(pop_range1)])
   
    t = t + 1;
end
theta = x';





