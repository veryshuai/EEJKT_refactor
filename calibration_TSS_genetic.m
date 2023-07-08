%This function initiates the genetic algorithm to solve for
%the parameters of the model.  Initial parameter guesses can be set here,
%as well as optional settings for the genetic algorithm.



%% relation between parameters and the X vector:
    
% mm.F_h       = exp(X(1));  % cost of maintaining a client- home 
% mm.scale_h   = X(2);       % Domestic profit function scale parameter
% mm.ah        = X(4)*X(3);  % Beta function, home (theta1) success parameter
% mm.bh        = X(4)*(1-X(3));% Beta function, home (theta1) failure parameter
% D_z          = X(5)/mm.pd_per_yr; % Jump size, match productivity shock
% mm.L_bF      = X(6)/mm.pd_per_yr; % Shipment order arrival hazard
% mm.gam       = X(7);       % Network effect parameter
% mm.cs_h      = exp(X(8));  % Cost scaling parameter, home market
% mm.sig_p     = X(9);       %standard deviation of productivity distribution
% mm.F_f       = exp(X(1)); % cost of maintaining a client- foreign
% mm.cs_f      = exp(X(10)); % Cost scaling parameter, foreign market
% mm.scale_f   = X(2);       % Export profit function scale parameter (same as home)
% mm.optimism  = 0; %parameter on prior distribution (positive means optimistic, negative pessamisitic)

%%

fprintf('\r\n STARTING A NEW RUN (restricted model): %s\n ', datestr(now,'yyyy_mmdd_HHMM'));

fileID1 = fopen('results/ga_running_output_restricted.txt','a');
 fprintf(fileID1,'\r\n STARTING A NEW RUN (restricted model) %s\n', datestr(now,'yyyy_mmdd_HHMM') );
fclose(fileID1);

fileID2 = fopen('results/ga_fitlog_restricted','a');
  fprintf(fileID2,'\r\n STARTING A NEW RUN %s\n', datestr(now,'yyyy_mmdd_HHMM') );
fclose(fileID2);
          
 
 theta = [-2.92217885480810,-21.7329307485543,0.163673081953874,0.142896261678048,...
         0.253189491294144,21.0101797166516,0.119991356279124,7.54692070316604,...
         3.14414508598657,15.3219457708342];
% theta for restricted model. Alternative fit metric: 12.5579
     
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





