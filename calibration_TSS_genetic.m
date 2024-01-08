%This function initiates the genetic algorithm to solve for
%the parameters of the model.  Initial parameter guesses can be set here,
%as well as optional settings for the genetic algorithm.



%% relation between parameters and the X vector:
    
% mm.F_h       = exp(X(1));  % cost of maintaining a client- home 
% mm.scale_h   = X(2);       % Domestic profit function scale parameter
% mm.scale_f   = X(2)+1;     % Export profit function scale parameter 
% mm.ah        = X(4)*X(3);  % Beta function, home (theta1) success parameter
% mm.bh        = X(4)*(1-X(3));% Beta function, home (theta1) failure parameter
% D_z          = X(5)/mm.pd_per_yr; % Jump size, match productivity shock
% mm.L_bF      = X(6)/mm.pd_per_yr; % Shipment order arrival hazard
% mm.gam       = X(7);       % Network effect parameter
% mm.cs_h      = exp(X(8));  % Cost scaling parameter, home market
% mm.sig_p     = X(9);       %standard deviation of productivity distribution
% mm.F_f       = exp(X(1));  % cost of maintaining a client- foreign
% mm.cs_f      = exp(X(10)); % Cost scaling parameter, foreign market
% mm.optimism  = 0;          %parameter on prior distribution 

%%

fprintf('\r\n STARTING A NEW RUN (known theta model): %s\r ', datestr(now,'yyyy_mmdd_HHMM'));
fprintf('\n Foreign profit scale set to twice home market profit scaler ');
fprintf('\n Degree distribution not targeted, minimizing original D \r');
fprintf('\n Updated hazard measure and hazard regression targeted \r');
fprintf('\n S = 50000 potential firms\r');


fileID1 = fopen('results/ga_running_output_NoLearning.txt','a');
  fprintf(fileID1,'\r\n STARTING A NEW RUN (known theta model) %s\n', datestr(now,'yyyy_mmdd_HHMM') );
  fprintf(fileID1,'\n Foreign profit scaler set to twice home market profit scaler \r');
  fprintf(fileID1,'\n Degree distribution not targeted, minimizing original D \r');
  fprintf(fileID1,'\n Updated hazard measure and hazard regression targeted \r');
  fprintf(fileID1,'\n S = 50000 potential firms \r');
fclose(fileID1);

fileID2 = fopen('results/ga_fitlog_NoLearning.txt','a');
  fprintf(fileID2,'\r\n STARTING A NEW RUN %s\n', datestr(now,'yyyy_mmdd_HHMM') );
  fprintf(fileID2,'\n Foreign profit scaler set to twice home market profit scaler \r');
  fprintf(fileID2,'\n Degree distribution not targeted, minimizing original D \r');
  fprintf(fileID2,'\n Updated hazard measure and hazard regression targeted \r');
  fprintf(fileID2,'\n S = 50000 potential firms \r');
fclose(fileID2);
 

% 
% theta = [-3.83253579377554	-19.6106680040131	0.137001306401593	0.282020343033900...
%          0.577076334796007	12.1194036685043	0.0492572810194975	4.88922911109957...
%          2.38047698896763	15.1377391760052]; % fit: 11.8453

% theta = [-3.76611438037509	-20.6927217338684	0.135477287531201	0.190633236901937...
%         0.579471364487730	11.3347872095937	0.0548954805115267	3.69372770061689...
%     	2.38225965393168	14.5381611202662];
    
theta = [-7.00428657247973 -23.0915400600166 0.129855165074604 0.147073857937515 0.637554188362434 12.5050138543584 0.0607676627106075 4.33149453421485 2.94198673283020 18.0202455396844];    

% X = theta;     

D0 = distance_nolearning(theta);

load('/storage/work/jxt32/EEJKT-codes/EEJKT_refactor_nolearn_11-2-23/Output/Optimization_2023_1108_0659_Run9of15.mat')
  
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
population(1,:) = theta;

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





