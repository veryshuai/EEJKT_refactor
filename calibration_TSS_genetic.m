%This function initiates the genetic algorithm to solve for
%the parameters of the model.  Initial parameter guesses can be set here,
%as well as optional settings for the genetic algorithm.


% relation between parameters and the X vector:

%% relation between parameters and the X vector:
    
% mm.F_h       = exp(X(1)); % cost of maintaining a client- home 
% mm.scale_h   = X(2);       % Profit function scale parameter, home
% mm.scale_f   = X(2);       % Profit function scale parameter, foreign
% mm.ah        = X(4)*X(3);  % Beta function, home and foreign success parameter
% mm.bh        = X(4)*(1-X(3));% Beta function, home and foreign failure parameter
% D_z          = X(5)/mm.pd_per_yr; % jump size, match profit shock
% mm.L_b       = X(6)/mm.pd_per_yr; % new shipment hazard
% mm.gam       = X(7);       % Network effect parameter
% mm.cs_h      = exp(X(8));  % Cost scaling parameter, home market
% mm.sig_p     = X(9);       %standard deviation of productivity distribution
% mm.F_f       = exp(X(10)); % cost of maintaining a client- foreign
% mm.cs_f      = exp(X(11)); % Cost scaling parameter, foreign market
% mm.optimism  = X(12);      %parameter on prior distribution (positive means optimistic, negative pessamisitic) 


%%
         
%  theta = [-10.80310  -2.61962   0.26410   0.30794   0.38029  11.15223...
%            0.06463   1.18974   1.03113  -0.17300   9.55341   1];

%  theta = [-10.80310  -2.61962   0.26410   0.30794   0.38029  11.15223...
%            0.06463   2.18974   1.03113  -0.17300   10.55341   1];

% theta = [-12.0234360850160,-8.86822062058644,0.216000078376388,0.702389472395804,...
%        0.374214942055985,11.6388730026429,0.164637004707622,0.907628494202198,...
%        0.915592893141143,-0.254901443835102,6.77853609204470,2]; 
   
 theta = [-12.41111  -8.54950   0.19560   0.70351   0.27889   9.26793...
           0.19288   0.69155   1.16032  -0.23543   7.72443 2.27114];

   D0 = distance(theta)

  K = length(theta);
        
% population size
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




