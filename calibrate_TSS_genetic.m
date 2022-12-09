%This function initiates the genetic algorithm to solve for
%the parameters of the model.  Initial parameter guesses can be set here,
%as well as optional settings for the genetic algorithm.


% relation between parameters and the X vector:

%% relation between parameters and the X vector:
    
% F_h      =  exp(X(1));     % home match fixed cost
% scale_h  =  X(2);          % log of home profit function scalar       
% delta    =  0.326;         % exogenous match death hazard (per year)                                 
% beta     = 1;              % cost function convexity parameter               
% ah       =  X(4)*X(3);     % X(3) is mean of beta dist ah / (ah + bh)
% bh       =  X(4)*(1-X(3)); % X(4) is (ah + bh)        
% L_z      =  4;             % buyer shock jump hazard (once per quarter)           
% D_z      =  X(5);          % buyer shock jump size            
% L_b      =  X(6);          % shipment hazard            
% gam      =  X(7);          % network effect
% cs_h     =  exp(X(8));     % cost function scalar, home market  
% sig_p    =  X(9);          % std. dev. of log productivity shock 
% F_f      =  exp(X(10));    % foreign match fixed cost
% cs_f     =  exp(X(11));    % cost function scalar, foreign market
% scale_f  =  X(12);         % log of foreign profit function scalar 

%%
% 
%   theta =  [-8.53630  -7.02531   0.23560   0.87633   0.40764  12.76228...
%             0.08514   1.51229   1.24776  -0.28887   8.37010  -9.03099]; 
%            % fit: 12.5248646781
%            
%    theta = [-10.68006  -7.48486   0.24552   1.13575   0.57243  11.71367 ...
%   0.06442   1.40703   1.37038  -0.25322   9.82913  -9.65594]; %fit = 11.2955538841          
%  

% theta = [-10.44911  -7.71754   0.24459   0.83398   0.55360  11.41532
%        0.07376   1.40956   1.35978  -0.22084   9.15662  -9.59909]; % fit = 11.2573441966 
   
 theta =  [-8.53630  -7.02531   0.23560   0.87633   0.40764  12.76228...
            0.08514   1.51229   1.24776  -0.28887   8.37010  -9.03099]; % fit = 11.6206


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

GN = 50; % number of generations per run
R =  5; % total number of runs

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





