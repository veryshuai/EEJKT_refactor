
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
% lowb = [-10; -10;   0.001;   0.1; 0.005;  0.001; -5;  -10; 0.001;  -10; -10; -10];
% upb =  [10;   20;   0.999; 100.0;   5.0;   100;   5;   20; 5.00;    10;  20;  20];

  theta = [ -6.69905  -6.88534   0.22734   1.15397   0.40905  12.99814 ...
             0.11426   1.88804   1.27705  -0.35441   7.44522  -9.96677]; 
           %  fit metric = 11.8476786561
           
  theta = [-10.68006  -7.48486   0.24552   1.13575   0.57243  11.71367 ...
  0.06442   1.40703   1.37038  -0.25322   9.82913  -9.65594]; %fit = 11.2955538841          
        
  Fit0 = distance(theta)      
        
% population size
PS = 2*length(theta);
disp(['Population size ' num2str(PS)])
H = rand(PS,length(theta));

bounds = [theta - 0.30*abs(theta); theta + 0.30*abs(theta)]';

lb = bounds(:,1);
ub = bounds(:,2);
IR = bounds';
lb = lb(:,ones(PS,1))';
ub = ub(:,ones(PS,1))';
rng = abs(ub-lb);
X0 = lb + rng.*H;
clear lb ub rng
lowb = bounds(:,1);
upb  = bounds(:,2);
population = X0;
    
    % Parallel setup
    clc
    
    % Put numbers in long format for printing
    format long;

    % random seed
    rng(80085);
    
% Options for ga only 
  
    gaoptions = gaoptimset('Display','iter','PopulationSize',PS,'Generations',1500,...
    'StallTimeLimit',86400,'TimeLimit',360000,'MutationFcn',@mutationadaptfeasible,...
    'FitnessScalingFcn',@fitscalingrank,'InitialPopulation',population,'UseParallel','always',...
    'PlotFcns',@gaplotbestf,'EliteCount',0); 

 % [X,fval,exitflag] = ga(@(X) distance(X),12,[],[],[],[],lowb, upb,[],gaoptions);
  [X,fval,exitflag] = ga(@(X) distance(X),12,[],[],[],[],lowb,upb,[],gaoptions);


    % Save results
    save est_results_genetic;

