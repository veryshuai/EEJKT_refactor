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

fprintf('\r\n STARTING A NEW RUN: %s\r ', datestr(now,'yyyy_mmdd_HHMM'));
fprintf('\n Unconstrained profit function scalars\r');
fprintf('\n Degree distribution not targeted, minimizing original D \r');
fprintf('\n Updated hazard measure and hazard regression targeted \r');
fprintf('\n S = 50000 potential firms\r');



fileID1 = fopen('results/ga_running_output_2sclNewHaz.txt','a');
  fprintf(fileID1,'\r\n STARTING A NEW RUN %s\n', datestr(now,'yyyy_mmdd_HHMM') );
  fprintf(fileID1,'\n Unconstrained profit function scalars\r');
  fprintf(fileID1,'\n Degree distribution not targeted, minimizing original D \r');
  fprintf(fileID1,'\n Updated hazard measure and hazard regression targeted \r');
  fprintf(fileID1,'\n Unconstrained profit function scalars\r');
  fprintf(fileID1,'\n S = 50000 potential firms \r');
fclose(fileID1);

fileID2 = fopen('results/ga_fitlog_2sclNewHaz.txt','a');
  fprintf(fileID2,'\r\n STARTING A NEW RUN %s\n', datestr(now,'yyyy_mmdd_HHMM') );
  fprintf(fileID2,'\n Unconstrained profit function scalars\r');
  fprintf(fileID2,'\n Degree distribution not targeted, minimizing original D \r');
  fprintf(fileID2,'\n Updated hazard measure and hazard regression targeted \r');
  fprintf(fileID2,'\n S = 50000 potential firms \r');
fclose(fileID2);
 

% theta =...
%     [ -6.52390 -19.91885   0.07778   0.22093   0.50050  11.43869...
%        0.05174   3.89969   2.63658  11.02215]; % fit = 11.3765
% 
% theta =[-5.73144888007375	-16.1844310171305    0.109738437680579...
%          0.242942844550223    0.605350029387283 11.0996895678257...
%          0.0467130185363216   3.07428217391970	 2.20454924542622...
%         11.1042769087445]; % fit = 10.9807
    
% theta = [-7.73724 -12.76932   0.11799   0.24781   0.46288   8.34401...
%           0.05868   4.17410   1.68781  12.07565 -11.76932];  % fit = 10.75460
      
% theta = [-6.13067795576060	-14.7878140561022	0.116104623271765...
%          0.287168108117947	0.551351332778652	9.70261839734432...
%          0.0646865540473156	3.90851546299691	1.93735060502043...
%          14.1444949193409	-12.9121281421765];  % fit = 10.4748
     
theta = [-5.44086676636325,-19.4382242330205,0.158825726539545,0.367383665973718,...
    0.483330444401690,10.6297324206023,0.0469201052969252,5.20915076877990,...
    2.18909944116805,11.7998950817073]; % fit = 10.4266     
      

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





