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
% mm.optimism  = 0 ; %parameter on prior distribution 

%%

fprintf('\r\n STARTING A NEW RUN (No learning): %s\n ', datestr(now,'yyyy_mmdd_HHMM'));

fileID1 = fopen('results/ga_running_output_nolearning.txt','a');
 fprintf(fileID1,'\r\n STARTING A NEW RUN (no learning) %s\n', datestr(now,'yyyy_mmdd_HHMM') );
fclose(fileID1);

fileID2 = fopen('results/ga_fitlog_nolearning.txt','a');
  fprintf(fileID2,'\r\n STARTING A NEW RUN %s\n', datestr(now,'yyyy_mmdd_HHMM') );
fclose(fileID2);
 
<<<<<<< HEAD
     
%theta = [-4.92444290760212,-27.9412276522655,0.173719817766680,...
%         0.149174564820941,0.435666463937585,14.6171833559373,...
%    0.131787405862920,8.81158882786645,3.37967374928709,10.6400808158221];
% BASELINE alt fit metric: 11.969925590604335

% theta = [-4.74060620625034	-46.4781761953731	0.0831810159195763	0.174061271943099	0.407684805092889	12.5014876279670 ...
%     0.200733646668880	9.28932144380271	4.68302354232828	15.2292368534740];
% % NO LEARNING alt fit metric: 12.2925255086572

theta = [-4.74336931061665	-39.3774942028084	0.102488335719205	0.124887701160301	0.365998047325849	13.0817977607945 ...
    0.180101053434558	9.89876415078752	4.00833277085712	13.6171487940690];
% NO LEARNING alt fit metric: 12.1572058735029

=======

% theta = [-4.92444290760212,-27.9412276522655,0.173719817766680,...
%          0.149174564820941,0.435666463937585,14.6171833559373,...
%     0.131787405862920,8.81158882786645,3.37967374928709,10.6400808158221,0.149174564820941];
% alt fit metric: 11.969925590604335 (with quadratic search costs and S=10000)

% theta = [-5.54845801486006,-27.9276029273286,0.183454179425287,0.134091994275606,...
%     0.312956304422016,13.0826080407669,0.115310320152066,9.69219215273992,...
%     3.45209285533265,11.7866318512518,0.118246569286443];
% alt fit metric: 11.650 with optimism parameter estimated

% theta = [-5.51055707725502,-22.9820628073306,0.152309883146146,0.150717626933370,...
%     0.197929034960503,11.9969904367609,0.0853512405334382,7.75034563017993,...
%      3.32790809640092,11.3087581117842,-26.9147671844486];
% alt fit metric: 11.8773 with no optimism parameter, separate profit scalars
 

% theta = [5,-22.9820628073306,0.152309883146146,0.150717626933370,...
%     0.197929034960503,11.9969904367609,0.0853512405334382,7.75034563017993,...
%      3.32790809640092,11.3087581117842,-26.9147671844486];
% alt fit metric: 11.8773 with no optimism parameter, separate profit scalars
 

% theta =[...
% -5.51055707725502,-22.9820628073306,0.152309883146146,0.150717626933370,...
%  0.197929034960503,11.9969904367609,0.0853512405334382,7.75034563017993,...
% 3.32790809640092,11.3087581117842,-26.9147671844486];
% baseline model, shipment hazard ratio 3.4, Marcela's new figures
% alt fit metric: 12.26
>>>>>>> TryAlt_refactor
% X = theta;     

theta = [...
  -6.48373445693280,-27.1116388656176,0.119431181295694,0.0976175374467408,...
   0.247110496961648,14.2241812754387,0.0630667418995650,8.80568820896566,...
   3.36769544702316,12.9980245788823,-24.9408986889662];
% alt fit metric: 12.42 with shipment hazard at home 3 X hazard abroad
% alt fit metric: 13.13 with shipment hazard at home 2 X hazard abroad

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
