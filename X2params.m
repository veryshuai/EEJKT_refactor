%% 
%  This file takes an X vector passed by the optimizaiton program
%  and turns it into parameters used in our simulation

%read in parameter names from the X vector


%% current mapping 

F_h      =  exp(X(1));     % home match fixed cost
scale_h  =  X(2);          % log of home profit function scalar       
delta    =  0.326;         % exogenous match death hazard (per year)                                 
beta     =  1;             % cost function convexity parameter               
ah       =  X(4)*X(3);     % X(3) is mean of beta dist ah / (ah + bh)
bh       =  X(4)*(1-X(3)); % X(4) is (ah + bh)        
L_z      =  4;             % annual buyer shock jump hazard (once per quarter)           
D_z      =  X(5);          % buyer shock jump size            
L_b      =  X(6);          % shipment hazard            
gam      =  X(7);          % network effect
cs_h     =  exp(X(8));     % cost function scalar, home market  
sig_p    =  X(9);          % std. dev. of log productivity shock 
F_f      =  exp(X(10));    % foreign match fixed cost
cs_f     =  exp(X(11));    % cost function scalar, foreign market
scale_f  =  X(12);         % log of foreign profit function scalar 

% delta    =  X(2);       % exogenous match death hazard (per year)   
% scale_f  =  X(3)+X(4);  % log of foreign profit function scalar 
% beta     =  1/(X(5)-1); % cost function convexity parameter 
% L_z      =  X(8);       % buyer shock jump hazard    
  

%% parameters no longer used that still require values

ag         =  .5;
bg         =  .5;
L_p        =  0;
D_p        =  0;
alp        =  0;
