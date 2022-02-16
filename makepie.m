% Called from solve_v1, this function calculates the present value of
% successful matches. Hazards are per period (not per year).

function [myPi,pi_z,c_val] = makepie(pie_scale,st,Z,Q0,Q0_d,Q_z,Q_z_d,erg_pz,F,mm)

rh              = mm.r;             % time pref. rate, per period
del             = mm.delta;         % exogenous match death rate per period
de              = mm.eta;           % demand elasticity
pi_tol          = mm.pi_tolerance;  % stopping rule for loop below 
n_z             = 2*mm.z_size+1;    % number of demand shocks states
L_b             = mm.L_b;           % shipment hazard.

h = rh + del + L_b + abs(Q0(1,1)) + abs(Q_z(1,1)); %hazard of "something" happening, plus discount rate:
% match death, new shipment, change in the macro state, or buyer shock.

% payoff to shipment before adjusting for buyer state (z)
payoff_except_z = (1/de) * exp(pie_scale) * exp((de-1)*st(1,:)+st(2,:)); 
%sf is estimated scalar, st(1,:) is productivity, st(2,:) is macro state
payoff = repmat(exp(Z)',size(payoff_except_z,2),1) .* repmat(payoff_except_z',1,size(Z,1)); 
% rows differ by payoff_except_z (i.e., macro state and phi combinations),
% columns differ by exp(Z) value

% get expected profits for each type of buyer in each macro/prod state (row macro/prod, column demand shock) 
pi_z = 1 / (del + rh) * payoff; %initial guess on expected lifetime profit stream
c_val = pi_z; % continuation value, that is the value of having a relationship of a particular type.
%We call it the continuation value because we think of it as the value of
%keeping the relationship around just after a shipment.  It is net of the
%cost of maintaining the relationship, F, paid just after receiving a shipment.

%Value function iteration (no choice here, just exogenous shocks)
myEps = 1e12;
it = 0;
max_iter = 50000;
while myEps > pi_tol && it<=max_iter   % Iterate on contraction
        
    it = it + 1;
    if it == max_iter
        display('WARNING: makepie.m did not converge!');
    end

    %new value function
    c_val_gross_new = (Q0_d*c_val + c_val * Q_z_d' + L_b * pi_z) / h ;
      %first term:  Q0_d is the hazard of a macro shock change multiplied by the probability of a move into a new macro state * continuation value, 
      %second term: Q_z_d is the hazard of a demand shock change multiplied by the probability of a move into a new demand shock * continuation value 
      %third term: value if the next event is a shipment          
    c_val_new = max(-F + c_val_gross_new,0);  %net out the cost of maintaining the relationship.
    pi_z_new = payoff + c_val_new;
    
    myEps = max(max(abs(pi_z-pi_z_new)));
    pi_z = pi_z_new; 
    c_val = c_val_new; 
end 

%expectation of profits over buyer states (matching is random, so this is the opject buyers care about). 
myPi = pi_z*erg_pz; 

end %end function
