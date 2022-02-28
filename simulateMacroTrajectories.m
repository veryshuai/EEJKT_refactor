function [macro_state_f, macro_state_h] = simulateMacroTrajectories(periods, policy)

macro_state_f    = zeros(periods,1);
macro_state_f(1) = 8; %  start macro trajectory at midpoint of distribution
macro_state_h    = zeros(periods,1);
macro_state_h(1) = 8; %  start macro trajectory at midpoint of distribution
for t = 2:periods
    macro_state_f(t) = find(policy.pmat_cum_msf(macro_state_f(t-1),:)>rand(1,1),1,'first'); % update macro state
    macro_state_h(t) = find(policy.pmat_cum_msh(macro_state_h(t-1),:)>rand(1,1),1,'first'); % update macro state
end

end