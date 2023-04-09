function simMoms = simulateMomentsMain(policy,mm)

rng(80085,'twister');

[macro_state_f, macro_state_h] = simulateMacroTrajectories(mm, policy);

sim_out = cell(mm.N_pt,1);

seeds = randi(1e6,size(mm.Phi,1),2);

parfor pt_ndx = 1:mm.N_pt
% for pt_ndx = 1:1:mm.N_pt 
%for pt_ndx = 90
% for pt_ndx = 106
% for pt_ndx = 108
% for pt_ndx = 101
% for pt_ndx = 113

    rng(seeds(mm.pt_type(pt_ndx,1),1),'twister');
    seed_crand(seeds(mm.pt_type(pt_ndx,1),2));

    if mm.sim_firm_num_by_prod_succ_type(pt_ndx)>0

        [sim_out{pt_ndx}] = simulateForeignMatches(pt_ndx,macro_state_f, mm, policy);

        [sim_out{pt_ndx}] = simulateHomeMatches(pt_ndx,macro_state_h, mm, policy,sim_out{pt_ndx});

        sim_out{pt_ndx} = splice_hf(sim_out{pt_ndx},policy,mm,pt_ndx);

    end

end

%% Uncomment commands below to generate data for spot checks
% check_type = mm.check_type;
% if sim_out{check_type}.domfirm_count > 0
%   check_cell_H = sim_out{mm.check_type}.iterH_check;
%   check_count_H = sim_out{mm.check_type}.stackH;
%    save 'iter_out_HomeChecks.mat' 'check_type' 'check_cell_H' 'check_count_H';  
% end
% if sim_out{check_type}.exptr_count >0
%   check_cell_F = sim_out{mm.check_type}.iterF_check;
%   check_count_F = sim_out{mm.check_type}.stackF;
%  save 'iter_out_XptrChecks.mat' 'check_type' 'check_cell_F' 'check_count_F';  
% end

% ALSO UNCOMMENT: lines 119-133 in simulateHomeMatchesInnerSim.m
%                 lines 74-86 in simulateForeignMatches.m
%%
sim_cum = aggregateSimulatedData(sim_out,mm);
simMoms = calculateSimulatedMoments(sim_cum,mm);
