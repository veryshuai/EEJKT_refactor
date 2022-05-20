function simMoms = simulateMomentsMain(policy,mm)

rng(80085,'twister');

[macro_state_f, macro_state_h] = simulateMacroTrajectories(mm, policy);

sim_out = cell(mm.N_pt,1);

seeds = randi(1e6,size(mm.Phi,1),2);

parfor pt_ndx = 1:1:mm.N_pt
% for pt_ndx = 1:1:mm.N_pt % use this for loop for debugging only
%     
%      if pt_ndx == 105
%           pause;
%      end

    rng(seeds(mm.pt_type(pt_ndx,1),1),'twister');
    seed_crand(seeds(mm.pt_type(pt_ndx,1),2));

    if mm.sim_firm_num_by_prod_succ_type(pt_ndx)>0

        sim_out{pt_ndx} = simulateForeignMatches(pt_ndx,macro_state_f, mm, policy);

        sim_out{pt_ndx} = simulateHomeMatches(pt_ndx,macro_state_h, mm, policy,sim_out{pt_ndx});

        sim_out{pt_ndx} = splice_hf(sim_out{pt_ndx},policy,mm);

    end

end

sim_cum = aggregateSimulatedData(sim_out,mm);

simMoms = calculateSimulatedMoments(sim_cum);
