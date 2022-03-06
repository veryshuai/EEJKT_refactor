%function moms = simulateMomentsMain(mm,policy)

rng(80085,'twister');
[macro_state_f, macro_state_h] = simulateMacroTrajectories(mm, policy);

sim_out_f = cell(mm.N_pt,1);
sim_out_h = cell(mm.N_pt,1);
sim_out_hf = cell(mm.N_pt,1);
sim_out = cell(mm.N_pt,1);

seeds = randi(1e6,size(mm.Phi,1),2);
parfor pt_ndx = 1:1:mm.N_pt
%for pt_ndx = 1:1:mm.N_pt % use this for loop for debugging only
    
    rng(seeds(mm.pt_type(pt_ndx,1),1),'twister');
    seed_crand(seeds(mm.pt_type(pt_ndx,1),2));

    if mm.sim_firm_num_by_prod_succ_type(pt_ndx)>0

        sim_out_f{pt_ndx} = matchdat_gen_f(pt_ndx,macro_state_f, mm, policy);

        sim_out_h{pt_ndx} = matchdat_gen_h(pt_ndx,macro_state_h, mm, policy);

        sim_out_hf{pt_ndx} = splice_hf(sim_out_h{pt_ndx},sim_out_f{pt_ndx},policy,mm);

        sim_out{pt_ndx} = cell2struct([struct2cell(sim_out_h{pt_ndx});struct2cell(sim_out_f{pt_ndx});struct2cell(sim_out_hf{pt_ndx})],...
            [fieldnames(sim_out_h{pt_ndx});fieldnames(sim_out_f{pt_ndx});fieldnames(sim_out_hf{pt_ndx})]);

    end

end

sim_cum = aggregateSimulatedData(sim_out,mm);

simMoms = calculateSimulatedMoments(sim_cum);
















