function [iter_out,transH] = simulateHomeMatches(pt_ndx,macro_state_h,mm,policy,iter_out)  

[iterH_in, iter_out] = simulateHomeInnerInitialize(mm, pt_ndx, macro_state_h, iter_out);

[iter_out,transH] = simulateHomeMatchesInnerSim(iter_out, mm, iterH_in, pt_ndx, policy);

end