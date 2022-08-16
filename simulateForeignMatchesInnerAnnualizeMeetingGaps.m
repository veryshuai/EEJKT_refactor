function [iter_in, iter_out] = simulateForeignMatchesInnerAnnualizeMeetingGaps(iter_in, mm, iter_out)
if iter_in.t>3*mm.pd_per_yr
    [time_gap,iter_in.mkt_exit] = time_gaps(iter_in,mm);
    %           time_gap: (1) firm_ID, (2) periods into interval, (3) time gap,
    %           (4) # new meetings, (5) iter_in.t, (6) cum. meetings, (7) cum succeses

    iter_out.time_gaps = [iter_out.time_gaps;time_gap];
end
end