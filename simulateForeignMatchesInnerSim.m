function iter_in = simulateForeignMatchesInnerSim(iter_in,mm,policy)

    [iter_in, drop_Zcut] = simulateForeignMatchesInnerSimUpdateClientCount(iter_in, mm, policy);
    iter_in = simulateForeignMatchesInnerSimUpdZHotel(mm, iter_in, policy);
    iter_in = simulateForeignMatchesInnerSimKickDormant(iter_in, mm);
    [iter_in, age] = simulateForeignMatchesInnerSimFirmAge(iter_in, mm);
    [iter_in, mat_tran] = simulateForeignMatchesInnerSimMatchLevelData(iter_in, mm, age, drop_Zcut);
    iter_in.N_match = iter_in.N_match + size(mat_tran,1); % cumulate match count within current year
   
end