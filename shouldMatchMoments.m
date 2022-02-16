function shouldMatchMoments(new_moments,new_fit,write_output)

    if write_output == "overwrite"
        old_moments = new_moments;
        old_fit = new_fit;
        save results/shouldMatchMomentsData old_moments old_fit
    elseif write_output == "test"
        load("results/shouldMatchMomentsData");
        match_all_moments = min(min(new_moments == old_moments));
        match_fit = new_fit == old_fit;
        assert(min(match_all_moments,match_fit)==1,"shouldMatchMoments: Benchmark moments do not match test moments");
    else
        error("shouldMatchOutput: Second argument is either overwrite or test.")
    end

end