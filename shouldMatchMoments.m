function shouldMatchMoments(new_moments,new_fit,write_output,filepath)

    if write_output == "overwrite"
        old_moments = new_moments;
        old_fit = new_fit;
        save(filepath,'old_moments','old_fit')
    elseif write_output == "test"
        load(filepath);
        match_all_moments = min(new_moments == old_moments,[],'all');
        match_fit = min(new_fit == old_fit,[],'all');
        assert(min(match_all_moments,match_fit)==1,"shouldMatchMoments: Benchmark moments do not match test moments");
    else
        error("shouldMatchOutput: Third argument must be either overwrite or test.")
    end

end