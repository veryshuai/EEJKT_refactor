function shouldMatchMoments(new_moment1,new_moment2,write_output,filepath)

    if write_output == "overwrite"
        old_moment1 = new_moment1;
        old_moment2 = new_moment2;
        save(filepath,'old_moment1','old_moment2')
    elseif write_output == "test"
        load(filepath);
        match_moment1 = min(new_moment1 == old_moment1,[],'all');
        match_moment2 = min(new_moment2 == old_moment2,[],'all');
        assert(min(match_moment1,match_moment2)==1,"shouldMatchMoments: Benchmark moments do not match test moments");
    else
        error("shouldMatchOutput: Third argument must be either overwrite or test.")
    end

end