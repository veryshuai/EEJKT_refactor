function policy = reshapePolicy(policy, mm)
    
    policy.value_f = permute(reshape(policy.value_f,[size(mm.X_f,1),size(mm.Phi,1),mm.n_size+1,mm.n_size+1,1,mm.net_size+1]),[4,3,5,6,2,1]);
    policy.lambda_f = permute(reshape(policy.lambda_f,[size(mm.X_f,1),size(mm.Phi,1),mm.n_size+1,mm.n_size+1,1,mm.net_size+1]),[4,3,5,6,2,1]);
    policy.c_val_f = permute(reshape(policy.c_val_f,[size(mm.X_f,1),size(mm.Phi,1),size(mm.Z,1)]),[3,2,1]);
    policy.value_h = permute(reshape(policy.value_h,[size(mm.X_h,1),size(mm.Phi,1),1,mm.dim1,mm.net_size+1]),[3,4,5,2,1]);
    policy.lambda_h = permute(reshape(policy.lambda_h,[size(mm.X_h,1),size(mm.Phi,1),1,mm.dim1,mm.net_size+1]),[3,4,5,2,1]);
    policy.c_val_h = permute(reshape(policy.c_val_h,[size(mm.X_h,1),size(mm.Phi,1),size(mm.Z,1)]),[3,2,1]);

end