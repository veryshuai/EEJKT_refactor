
function [lambda_f,lambda_h,pi_h,pi_f,c_val_h,c_val_f,my_flag] = solve(mm)
% takes current parameters and calculates value and policy functions

    % Read in parameters
    scale_f     = mm.scale_f;       % Export profit function scale parameter
    scale_h     = mm.scale_h;       % Domestic profit function scale parameter
    net_size    = mm.net_size;      % max number of network effects
    n_size      = mm.n_size;        % Maximum number of informative signals per firm
    phi_size    = 2*mm.phi_size+1;  % number of different discretized seller productivities
    x_size      = 2*mm.x_size+1;    % Number of different discretized macro shocks; same for home and foreign
    z_size      = 2*mm.z_size+1;    % Number of discretized buyer states
    dim0        = mm.dim0;          % Number of possible theta0 values (common to both markets)
    dim1        = mm.dim1;          % Number of possible theta1 values (specific to home market);
    th_g        = mm.theta0;        % common theta grid
    th_h        = mm.theta1;        % home theta grid
    th_f        = mm.theta2;        % foreign theta grid
    af          = mm.af;            % beta function parameter for foreign theta draws
    bf          = mm.bf;            % beta function parameter for foreign theta draws
    alp         = mm.alpha;         % weight of "common" theta in determining match probabilities (set to 0)

    Z           = mm.Z;        %buyer productivities 
    st_h        = mm.st_h;     %home state (for use with intensity matrix)
    st_f        = mm.st_f;     %foreign state (for use with intensity matrix)
    Q_z         = mm.Q_z;      %intensity matrix for buyer productivities 
    Q_z_d       = mm.Q_z_d;    %with zeros on the diagonal
    erg_pz      = mm.erg_pz;   %ergodic distribution of buyer productivities 
    Q0_h        = mm.Q0_h;     %intensity matrix for home state
    Q0_f        = mm.Q0_f;     %intensity matrix for foreign state
    Q0_h_d      = mm.Q0_h_d;   %home with zeros on diagonal
    Q0_f_d      = mm.Q0_f_d;   %foreign with zeros on diagonal
    Phi         = mm.Phi;      %vector of seller productivites
    
    
    F_f         = mm.F_f;      %fixed costs of maintaining a relationship foreign
    F_h         = mm.F_h;      %fixed costs of maintaining a relationship home


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%=====Solution=====%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %% CALCULATE PROFIT STREAM OF SUCCESSFUL MATCH
    
    [pi_f,~,c_val_f] = makepie(scale_f,st_f,Z,Q0_f,Q0_f_d,Q_z,Q_z_d,erg_pz,F_f,mm);
    [pi_h,~,c_val_h] = makepie(scale_h,st_h,Z,Q0_h,Q0_h_d,Q_z,Q_z_d,erg_pz,F_h,mm);
    
    %% CALCULATE POSTERIOR MATCH PROBABILITIES
    
    % foreign
    a_f = makepost(n_size,th_g,af,bf,alp); %[trials,successes,th_g]

    %home
    a_h = ones(size(th_g,2),size(th_h,2));
    for m = 1:size(th_g,2)
        for j = 1:size(th_h,2)
            a_h(m,j) = th_g(m) * alp + th_h(j)*(1-alp);
        end
    end
    
    %% SOLVE FOR SEARCH EFFORTS AT HOME AND AT FOREIGN
    
    [~,lh,flag_h] = val_loop_h(Q0_h,Q0_h_d,a_h,pi_h,mm); %lh = (exogenous shock vector,theta_g,theta_h,client no)
    [~,lf,flag_f] = val_loop_f(Q0_f,Q0_f_d,a_f,pi_f,mm); %lf = (exogenous shock vector,successes,trials,theta_g,client no) 
    my_flag = max(flag_h,flag_f); %flag = 1 means kill, flag = 0 not kill

    %For use in "no learning" counterfactual
    %punishment = 0;
    %lf = zeros(size(Q0_f,1),n_size+1,n_size+1,dim0,net_size+1);
    
    %% RESHAPE 
    % For now, put matrix result into arrays to match old estimation routine format 
    if my_flag == 0
    
        % reshape policy functions
        lh_new = reshape(lh,x_size,phi_size,dim0,dim1,net_size+1); %splitting up the exogenous shock vector to macro and productivity shocks
        lf_new = reshape(lf,x_size,phi_size,n_size+1,n_size+1,dim0,net_size+1); %same
        
        % read policy function into cell arrays 
        lambda_h = cell(dim0,dim1,net_size); %(theta_g,theta_h,net effect)
        lambda_f = cell(n_size+1,n_size+1,dim0,net_size+1); %(successes,trials,theta_g,net effect)
        for j = 1:n_size+1
            for k = 1:n_size+1
                for m = 1:dim0
                    for n = 1:net_size+1
                    lambda_f{k,j,m,n} = lf_new(:,:,j,k,m,n)';
                    end
                end
            end
        end
        for j = 1:dim0
            for k = 1:dim1
                for l = 1:net_size+1
                lambda_h{j,k,l} = lh_new(:,:,j,k,l)'; %transpose the matrices so the productivity phi is the rows and macro shock x is the columns (as in old program)
                end
            end
        end
        
        %also reshape profits
        pi_f = reshape(pi_f,x_size,phi_size)';
        pi_h = reshape(pi_h,x_size,phi_size)';
        
        %reshape value functions
        c_val_f_cell = cell(z_size,1);
        c_val_h_cell = cell(z_size,1);
        for k = 1:z_size
           c_val_f_cell{k,1} = reshape(c_val_f(:,k),x_size,phi_size)';
           c_val_h_cell{k,1} = reshape(c_val_h(:,k),x_size,phi_size)';
        end
        c_val_f = c_val_f_cell;
        c_val_h = c_val_h_cell;

    else %IN CASE OF ERROR
        lambda_f = 1;lambda_h = 1;pi_h = 1;pi_f = 1;c_val_h = 1;c_val_f = 1;
    end

end
