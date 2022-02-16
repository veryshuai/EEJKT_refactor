% This script switches all home specific variables to foreign to calculate full information foreign search behavior

mm.scale_h   = scale_f;    % Domestic profit function scale parameter

[Q0_h,Q0_h_d,st_h] = makebigq(Q_p,Q_f,mm.phi_size,mm.x_size,Phi,X_f);

mm.actual_h     = actual_f; %actual indexes of home macro shocks for available years
mm.X_h          = X_f;      %home macro shocks
mm.st_h         = st_f;     %home state (for use with intensity matrix)

mm.L_h          = L_f;      %arrival rate for jumps in home macro shock
mm.D_h          = D_f;      %size of jump in home macro shock
mm.Q_h          = Q_f;      %intensity matrix for home macro shock

mm.Q0_h         = Q0_f;     %intensity matrix for home state
mm.Q0_h_d       = Q0_f_d;   %home with zeros on diagonal
