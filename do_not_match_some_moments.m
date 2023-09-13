function [Data, error, Model] = do_not_match_some_moments(Data, error, Model)
% To exclude log-log regression, allow next block

    Data(11:13) = 0;
%   Model(11) = 0;
    error(11:13) = 0;
    

% To exclude avg. log exports from DANE-based X-D regression, allow next block

    Data(15) = 0;
%   Model(15) = 0;
    error(15) = 0;
    
% To exclude variance of DANE-based X-D regression, allow next block

   Data(17) = 0;
%   Model(17) = 0;
   error(17) = 0;
    
% To exclude avg. log domestic sales from DANE-based AR1 regression, allow next block

    Data(18) = 0;
%   Model(15) = 0;
    error(18) = 0;
        
% To exclude variance of error in domestic AR1, allow next block  
    
  Data(20) = 0;
%   Model(20) = 0;
  error(20) = 0;
%
end