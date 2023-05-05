function [Data, error, Model] = do_not_match_some_moments(Data, error, Model)
  

% To exclude avg. log exports from DANE-based X-D regression, allow next block

    Data(12) = 0;
%   Model(12) = 0;
    error(12) = 0;
    
% To exclude variance of DANE-based X-D regression, allow next block

   Data(14) = 0;
%   Model(14) = 0;
   error(14) = 0;
    
% To exclude avg. log domestic sales from DANE-based AR1 regression, allow next block

    Data(15) = 0;
%   Model(15) = 0;
    error(15) = 0;
        
% To exclude variance of error in domestic AR1, allow next block  
    
  Data(17) = 0;
%   Model(17) = 0;
  error(17) = 0;
%
end