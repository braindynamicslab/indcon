function y = f_nan2ave(x)
% - replace missing values (NaN's) with average value
%
% USAGE: y = f_nan2ave(x)
% 
% x = input column vector
% y = output vector

% -----Author:-----
% by David L. Jones, Feb-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input:-----
if size(x,2)>1
   error ('X should be a column vector!')
end

y           = x; % make a copy
y(isnan(x)) = nanmean(x); 