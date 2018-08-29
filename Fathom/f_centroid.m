function [centroid,SE] = f_centroid(x,grps,method)
% - returns coordinates of the centroid of X, optionally partitioned into groups
%
% USAGE: [centroid,SE] = f_centroid(x,grps,method);
%
% x    = n-dimensional coordinates (rows = observations,
%        (cols = dimensions
% grps = optional vector of integers specifying group membership
%        e.g., grps = [1 1 2 2 3 3 3];
% method = use mean (= 1, default) or median (= 2);
% 
% centroid = spatial mean or median
% SE       = corresponding standard error

% -----Author:-----
% by David L. Jones, Apr-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% 08-Dec-2002: made grouping vector optional, error checking
% Nov-2010: replaced 'unique' with 'f_unique'
% Dec-2012: added support for SE

% -----Check input & set defaults:-----
if (nargin < 3), method = 1; end; % default use mean

if (size(x,1)<2)
   error('You need at least 2 points to compute the centroid!');
end
% -------------------------------------

if (nargin < 2)
   switch method
      case 1
         centroid = mean(x);
      case 2
         centroid = median(x);
      otherwise
         error('Invalid METHOD!')
   end
   
else
   grps   = grps(:);       % make sure it's a column vector
   uGrps  = f_unique(grps);  % f_unique grps, unsorted
   noGrps = length(uGrps); % number of unique groups
   
   centroid(noGrps,size(x,2)) = 0; % preallocate results array
   SE(noGrps,size(x,2))       = 0; % 
   
   for i=1:noGrps
      subsetRows = find(grps==uGrps(i)); % get indices of rows to extract
      subsetX    = x(subsetRows,:);      % extract each group separately
      
      switch method
         case 1
            centroid(i,:) = mean(subsetX,1);         % compute centroid for this group
            SE(i,:)       = f_stdErr(subsetX);       % standard error of the mean
         case 2
            centroid(i,:) = median(subsetX,1);       % compute centroid for this group
            SE(i,:)       = f_stdErr(subsetX)*1.253; % standard error of the median
         otherwise
            error('Invalid METHOD!');
      end
   end
end
