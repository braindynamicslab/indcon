%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Cluster Analysis:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% by David L. Jones, 14-Apr-2010
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Load the file 'ekofisk.mat', which consists of North Sea Ekofisk oil field
% data from Gray et al., (1988) and has the following variables:
% 
% biotic     = abundance of 174 species of soft-bottom benthic macrofauna from 39 sites
% biotic_txt = cell array of corresponding column labels
% env        = values of 6 environmental variables from the same 39 sites
% env_txt    = cell array of corresponding column labels:
%       'Dist': distance (km) from center of drilling activity
%        'THC': total hydrocarbon concentration
%      'Redox': redox potential (~ organic matter)
%       '%Mud': percent of sediment that was mud
%    Phi_mean': mean sediment particle size
%         'Ba': sediment Barium concentration (nontoxic, tracer for drilling muds)
%         'Sr': sediment Strontium concentration
%         'Cu': sediment Copper concentration
%         'Pb': sediment Lead concentration
%         'Ni': sediment Nickel concentration
% sites_txt  = cell array of site labels
% 

% Perform an UPGMA-based cluster analysis of sites based on species abundance
% and composition. Examine the resulting dendrogram and identify aggregations of
% sites that appear to form natural groupings. Create a new variable that
% specifies membership of sites within these new groups.
% 
% HINT: use the 'f_cluster' function to perform the cluster analysis.

% ANSWER:
load ekofisk.mat; % load data

biotic_2 = f_normal(biotic,'2');         % square-root transform
disBC    = f_dis(biotic_2,'bc');         % Bray-Curtis
idx      = f_cluster(disBC,1,sites_txt); % UPGMA cluster analysis



% Examine the resulting dendrogram and identify aggregations of sites that
% appear to form natural groupings. %  Groups can be defined depending on
% how you arbirtarily decide on a cutoff. Show sites ordered by dendrogram:

sites_txt(idx)

% ans = 
% 
%     'S37'   GROUP 1
%     'S30'
%     'S1'    GROUP 2
%     'S5'
%     'S10'
%     'S11'
%     'S16'
%     'S17'
%     'S20'
%     'S6'
%     'S7'
%     'S21'  
%     'S25'   GROUP 3
%     'S35'
%     'S38'
%     'S3'
%     'S31'
%     'S36'
%     'S29'
%     'S12'
%     'S14'
%     'S15'
%     'S39'
%     'S13'
%     'S18'
%     'S8'
%     'S19'
%     'S2'
%     'S32'
%     'S9'
%     'S4'
%     'S24'
%     'S33'
%     'S28'
%     'S34'
%     'S22'
%     'S23'
%     'S26'
%     'S27'
%


% Display the concentration of THC in the same order as plotted on the dendrogram.

env(idx,2)

% ans =
%         1896
%          232
%            8
%            8
%           11 
%           13
%            9
%            8
%           13
%           12
%           10
%           14
%           57
%           56
%          127
%           17
%          251
%          443
%           61
%           13
%           11
%           13
%           25
%           14
%           10
%            8
%           17
%           17
%          251
%            8
%           18
%           13
%           38
%           15
%           19
%           11
%            9
%           39
%           32

