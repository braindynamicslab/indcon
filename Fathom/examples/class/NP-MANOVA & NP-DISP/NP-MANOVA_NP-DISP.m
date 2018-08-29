%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           NP-MANOVA & NP-DISP                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% by David L. Jones, 10-Nov-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.
 
% The file 'tikus.mat' contains data on coral assemblages from Tikus
% Island, Indonesia taken along 10 replicate 30 m line transects during
% each of 6 years (Warwick et al., 1990). The 1982 El Niño triggered a
% disturbance event which resulted in widespread coral bleaching in this
% region. As such, the 1981, 1983, and 1985 samples were collected before,
% during, and after the El Niño disturbance event, repectively.
% 
% Warwick, R. M., K. R. Clarke, and Suharsono. 1990. A statistical analysis
% of coral community responses to the 1982-1983 El Niño in the Thousand
% Islands, Indonesia. Coral Reefs 8:171-179.
% 
% -----Variables:----
% coral     = percent cover of 75 species in each transect
% coral_txt = cell array of corresponding column labels
% yr        = year (1 = 1981, 3 = 1983, 4 = 1984, 5 = 1985, 7 = 1987, 8 = 1988)
% yr_txt    = cell array of corresponding row labels

% 1) Examine the effect of this disturbance event on the coral assemblages
% by determining whether there were differences in the variation in the
% abundance and composition of corals before (1981), during (1983), and
% after (1985) the 1982 El Niño Event. Interpret your results in terms of a
% null hypothesis.

% ANSWER:
% 
% Load data:
load tikus.mat

% Apply a mild data transformation to down-weight the influence of the most
% abundant species:
coral_2 = f_normal(coral,'2');

% Create a Bray-Curtis dissimilarity matrix from the transformed data:
disBC = f_dis(coral_2,'bc'); % Bray-Curtis

% 1981 = before 
% 1983 = during
% 1985 = after disturbance event
% Create index to just these years:
idx = (yr==1 | yr==3 | yr==5);

% Quantify multivariate dispersion for each of these years:
md = f_npDisp(disBC(idx,idx),yr(idx),1,1);
% ==================================================
%       Quantify Multivariate Dispersions:
% ==================================================
% 
% # Pos Eigenvalues = 20
% # Neg Eigenvalues = 9
% 
% Average distance to spatial median:
% Group 1 = 0.4344
% Group 3 = 0.6407
% Group 5 = 0.3787
 
% Test for significant differences in dispersion among years:
f_npManova(f_dis(md.z,'euc'),yr(idx),1000,1);
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'      'df'    'SS'         'MS'          'F'         'p'    
%     'factor 1'    [ 2]    [0.38096]    [ 0.19048]    [12.781]    [0.001]
%     'residual'    [27]    [0.40239]    [0.014903]    [   NaN]    [  NaN]
%     'total'       [29]    [0.78335]    [     NaN]    [   NaN]    [  NaN]
% 
% 
%       # iterations =        1000 
% --------------------------------------------------
% 
% -> Reject the null hypothesis of no difference in multivariate
% dispersions among years (at alpha = 0.05).

% Determine which years differed?
f_npManovaPW(f_dis(md.z,'euc'),yr(idx),1000,1);
% 
% Results of pair-wise comparisons between each factor level:
% ----------------------------------------------------------
% 1 vs. 2: t = 4.3706 p = 0.003
% 1 vs. 3: t = 0.9615 p = 0.339
% 2 vs. 3: t = 4.5156 p = 0.002
% 
% -> The dispersion of 1983 differs from 1981 and 1985, while 1981/1985 are
% not significantly different.

% INTERPRETATION: There was MORE variation in the abundance and composition
% (beta diversity) of coral species during the disturbance event (1983)
% compared to that observed either before (1981) or after (1985) the event.
% This suggests that some type of environmental stress imposed by the
% disturbance event, in the form of coral bleaching, was manifested in
% these assemblages of corals as increased multivariate dispersion.

% 2) How could you examine the same question, but focus entirely on
% differences in the variation in species composition?

% ANSWER: Use, say, the Jaccard index, which only considers species
% presence/absence (NOT abundance):
disJAC = f_dis(coral_2,'jac');

% Or, just convert your data to presence/absence and use the Bray-Curtis
% metric:
disBIN  = f_dis(f_normal(coral_2,'01'),'bc');
disBIN2 = f_dis((coral_2>0),'bc');

% 3) Were there differences in the number of coral species present before,
% during, and after the disturbance event?

% Here's one of many ways to accomplish this:
noSppYr_1 = sum(sum(coral(yr==1,:))>0)
noSppYr_3 = sum(sum(coral(yr==3,:))>0)
noSppYr_5 = sum(sum(coral(yr==5,:))>0)

% Yes, there were 54, 21, and 33 species present before, during, and after
% the disturbance event.

% 4) Tuomisto et al. (2003) examined the distribution of plant species from
% 163 sites across 4 regions of western Amazonian forests. Their study
% focused on 286 species of ferns (Pteridophytes) and 265 species of shrubs
% and small trees (Melastomataceae) inventoried using 500 m by 5 m line
% transects. Tuomisto, H., K. Ruokolainen, and M. Yli-Halla. 2003.
% Dispersal, environment, and floristic variation of western Amazonian
% forests. Science  299:241-244.
% 
% Their data is contained in file 'amazon.mat' as follows:
% 
% site = structure of SITE data with the following fields:
%  .name = site name
%  .txt  = site code
%  .reg  = site region (Colombia, Ecuador, N_Peru, S_Peru)
%  .lat  = latitude  (degrees S)
%  .lon  = longitude (degrees W)
%  .inu  = terrain is seasonally inundated (1 = yes, 0 = no)
%  .swa  = terrain is a swamp forest       (1 = yes, 0 = no)
%  .ter  = terrain is tierra firme         (1 = yes, 0 = no)
%  .san  = terrain is white sand           (1 = yes, 0 = no)
%  .sea  = climate has a dry season        (1 = yes, 0 = no)
% 
% bio = structure of corresponding BIOTIC data with the following fields:
%  .txt  = site code
%  .Pter = Jaccard index based dissimilarity matrix for the fern taxa
%  .Mel  = Jaccard index based dissimilarity matrix for the schrub/tree taxa
% 
% soil = structure of corresponding ENVIRONMENTAL data with the following fields:
%  .txt  = site code
%  .dis = Euclidean distance matrix created from the following 8 soil
%          variables: pH, loss on ignition, percentage of clay and silt,
%          and natural log transformed concentration of Ca, K, Mg, Na, and
%          Al).
% 
% Re-examine these data from Class 10 and determine if there were significant
% differences in the composition of fern species among the 4 regions of western
% Amazonian forests. Interpret your results.
% 
% HINT: You might need to create a vector of integers based on site.reg
% that specifies the region each sample site originated. There are MANY
% ways to accomplish this, including using the 'ismember' and/or
% 'f_dummy2cat' commands.

% ANSWER:
% Load the data:
load amazon.mat

% Create grouping vector indicating the region each site came from
% 1 = Colombia
% 2 = Ecuador
% 3 = N_Peru
% 4 = S_Peru
% 
grp = sum([ismember(site.reg,'Colombia') ...
ismember(site.reg,'Ecuador')*2 ...
ismember(site.reg,'N_Peru')*3 ...
ismember(site.reg,'S_Peru')*4],2);

% Alternative method:
grp2 = f_dummy2cat([ismember(site.reg,'Colombia') ...
ismember(site.reg,'Ecuador') ...
ismember(site.reg,'N_Peru') ...
ismember(site.reg,'S_Peru')]);

% Another method (logical indices):
grp3 = zeros(size(bio.Pter,1),1); % initialize
grp3(ismember(site.reg,'Colombia')) = 1;
grp3(ismember(site.reg,'Ecuador'))  = 2;
grp3(ismember(site.reg,'N_Peru'))   = 3;
grp3(ismember(site.reg,'S_Peru'))   = 4;

% Another method (find command):
grp4 = zeros(size(bio.Pter,1),1); % initialize
idx_1 = find(ismember(site.reg,'Colombia')==1);
idx_2 = find(ismember(site.reg,'Ecuador') ==1);
idx_3 = find(ismember(site.reg,'N_Peru') ==1);
idx_4 = find(ismember(site.reg,'S_Peru') ==1);
grp4(idx_1) = 1;
grp4(idx_2) = 2;
grp4(idx_3) = 3;
grp4(idx_4) = 4;

% Generic method (FOR loop):
txt  = site.reg;          % assign cell array of strings to generic variable txt
uGrp = unique(txt);       % get sorted list of unique groups
nGrp = numel(uGrp);       % get # of unique groups
g    = nan(numel(txt),1); % initialize
for i=1:nGrp
   g(ismember(txt,uGrp(i)),1) = i;
end

% Test the null hypothesis that there is no difference in the composition
% of species among the 4 regions:
f_npManova(bio.Pter,grp,1000,1);
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'      'df'     'SS'        'MS'        'F'         'p'    
%     'factor 1'    [  3]    [11.423]    [3.8078]    [13.767]    [0.001]
%     'residual'    [159]    [43.979]    [0.2766]    [   NaN]    [  NaN]
%     'total'       [162]    [55.403]    [   NaN]    [   NaN]    [  NaN]
% 
% 
%       # iterations =        1000 
% --------------------------------------------------
% 
% -> Reject the null hypothesis of no significant difference in species
% composition among regions (at alpha = 0.05).

% Check that we met the assumptions of MANOVA regarding Homogeneity of
% Multivariate Dispersions:
% 
% Quantify multivariate dispersion for each region:
md_Pter = f_npDisp(bio.Pter,grp,1,1);
% 
% ==================================================
%       Homogeneity of Multivariate Dispersions:
% ==================================================
% 
% # Pos Eigenvalues = 116
% # Neg Eigenvalues = 46
% 
% Average distance to spatial median:
% Group 3 = 0.5872
% Group 2 = 0.4273
% Group 1 = 0.5192
% Group 4 = 0.5552

% Test for significant differences in dispersion among regions:
f_npManova(f_dis(md_Pter.z,'euc'),grp,1000,1);
% 
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'      'df'     'SS'         'MS'          'F'         'p'    
%     'factor 1'    [  3]    [0.73106]    [ 0.24369]    [24.163]    [0.001]
%     'residual'    [159]    [ 1.6035]    [0.010085]    [   NaN]    [  NaN]
%     'total'       [162]    [ 2.3346]    [     NaN]    [   NaN]    [  NaN]
% 
% 
%       # iterations =        1000 
% --------------------------------------------------
% 
% -> Reject the null hypothesis of no difference in multivariate
% dispersions among regions (at alpha = 0.05).

% INTERPRETATION: Results from the MANOVA indicate there were significant
% differences in the composition of species among regions, however, since
% we failed to meet the assumpton of homogeneous within-group variances
% required by MANOVA, we may have suffered a TYPE I error (false positive).
% This is because the NP-DISP analysis indicated there were significant
% differences in the variation of species composition among regions.
% However, simply knowing there are differences in beta diversity among regions
% may itself provide valuable insights into the ecology of the fern species
% in this area.
