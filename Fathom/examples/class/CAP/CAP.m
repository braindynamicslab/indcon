%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Canonical Analysis of Principal Coordinates               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% by David L. Jones, 17-Nov-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% The file 'ekofisk.mat' consists of North Sea Ekofisk oil field data from
% Gray et al., (1988) and has the following variables: 
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

% 1) You are interested in determining whether the community composition and
% abundance of the macrobenthic fauna in the Ekofisk oil field is distinct
% among geographic zones. Geographic zones for each site are defined
% according to their distance from the center of drilling activity as
% follows:
% 
% zone 1: 0 to 250 m
% zone 2: greater than 250m, but not exceeding 1 km
% zone 3: greater than 1 km, but not exceeding 3.5 km
% zone 4: greater than 3.5 km
% 
% Perform a formal statistical analysis using MANOVA to determine whether
% differences in macrobenthos species composition and abundance exist among
% the Ekofisk oil field zones. Interpret your results.

% -----ANSWER:-----
% Load data
load ekofisk.mat;

% Apply a mild data transformation to down-weight the influence of the most
% abundant species:
biotic_2 = f_normal(biotic,'2'); % square-root transform

% Create a Bray-Curtis dissimilarity matrix from the transformed data:
disBC = f_dis(biotic_2,'bc');

% Create grouping vector indicating the zone each site came from:
% zone 1: 0 to 250 m
% zone 2: greater than 250m, but not exceeding 1 km
% zone 3: greater than 1 km, but not exceeding 3.5 km
% zone 4: greater than 3.5 km

% Create grouping variable using logical indices:
geo = zeros(size(env,1),1) + 4; % initialize
geo(env(:,1)<=3.5)  = 3;
geo(env(:,1)<=1)    = 2;
geo(env(:,1)<=0.25) = 1;

% Alternatively, create grouping variable using the find command:
idx1 = find(env(:,1)<=0.250);
idx2 = find(env(:,1)>0.250 & env(:,1)<=1);
idx3 = find(env(:,1)>1 & env(:,1)<=3.5);
idx4 = find(env(:,1)>3.5);
geo2 = repmat(NaN,size(env(:,1))); % initialize
geo2(idx1) = 1;
geo2(idx2) = 2;
geo2(idx3) = 3;
geo2(idx4) = 4;

% Alternatively, use a FOR loop:
geo3 = zeros(size(env,1),1); % initialize
for i = 1:size(env,1)
   if env(i,1)<=0.250, geo3(i) = 1;
   elseif env(i,1)>0.250 && env(i,1)<=1, geo3(i) = 2;
   elseif env(i,1)>1 && env(i,1)<=3.5, geo3(i) = 3;
   elseif env(i,1)>3.5, geo3(i) = 4;
   end
end

% Test the null hypothesis that there is no significant difference in the
% composition and abundance of macrobenthic fauna species among the 4 zones:
resultBC = f_npManova(disBC,geo,1000,1);
% 
% Permuting the data 999 times...
% 
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'      'df'    'SS'        'MS'          'F'         'p'    
%     'factor 1'    [ 3]    [1.2536]    [ 0.41785]    [6.3634]    [0.001]
%     'residual'    [35]    [2.2983]    [0.065665]    [   NaN]    [  NaN]
%     'total'       [38]    [3.5518]    [     NaN]    [   NaN]    [  NaN]
% 
%       # iterations =        1000 
% --------------------------------------------------
% -> Reject the null hypothesis of no significant difference in species
% composition and abundance among zones (at alpha = 0.05).

% Determine which zones differed:
f_npManovaPW(disBC,geo,1000,1);
% 
% Results of pair-wise comparisons between each factor level:
% ----------------------------------------------------------
% 1 vs. 2: t = 2.2673 p = 0.001
% 1 vs. 3: t = 2.8716 p = 0.001
% 1 vs. 4: t = 3.1752 p = 0.001
% 2 vs. 3: t = 1.4337 p = 0.001
% 2 vs. 4: t = 2.7311 p = 0.001
% 3 vs. 4: t = 2.1497 p = 0.001
% 
% -> Reject the null hypothesis of no significant difference in species
% composition and abundance among all possible pairs of zones (at alpha = 0.05).

% Check that we met the assumptions of MANOVA regarding Homogeneity of
% Multivariate Dispersions:
% 
% Quantify multivariate dispersion for each zone:
md = f_npDisp(disBC,geo,1,1);
% ==================================================
%       Quantify Multivariate Dispersions:
% ==================================================
% 
% # Pos Eigenvalues = 33
% # Neg Eigenvalues = 5
% 
% Average distance to spatial median:
% Group 1 = 0.2779
% Group 2 = 0.2180
% Group 3 = 0.2260
% Group 4 = 0.2596
% --------------------------------------------------

% Test for significant differences in dispersion among zones:
f_npManova(f_dis(md.z,'euc'),geo,1000,1);
% 
% Permuting the data 999 times...
% 
% ==================================================
%     Nonparametric (Permutation-based) MANOVA:
% --------------------------------------------------
%     'Source'      'df'    'SS'          'MS'           'F'         'p'    
%     'factor 1'    [ 3]    [0.019944]    [0.0066479]    [2.3781]    [0.093]
%     'residual'    [35]    [ 0.09784]    [0.0027954]    [   NaN]    [  NaN]
%     'total'       [38]    [ 0.11778]    [      NaN]    [   NaN]    [  NaN]
% 
%       # iterations =        1000 
% --------------------------------------------------
% -> Accept the null hypothesis of no significant difference in multivariate
% dispersions among zones (at alpha = 0.05).

% INTERPRETATION: Results from the MANOVA indicate there were significant
% differences in macrobenthos species composition and abundance among the
% Ekofisk oil field zones at an alpha=0.05. The post hoc multiple
% comparison tests revealed there were significant differences among ALL 4
% zones, which indicates there were distinct benthic macrofaunal
% communities in all 4 zones. 

% 2) Follow up your MANOVA with a CAP-based Canonical Discriminant
% Analysis. Please interpret ALL of your results. Provide a multivariate
% visualization depicting the differences among zones according to species
% abundance and composition. Provide an estimate of the discrimant
% functions' ability to accurately classify a benthic macrofauna sample of
% unknown origin to the correct zone.

% Find optimal value of m:
f_capOptimal(biotic_2,'bc',geo,0,1);
% 
% ==================================================
%  Diagnostics for CAP - db-CDA:
% --------------------------------------------------
% m:   propG:   RSS: d_1^2: d_2^2: d_3^2:  Correct:
% 5  57.9581  2.7302 0.8822 0.6490 0.0567   71.79% 
% 6  62.0871  2.7406 0.8822 0.6520 0.0594   71.79% 
% 7  65.9613  2.8544 0.8840 0.7020 0.0803   66.67% 
% 8  69.1493  2.8887 0.8845 0.7021 0.0979   64.10% 
% 9  72.1367  2.8832 0.8894 0.7057 0.1752   69.23% 
% 10  74.7182  2.9072 0.8903 0.7418 0.2138   66.67% 
% 11  77.2580  2.9391 0.8914 0.7452 0.3359   71.79% 
% 12  79.7368  2.9683 0.9050 0.7517 0.4058   76.92% <- optimal m
% 13  81.9887  3.0097 0.9083 0.7558 0.4820   74.36% 
% 14  84.1416  3.2480 0.9093 0.7611 0.4882   71.79% 
% 15  86.0304  3.3015 0.9181 0.7626 0.4882   71.79% 
% 16  87.7650  3.3439 0.9213 0.7696 0.4892   74.36% 
% 17  89.3458  3.3594 0.9218 0.7980 0.4893   74.36% 
% 18  90.7972  3.3937 0.9221 0.7983 0.5143   74.36% 
% 19  92.1669  3.5405 0.9241 0.8001 0.5169   71.79% 
% 20  93.3866  3.5580 0.9342 0.8293 0.5186   71.79% 
% 21  94.4794  3.5576 0.9369 0.8351 0.5525   71.79% 
% 22  95.5353  3.6908 0.9369 0.8351 0.5950   64.10% 
% 23  96.4710  3.6933 0.9392 0.8404 0.6761   64.10% 
% 24  97.3932  4.1314 0.9514 0.8415 0.6952   71.79% 
% 25  98.0996  4.2977 0.9551 0.8507 0.7176   64.10% 
% 26  98.7213  4.3605 0.9685 0.8738 0.7178   69.23% 
% 27  99.2642  4.3477 0.9723 0.8969 0.7334   66.67% 
% 28  99.7344  4.6580 0.9734 0.9069 0.7378   64.10% 
% --------------------------------------------------
% m       = # PCoA axes retained
% propG   = proportion of yDis explained by m PCoA axes
% RSS     = leave-one-out residual sums-of-squares
% d_1^2   = squared canonical correlation for axis 1
% Correct = leave-one-out classification success
% 
% Central Tendency = centroid 
% Optimal value of m for CDA may be 12
% ==================================================
% 
% -> the optimal value of m is 12 because using the first 12 PCoA axes
% accounts for more than 60% of the total variation in the original
% dissimilarity matrix (i.e., it accounts for nearly 80%) and it provides
% the highest classification success rate (i.e., 77% of all observations 
% were correctly classified to their correct geographic zone).

% Perform CAP using optimal m:
cap = f_cap(biotic_2,'bc',geo,[],0,1000,1,12,1);
% ==================================================
%  CAP - Canonical Discriminant Analysis:
% --------------------------------------------------
% Trace Stat    = 2.0624  p =  0.00100 
% Greatest Root = 0.9050  p =  0.00100 
% No. of permutations = 1000 
% --------------------------------------------------
% No. of axes of Q used (m)     = 12 
% Variability of yDis explained = 79.74 % 
% Canonical Correlations:
%   0.9513  0.8670  0.6370
% Squared Canonical Correlations (= delta^2):
%   0.9050  0.7517  0.4058
% ==================================================
% 
% -> We reject the null hypothesis of no significant difference in species
% composition and abundance among zones, based on obtaining statitically
% significant Trace and the Greatest Root statistics (at alpha = 0.05). The
% 3 canonical correlations provide a relative measure of how well canonical
% axes I, II, and III separate the sites into distinct geographic zones.
% Obviously, axis I separates sites according to zone better than axis II
% and III, etc.

% ==================================================
%             LOO CROSS-VALIDATION
%             Classification Success: 
% --------------------------------------------------
% Group        Correct  
%    1             66.7 % 
%    2             70.0 % 
%    3             83.3 % 
%    4             81.8 % 
% 
% 
% Total Correct  = 76.92 % 
% Total Error    = 23.08 % 
% 
% -> Overall, the discriminant functions are able to correctly classify an
% unknown observation 77% of the time, which indicates the benthic
% macrofaunal assemblages associated with each of the 4 zones were quite
% distinct.

% --------------------------------------------------
%      Confusion Matrix (%): 
% group: 1      2      3      4      
%      1   66.7   33.3    0.0    0.0 
%      2    0.0   70.0   30.0    0.0 
%      3    0.0    0.0   83.3   16.7 
%      4    0.0    0.0   18.2   81.8 
% 
% ================================================== 
% 
% -> It is clear from the confusion matrix that mis-classification of sites
% only occurred between contiguous zones. This is not surprising, since we
% lumped sites into distance classes using arbitrary boundaries. Also
% (except for zone 4), mis-classifications were always the result of sites
% being assigned to the zone one step further from the center of drilling
% activity than they actually occurred (rather than closer). 

% Create group labels for plot legend:
gLabels  = {'zone 1' 'zone 2' 'zone 3' 'zone 4'}; % group labels for legend
% gLabels = cellstr(num2str(f_unique(geo)))';      % alternative method

% Create canonical plot:
f_capPlot(cap,gLabels,[],biotic_2,biotic_txt,0.05,'none');

% ----Optional:-----
% Plot vectors for the 5 most important taxa for each axis:
[h,yCor] = f_capPlot(cap,gLabels,[],biotic_2,biotic_txt,0.05,'none');
[null,idx_I]  = sortrows(abs(yCor),-1); % sort rows descending by column 1
[null,idx_II] = sortrows(abs(yCor),-2); % sort rows descending by column 2
idx = f_unique([idx_I(1:5);idx_II(1:5)]);
f_capPlot(cap,gLabels,[],biotic_2(:,idx),biotic_txt(idx),0.05,'none');

% 3) Which taxon displayed the greatest difference in abundance
% between zones 3 and 4?
% 
% ANSWER: Amphictene_auricoma, which had the longest vector associated with
% canonical axis II

% 4) Which taxon displayed the greatest difference in abundance
% between zones 1 and 4?
% 
% ANSWER: Chaetozone_setosa, which had the longest correlation vector associated
% with canonical axis I.

% 5) Provide a brief, overall ecological interpretation of the results you
% obtained from the Ekofisk data analysis.
% 
% ANSWER: There is a gradient (gradation of change) associated with the
% distance from the center of oil drilling activity that imposed a
% significant effect on the distribution and abundance of the macrobenthos
% data.

% 6) What do the trace statistic, greatest root statistic, and their
% associated p-values generated by the CAP procedure indicate?
% 
% ANSWER: the trace statistic is derived from the sum of the canonical
% eigenvalues and, along with its associated p-value, allows one to test
% the null hypothesis that there is no significant difference in the
% abundance and composition of species among geographic zones. The greatest
% root statistic is similar, but is based on the first canonical eigenvalue
% rather than the sum of all of them.

% 7) How does Leave-One-Out cross-validation work?
% 
% ANSWER: After all observations are assigned membership to one of the _a
% prior_ specified groups, one observation is placed in a test set while
% the remaining observations are placed in a training set. The training set
% is used to create a discriminant function and the resulting canonical
% eigenvectors are used to define the canonical space which best separates
% observations by their group memberships. The observations in the training
% set are then projected into this canonical space and the coordinates of
% each group's central tendency (centroid or spatial median) in this
% canonical space is calculated. The observation left out of this process
% (the test set) is then projected into the canonical space and is assigned
% membership to the group whose centroid (or spatial median) it is in
% closest proximity to. If this assignment (predicted group membership) is
% the same as the actual group membership, then the classification was
% successful. This entire procedure is repeated over and over until each
% observation have served once as a test set. The successful and unsuccessful
% classifications are then tabulated to yield an overall classification
% success rate and also broken down for each group in a confusion matrix.

% 8) How do you interpret the confusion matrix generated by the CAP procedure?
% 
% ANSWER: The confusion matrix breaks down the classification success rates
% for each group, with the diagonal elements indicating the proportion of
% sites that were correctly classified. The non-diagonal elements
% indicate where the confusion occurred when there were observations that were
% mis-classified.

% 9) What does it mean for a CDA to be 'overparameterized'?
% 
% A discriminant analysis is said to be 'overparameterized' if the number
% of variables in the model approaches or exceedes the number of
% observations. One of the disadvantages of this is that (either by chance 
% or by simply adding enough variables) groups may appear very distinct in
% a canonical plot, but the model has a low classification success rate.  

% NOTE: For questins 6-9, consult the Anderson & Willis (2003) as necessary.
