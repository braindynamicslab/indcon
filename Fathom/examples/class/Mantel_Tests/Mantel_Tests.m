%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                Mantel Tests:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% by David L. Jones, 2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Tuomisto et al. (2003) examined the distribution of plant species from
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        ANALYSIS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% It has been hypothesized that the plant species composition of Amazonian
% tierra firme forests are uniform over large areas, which suggests forests
% are dominated by a limited suite of competitively superior tree species.
% If this hypothesis is true, then differences in the composition of plant
% species among sites (= beta diversity) should be relatively uniform and
% not vary with differences in geographic distance or environmental
% conditions among sites. Use a series Mantel Tests to test address the
% following questions:

% 1) Does the variation in beta diversity of Pteridophytes from tierra
% firme forests depend on the geographic distance among sites?
% 
% 2) What about the Melastomataceae?
% 
% 3) Does the variation in beta diversity of Pteridophytes from tierra
% firme forests depend on differences in the soil condition among sites?
% 
% 4) What about the Melastomataceae?
% 
% 5) What were your null hypotheses?
% 
% 6) What can you conclude from your results?

% HINTS: The 'f_latlong' function may be useful for calculating the
% distances among sites. Since the focus of your analysis is on the tierra
% firme sites, you'll need to use the 'find' (or 'logical' or 'ismember')
% command to query (i.e., index) the apprproiate subset of the data.

% -----Answer:-----
% Load the data:
load amazon.mat
 
% Create a matrix of geographic distances among sites:
geo = f_latlong([site.lon site.lat]);

% Get index to tierra firme sites:
idx = find(site.ter==1);
% idx = logical(site.ter);    % alternative form
% idx = ismember(site.ter,1); % alternative form

% Correlation betweeen differences in beta diversity and geographic distance:
[r,p] = f_mantel(bio.Pter(idx,idx),geo(idx,idx),1,1000)
% r = 0.40045
% p = 0.001
% 
[r,p] = f_mantel(bio.Mel(idx,idx),geo(idx,idx),1,1000)
% r = 0.47402
% p = 0.001

% Correlation betweeen beta diversity and environmental distance:
[r,p] = f_mantel(bio.Pter(idx,idx),soil.dis(idx,idx),1,1000)
% r = 0.46213
% p = 0.001
% 
[r,p] = f_mantel(bio.Mel(idx,idx),soil.dis(idx,idx),1,1000)
% r = 0.44841
% p = 0.001


% Conclusions: I tested the following 4 null hypotheses: there was no
% significant correlation between variation in beta divesity of
% Pteridophytes (or Melastomataceae) with geographic distance (or soil
% conditions). On all accounts, the null hypothesis was rejected. The
% tierra firme forests are NOT uniformly inhabitated by tree species. Sites
% that were further apart (or with more dissimilar soil conditions) tended
% to be more different, in terms of trees species composition, than sites
% that were closer (or with more similar soil conditions).


% 7) How does the Jaccard index differ from the Bray-Curtis metric?
% 
% ANSWER: It is similar, but the Jaccard index is used for binary
% (presence/absence data). So, the Jaccard index is a measure of beta
% diversity based on species composition while the Bray-Curtis metric is
% based on species composition and abundance.

% 8) Create a PCoA ordination diagram of the shrub data. Visually inspect
% the plot to determine whether there appear to be any regions patterns in
% these data.

% ANSWER:

% Create a tierra firme subset:
tf.reg  = site.reg(idx);
tf.Pter = bio.Pter(idx,idx);

% PCoA:
pcoa = f_pcoa(tf.Pter,0,1,1);

idx_1 = ismember(tf.reg,'N_Peru');
idx_2 = ismember(tf.reg,'Ecuador');
idx_3 = ismember(tf.reg,'Colombia');
idx_4 = ismember(tf.reg,'S_Peru');

figure; hold on;
h(1) = plot(pcoa.scores(idx_1,1),pcoa.scores(idx_1,2),'bo');
h(2) = plot(pcoa.scores(idx_2,1),pcoa.scores(idx_2,2),'r.');
h(3) = plot(pcoa.scores(idx_3,1),pcoa.scores(idx_3,2),'g*');
h(4) = plot(pcoa.scores(idx_4,1),pcoa.scores(idx_4,2),'k^');

legend(h,f_unique(tf.reg),'Interpreter','none');

box on;
xlabel('PCoA Axis I');
ylabel('PCoA Axis II');

% ALTERNATIVE (not as pretty, but still works):
f_pcoaPlot(pcoa,tf.reg,[],[],0,'none')

% CONCLUSIONS:
% There are clear regional pattern in the shrub data. Most notably the
% northern sites (N_Peru, Ecuador, and Colombia) appear distinctly
% separate from S_Peru. This suggest climatic differences may influence
% species composition.


% 9) It has been suggested that distance-based methods, such as ANOSIM and
% the Mantel Test, are not necessarily your best choice for analyzing beta
% diversity. What is beta diversity and why would these methods not be your
% method of choice for analyzing it? Keep your answers brief (1-2
% paragraphs at the most) and to the point, yet thorough and grammatically
% correct. Consult one or more of the following references (available on
% Blackboard) as necessary.
% 
% Legendre, P., D. Borcard, and P. R. Peres-Neto. 2005. Analyzing beta
% diversity: partitioning the spatial variation of community composition
% data. Ecological Monographs 75: 435-450.
% 
% Tuomisto, H. and K. Ruokolainen. 2006. Analyzing or explaining beta
% diversity? Understanding the targest of different methods of analysis.
% Ecology 87(11): 2697-2708
% 
% Pélissier, R., P. Couteron, and S. Dray. 2008. Analyzing or explaining
% beta diversity: Comment. Ecology 89: 3227-3232.
% 
% Laliberté, É. 2008. Analyzing or explaining beta diversity: Comment.
% Ecology 89: 3232-3237.
% 
% Legendre, P., D. Borcard, and P. R. Peres-Neto. 2008. Analyzing or
% explaining beta diversity: Comment. Ecology 89: 3238-3244.
% 
% Tuomisto, H., and K. Ruokolainen. 2008. Analyzing or explaining beta
% diversity: Reply. Ecology 89: 3244-3256.

% ANSWER: Beta diversity is the variation in the abundance and composition
% of species among sites. It can be measured as the total sum-of-squares (=
% variance) of a data set consisting of the abundance/presence of species
% by sites. Alternatively, it can be measured as the mean of the
% dissimilarities computed from the same data set, but not by the variance
% of the dissimilarity matrix.
% 
% Using canonical analysis (e.g., RDA) to examine raw data (species
% abundance/composition among sites) allows one to analyze Beta diversity.
% Using distance-based methods (e.g., Mantel Tests or ANOSIM) to examine
% dissimilarity matrices based on species abundance/composition data among
% sites from different regions allows one to explain differences in Beta
% diversity among those regions. Canonical analysis correctly partitions
% the variation in species abundance/composition data, while distance-based
% methods do not because they underestimate the percent variation
% explained. The latter methods should only be used to examine the
% variation in Beta diversity, not Beta diversity itself. Unlike
% distance-based methods, canonical analyses provide more powerful
% significance tests, estimates of the contribution (importance) of
% individual variables in the overall response-predictor relationship, and
% ordination plots that summarize (visualize) these relationships.
