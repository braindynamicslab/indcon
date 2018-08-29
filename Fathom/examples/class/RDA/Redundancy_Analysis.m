%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Redundancy Analysis (RDA):                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% by David L. Jones, 29-Nov-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% The file 'tankers.mat' contains data from Aguirre-Macedo et al. (2008)
% concerning abundances of phytoplankton and pathogenic bacteria inhabiting the
% ballast water of 30 oil tankers docked at a PEMEX crude oil terminal in
% Mexico's Campeche Sound, along with associated environmental/nutrient data.
% Also included are similar data for coral reefs within the nearby Cayo Arcas
% reef system. The variables included in this data set are as follows:

% -----Variables:-----
% Y     = biotic and associated variables
% Y_txt = cell array of corresponding column labels as follows:
%   'resi'  = time of water residence (days)
%   'orig'  = ballast water origin (1=river, 2=coastal, 3=oceanic, 4=reef)
%   'totc'  = Total Coliforms
%   'feca'  = Fecal Coliforms
%   'coli'  = E. coli 0157
%   'coli2' = Other E. coli
%   'ente'  = Enterococcus sp.
%   'vibr'  = Vibrio cholerae
%   'serr'  = Serratia marcescens
%   'begg'  = Beggiatoa spp.
%   'desu'  = Desulfovibrio spp.
%   'phor'  = Phormidium corallyticum
%   'sphi'  = Sphingomonas spp.
%   'chla'  = Chlorophyll a
%   'chlb'  = Chlorophyll b
%   'chlc'  = Chlorophyll c
%   'caro'  = Carotenes
%   'tri'   = Trix index of eutrophic condition
% 
% X    = environmental/nutrient variables
% X_txt = cell array of corresponding column labels as follows:
%   'T'    = temperature
%   'S'    = salinity
%   'con'  = conductivity
%   'pH'   = pH
%   'alk'  = alkalinity
%   'CO2'  = total carbon dioxide
%   'tur'  = turbidity
%   'DO2'  = dissolved oxygen
%   'O2.S' = oxygen saturation
%   'NH4'  = ammonium
%   'NO2'  = nitrite
%   'NO3'  = nitrate
%   'DIN'  = dissolved inorganic nitrogen
%   'P'    = phosphates
%   'Si'   = silicate
%   'TSS'  = total suspended solids
%   'OSS'  = organic suspended solids
%   'ISS'  = inorganic suspended solids
%   'WQI'  = water quality index
% 
% sites     = codes for each of the 30 tankers + 2 reef sites
% sites_txt = cell array of corresponding text labels

% Aguirre-Macedo et al. (2008) used a classical RDA approach to analyze these
% data. Your goal is to repeat their analysis, but use a method more appropriate
% for ecological data. Therefore, use distance-based redundancy analysis to
% create a canonical ordination (constrained PCoA) of the tanker biotic data
% that depicts the variation among sites in terms of composition and abundance
% of microbes + chlorophyll-a that is explained by the environmental/nutrient
% data. However, limit your analysis to those predictor variables that
% Aguirre-Macedo et al. (2008) considered relevant in their analysis. Create a
% 2-d db-RDA ordination distance biplot with vectors for both the response and
% predictor variables. Provide a statistical interpretation of the ordination
% diagram in terms of the relationship(s) between the response and predictor
% variables. Do your results differ from those of Aguirre-Macedo et al. (2008),
% especially their Fig. 2? If so, how and why?

% HINT: Use the 'f_rdaDB' and 'f_rdaPlot' functions. For 'f_rdaDB', set W=0
% because you are not working with covariables. For 'f_rdaPlot' set wascores=0,
% offset=0, fmt='none', iter=0. Setting scale=1 won't affect the size of the
% biplot vectors, but you may want to experiment with this parameter, setting it
% to smaller or larger values so your plots are easier to read. Once you've
% decided on a value for 'scale', try adding a small value to 'offset' so the
% endpoints of the biplot vectors don't overlap with their corresponding text
% labels. Also, you might try examining the data before proceeding and make use
% of summary statistics such as mean, range, etc. to create and trim new
% variables as needed.

% Load the data:
load tankers.mat;

% Create new variable + text labels for biotic data:
bio     = Y(:,3:14);
bio_txt = Y_txt(3:14);

% Examine biotic data:
[min(bio); max(bio); mean(bio)]'

% Remove those taxa that are always 0:
bio(:,8:10)   = [];
bio_txt(8:10) = [];

% Create a symmetric dissimilarity matrix using an ecologically relevant metric,
% transform the data to down weight the influence of the more abundant taxa:
disBC = f_dis(f_normal(bio,'4'),'bc');
% 
% The biotic response data consists of 9 variables (8 microbes + chlorophyll-a),
% some provided as abundance and others as presence/absence. Thus I used a
% relatively strong (4th root) transformation to down weight the importance of
% the more abundant taxa. Since 5 of the 9 biotic response variables are
% provided as presence/absence, normally I would have converted all the data to
% presence/absence before calculating a Bray-Curtis dissimilarity (or the
% equivalent would be just using the Jaccard coefficient). However, since we're
% more/less trying to follow the approach of Aguirre-Macedo et al. (2008) I
% retained the abundance data.

% Create new variable + text labels predictors:
idxX = [10 8 14 2 15 1 6 7]'; % get index to 'relevant' variables
env     = X(:,idxX);
env_txt = X_txt(idxX);

% Perform db-RDA:
rda_DB = f_rdaDB(disBC,size(bio,2),env,0,1000,1);
% 
% ==================================================
% REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F = 2.0857    p    =  0.01400 
% R2 = 0.4204   R2adj =  0.21887 
% No. of permutations = 1000 
% --------------------------------------------------
% 
% Canonical Eigenvalues:
%   0.9779  0.2561  0.0459  0.0184  0.0093  0.0000  0.0000  0.0000
% Residual Eigenvalues:
%   0.8588  0.4432  0.2216  0.1564  0.1167  0.0515  0.0400
% 
% Species-Environment Correlations (r):
%   0.8706  0.4959  0.3281  0.2976  0.1049  0.7958  0.6815  0.6050
% 
% Fraction of variance explained:
% ------------------------------
% Canonical axes (total = 0.4091): 
%   0.3060  0.0801  0.0144  0.0058  0.0029  0.0000  0.0000  0.0000
% Cumulative: 
%   0.3060  0.3861  0.4005  0.4062  0.4091  0.4091  0.4091  0.4091
% 
% Residual axes  (total = 0.5909):
%   0.2687  0.1387  0.0693  0.0490  0.0365  0.0161  0.0125
% Cumulative: 
%   0.2687  0.4074  0.4767  0.5257  0.5622  0.5783  0.5909
% ==================================================

% Create RDA distance biplot ordination:
f_rdaPlot(rda_DB,f_normal(bio,'4'),0,[1 1/50],[0.025],bio_txt,env_txt,...
   'none',0,sites_txt);

% INTERPRETATION:
% -> Based on these results, I would reject the null hypothesis of no
% relationship between the response and predictor variables at alpha = 0.05.
% Thus, there is a significant relationship between the response and predictor
% variables, with this subset of the environmental and nutrient variables
% explaining a fair portion (42%) of the total variation in species abundance
% and composition among sites. Note that a large portion of the explained
% variation (30%) is accounted for by just one axis (Canonical Axis I).
% 
% Overall, the most important conclusions to be derived from the ordination
% diagram are:
% 1) Salinity (S), silicates (Si), turbidity (tur), and dissolved oxygen (DO2)
% are the 4 most important predictors (in decreasing order of importance) in
% terms of explaining the differences in species abundance and composition among
% the 32 sites.
% 
% 2) 'feca' (fecal coliforms), 'totc' (total coliforms), and 'ente'
% (Enterococcus sp.) are the 3 taxa displaying the greatest variability, in
% terms of abundance and composition, among sites. Both 'feca' and 'totc' are
% assoicated with lower salinities, higher silicates, higher turbidities, and
% lower dissolved oxyen values, while 'ente' are associated with the opposite
% trends.
% 
% Sites 13,26,27,30 were placed on opposite ends of the S/DO2/Si/tur gradient
% compared to sites 7,9.
% 
% My results are similar to those provided by Aguirre-Macedo et al. (A-M), but
% with some important distinctions. First, their ordination diagram shows a
% similar overall pattern with their plot flipped along the horizontal plane.
% A-M used a classical RDA with implicitly uses the Euclidean distance measure
% instead of my use of a Bray-Curtis based db-RDA. Thus, my plot provides a more
% realistic picture of differences among sites in terms of species abundance and
% composition, while their method is overly sensitive to total abundances
% (across taxa) and double-zeros. For example, my plot shows sites 1 and 8
% separated from most other sites, while A-M's plot indicates these sites are
% similar to sites 2,4,5 as well. Looking at the raw data tables indicates, for
% example, sites 1,2,5 have few species in common but share similar absences
% of species.

% Also fecal and total coliforms are supposed to be the two most important taxa
% driving the variation among sites in A-M's ordination biplot, since it
% displays the longest vector and is aligned parallel to Axis I. However, sites
% 9 and 10 are separated by an extreme distance in their plot despite having
% IDENTICAL values of these taxa in the raw data tables. In contrast, the db-RDA
% ordination plot offers a better representation of the similarity among these
% sites. This is most likely due the classical RDA's reliance on a distance metric
% that treats sites displaying similar absences of taxa as similar.
