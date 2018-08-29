%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Ecological Dissimilarity:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% by David L. Jones, 29-Sep-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% 1) The location of 6 sampling sites are depicted in the figure in the
% file 'plot.pdf'. Print out a hard copy of this figure and bring it to
% class.
% 
% 2) Using a ruler, measure the distances (in Cartesian space) between all
% possible pairs of sites. Record them on this page in decimal inches or metric units:
% 
% 1 & 2 =  46 mm
% 1 & 3 = 129 mm
% 1 & 4 = 130 mm
% 1 & 5 = 147 mm
% 1 & 6 = 107 mm
% 2 & 3 =  91 mm
% 2 & 4 =  86 mm
% 2 & 5 = 102 mm
% 2 & 6 =  67 mm
% 3 & 4 =  37 mm
% 3 & 5 =  69 mm
% 3 & 6 =  96 mm
% 4 & 5 =  28 mm
% 4 & 6 =  62 mm
% 5 & 6 =  58 mm
% 
% ANSWERS will vary according to the units/precision that interpoint
% distances were measured.
 
% 3) Import these data into Matlab as a column vector named 'vec' along with a
% cell array of site labels:
% >> vec     = [value1 value2 value3...]';
% >> vec_txt = {'1' '2' '3' '4' '5' '6'}; 

vec     = [46 129 130 147 107 91 86 102 67 37 69 96 28 62 58]';
vec_txt = {'1' '2' '3' '4' '5' '6'}'; 

% 4) Wrap these data up into a distance matrix named 'dis' using either the
%    'f_rewrap' or 'squareform' function.

dis = f_rewrap(vec);

% 5) Examine matrix 'dis'. What do you notice about the placement of the values
% from 'vec', its overall shape, and how it compares to its transpose?
% 
% ANSWER: It is a square, symmetric dissimilarity (distance) matrix with all
% zeros along the diagonal indicating an object is 0 units from itself, the upper
% tridiagonal is a mirror image of the lower tridiagonal, and it is identical to
% its transpose.

% 6) Load the data in the file 'crds.dat' and use the 'plot' command to
% visualize these coordinates; use the 'text' command to ADD the point
% labels (vec_txt) to your plot.

load crds.dat
figure;
plot(crds(:,1),crds(:,2),'bo');
text(crds(:,1)+0.35,crds(:,2),vec_txt); % labels are offset by 0.3 units along X axis
axis([-100 -65 24 45]);                 % increase bounds of axis

% 7) How does your Matlab plot compare to the figure in 'plot.pdf'?
% 
% ANSWER: they should be nearly identical.

% 8) Create a Euclidean distance matrix from the data imported from
% 'crds.dat' using the 'f_dis' function. Call it 'edis'.

edis = f_dis(crds,'euc');

% 9) Calculate a simple statistic to examine the congruence between the
% interpoint distances you measured for 'dis' and those calculated for
% 'edis'. Use the 'f_unwrap' command to unwrap the lower tri-diagonals of
% each symmetric dissimilarity matrix into a column vector named 'vec_dis'
% and 'vec_edis', respectively. Use the 'corrcoef' function to calculate the
% Pearson correlation coefficient between these two vectors. Note the
% 'corrcoef' command generates a symmetric matrix of correlations, we are
% only interested in retaining the 2nd value that is returned in 'r' (i.e.,
% r(2)). Save that value in a variable named 'myR'.

% Unwrap lower tri-diagonals as a vector:
vec_dis  = f_unwrap(dis);
vec_edis = f_unwrap(edis);
% 
% Get observed value of r:
r    = corrcoef(vec_dis,vec_edis);
myR  = r(2);
% myR = 0.98537

% 10) What does the statistic captured in the variable 'myR' tell you about
% the congruence between 'dis' and 'edis'?
% 
% ANSWER: 'myR' is used to measure the congruence between the two distance
% matrices 'dis' and 'edis' by providing a measure of the linear
% correlation between the two sets of interpoint distances. Since we're
% interested in comparing the relative distances among one distance matrix to
% those of the other, the ABSOLUTE values of distance are not important.
% Values of Pearson's 'r' range from 0 to 1 indicating no linear
% correlation or perfect linear correlation, respectively. In this example,
% my value was 0.9853, which indicates the two distance matrices were
% nearly perfectly identical in terms of relative inter-point distance
% among sites. Pearson's 'r' could have equaled 1 if there were no
% measurement errors when using a ruler to gauge distances (and/or there
% were no round-off errors).

% 11) Load the file 'spiders_env.mat', which contains data from van der
% Aart & Smeenk-Enserink (1975) on the abundances of Hunting spiders along
% with associated environmental data. The variable 'spiders' contains the
% counts of 12 species of spiders from 28 sampling sites, with
% corresponding column labels in 'spiders_txt'. The variable 'env' provides
% the values of 6 environmental variables (each measured on a different
% scale) collected from the same 28 sampling sites, with corresponding
% column labels in env_txt.
% 
% van der Aart, P. J. M. & N. Smeenk-Enserink. 1975. Correlations between
% distributions of hunting spiders (Lycosidae, Ctenidae) and environmental
% characteristics in a dune area. Netherlands Journal of Zoology 25: 1-45.

load spiders_env.mat;

% 12) Create a Bray-Curtis dissimilarity matrix for the spider data named BC
% and a Euclidean distance matrix for the environmental data named ED using the
% 'f_dis' function. Note the biotic data need to be transformed in order to
% down weight the influence of the more abundant taxa. Also, the abiotic
% data need to be standardized because these variables were measured using
% different units. See the 'f_transform' and 'f_stnd' functions. Examine
% the dimensions of the resulting distance matrices to make sure they were
% created correctly.

BC = f_dis(f_transform(spiders,2),'bc'); % transform biotic data via fourth-root
ED = f_dis(f_stnd(env),'euc');           % standardize abiotic data

% Check that dimensions are compatible:
size(spiders,1)
size(BC)
size(env,1)
size(ED)

% -> both 'spiders' and 'env' have 28 observations, so the 2 distance
% matrices should be 28 x 28.

% 13) Perform a simple test to determine whether the environmental factors
% measured at the sampling sites affect the distribuion/abundance/species
% composition of the spiders observed at the same sites. See step 9 above.

% Unwrap lower tri-diagonals as a vector:
vec_BC = f_unwrap(BC);
vec_ED = f_unwrap(ED);

% Get observed value of r:
r    = corrcoef(vec_BC,vec_ED);
obsR = r(2);

% 14) Is there evidence of a relationship between the biotic data and the
% measured abiotic factors in this dataset? Would kind of underlying assumptions
% would you have to make in order to arrive at this conclusion? In your opinion,
% would you consider any such effects to be symmetrical or asymmetrical?
% 
% ANSWER:
% Using the following command a parametric p-value associated with the
% correlation coefficient provides evidence of a significant relationship,
% assuming the data are independent-identically-distributed, normally
% distributed, randomly sampled, etc.
[r_par,p_par] = corrcoef(vec_BC,vec_ED); % parametric test

% It is hoped most people would consider such relationships to be asymmetrical;
% i.e, in this case the environmental factors are influencing the
% distibution/abundance of the spiders but NOT necessarily vice versa. This is
% an important distincion to make and will determine, for example, whether you'd
% subsequently perform Canonical Correlation Analysis (symmetrical) vs. RDA/CCA,
% etc.

% KUDOS to anyone who figures out this is a MANTEL test and/or
% creates a randomization test to assess the significance of r.
% 
[r_man,p_man] = f_mantel(BC,ED,0,1000);       % Mantel Test
[r_rnd,p_rnd] = f_corr(vec_BC,vec_ED,0,1000); % Randomization Test
