%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                Matlab Basics:                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% by David L. Jones, Sep-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Load data:
load biscayne.mat

% Transpose matrices so ROWS = observations, COLS = variables:
X = X';
Y = Y';

% Inspect the data:
size(X) % determine # rows and cols
size(Y) 

% Display the data:
X
Y

% Use a cell array of strings to create column labels:
X_txt = {'tra' 'str' 'sea' 'yr' 'tem' 'do' 'sal' 'dep' 'vis' 'E' 'N'  'jdate' 'FW'};
Y_txt = {'Lut_gris'}'

% Use comments to define abbreviations:
% Variables:
% Y = # gray snapper per 60 m^2
% X = matrix of predictors
% 
% X_txt  = cell array of predictor labels as follows:
% tra    = transect
% str    = strata (1=Key, 2=Mainland)
% sea    = season (1=dry, 2=wet)
% yr     = YYYY
% tem    = water temperature (degrees C)
% do     = dissolved O2 (mg/L)
% sal    = salinity (PPT)
% dep    = average depth (cm)
% E      = UTM Easting
% N      = UTM Northing
% jdate  = julian date
% FW     = distance (km) from nearest FW canal


% Calculate some basic statistics, columnwise:
mean(X)
min(X)
max(X)
% -> note: transpose the output to improve readability

% Use a FOR loop to plot each variable separately (press any key to continue):
for i=1:size(X,2)
   plot(X(:,i))
   title(X_txt(i))
   disp('Press any key to continue...');
   pause
end

% Close all figure windows:
close all;

% Extract 'tem' and 'do' as separate variables:
W = X(:,5); % temperature
Z = X(:,6); % dissolved oxygen

% Plot each variable (are either normally distributed?):
plot(W);
plot(Z);

% Cleanup the workspace by closing all previous figure windows:
close all;

% Plot each variable in the SAME window:
figure; hold on;
plot(W,'b-');    % blue lines
plot(Z,'k:');    % black dots
plot(Z*10,'r-'); % red lines, multiplied by 10

% Create a new figure & plot the distribution of each:
figure; hist(W,25); title('Temperature');
figure; hist(Z,25); title('Dissolved Oxygen');

% Create a new variable by concatenating W & Z:
A = [W Z]; % horizontally
B = [W;Z]; % vertically

% Examine their sizes:
size(A)
size(B)


% Clean up the workspace:
close all;
clear A B;

% Transform the data:
W_log = log10(W); % log base 10 transformation
Z_sq  = Z.^0.5;   % square root transform

figure; hold on;
plot(W,'k-');
plot(W_log,'r-');

figure; hold on;
plot(Z,'k-');
plot(Z_sq,'r-');

% How would you transform the data via a natural log?
% Is there another way to perform a square-root transform?

% -----Fun with Queries:-----
% Create data:
A = [1 1 2 2 3 3 3 5 5 9 9 7 7 7 4 4 4 4 6 8 2 5 9 6 5 3 5 6 7]';

% Caluclate the mean:
mean(A)

% Calculate the mean for those values less than 6:
idx = find( (A<6) == 1 ); % get indices to elements less than 6
A(idx)                    % examine the results of your query
mean(A(idx))

% Calculate the mean for values > 2 and < 9 
idx = find( ((A>2) & (A<9))  == 1 ); % get indices to elements matching query
A(idx)                               % examine the results of your query
mean(A(idx))

% -----Determine the average temperature in the dry vs. wet season:-----
idx_D = X(:,3)==1; % logical indices to DRY season
idx_W = X(:,3)==2; % logical indices to WET season

% Compare logical indices to original data:
[X(:,3) idx_D idx_W]

mean(X(idx_D,5)) % mean temperature of DRY season
mean(X(idx_W,5)) % mean temperature of WET season

% Test the null hypothesis that there is no significant difference in
% TEMPERATURE among SEASON:
H = ttest2(X(idx_D,5),X(idx_W,5))

% What do the results indicate?
% Repeat the analysis, but perform an appropriate 2-tailed test.
% Perform a similar analysis, but examine SALINITY in terms of STRATA.

% Is there a relationship between TEMPERATURE and SALINITY?
T = X(:,5); % extract temperature
S = X(:,7); % extract salinitys

% Examaine the correlation between T and S:
corr(T,S)
% What is R?
% What is R-squared?

% Is there a significant correlation?
[rho,pval] = corr(T,S,'tail','both')

% Plot T vs. S:
figure;
plot(T,S,'k.');
xlabel('Temperature');
ylabel('Salinity');

% Perform a simple regression (requires the FATHOM Toolbox for Matlab):
model = f_mregress(T,S,0,0,1);
figure; hold on; plot(T,S,'k.')
plot(T,model.yfit,'r.')

% Examine the relationship between the output of 'corr' and that of
% f_mregress:
sqrt(model.R2)
abs(corr(T,S))
