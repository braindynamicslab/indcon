%WGGG
%
%
% aim - extract individual contribution from group-derived SCM using two
% approaches: [1] Leave One out (LOO) and Add One patient (AOP). See the 
% article for more info - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4589508/
% 
% if you use this code in your work, please cite the following paper:
%
%   Saggar, M. et al. Estimating individual contribution from group-based 
%   structural correlation networks. NeuroImage 120, 274?284 (2015).
% 
%
% author - saggar@stanford.edu
% date - 5.27.2013

%% load data

fxs = load('fxs_freesurferData.mat'); % size = subjects x regions for patients group
td = load('td_freesurferData.mat'); % size = subjects x regions for typically developing group
region = load('region_names'); % cell containing labels for each region

%% extract individual contribution using two methods: LOO and AOP
% approach LOO: Leaving one participant out to extract the individual participant's contribution
% Note: We use f_mantel function function from The Fathom Toolbox for Matlab: multivariate ecological and oceanographic data analysis. College of Marine Science, University of South Florida, St. Petersburg, Florida, USA. Available from: http://www.marine.usf.edu/user/djones/
nsub_fxs = size(fxs,1);
nsub_td = size(td,1);
ic_fxs_mantel = [];
for i = 1:1:size(fxs,1)
   ic_fxs_mantel(i) = f_mantel(putDiagZeros(corr(fxs)),putDiagZeros(corr(fxs(setdiff(1:nsub_fxs,i),:))));
end

ic_td_mantel = [];
for i = 1:1:size(td,1)
   ic_td_mantel(i) = f_mantel(putDiagZeros(corr(td)),putDiagZeros(corr(td(setdiff(1:nsub_td,i),:))));
end

% approach AOP: Adding one patient to the typical developing group to extract AOP-based ind. contribution
% Note: We use f_mantel function function from The Fathom Toolbox for Matlab: multivariate ecological and oceanographic data analysis. College of Marine Science, University of South Florida, St. Petersburg, Florida, USA. Available from: http://www.marine.usf.edu/user/djones/
ic_fxs_td_mantel = [];
for i = 1:1:size(fxs,1)
   ic_fxs_td_mantel(i) = f_mantel(putDiagZeros(corr([td; fxs(i,:)])),putDiagZeros(corr(td)));
end

% calculating distance & normalizing to get Z-value
ic_fxs_mantel = 1 - ic_fxs_mantel;
ic_fxs_mantel_z = (ic_fxs_mantel-mean(ic_fxs_mantel))./std(ic_fxs_mantel); 

ic_td_mantel = 1 - ic_td_mantel;
ic_td_mantel_z = (ic_td_mantel-mean(ic_td_mantel))./std(ic_td_mantel);

ic_fxs_td_mantel = 1 - ic_fxs_td_mantel;
ic_fxs_td_mantel_z = (ic_fxs_td_mantel-mean(ic_fxs_td_mantel))./std(ic_fxs_td_mantel);


%% extract regional differences and plot them.
tmp=[];
tmp1=[];
tmp2=[];
n_regions = length(region);
for i = 1:1:nsub_td 
    tmp(i,:)=mean(abs(corr([td; fxs(i,:)]) - corr(td))); 
    tmp1(i,:)=mean(abs(corr(fxs)-corr(fxs(setdiff(1:size(fxs,1),i),:))));
    tmp2(i,:)=mean(abs(corr(td)-corr(td(setdiff(1:size(td,1),i),:)))); 
end
figure; 

regionnew={};
for i=1:1:n_regions
    regionnew{i}=strrep(region{i},'_',' '); 
end

errorbar(1:n_regions, mean(tmp), std(tmp)./sqrt(nsub_td),'ro-.');
set(gca,'LineWidth',2);axis('tight');
xticks(1:86);
xticklabels(regionnew);
xtickangle(90);

%% stability analysis

for p = 1:1:10
    k = 1;
    for j = 2:1:size(td,1)
        td_subs = randperm(size(td,1));
        td1 = td(td_subs(1:j),:);
        ic_fxs_td_aop = [];
        for i = 1:1:size(fxs,1)
           ic_fxs_td_aop(i) = f_mantel(putDiagZeros(corr([td1; fxs(i,:)])),putDiagZeros(corr(td1)));
           
        end
        corr_stability(p,k) = corr(ic_fxs_td_mantel', ic_fxs_td_aop');
        k=k+1;
    end
end
