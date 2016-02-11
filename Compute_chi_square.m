%% IRM classes member
n_insect = size(Zinsect,1);
n_plant = size(Zplant,1);

for i = 1: n_insect
    Insect_IRM_class_cell{i} = UniC1(Zinsect(i,:)==1,:);
end;

for i = 1: n_plant
    Plant_IRM_class_cell {i} = UniC2(Zinsect(i,:)==1,:);
end;
%%


%% Insect_family Chi_sqared_test

%  Insect_original_family class
Insect_original_family = Insect(:,4);
Insect_original_family = Insect_original_family(ia1,:); 

%  Unique Insect family 
[Insect_family_unique_class,i_f_ia,i_f_ic] = unique(Insect_original_family); 

%  Insect_family_class_cell
for i=1:size(Insect_family_unique_class,1)
    Insect_family_class_cell{i} = UniC1(i_f_ic == i,:); 
end;

% The joint relational metrix of Insect_family and RIM
for i=1:size(Insect_IRM_class_cell,2) 
    for j=1:size(Insect_family_unique_class,1)
        R_IRM_and_insect_f(i,j) = numel(intersect(Insect_IRM_class_cell{i},Insect_family_class_cell{j})); 
    end; 
end; 

% Chi-squared test for Inscet_family_IRM
P_inscet_family_IRM = chi2Tests(R_IRM_and_insect_f,'Pe');
%%


%% Insect_super_family Chi_sqared_test

%  Insect_original_super_family class
Insect_original_super_family = Insect(:,3);
Insect_original_super_family = Insect_original_super_family(ia1,:); 

%  Unique Insect super family  
[Insect_super_family_unique_class,i_s_f_ia,i_s_f_ic] = unique(Insect_original_super_family); 

%  Insect_super_family_class_cell
for i=1:size(Insect_super_family_unique_class,1)
    Insect_super_family_class_cell{i} = UniC1(i_s_f_ic == i,:); 
end;

% The joint relational metrix of Insect-Super_family and RIM
for i=1:size(Insect_IRM_class_cell,2) 
    for j=1:size(Insect_super_family_unique_class,1)
        R_IRM_and_insect_s_f(i,j) = numel(intersect(Insect_IRM_class_cell{i},Insect_super_family_class_cell{j})); 
    end; 
end; 

% Chi-squared test for Inscet_super_family_IRM

P_inscet_super_family_IRM = chi2Tests(R_IRM_and_insect_s_f,'Pe');
%%



%% Insect_order Chi_sqared_test

%  Insect_original_order class
Insect_original_order = Insect(:,1);
Insect_original_order = Insect_original_order(ia1,:); 

%  Unique Insect order 
[Insect_order_unique_class,i_s_o_ia,i_s_o_ic] = unique(Insect_original_order); 

%  Insect_order_cell
for i=1:size(Insect_order_unique_class,1)
    Insect_order_class_cell{i} = UniC1(i_s_o_ic == i,:); 
end;

% The joint relational metrix of Insect_order and RIM
for i=1:size(Insect_IRM_class_cell,2) 
    for j=1:size(Insect_order_unique_class,1)
        R_IRM_and_insect_o(i,j) = numel(intersect(Insect_IRM_class_cell{i},Insect_order_class_cell{j})); 
    end; 
end; 

% Chi-squared test for Inscet_order_IRM

P_inscet_order_IRM = chi2Tests(R_IRM_and_insect_o,'Pe');
%%


%% Plant_family Chi_sqared_test

%  Plant_original_family class
Plant_original_family = Plant(:,4);
Plant_original_family = Plant_original_family(ia2,:); 

%  Unique Plant family 
[Plant_family_unique_class,p_f_ia,p_f_ic] = unique(Plant_original_family); 

%  Plant_family_class_cell
for i=1:size(Plant_family_unique_class,1)
    Plant_family_class_cell{i} = UniC2(p_f_ic == i,:); 
end;

% The joint relational metrix of Plant_family and RIM
for i=1:size(Plant_IRM_class_cell,2) 
    for j=1:size(Plant_family_unique_class,1)
        R_IRM_and_Plant_f(i,j) = numel(intersect(Plant_IRM_class_cell{i},Plant_family_class_cell{j})); 
    end; 
end; 

% Chi-squared test for Plant_family_IRM
P_plant_family_IRM = chi2Tests(R_IRM_and_Plant_f,'Pe');
%%


%%  Plant_original_order class
Plant_original_order = Plant(:,3);
Plant_original_order = Plant_original_order(ia2,:); 

%  Unique Plant order 
[Plant_order_unique_class,p_o_ia,p_o_ic] = unique(Plant_original_order); 

%  Plant_order_class_cell
for i=1:size(Plant_order_unique_class,1)
    Plant_order_class_cell{i} = UniC2(p_o_ic == i,:); 
end;

% The joint relational metrix of Plant_order and RIM
for i=1:size(Plant_IRM_class_cell,2) 
    for j=1:size(Plant_order_unique_class,1)
        R_IRM_and_Plant_o(i,j) = numel(intersect(Plant_IRM_class_cell{i},Plant_order_class_cell{j})); 
    end; 
end; 

% Chi-squared test for Plant_order_IRM
P_plant_order_IRM = chi2Tests(R_IRM_and_Plant_o,'Pe');

%%


%% Plant_clade2 Chi_sqared_test

%  Plant_original_clade2 class
Plant_original_clade2 = Plant(:,2);
Plant_original_clade2 = Plant_original_clade2(ia2,:); 

%  Unique Plant clade2 
[Plant_clade2_unique_class,p_c2_ia,p_c2_ic] = unique(Plant_original_clade2); 

%  Plant_clade2_class_cell
for i=1:size(Plant_clade2_unique_class,1)
    Plant_clade2_class_cell{i} = UniC2(p_c2_ic == i,:); 
end;

% The joint relational metrix of Plant_clade2 and IRM
for i=1:size(Plant_IRM_class_cell,2) 
    for j=1:size(Plant_clade2_unique_class,1)
        R_IRM_and_Plant_c2(i,j) = numel(intersect(Plant_IRM_class_cell{i},Plant_clade2_class_cell{j})); 
    end; 
end; 

% Chi-squared test for Plant_clade2_IRM
P_plant_clade2_IRM = chi2Tests(R_IRM_and_Plant_c2,'Pe');

%%

%%
%  Plant_original_clade1 class
Plant_original_clade1 = Plant(:,1);
Plant_original_clade1 = Plant_original_clade1(ia2,:); 

%  Unique Plant clade1 
[Plant_clade1_unique_class,p_c1_ia,p_c1_ic] = unique(Plant_original_clade1); 

%  Plant_clade1_class_cell
for i=1:size(Plant_clade1_unique_class,1)
    Plant_clade1_class_cell{i} = UniC2(p_c1_ic == i,:); 
end;

% The joint relational metrix of Plant_clade1 and IRM
for i=1:size(Plant_IRM_class_cell,2) 
    for j=1:size(Plant_clade1_unique_class,1)
        R_IRM_and_Plant_c1(i,j) = numel(intersect(Plant_IRM_class_cell{i},Plant_clade1_class_cell{j})); 
    end; 
end; 

% Chi-squared test for Plant_clade1_IRM
P_plant_clade1_IRM = chi2Tests(R_IRM_and_Plant_c1,'Pe');


%%
function P = chi2Tests(X, method)
% Chi-square tests of homogeneity and independence
%   Computes the P-value for IxJ-table row/col independence.
%   Program by Steinar Thorvaldsen, steinar.thorvaldsen@uit.no, Dec. 2004. 
%   Ref.: DeltaProt toolbox at http://services.cbu.uib.no/software/deltaprot/
%   Last changes 22. Dec 2010.
%   Requires Matlab 7.1 or newer, and Matlab Statistics toolbox.
%   Input:
%        X:      data matrix (IxJ-table) of the observed frequency cells
%        method: 'RC': Read-Cressie power divergence statistics (default), lambda= 2/3 
%                'Pe': Standard Pearson chi2-distance, lambda= 1
%                'LL': Log Likelihood ratio distance, lambda= 0   
%   Output:
%        P-value
%
%   Use: P = chi2Tests(Observed,'RC') 

%	The P-value is computed through approximation with chi-2 distribution
%	under the null hypothesis for all methods.
%   The 'RC'-method is sligtly better than the 'Pe'-method in small tables 
%   with unbalanced column margins

%   Please, use the following reference:
%   Thorvaldsen, S. , Fl�, T. and Willassen, N.P. (2010) DeltaProt: a software toolbox 
%       for comparative genomics. BMC Bioinformatics 2010, Vol 11:573.
%       See http://www.biomedcentral.com/1471-2105/11/573

%   Other reference:
%   Rudas, T. (1986): A Monte Carlo Comparision of Small Sample Behaviour
%       of The Pearson, the Likelihood Ratio and the Cressie-Read Statistics.
%       J.Statist. Comput. Simul, vol 24, pp 107-120.
%   Read, TRC and Cressie, NAC (1988): Goodness of Fit Statistics for
%       Discrete Multivariate Data. Springer Verlag.
%   Ewens, WJ and Grant, GR (2001): Statistical Methods in Bioinformatics.
%       Springer Verlag.

if nargin < 2
    method = 'RC'; %default value
else
    method=upper(method);
    switch method %validate the method parameter:
    case {'RC' 'PE' 'LL'}
        % these are ok
    otherwise
        error('Chi2Test:UnknownMethod', ...
          'The ''method'' parameter value must be ''RC'', ''Pe'', or ''LL''.');
    end %switch
end %if

% lambda - power for statistic:
if method == 'RC'
    lambda=2/3;
elseif method == 'PE'
    lambda=1;
elseif method == 'LL'
    lambda=0;
end

if any(any(X < 0))
   error('chi2Test expects counts that are nonnegative values');
   return;
end

[I J] = size(X);
if I < 2 | J < 2,
    error('Matrix of observation must at least be of size 2x2');
    return;
end

rs = sum(X')';
cs = sum(X);
Ns = sum(rs);
if Ns <= 0 % no counts found in table
    P=NaN;
    return;
end
Eo = (rs*cs)./Ns;  %matrix of null means

% Make sure expected values are not too small
Emin=1.5;
if any(any(Eo < Emin))
    %disp ('Note: Some expected values in chi2Test are small (<1.5)...');
end

if any(any(Eo <= 0))
   t=NaN;
elseif lambda == 0 % logL
   i = find(X>0);
   t = 2 * sum(X(i).*(log(X(i)./Eo(i)))); % the test statistic
elseif lambda == 1 % Pearson chi2
   t = sum(sum((Eo-X).^2 ./ Eo)); % Pearson test statistic
else % Read-Cressie test statistic
   t = (2/(lambda*(lambda+1)))*sum(sum(X.*((X./Eo).^lambda -1))); 
end

df=prod(size(X)-[1,1]); % degree of freedom
P = 1-chi2cdf(t,df);
