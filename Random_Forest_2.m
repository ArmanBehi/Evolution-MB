%% Preparing data for Feed
% Specify the name and path of the Excel file to be imported
filename = '0_30_selected_variable';
filepath = 'D:\New folder';

% Construct the full file path
full_path = fullfile(filepath, filename);

% Read the data from the Excel file into a MATLAB table
d0_30 = readtable(full_path);

d0_30 = Interested_attribute(d0_30);

Group_splitTables_0_30 = GroupwiseSplitBinnedData(d0_30);

for i=1:length(Group_splitTables_0_30)
    Group_splitTables_0_30{i,1} = convertToDoubleAndNaN(Group_splitTables_0_30{i,1});
end
% Assuming you have table1, table2, table3, table4, and table5 already defined
X = vertcat(Group_splitTables_0_30{1,1},Group_splitTables_0_30{2,1}...
    ,Group_splitTables_0_30{3,1},Group_splitTables_0_30{4,1});

% Replace this with your actual class labels
classLabels = [ones(size(Group_splitTables_0_30{1,1}, 1), 1);  % Assuming table1 has 'n' rows
               2*ones(size(Group_splitTables_0_30{2,1}, 1), 1);  % Assuming table2 has 'm' rows
               3*ones(size(Group_splitTables_0_30{3,1}, 1), 1);  % Assuming table3 has 'p' rows
               4*ones(size(Group_splitTables_0_30{4,1}, 1), 1);  % Assuming table4 has 'q' rows
              ]; % Assuming table5 has 'r' rows
%% Train Bagged Ensemble of Regression Trees
%'NumVariablesToSample','all' — Use all predictor variables at each node ...
% to ensure that each tree uses all predictor variables.


%'PredictorSelection','interaction-curvature' — Specify usage...
% of the interaction test to select split predictors.

%'Surrogate','on' — Specify usage of surrogate splits to increase...
% accuracy because the data set includes missing values.

t = templateTree('NumVariablesToSample','all',...
    'PredictorSelection','interaction-curvature','Surrogate','on');
rng(1); % For reproducibility
Mdl = fitrensemble(X,classLabels,'Method','Bag','NumLearningCycles',500, ...
    'Learners',t);

%Estimate the model R2 using out-of-bag predictions:
yHat = oobPredict(Mdl);
R2 = corr(Mdl.Y,yHat)^2; %Mdl explains 87% of the variability around the mean.
%% Predictor Importance Estimation
%Estimate predictor importance values by permuting out-of-bag observations
%among the trees:
impOOB = oobPermutedPredictorImportance(Mdl);  %The estimates are not biased toward predictors containing many levels.

% Sort the numeric values and get the sorting indices
[sortedValues, sortingIndices] = sort(impOOB,'descend');
VarNames = Group_splitTables_0_30{1,1}.Properties.VariableNames;

for i=1:numel(VarNames)
    VarNames_sorted(1,i) = VarNames(sortingIndices(i));
end

%Plot the Compare the predictor importance estimates.
n = 18;
figure(2)
barh(sortedValues(1:n))
title('Attribute Importances from Random Forest')
xlabel('Predictor Attributes')
ylabel('Attribute Importance Score')
h = gca;
h.XTick = 1:n;
h.XTickLabel = VarNames_sorted(1:n);
h.XTickLabelRotation = 45;
h.TickLabelInterpreter = 'none';
h.FontSize = 15;

%The Predictive Measure of Association is a value that indicates the similarity...
% between decision rules that split observations. Larger values indicate more highly...
% correlated pairs of predictors.
[impGain,predAssociation] = predictorImportance(Mdl);

figure(3)
n = 10;
imagesc(predAssociation(1:n,1:n))
title('Predictor Association Estimates')
colorbar
h = gca;
h.XTickLabel = Mdl.PredictorNames;
h.XTickLabelRotation = 45;
h.TickLabelInterpreter = 'none';
%h.YTickLabel = Mdl.PredictorNames;
h.XTickLabel = VarNames_sorted(1:n);
h.YTickLabel = VarNames_sorted(1:n);
h.XTick = 1:n;
h.YTick = 1:n;




