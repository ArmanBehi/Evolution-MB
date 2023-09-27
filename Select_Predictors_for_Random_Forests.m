%% Load and Preprocess Data
%a model that predicts the fuel economy of a car given its number of cylinders, engine displacement, horsepower, ...
% weight, acceleration, model year, and country of origin. 
load carbig
Cylinders = categorical(Cylinders);
Model_Year = categorical(Model_Year);
Origin = categorical(cellstr(Origin));
X = table(Cylinders,Displacement,Horsepower,Weight,Acceleration,Model_Year,Origin);

%% Determine Levels in Predictors
%For each predictor, determine the number of levels in the data
countLevels = @(x)numel(categories(categorical(x)));
numLevels = varfun(countLevels,X,'OutputFormat','uniform');

%Plot the levels
figure(1)
bar(numLevels)
title('Number of Levels Among Predictors')
xlabel('Predictor variable')
ylabel('Number of levels')
h = gca;
h.XTickLabel = X.Properties.VariableNames(1:end-1);
h.XTickLabelRotation = 45;
h.TickLabelInterpreter = 'none';

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
Mdl = fitrensemble(X,MPG,'Method','Bag','NumLearningCycles',500, ...
    'Learners',t);

%Estimate the model R2 using out-of-bag predictions:
yHat = oobPredict(Mdl);
R2 = corr(Mdl.Y,yHat)^2; %Mdl explains 87% of the variability around the mean.

%% Predictor Importance Estimation
%Estimate predictor importance values by permuting out-of-bag observations
%among the trees:
impOOB = oobPermutedPredictorImportance(Mdl);  %The estimates are not biased toward predictors containing many levels.

%Plot the Compare the predictor importance estimates.
figure(2)
bar(impOOB)
title('Unbiased Predictor Importance Estimates')
xlabel('Predictor variable')
ylabel('Importance')
h = gca;
h.XTickLabel = Mdl.PredictorNames;
h.XTickLabelRotation = 45;
h.TickLabelInterpreter = 'none';

%The Predictive Measure of Association is a value that indicates the similarity...
% between decision rules that split observations. Larger values indicate more highly...
% correlated pairs of predictors.
[impGain,predAssociation] = predictorImportance(Mdl);

figure(3)
imagesc(predAssociation)
title('Predictor Association Estimates')
colorbar
h = gca;
h.XTickLabel = Mdl.PredictorNames;
h.XTickLabelRotation = 45;
h.TickLabelInterpreter = 'none';
h.YTickLabel = Mdl.PredictorNames;




