%% ====================================================
%  Load Data Relative to Kinematics
%% ====================================================
% Load in data file for  arm raising task
allData = readtable('AllExtractedMetrics2.xlsx');
allData.TIBIALIS_ANTERIOR = str2double(allData.TIBIALIS_ANTERIOR);

muscleVars = {'BICEPS','TRICEPS','RECTUS_ABDOMINUS','OBLIQUES', ...
              'LATISSIMUS','ERECTOR','TIBIALIS_ANTERIOR','GASTROCNEMIUS'};

% Convert EMG data to milliseconds (multiply by 1000)
for m = 1:length(muscleVars)
    muscle = muscleVars{m};
    allData.(muscle) = allData.(muscle) * 1000;
end

% Ensure proper categorical coding
allData.GROUP = categorical(allData.GROUP);
allData.SUBJ_NUM = categorical(allData.SUBJ_NUM);
allData = allData(:,[1:4,13:20]);

% Convert to long format
longData = stack(allData, muscleVars, ...
    'NewDataVariableName','Timing', ...
    'IndexVariableName','MUSCLE');

longData.MUSCLE = categorical(longData.MUSCLE);
longData.GROUP = categorical(longData.GROUP);
longData.SUBJ_NUM = categorical(longData.SUBJ_NUM);
longData = longData(~isnan(longData.Timing), :);

%% Model 1: Fixed Effects only
% Timing ~ GROUP*MUSCLE
M1 = fitlme(longData, 'Timing ~ GROUP*MUSCLE', 'FitMethod', 'ML');

%% Model 2: Add random intercept for SUBJECT
% Timing ~ GROUP*MUSCLE + (1|SUBJ)
M2 = fitlme(longData, 'Timing ~ GROUP*MUSCLE + (1|SUBJ_NUM)', 'FitMethod', 'ML');

% Compare M1 and M2
compare(M1, M2)

%% Model 2: Add random intercept for SUBJECT
% Timing ~ GROUP*MUSCLE + (1|SUBJ)
M5 = fitlme(longData, 'Timing ~ GROUP + GROUP:MUSCLE + (1|SUBJ_NUM)', 'FitMethod', 'ML');

% Compare M1 and M2
compare(M5,M1)
%% Model 3: Allow muscle to vary within subject
% Timing ~ GROUP*MUSCLE + (MUSCLE|SUBJ)
M3 = fitlme(longData, ...
    'Timing ~ GROUP*MUSCLE + (MUSCLE | SUBJ_NUM)', ...
    'FitMethod','REML');

% compare(M1, M3)

%% Using M3 for Analysis
[beta,~,stats] = fixedEffects(M3,'DFMethod','satterthwaite');
coefNames = M3.Coefficients.Name;  % Names of coefficients

muscles = categories(longData.MUSCLE);
groups = categories(longData.GROUP);

fprintf('Group Differences within Each Muscle (Satterthwaite-corrected DF)\n');
fprintf('---------------------------------------------------------------\n');

for i = 1:length(muscles)
    
    % Build contrast vector L
    L = zeros(1,length(coefNames));
    
    % Always include GROUP main effect
    idx_group = strcmp(coefNames,'GROUP_NT');
    L(idx_group) = 1; % Reference is first level of muscle
    
    % Add interaction term if not reference muscle
    interactionName = ['GROUP_NT:MUSCLE_' char(muscles{i})];
    idx_inter = strcmp(coefNames,interactionName);
    if any(idx_inter)
        L(idx_inter) = 1;
    end
    
    % Contrast estimate (linear combo of coeffs -- estimated group diff for
    % muscle)
    contrastEstimate = L * beta;
    
    % Standard error of contrast
    Covb = M3.CoefficientCovariance;
    contrastSE = sqrt(L * Covb * L'); % To get standard error. Accounts for correlation btwn coeffs (group and group x muscle)
    
    % t-statistic
    tVal = contrastEstimate / contrastSE;
    
    % DF from Satterthwaite
    df = stats.DF;  % satterthwaite-corrected DF
    idx_nonzero = find(L~=0);
    df_contrast = min(stats.DF(idx_nonzero)); % Of coeff used, use smallest df (more conservative)
    
    % Two-tailed p-value
    pVal = 2*(1 - tcdf(abs(tVal), df_contrast));
    
    fprintf('%-20s  b = %.3f, SE = %.3f, t(%.1f) = %.3f, p = %.4f\n', ...
    muscles{i}, contrastEstimate, contrastSE, df_contrast, tVal, pVal);
    
end

%% Holm Corrected version
[beta,~,stats] = fixedEffects(M3,'DFMethod','satterthwaite');
coefNames = M3.Coefficients.Name;

muscles = categories(longData.MUSCLE);

fprintf('Group Differences within Each Muscle (Satterthwaite-corrected DF)\n');
fprintf('---------------------------------------------------------------\n');

% Preallocate
nMuscles = length(muscles);
contrastEstimate = zeros(nMuscles,1);
contrastSE = zeros(nMuscles,1);
tVal = zeros(nMuscles,1);
df_contrast = zeros(nMuscles,1);
pVal = zeros(nMuscles,1);

Covb = M3.CoefficientCovariance;

for i = 1:nMuscles

    % Build contrast vector
    L = zeros(1,length(coefNames));

    % Main effect of group
    idx_group = strcmp(coefNames,'GROUP_NT');
    L(idx_group) = 1;

    % Interaction term
    interactionName = ['GROUP_NT:MUSCLE_' char(muscles{i})];
    idx_inter = strcmp(coefNames,interactionName);
    if any(idx_inter)
        L(idx_inter) = 1;
    end

    % Contrast estimate
    contrastEstimate(i) = L * beta;

    % Standard error
    contrastSE(i) = sqrt(L * Covb * L');

    % t statistic
    tVal(i) = contrastEstimate(i) / contrastSE(i);

    % Satterthwaite DF
    idx_nonzero = find(L~=0);
    df_contrast(i) = min(stats.DF(idx_nonzero));

    % Raw p-value
    pVal(i) = 2*(1 - tcdf(abs(tVal(i)), df_contrast(i)));

end

% Holm correction
[pSorted,sortIdx] = sort(pVal);
m = length(pVal);

pHolmSorted = (m-(1:m)+1)'.*pSorted;

% Enforce monotonicity
for k = 2:m
    pHolmSorted(k) = max(pHolmSorted(k), pHolmSorted(k-1));
end

pHolmSorted = min(pHolmSorted,1);

% Return to original order
pHolm = zeros(size(pVal));
pHolm(sortIdx) = pHolmSorted;

% Display results
for i = 1:length(muscles)
    fprintf('%-20s  b = %.3f, SE = %.3f, t(%.1f) = %.3f, p = %.4f, Holm p = %.4f\n',...
        muscles{i}, contrastEstimate(i), contrastSE(i), ...
        df_contrast(i), tVal(i), pVal(i), pHolm(i));
end

%%
M_reduced = fitlme(longData, ...
    'Timing ~ MUSCLE + (MUSCLE | SUBJ_NUM)', ...
    'FitMethod','ML');
compare(M_reduced, M3)


%%
biceps = allData.BICEPS;
triceps = allData.TRICEPS - biceps;
rectusAbd = allData.RECTUS_ABDOMINUS - biceps;
obliques = allData.OBLIQUES-biceps;
latissimus = allData.LATISSIMUS-biceps;
erector = allData.ERECTOR-biceps;
tibialisAnt = allData.TIBIALIS_ANTERIOR-biceps;
gastrocnemius = allData.GASTROCNEMIUS-biceps;

allDataToBiceps = table(allData.SUBJ_NUM, allData.SUBJ_ID, allData.GROUP, allData.TRIAL, ...
    [triceps], [rectusAbd], [obliques], ...
    [latissimus], [erector], [tibialisAnt], ...
    [gastrocnemius], ...
    'VariableNames', {'SUBJ_NUM','SUBJ_ID','GROUP','TRIAL','TRICEPS','RECTUS_ABDOMINUS','OBLIQUES', ...
              'LATISSIMUS','ERECTOR','TIBIALIS_ANTERIOR','GASTROCNEMIUS'});

% Ensure proper categorical coding
allDataToBiceps.GROUP = categorical(allDataToBiceps.GROUP);
allDataToBiceps.SUBJ_NUM = categorical(allDataToBiceps.SUBJ_NUM);

muscleVarsToBiceps = {'TRICEPS','RECTUS_ABDOMINUS','OBLIQUES', ...
              'LATISSIMUS','ERECTOR','TIBIALIS_ANTERIOR','GASTROCNEMIUS'};


% Convert to long format
longDataToBiceps = stack(allDataToBiceps, muscleVarsToBiceps, ...
    'NewDataVariableName','Timing', ...
    'IndexVariableName','MUSCLE');

longDataToBiceps.MUSCLE = categorical(longDataToBiceps.MUSCLE);
longDataToBiceps.GROUP = categorical(longDataToBiceps.GROUP);
longDataToBiceps.SUBJ_NUM = categorical(longDataToBiceps.SUBJ_NUM);
longDataToBiceps = longDataToBiceps(~isnan(longDataToBiceps.Timing), :);

% Timing ~ GROUP*MUSCLE + (MUSCLE|SUBJ)
M4 = fitlme(longDataToBiceps, ...
    'Timing ~ GROUP*MUSCLE + (MUSCLE | SUBJ_NUM)', ...
    'FitMethod','ML');

[beta,~,stats] = fixedEffects(M4,'DFMethod','satterthwaite');
coefNames = M4.Coefficients.Name;  % Names of coefficients

muscles = categories(longDataToBiceps.MUSCLE);
groups = categories(longDataToBiceps.GROUP);

fprintf('Group Differences within Each Muscle (Satterthwaite-corrected DF)\n');
fprintf('---------------------------------------------------------------\n');

for i = 1:length(muscles)
    
    % Build contrast vector L
    L = zeros(1,length(coefNames));
    
    % Always include GROUP main effect
    idx_group = strcmp(coefNames,'GROUP_NT'); % adjust reference level if needed
    L(idx_group) = 1;
    
    % Add interaction term if not reference muscle
    interactionName = ['GROUP_NT:MUSCLE_' char(muscles{i})];
    idx_inter = strcmp(coefNames,interactionName);
    if any(idx_inter)
        L(idx_inter) = 1;
    end
    
    % Contrast estimate
    contrastEstimate = L * beta;
    
    % Standard error of contrast
    Covb = M4.CoefficientCovariance;
    contrastSE = sqrt(L * Covb * L');
    
    % t-statistic
    tVal = contrastEstimate / contrastSE;
    
    % DF from Satterthwaite
    df = stats.DF;  % vector of DF per coefficient, but using weighted approximation works
    % For simple contrasts, approximate with minimum DF among non-zero L entries
    idx_nonzero = find(L~=0);
    df_contrast = min(stats.DF(idx_nonzero));
    
    % Two-tailed p-value
    pVal = 2*(1 - tcdf(abs(tVal), df_contrast));
    
    fprintf('%-20s  b = %.3f, SE = %.3f, t(%.1f) = %.3f, p = %.4f\n', ...
    muscles{i}, contrastEstimate, contrastSE, df_contrast, tVal, pVal);
    
end
%%
%% Holm Corrected version
[beta,~,stats] = fixedEffects(M4,'DFMethod','satterthwaite');
coefNames = M4.Coefficients.Name;

muscles = categories(longDataToBiceps.MUSCLE);

fprintf('Group Differences within Each Muscle (Satterthwaite-corrected DF)\n');
fprintf('---------------------------------------------------------------\n');

% Preallocate
nMuscles = length(muscles);
contrastEstimate = zeros(nMuscles,1);
contrastSE = zeros(nMuscles,1);
tVal = zeros(nMuscles,1);
df_contrast = zeros(nMuscles,1);
pVal = zeros(nMuscles,1);

Covb = M4.CoefficientCovariance;

for i = 1:nMuscles

    % Build contrast vector
    L = zeros(1,length(coefNames));

    % Main effect of group
    idx_group = strcmp(coefNames,'GROUP_NT');
    L(idx_group) = 1;

    % Interaction term
    interactionName = ['GROUP_NT:MUSCLE_' char(muscles{i})];
    idx_inter = strcmp(coefNames,interactionName);
    if any(idx_inter)
        L(idx_inter) = 1;
    end

    % Contrast estimate
    contrastEstimate(i) = L * beta;

    % Standard error
    contrastSE(i) = sqrt(L * Covb * L');

    % t statistic
    tVal(i) = contrastEstimate(i) / contrastSE(i);

    % Satterthwaite DF
    idx_nonzero = find(L~=0);
    df_contrast(i) = min(stats.DF(idx_nonzero));

    % Raw p-value
    pVal(i) = 2*(1 - tcdf(abs(tVal(i)), df_contrast(i)));

end

% Holm correction
[pSorted,sortIdx] = sort(pVal);
m = length(pVal);

pHolmSorted = (m-(1:m)+1)'.*pSorted;

% Enforce monotonicity
for k = 2:m
    pHolmSorted(k) = max(pHolmSorted(k), pHolmSorted(k-1));
end

pHolmSorted = min(pHolmSorted,1);

% Return to original order
pHolm = zeros(size(pVal));
pHolm(sortIdx) = pHolmSorted;

% Display results
for i = 1:length(muscles)
    fprintf('%-20s  b = %.3f, SE = %.3f, t(%.1f) = %.3f, p = %.4f, Holm p = %.4f\n',...
        muscles{i}, contrastEstimate(i), contrastSE(i), ...
        df_contrast(i), tVal(i), pVal(i), pHolm(i));
end


M_reduced = fitlme(longDataToBiceps, ...
    'Timing ~ MUSCLE + (MUSCLE | SUBJ_NUM)', ...
    'FitMethod','ML');
compare(M_reduced, M4)