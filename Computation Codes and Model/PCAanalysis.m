% ==========================
% PCA Analysis of 12D Vector
% ==========================

% This script first compiles the data and then runs a pca analysis on it.
% Feature vector (12D):
    % (1) Magnitude of peak velocity (mm/sec)
    % (2) Movement Time (s)
    % (3) Entropy during prep phase
    % (4) Entropy during acceleration
    % (5) Entropy deuring deceleration
    % (6) Triceps onset relative to biceps
    % (7) Rectus abdominus onset relative to biceps
    % (8) Obliques onset relative to biceps
    % (9) Latissimus onset relative to biceps
    % (10) Erector onset relative to biceps
    % (11) Tibialis anterior onset relative to biceps
    % (12) Gastrocnemius onset relative to biceps

% ====================
% Load and split data
% ====================
allData = readtable('AllExtractedMetrics2');
allData.TIBIALIS_ANTERIOR = str2double(allData.TIBIALIS_ANTERIOR);

% Split by group
allNTData = allData(strcmp(allData.GROUP, 'NT'),:);
allASDData = allData(strcmp(allData.GROUP, 'ASD'),:);
subjNTid = unique(allNTData.SUBJ_NUM);
subjASDid = unique(allASDData.SUBJ_NUM);

% ======================
% Subject-level averages
% ======================

% NT Averages
[peakVelNT] = subjectAverageCalculator([allNTData.SUBJ_NUM allNTData.PEAK_MAG])/1000; % m/s
[mtNT] = subjectAverageCalculator([allNTData.SUBJ_NUM allNTData.MT])*1000; % ms
[entropyPrepNT] = subjectAverageCalculator([allNTData.SUBJ_NUM allNTData.PREP_ENTROPY]);
[entropyAccelNT] = subjectAverageCalculator([allNTData.SUBJ_NUM allNTData.LIFT1_ENTROPY]);
[entropyDecelNT] = subjectAverageCalculator([allNTData.SUBJ_NUM allNTData.LIFT2_ENTROPY]);
[tricepsNT] = subjectAverageCalculator([allNTData.SUBJ_NUM (allNTData.TRICEPS - allNTData.BICEPS)]);
[rectusAbdNT] = subjectAverageCalculator([allNTData.SUBJ_NUM (allNTData.RECTUS_ABDOMINUS - allNTData.BICEPS)]);
[obliquesNT] = subjectAverageCalculator([allNTData.SUBJ_NUM (allNTData.OBLIQUES - allNTData.BICEPS)]);
[latissimusNT] = subjectAverageCalculator([allNTData.SUBJ_NUM (allNTData.LATISSIMUS - allNTData.BICEPS)]);
[erectorNT] = subjectAverageCalculator([allNTData.SUBJ_NUM (allNTData.ERECTOR - allNTData.BICEPS)]);
[tibialisAntNT] = subjectAverageCalculator([allNTData.SUBJ_NUM (allNTData.TIBIALIS_ANTERIOR - allNTData.BICEPS)]);
[gastrocnemiusNT] = subjectAverageCalculator([allNTData.SUBJ_NUM (allNTData.GASTROCNEMIUS - allNTData.BICEPS)]);

% ASD Averages
[peakVelASD] = subjectAverageCalculator([allASDData.SUBJ_NUM allASDData.PEAK_MAG])/1000; % m/s
[mtASD] = subjectAverageCalculator([allASDData.SUBJ_NUM allASDData.MT])*1000; % ms
[entropyPrepASD] = subjectAverageCalculator([allASDData.SUBJ_NUM allASDData.PREP_ENTROPY]);
[entropyAccelASD] = subjectAverageCalculator([allASDData.SUBJ_NUM allASDData.LIFT1_ENTROPY]);
[entropyDecelASD] = subjectAverageCalculator([allASDData.SUBJ_NUM allASDData.LIFT2_ENTROPY]);
[tricepsASD] = subjectAverageCalculator([allASDData.SUBJ_NUM (allASDData.TRICEPS - allASDData.BICEPS)]);
[rectusAbdASD] = subjectAverageCalculator([allASDData.SUBJ_NUM (allASDData.RECTUS_ABDOMINUS - allASDData.BICEPS)]);
[obliquesASD] = subjectAverageCalculator([allASDData.SUBJ_NUM (allASDData.OBLIQUES - allASDData.BICEPS)]);
[latissimusASD] = subjectAverageCalculator([allASDData.SUBJ_NUM (allASDData.LATISSIMUS - allASDData.BICEPS)]);
[erectorASD] = subjectAverageCalculator([allASDData.SUBJ_NUM (allASDData.ERECTOR - allASDData.BICEPS)]);
[tibialisAntASD] = subjectAverageCalculator([allASDData.SUBJ_NUM (allASDData.TIBIALIS_ANTERIOR - allASDData.BICEPS)]);
[gastrocnemiusASD] = subjectAverageCalculator([allASDData.SUBJ_NUM (allASDData.GASTROCNEMIUS - allASDData.BICEPS)]);

% ======================
% Build Feature Matrices
% ======================

allNTaverageData = [peakVelNT mtNT entropyPrepNT entropyAccelNT entropyDecelNT...
    tricepsNT rectusAbdNT obliquesNT latissimusNT erectorNT tibialisAntNT gastrocnemiusNT];

allASDaverageData = [peakVelASD mtASD entropyPrepASD entropyAccelASD entropyDecelASD...
    tricepsASD rectusAbdASD obliquesASD latissimusASD erectorASD tibialisAntASD gastrocnemiusASD];

% Eliminate subjects with most EMG data missing
allaverages = [allNTaverageData([1:17,20:22],:); allASDaverageData([1:4,6:9,11,12],:)];
groupLabels = [1*ones(20,1); 2*ones(10,1)];

% Remove mean
data_centered = allaverages - mean(allaverages, 1, 'omitnan');

% Find min and max for each column 
data_min = min(data_centered, [], 1);
data_max = max(data_centered, [], 1);

% Normalize to [-1, 1]
data_normalized = 2 * (data_centered - data_min) ./ (data_max - data_min) - 1;

allaverages = data_normalized;
[nObs, nFeat] = size(allaverages);

% =============
% Perform PPCA
% ============= 
numPCs = 3;
tol = 1e-6;
maxIter = 1000;

if any(isnan(allaverages(:)))
    warning('PPCA may not handle missing values well. Consider imputation.');
end

[coeff_ppca, score_ppca, latent_ppca, ~, explained_ppca, mu_ppca] = ppca(allaverages, numPCs);
explained_ppca = 100 * latent_ppca / sum(latent_ppca); 
disp('Variance explained by PPCA components:'); 
disp(explained_ppca);

% ===============
% Plot Loadings^2
% =============== 
figure; hold on;
bar(coeff_ppca(:,1:2).^2);
title('PPCA Loadings^2 (first 3 PCs)');
xlabel('Features'); ylabel('Loading Value');
legend('PC1', 'PC2', 'PC3')

% =======================================================
% Compute and plot euclidean distance from global NT mean
% =======================================================

X = score_ppca(:, 1:2);
Y = groupLabels;

mu_eucPCA = mean(X(groupLabels == 1, :));
eucPCA = sqrt(sum((X - mu_eucPCA).^2, 2));

% Store distances by group
euc_group1PCA = eucPCA(Y == 1);
euc_group2PCA = eucPCA(Y == 2);

% Box plot of euclidean distance from global mean by group
figure; hold on;
BoxPlotGenerator(euc_group1PCA, euc_group2PCA, 'Euclidean Distance in PCA Space', 'Euclidean Distance')

% Scatter plot of subjects in PC space
figure; hold on;
for i = 1:20
    plot([score_ppca(i,1), mu_eucPCA(1)], [score_ppca(i,2), mu_eucPCA(2)],'--','color',[0.553, 0.612, 1]);
end
for i = 21:30
    plot([score_ppca(i,1), mu_eucPCA(1)], [score_ppca(i,2), mu_eucPCA(2)], '--', 'color', [1, 0.545, 0.627]);
end
scatter(score_ppca(1:20,1), score_ppca(1:20,2), 60, 'b', 'filled');
scatter(score_ppca(21:end,1), score_ppca(21:end,2), 60, 'r', 'filled');

scatter(mu_eucPCA(1), mu_eucPCA(2), 200, 'k', 'filled',"pentagram");
xlabel(['PC1 (' num2str(explained_ppca(1), '%.1f') '%)']);
ylabel(['PC2 (' num2str(explained_ppca(2), '%.1f') '%)']);
title('3D PCA of 13D Feature Space');
legend('NT','ASD');
grid on;

% Run statistics on euclidean distance
allEucData = [euc_group1PCA; euc_group2PCA];
subjects = length(euc_group1PCA) + length(euc_group2PCA);

eucData = table();
eucData.SUBJ_NUM = [1:1:30]';
eucData.GROUP = [allData.GROUP(1:length(euc_group1PCA)); allData.GROUP(300:length(euc_group2PCA)+299)];
eucData.EUC = [euc_group1PCA; euc_group2PCA];

eucData.GROUP = categorical(eucData.GROUP);

lme = fitlme(eucData, 'EUC ~ GROUP + (1|SUBJ_NUM)');

% Corrected values using Satterthwaite
[betasFixRaw,~,FixEst] = fixedEffects(lme,'DFMethod','satterthwaite');
disp(FixEst)

L = zeros(1, height(FixEst));
L(strcmp(FixEst.Name, 'GROUP_NT')) = 1;
[p, F, DF1, DF2] = coefTest(lme, L);

F_values = F;
p_values = p;
DF1_values = DF1;
DF2_values = DF2;

disp(p_values)