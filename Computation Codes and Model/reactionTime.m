
% Extra kinematic analyses

% Starting with Reaction Time
% load('UpdatedMovementOnsetInfo.mat');

% for each subject need to load in movement onset

allParticipants = [222	223	226	230	234	231 239	250	244 253	232	263	259	258	266	229	264	265	271	280	279	276	282	277	218	20	21 240	233	267	269	215	270	274 284	273 319 326];
Path = 'C:\Users\Action Lab\OneDrive - Northeastern University\Desktop\OneDrive - Northeastern University\Action Lab\Arm Raising\Armlifting\MATLAB\My Codes\QTM Files\';
Path = 'C:\Users\Action Lab\OneDrive - Northeastern University\Desktop\OneDrive - Northeastern University\Action Lab\Arm Raising\Armlifting\MATLAB\FinalFiles\SubmissionFiles\Data\';
% Path = 'C:\Users\xvssc\OneDrive - Northeastern University\Action Lab\Arm Raising\Armlifting\MATLAB\My Codes\QTM Files\';


for subjSeq = 1:length(allParticipants)

    if allParticipants(subjSeq) == 266 || allParticipants(subjSeq) == 229
        continue;
    end
    if subjSeq < 25
        load('initiationNT2.mat');
        movementOnset = initiationNT{subjSeq};
    else
        load('initiationASD2.mat');
        movementOnset = initiationASD{subjSeq};
    end
    subj = allParticipants(subjSeq);  
    fileName = sprintf('%d_Control_Lift',subj); 
    FILE = load([Path fileName '.mat']);
    QTM = getfield(FILE,sprintf('qtm_%s',fileName));
    fileName2 = sprintf('qtm_%d_Control_Lift',subj);
    
    load([Path fileName '_qtm.mat']);
    load([Path fileName '.mat']);
    
    force = FILE.(fileName2).Force(1).Force(3,:);

    if subj == 253
        force = FILE.(fileName2).Force(2).Force(3,:);
    end
    
    cueSignal = FILE.(fileName2).Analog.Data(2,:);

    % movementOnset = subject(subjSeq).movementOnset;
    
    

    % if subjSeq == 7
    %     movementOnset = movementOnset(3:end);
    % end
    movementOnsetInCueFrames = movementOnset*10;
    cueTimes = []; drumrollTimes = [];
    for i = 1:length(movementOnsetInCueFrames)
        trialCueSignal = cueSignal(movementOnsetInCueFrames(i)-3000:movementOnsetInCueFrames(i));
        [row, col] = find(trialCueSignal > 3);
        if ~isempty(col)
                cueTimes = [cueTimes col(1)+(movementOnsetInCueFrames(i)-3000)];
        else
            cueTimes = [cueTimes NaN];
        end
        [row, col] = find(trialCueSignal > 1.98 & trialCueSignal < 3);
        if ~isempty(col)
            drumrollTimes = [drumrollTimes col(1)+(movementOnsetInCueFrames(i)-3000)];
        else
            drumrollTimes = [drumrollTimes NaN];
        end
        % for(x = 1:length(trialCueSignal))
        %     if trialCueSignal(x) > 3 & trialCueSignal(x-1) < 0.5
        %         cueTimes = [cueTimes x];
        %     end
        % end
    end
    % if ~isempty(cueTimes)
    %     figure; plot(cueSignal); hold on; xline(movementOnsetInCueFrames, 'k'); hold on; xline(cueTimes, 'r');
    % end
    
    onsetData(subjSeq).movementOnsetSec = movementOnset/200;
    onsetData(subjSeq).cueOnsetSec = cueTimes/2000;
    onsetData(subjSeq).reactionTime = (movementOnset/200) - (cueTimes/2000);
    onsetData(subjSeq).drumrollTimeSec = drumrollTimes/2000;
    onsetData(subjSeq).startFromDrumPercent = ((movementOnset/200)-(drumrollTimes/2000))/((cueTimes/2000)-(drumrollTimes/2000));
    onsetData(subjSeq).startFromDrumSec = (movementOnset/200) - (drumrollTimes/2000);
    onsetData(subjSeq).drumrollDuration = cueTimes/2000 - drumrollTimes/2000;
end

%%
allMovementData = [];
for i = 1:38
    if i == 15 || i == 16
        continue;
    else
        subjectMovement = onsetData(i).movementOnsetSec;
        subjectCue = onsetData(i).cueOnsetSec;
        subjectID = ones(size(subjectCue))*i;
        subjectData = [subjectID' subjectMovement' subjectCue'];
        allMovementData = [allMovementData; subjectData];
    end
end
%%
% NT averages
NTrxns = [];
NTrxnsAll = [];
for i = 1:24
    if i == 15 || i == 16
        NTrxnsAll = NTrxnsAll;
    else
        subjrxn = onsetData(i).reactionTime;
        subjrxn_average = mean(subjrxn, 'omitnan');
        NTrxns = [NTrxns subjrxn_average];
        NTrxnsAll = [NTrxnsAll; onsetData(i).reactionTime'];
    end
end

ASDrxns = [];
ASDrxnsAll = [];
for i = 25:38
    subjrxn = onsetData(i).reactionTime;
    subjrxn_average = mean(subjrxn, 'omitnan');
    ASDrxns = [ASDrxns subjrxn_average];
    ASDrxnsAll = [ASDrxnsAll; onsetData(i).reactionTime'];
end

figure; hold on; tiledlayout(1,5);
nexttile; hold on; BoxPlotGenerator(NTrxns'*1000, ASDrxns'*1000, 'Reaction Times', 'Timing (ms)');

[h, p] = ttest2(NTrxns', ASDrxns')

%%
allData = readtable('AllExtractedMetrics2.xlsx');
allData.TIBIALIS_ANTERIOR = str2double(allData.TIBIALIS_ANTERIOR);


allNTData = allData(strcmp(allData.GROUP, 'NT'),:); % Parse out NT data
allASDData = allData(strcmp(allData.GROUP, 'ASD'),:); % Parse out ASD data

allsubjID = unique(allData.SUBJ_NUM);
subjNTid = unique(allNTData.SUBJ_NUM); % Parse out NT subject numbers
subjASDid = unique(allASDData.SUBJ_NUM); % Parse out ASD subject numbers

% Statistics for peak tangential arm velocity
GROUP = categorical(allData.GROUP);
SUBJ_NUM = allData.SUBJ_NUM;
RXN_TIME = [NTrxnsAll; ASDrxnsAll];

for i = 1:length(RXN_TIME)
    if RXN_TIME(i) < 0.15 || RXN_TIME(i) > 0.5
        RXN_TIME(i) = nan;
    end
end

reactionTable = table(GROUP, SUBJ_NUM, RXN_TIME);

lme_peakvel = fitlme(reactionTable, 'RXN_TIME ~ GROUP + (1|SUBJ_NUM)');

if isa(lme_peakvel, 'GeneralizedLinearMixedModel')
    [betasFixRaw,~,FixEst] = fixedEffects(lme_peakvel, 'DFMethod','residual');
else
    [betasFixRaw,~,FixEst] = fixedEffects(lme_peakvel,'DFMethod','satterthwaite');
end
disp(FixEst)

%% sort by group
NTrxnsAll = NTrxnsAll / -1;
ASDrxnsAll = ASDrxnsAll / -1;


%% Load in data file for  arm raising task
allData = readtable('AllExtractedMetrics2.xlsx');
allData.TIBIALIS_ANTERIOR = str2double(allData.TIBIALIS_ANTERIOR);


allNTData = allData(strcmp(allData.GROUP, 'NT'),:); % Parse out NT data
allASDData = allData(strcmp(allData.GROUP, 'ASD'),:); % Parse out ASD data

allsubjID = unique(allData.SUBJ_NUM);
subjNTid = unique(allNTData.SUBJ_NUM); % Parse out NT subject numbers
subjASDid = unique(allASDData.SUBJ_NUM); % Parse out ASD subject numbers

%%

tricepsNT = allNTData.TRICEPS - allNTData.BICEPS;
rectusAbdNT = allNTData.RECTUS_ABDOMINUS - allNTData.BICEPS;
obliquesNT = allNTData.OBLIQUES - allNTData.BICEPS;
latissimusNT = allNTData.LATISSIMUS - allNTData.BICEPS;
erectorNT = allNTData.ERECTOR - allNTData.BICEPS;
tibialisAntNT = allNTData.TIBIALIS_ANTERIOR - allNTData.BICEPS;
gastrocnemiusNT = allNTData.GASTROCNEMIUS - allNTData.BICEPS;

tricepsASD = allASDData.TRICEPS - allASDData.BICEPS;
rectusAbdASD = allASDData.RECTUS_ABDOMINUS - allASDData.BICEPS;
obliquesASD = allASDData.OBLIQUES - allASDData.BICEPS;
latissimusASD = allASDData.LATISSIMUS - allASDData.BICEPS;
erectorASD = allASDData.ERECTOR - allASDData.BICEPS;
tibialisAntASD = allASDData.TIBIALIS_ANTERIOR - allASDData.BICEPS;
gastrocnemiusASD = allASDData.GASTROCNEMIUS - allASDData.BICEPS;

%%
% NT Data
bicepsNT = allNTData.BICEPS;
tricepsNT = allNTData.TRICEPS;
rectusAbdNT = allNTData.RECTUS_ABDOMINUS;
obliquesNT = allNTData.OBLIQUES;
latissimusNT = allNTData.LATISSIMUS;
erectorNT = allNTData.ERECTOR;
tibialisAntNT = allNTData.TIBIALIS_ANTERIOR;
gastrocnemiusNT = allNTData.GASTROCNEMIUS;
% ASD Data
bicepsASD = allASDData.BICEPS;
tricepsASD = allASDData.TRICEPS;
rectusAbdASD = allASDData.RECTUS_ABDOMINUS;
obliquesASD = allASDData.OBLIQUES;
latissimusASD = allASDData.LATISSIMUS;
erectorASD = allASDData.ERECTOR;
tibialisAntASD = allASDData.TIBIALIS_ANTERIOR;
gastrocnemiusASD = allASDData.GASTROCNEMIUS;

NTmuscles = [bicepsNT tricepsNT rectusAbdNT obliquesNT latissimusNT erectorNT tibialisAntNT gastrocnemiusNT];
NTspread = [];
for i = 1:size(NTmuscles,1)
    firstMuscle = min(NTmuscles(i,:));
    lastMuscle = max(NTmuscles(i,:));
    spread = lastMuscle - firstMuscle;
    NTspread = [NTspread; spread];
end

ASDmuscles = [bicepsASD tricepsASD rectusAbdASD obliquesASD latissimusASD erectorASD tibialisAntASD gastrocnemiusASD];
ASDspread = [];
for i = 1:size(ASDmuscles,1)
    firstMuscle = min(ASDmuscles(i,:));
    lastMuscle = max(ASDmuscles(i,:));
    spread = lastMuscle - firstMuscle;
    ASDspread = [ASDspread; spread];
end
[NTspread_averages] = subjectAverageCalculator([allNTData.SUBJ_NUM NTspread])*1000;
[ASDspread_averages] = subjectAverageCalculator([allASDData.SUBJ_NUM ASDspread])*1000;

figure; hold on; BoxPlotGenerator(NTspread_averages, ASDspread_averages, 'Timing Spread', 'Timing (ms)');
[h, p] = ttest2(NTspread_averages', ASDspread_averages')



%% EMG Data relative to cue

% NT Data
bicepsNT = allNTData.BICEPS - NTrxnsAll;
tricepsNT = allNTData.TRICEPS - NTrxnsAll;
rectusAbdNT = allNTData.RECTUS_ABDOMINUS - NTrxnsAll;
obliquesNT = allNTData.OBLIQUES - NTrxnsAll;
latissimusNT = allNTData.LATISSIMUS - NTrxnsAll;
erectorNT = allNTData.ERECTOR - NTrxnsAll;
tibialisAntNT = allNTData.TIBIALIS_ANTERIOR - NTrxnsAll;
gastrocnemiusNT = allNTData.GASTROCNEMIUS - NTrxnsAll;
% ASD Data
bicepsASD = allASDData.BICEPS - ASDrxnsAll;
tricepsASD = allASDData.TRICEPS - ASDrxnsAll;
rectusAbdASD = allASDData.RECTUS_ABDOMINUS - ASDrxnsAll;
obliquesASD = allASDData.OBLIQUES - ASDrxnsAll;
latissimusASD = allASDData.LATISSIMUS - ASDrxnsAll;
erectorASD = allASDData.ERECTOR - ASDrxnsAll;
tibialisAntASD = allASDData.TIBIALIS_ANTERIOR - ASDrxnsAll;
gastrocnemiusASD = allASDData.GASTROCNEMIUS - ASDrxnsAll;

% Find average of activation onset across trials for each subject and each
% muscle (in ms)
[bicepsNT_averages] = subjectAverageCalculator([allNTData.SUBJ_NUM bicepsNT])*1000;
[tricepsNT_averages] = subjectAverageCalculator([allNTData.SUBJ_NUM tricepsNT])*1000;
[rectusAbdNT_averages] = subjectAverageCalculator([allNTData.SUBJ_NUM rectusAbdNT])*1000;
[obliquesNT_averages] = subjectAverageCalculator([allNTData.SUBJ_NUM obliquesNT])*1000;
[latissimusNT_averages] = subjectAverageCalculator([allNTData.SUBJ_NUM latissimusNT])*1000;
[erectorNT_averages] = subjectAverageCalculator([allNTData.SUBJ_NUM erectorNT])*1000;
[tibialisAntNT_averages] = subjectAverageCalculator([allNTData.SUBJ_NUM tibialisAntNT])*1000;
[gastrocnemiusNT_averages] = subjectAverageCalculator([allNTData.SUBJ_NUM gastrocnemiusNT])*1000;

[bicepsASD_averages] = subjectAverageCalculator([allASDData.SUBJ_NUM bicepsASD])*1000;
[tricepsASD_averages] = subjectAverageCalculator([allASDData.SUBJ_NUM tricepsASD])*1000;
[rectusAbdASD_averages] = subjectAverageCalculator([allASDData.SUBJ_NUM rectusAbdASD])*1000;
[obliquesASD_averages] = subjectAverageCalculator([allASDData.SUBJ_NUM obliquesASD])*1000;
[latissimusASD_averages] = subjectAverageCalculator([allASDData.SUBJ_NUM latissimusASD])*1000;
[erectorASD_averages] = subjectAverageCalculator([allASDData.SUBJ_NUM erectorASD])*1000;
[tibialisAntASD_averages] = subjectAverageCalculator([allASDData.SUBJ_NUM tibialisAntASD])*1000;
[gastrocnemiusASD_averages] = subjectAverageCalculator([allASDData.SUBJ_NUM gastrocnemiusASD])*1000;

% Statistics for muscle activation results
muscleLabels = {'Biceps', 'Triceps', 'Rectus Abdominus', 'Obliques', 'Latissimus', ...
    'Erector', 'Tibialis Anterior', 'Gastrocnemius'};
muscles = {'BICEPS', 'TRICEPS', 'RECTUS_ABDOMINUS', 'OBLIQUES', 'LATISSIMUS', ...
    'ERECTOR', 'TIBIALIS_ANTERIOR', 'GASTROCNEMIUS'};

BICEPS = [bicepsNT; bicepsASD];
TRICEPS = [tricepsNT; tricepsASD];
RECTUS_ABDOMINUS = [rectusAbdNT; rectusAbdASD];
OBLIQUES = [obliquesNT; obliquesASD];
LATISSIMUS = [latissimusNT; latissimusASD];
ERECTOR = [erectorNT; erectorASD];
TIBIALIS_ANTERIOR = [tibialisAntNT; tibialisAntASD];
GASTROCNEMIUS = [gastrocnemiusNT; gastrocnemiusASD];

updatedTable = table();
updatedTable.SUBJ_NUM = allData.SUBJ_NUM;
updatedTable.GROUP = allData.GROUP;
updatedTable.BICEPS = [bicepsNT; bicepsASD];
updatedTable.TRICEPS = [tricepsNT; tricepsASD];
updatedTable.RECTUS_ABDOMINUS = [rectusAbdNT; rectusAbdASD];
updatedTable.OBLIQUES = [obliquesNT; obliquesASD];
updatedTable.LATISSIMUS = [latissimusNT; latissimusASD];
updatedTable.ERECTOR = [erectorNT; erectorASD];
updatedTable.TIBIALIS_ANTERIOR = [tibialisAntNT; tibialisAntASD];
updatedTable.GASTROCNEMIUS = [gastrocnemiusNT; gastrocnemiusASD];

muscleData = table();
for i = 1:length(muscles)
    temp = updatedTable(:, {'SUBJ_NUM', 'GROUP'});
    temp.MUSCLE = repmat(muscleLabels(i), height(updatedTable), 1);
    temp.TIMING = updatedTable.(muscles{i});
    muscleData = [muscleData; temp];
end

% Make sure categorical variables are properly defined
muscleData.GROUP = categorical(muscleData.GROUP);
muscleData.MUSCLE = categorical(muscleData.MUSCLE);

lme_muscle = fitlme(muscleData, 'TIMING ~ GROUP*MUSCLE + (1|SUBJ_NUM)');

% Corrected values using Satterthwaite
if isa(lme_muscle,'GeneralizedLinearMixedModel')
    [betasFixRaw,~,FixEst] = fixedEffects(lme_muscle,'DFMethod','residual');
else
    [betasFixRaw,~,FixEst] = fixedEffects(lme_muscle,'DFMethod','satterthwaite');
end
disp(FixEst)

coefNames = lme_muscle.CoefficientNames;
nCoef = numel(coefNames);

% muscles to test (include baseline Erector explicitly)
musclesTest = {'Biceps','Erector','Gastrocnemius','Latissimus','Obliques', ...
               'Rectus Abdominus','Tibialis Anterior','Triceps'};

pvals = nan(numel(musclesTest),1);

for m = 1:numel(musclesTest)
    thisMuscle = musclesTest{m};
    
    % start with zeros
    H = zeros(1,nCoef);
    
    % group main effect
    idxGroup = strcmp(coefNames,'GROUP_NT');
    H(idxGroup) = 1;
    
    % add interaction if this muscle ≠ baseline (Erector)
    intName = ['GROUP_NT:MUSCLE_' thisMuscle];
    if any(strcmp(coefNames,intName))
        idxInt = strcmp(coefNames,intName);
        H(idxInt) = 1;
    end
    
    % test group difference at this muscle
    pvals(m) = coefTest(lme_muscle,H);
end

% Bonferroni correction
m = numel(pvals);
pvals_bonf = min(pvals * m, 1);

% put into table
resultsTable = table(musclesTest',pvals_bonf, ...
    'VariableNames',{'Muscle','GroupDiff_pValue'});
disp(resultsTable)


%%

% Generate boxplots for muscle activation
boxplot_data = [bicepsNT_averages; bicepsASD_averages; ...
                tricepsNT_averages; tricepsASD_averages; ...
                rectusAbdNT_averages; rectusAbdASD_averages; ...
                obliquesNT_averages; obliquesASD_averages; ...
                latissimusNT_averages; latissimusASD_averages; ...
                erectorNT_averages; erectorASD_averages; ...
                tibialisAntNT_averages; tibialisAntASD_averages; ...
                gastrocnemiusNT_averages; gastrocnemiusASD_averages];

NTsubjSize = length(tricepsNT_averages);
ASDsubjSize = length(tricepsASD_averages);

group_labels = [ones(NTsubjSize,1); 2*ones(ASDsubjSize,1); ...
                3*ones(NTsubjSize,1); 4*ones(ASDsubjSize,1); ...
                5*ones(NTsubjSize,1); 6*ones(ASDsubjSize,1); ...
                7*ones(NTsubjSize,1); 8*ones(ASDsubjSize,1); ...
                9*ones(NTsubjSize,1); 10*ones(ASDsubjSize,1); ...
                11*ones(NTsubjSize,1); 12*ones(ASDsubjSize,1); ...
                13*ones(NTsubjSize,1); 14*ones(ASDsubjSize,1); ...
                15*ones(NTsubjSize,1); 16*ones(ASDsubjSize,1)];

figure; hold on;
nColors = 8;
blue = [0, 0, 1];       green = [0, 0.6, 0.3]; 
red = [1, 0, 0];        softOrange = [1, 0.6, 0.2];
colorsNT = zeros(nColors, 3);
colorsASD = zeros(nColors, 3);
for i = 1:3
    colorsNT(:, i) = linspace(blue(i), green(i), nColors);
    colorsASD(:, i) = linspace(red(i), softOrange(i), nColors);
end

% Highlight significant muscles
x_coords = [-300 600 600 -300]; y_coords = [12.3 12.3 11.2 11.2]; % Highlight biceps (significant muscle)
fill(x_coords, y_coords,[0.5, 0.5, 0.5],'FaceAlpha', 0.25, 'EdgeColor', 'none');

x_coords = [-300 600 600 -300]; y_coords = [6.3 6.3 5.2 5.2]; % Highlight Latissimus (significant muscle)
fill(x_coords, y_coords,[0.5, 0.5, 0.5],'FaceAlpha', 0.25, 'EdgeColor', 'none');

x_coords = [-300 600 600 -300]; y_coords = [10.8 10.8 9.7 9.7]; % Highlight Latissimus (significant muscle)
fill(x_coords, y_coords,[0.5, 0.5, 0.5],'FaceAlpha', 0.25, 'EdgeColor', 'none');


% Add vertical jitter manually to y-coordinates
jitter_amount = 0.1;
scatter(bicepsNT_averages, 12*ones(NTsubjSize,1) + (rand(NTsubjSize,1)-0.5)*jitter_amount, 50, colorsNT(1,:), 'filled', 'MarkerFaceAlpha', 0.7);
scatter(bicepsASD_averages, 11.5*ones(ASDsubjSize,1) + (rand(ASDsubjSize,1)-0.5)*jitter_amount, 50, colorsASD(1,:), 'filled', 'MarkerFaceAlpha', 0.7);
scatter(tricepsNT_averages, 10.5*ones(NTsubjSize,1) + (rand(NTsubjSize,1)-0.5)*jitter_amount, 50, colorsNT(2,:), 'filled', 'MarkerFaceAlpha', 0.7);
scatter(tricepsASD_averages, 10*ones(ASDsubjSize,1) + (rand(ASDsubjSize,1)-0.5)*jitter_amount, 50, colorsASD(2,:), 'filled', 'MarkerFaceAlpha', 0.7);
scatter(rectusAbdNT_averages, 9*ones(NTsubjSize,1) + (rand(NTsubjSize,1)-0.5)*jitter_amount, 50, colorsNT(3,:), 'filled', 'MarkerFaceAlpha', 0.7);
scatter(rectusAbdASD_averages, 8.5*ones(ASDsubjSize,1) + (rand(ASDsubjSize,1)-0.5)*jitter_amount, 50, colorsASD(3,:), 'filled', 'MarkerFaceAlpha', 0.7);
scatter(obliquesNT_averages, 7.5*ones(NTsubjSize,1) + (rand(NTsubjSize,1)-0.5)*jitter_amount, 50, colorsNT(4,:), 'filled', 'MarkerFaceAlpha', 0.7);
scatter(obliquesASD_averages, 7*ones(ASDsubjSize,1) + (rand(ASDsubjSize,1)-0.5)*jitter_amount, 50, colorsASD(4,:), 'filled', 'MarkerFaceAlpha', 0.7);
scatter(latissimusNT_averages, 6*ones(NTsubjSize,1) + (rand(NTsubjSize,1)-0.5)*jitter_amount, 50, colorsNT(5,:), 'filled', 'MarkerFaceAlpha', 0.7);
scatter(latissimusASD_averages, 5.5*ones(ASDsubjSize,1) + (rand(ASDsubjSize,1)-0.5)*jitter_amount, 50, colorsASD(5,:), 'filled', 'MarkerFaceAlpha', 0.7);
scatter(erectorNT_averages, 4.5*ones(NTsubjSize,1) + (rand(NTsubjSize,1)-0.5)*jitter_amount, 50, colorsNT(6,:), 'filled', 'MarkerFaceAlpha', 0.7);
scatter(erectorASD_averages, 4*ones(ASDsubjSize,1) + (rand(ASDsubjSize,1)-0.5)*jitter_amount, 50, colorsASD(6,:), 'filled', 'MarkerFaceAlpha', 0.7);
scatter(tibialisAntNT_averages, 3*ones(NTsubjSize,1) + (rand(NTsubjSize,1)-0.5)*jitter_amount, 50, colorsNT(7,:), 'filled', 'MarkerFaceAlpha', 0.7);
scatter(tibialisAntASD_averages, 2.5*ones(ASDsubjSize,1) + (rand(ASDsubjSize,1)-0.5)*jitter_amount, 50, colorsASD(7,:), 'filled', 'MarkerFaceAlpha', 0.7);
scatter(gastrocnemiusNT_averages, 1.5*ones(NTsubjSize,1) + (rand(NTsubjSize,1)-0.5)*jitter_amount, 50, colorsNT(8,:), 'filled', 'MarkerFaceAlpha', 0.7);
scatter(gastrocnemiusASD_averages, 1*ones(ASDsubjSize,1) + (rand(ASDsubjSize,1)-0.5)*jitter_amount, 50, colorsASD(8,:), 'filled', 'MarkerFaceAlpha', 0.7);

positions = [12 11.5 10.5 10 9 8.5 7.5 7 6 5.5 4.5 4 3 2.5 1.5 1];
h1 = boxplot(boxplot_data, group_labels, 'Symbol', '', 'Orientation', 'Horizontal', 'Positions', positions);

set(h1(:,1), 'Color', colorsASD(8,:), 'LineWidth', 1.5);
set(h1(:,2), 'Color', colorsNT(8,:), 'LineWidth', 1.5);
set(h1(:,3), 'Color', colorsASD(7,:), 'LineWidth', 1.5); 
set(h1(:,4), 'Color', colorsNT(7,:), 'LineWidth', 1.5);
set(h1(:,5), 'Color', colorsASD(6,:), 'LineWidth', 1.5); 
set(h1(:,6), 'Color', colorsNT(6,:), 'LineWidth', 1.5);
set(h1(:,7), 'Color', colorsASD(5,:), 'LineWidth', 1.5); 
set(h1(:,8), 'Color', colorsNT(5,:), 'LineWidth', 1.5);
set(h1(:,9), 'Color', colorsASD(4,:), 'LineWidth', 1.5); 
set(h1(:,10), 'Color', colorsNT(4,:), 'LineWidth', 1.5);
set(h1(:,11), 'Color', colorsASD(3,:), 'LineWidth', 1.5); 
set(h1(:,12), 'Color', colorsNT(3,:), 'LineWidth', 1.5);
set(h1(:,13), 'Color', colorsASD(2,:), 'LineWidth', 1.5); 
set(h1(:,14), 'Color', colorsNT(2,:), 'LineWidth', 1.5);
set(h1(:,15), 'Color', colorsASD(1,:), 'LineWidth', 1.5); 
set(h1(:,16), 'Color', colorsNT(1,:), 'LineWidth', 1.5);

set(gca, ...
    'YTickLabel', {'Gastrocnemius ASD', 'Gastrocnemius NT', ...
                    'Tibialis Ant ASD', 'Tibialis Ant NT', ...
                    'Erector ASD', 'Erector NT', ...
                    'Latissimus ASD', 'Latissimus NT', ...
                    'Obliques ASD', 'Obliques NT', ...
                    'RectusAbd ASD', 'RectusAbd NT', ...
                    'Triceps ASD', 'Triceps NT', ...
                    'Biceps ASD', 'Biceps NT'},'fontsize',12);
xlabel('Time (ms)'); title('Activation Onset Time Relative to Biceps Activation')
xline(0, 'k', 'linewidth', 3, 'Alpha', 1);
h2=findobj('LineStyle','--'); set(h2, 'LineStyle','-');

ylim([0.5 12.5]);
%xlim([-230, 50])



%%

% Number of False starts between NT and ASD
falseStartNT = 0;
for i = 1:24
    if ~isempty(onsetData(i).FalseStart)
        falseStartNT = falseStartNT + onsetData(i).FalseStart;
    end
end

falseStartASD = 0;
for i = 25:38
    falseStartASD = falseStartASD + onsetData(i).FalseStart;
end

% Total number of good trials
trialGoodNT = 0;
for i = 1:24
    trialGoodNT = trialGoodNT + length(onsetData(i).movementOnsetSec);
end

trialGoodASD = 0;
for i = 25:38
    trialGoodASD = trialGoodASD + length(onsetData(i).movementOnsetSec);
end



