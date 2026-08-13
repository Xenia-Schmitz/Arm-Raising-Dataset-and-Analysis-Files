% ============================
% Kinematic Metrics Extraction
% ============================
% This code extracts the following kinematic metrics:
%   (1) Initiation time based on 10% of peak tangential velocity
%   (2) Termination time based on 10% of peak tangential velocity
%   (3) Magnitude of peak tangential velocity
%   (4) Timing of peak tangential velocity

% ===================
% Subject Information
% ===================
allParticipants = [222 223	226	230	234 231	239	250	244	253	232	263	259	258	266	229	264	265	271	280	279	276	282	277	218	20	21 240	233	267	269	215	270	274	284 273 319 326];

orderedNT = [222 223 226 230 231 232 234 239 244 250 253 258 259 263 264 265 271 276 277 279 280 282];
orderedASD = [20 21 215 218 233 240 267 269 270 273 274 284 319 326];

excluded = [266 229];

% ================================
% Paths, Constants, Initialization
% ================================
PathRawData = 'C:\Users\Action Lab\OneDrive - Northeastern University\Desktop\OneDrive - Northeastern University\Action Lab\Arm Raising\Armlifting\MATLAB\Sabrina Codes\';
fs = 2000;
dt = 1/fs;

sgolayOrder = 4;
sgolayWindow = 11;

ControlLiftData = load('UpdatedMovementOnsetInfo.mat');
goodTrialNT = load('ParsedTrialsNT.mat');
goodTrialsASD = load('ParsedTrialsASD.mat');

velocityBARnt = {};         velocityBARasd = {};
terminationNT = {};         terminationASD = {};
initiationNT = {};          initiationASD = {};
peakNT = {};                peakASD = {};
startNT = 0;                startASD = 0;


for subjSeq = 1:length(allParticipants)

    subj = allParticipants(subjSeq); 
    if ismember(subj, excluded)
        continue;
    end

    % Determine Group
    if subjSeq <= 24
        subjID = find(orderedNT == subj);
        goodTrials = goodTrialNT.ParsedTrialsNT{subjID};
    else
        subjID = find(orderedASD == subj);
        goodTrials = goodTrialASD.ParsedTrialsASD{subjID};   
    end

    % Load Movement Onset
    movementOnset = round(ControlLiftData.subject(subjSeq).movementOnset)';

    % Load Kinematics
    fileName = sprintf('%d_Control_Lift',subj);
    fileName2 = sprintf('qtm_%d_Control_Lift',subj);
    
    FILE = load([PathRawData fileName '.mat']);
    load([PathRawData fileName '_qtm.mat']);

    labels = FILE.(filename2).Trajectories.Labeled.Labels;

    % Identify and extract relevant markers
    for x = 1:length(labels)
        currLabel = FILE.(fileName2).Trajectories.Labeled.Labels(x);
        if any(strcmp(currLabel, 'Shoulder'))
            shoulderMarker = x;
        end
        if any(strcmp(currLabel, 'Bar2'))
            barMarker = x;
        end
    end

    shoulderDataY = FILE.(fileName2).Trajectories.Labeled.Data(shoulderMarker,:,:); shoulderDataY = -1*reshape(cell2mat(num2cell(shoulderDataY(1,2,:))),1,length(shoulderDataY));
    shoulderDataZ = FILE.(fileName2).Trajectories.Labeled.Data(shoulderMarker,:,:); shoulderDataZ = -1*reshape(cell2mat(num2cell(shoulderDataZ(1,3,:))),1,length(shoulderDataZ));

    barDataY = FILE.(fileName2).Trajectories.Labeled.Data(barMarker,:,:); barDataY = -1*reshape(cell2mat(num2cell(barDataY(1,2,:))),1,length(barDataY));
    barDataZ = FILE.(fileName2).Trajectories.Labeled.Data(barMarker,:,:); barDataZ = -1*reshape(cell2mat(num2cell(barDataZ(1,3,:))),1,length(barDataZ));
    
    % Filter Bar Data
    y = sgolayfilt(barDataY, sgolayOrder, sgolayWindow);
    z = sgolayfilt(-barDataZ, sgolayOrder, sgolayWindow);

    % Compute tangential velocity
    dp = diff([y', z']);
    speedBAR = sqrt(sum(dp.^2, 2));
    BARtan = speedBAR./dt;

    % Peak velocity detection
    velocityBAR = {}; timeBAR = {};
    terminationValues = []; peakValues = []; initiationValues = [];

    for i = 1:length(movementOnset)

        if subjSeq <= 24
            startNT = startNT + 1;
            movementLength = mvmntNT(startNT);
        else
            startASD = startASD + 1;
            movementLength = mvmntASD(startASD);
        end
        
        velocityTemp = BARtan(movementOnset(i)-120:movementOnset(i)+mvmntNT(startNT)+80);
        [val peakIDX] = max(velocityTemp(1:end-100));

        firsthalf = velocityTemp(1:peakIDX);
        secondhalf = velocityTemp(peakIDX:end);

        initiation = find(firsthalf > 0.1*val, 1, 'first');
        termination = find(secondhalf < 0.1*val, 1, 'first');

        velocityBAR{i} = velocityTemp(initiation:peakIDX+termination);

        if isempty(termination)
            terminationValues = [terminationValues NaN];
            peakValues = [peakValues NaN];
        else
            terminationValues = [terminationValues movementOnset(i)-120+peakIDX+termination];
            peakValues = [peakValues movementOnset(i)-120+peakIDX];
        end

        if isempty(initiation)
            initiationValues = [initiationValues NaN];
        else
            initiationValues = [initiationValues movementOnset(i)-120+initiation];
        end
    end

    % Store values
    if subjSeq <= 24
        velocityBARnt{subjSeq} = velocityBAR;
        terminationNT{subjSeq} = terminationValues;
        initiationNT{subjSeq} = initiationValues;
        peakNT{subjSeq} = peakValues;
    else
        velocityBARasd{subjSeq} = velocityBAR;
        terminationASD{subjSeq} = terminationValues;
        initiationASD{subjSeq} = initiationValues;
        peakASD{subjSeq} = peakValues;
    end
end

%%
%% Load in data file for  arm raising task
allData = readtable('AllExtractedMetrics2.xlsx');
allData.TIBIALIS_ANTERIOR = str2double(allData.TIBIALIS_ANTERIOR);


allNTData = allData(strcmp(allData.GROUP, 'NT'),:); % Parse out NT data
allASDData = allData(strcmp(allData.GROUP, 'ASD'),:); % Parse out ASD data

allsubjID = unique(allData.SUBJ_NUM);
subjNTid = unique(allNTData.SUBJ_NUM); % Parse out NT subject numbers
subjASDid = unique(allASDData.SUBJ_NUM); % Parse out ASD subject numbers


%% KINEMATIC PLOTTING

%%%%%%%%%%%%%%%%%%%%%% PLOTTING OF KINEMATIC METRICS %%%%%%%%%%%%%%%%%%%%%%

% Find average of PEAK_MAG, MT, LIFT1, LIFT2 for each subject
[peakVelNT] = subjectAverageCalculator([allNTData.SUBJ_NUM allNTData.PEAK_MAG]);
[mtNT] = subjectAverageCalculator([allNTData.SUBJ_NUM allNTData.MT]);
[lift1TimeNT] = subjectAverageCalculator([allNTData.SUBJ_NUM allNTData.LIFT1]);
[lift2TimeNT] = subjectAverageCalculator([allNTData.SUBJ_NUM allNTData.LIFT2]);

[peakVelASD] = subjectAverageCalculator([allASDData.SUBJ_NUM allASDData.PEAK_MAG]);
[mtASD] = subjectAverageCalculator([allASDData.SUBJ_NUM allASDData.MT]);
[lift1TimeASD] = subjectAverageCalculator([allASDData.SUBJ_NUM allASDData.LIFT1]);
[lift2TimeASD] = subjectAverageCalculator([allASDData.SUBJ_NUM allASDData.LIFT2]);

% Plot boxplots for each of the metrics overlayed with subject averages as scatter
%figure; hold on; % tiledlayout("horizontal");
nexttile; hold on; BoxPlotGenerator(peakVelNT/1000, peakVelASD/1000, 'Peak Tangential Bar Velocity', 'Velocity (m/s)');
nexttile; hold on; BoxPlotGenerator(mtNT*1000, mtASD*1000, 'Movement Time', 'Time (ms)');
nexttile; hold on; BoxPlotGenerator(lift1TimeNT*1000, lift1TimeASD*1000, 'Acceleration Phase Duration', 'Time (ms)');
nexttile; hold on; BoxPlotGenerator(lift2TimeNT*1000, lift2TimeASD*1000, 'Deceleration Phase Duration', 'Time (ms)');
%%
% Statistics for peak tangential arm velocity
kinematicsTable = allData(:,1:8);
kinematicsTable.GROUP = categorical(kinematicsTable.GROUP);
lme_peakvel = fitlme(kinematicsTable, 'PEAK_MAG ~ GROUP + (1|SUBJ_NUM)');

[beta,~,stats] = fixedEffects(lme_peakvel,'DFMethod','satterthwaite');
disp(stats)

%%
% Statistics for movement time
lme_MT = fitlme(kinematicsTable, 'MT ~ GROUP + (1|SUBJ_NUM)');
[beta,~,stats] = fixedEffects(lme_MT,'DFMethod','satterthwaite');
disp(stats)

%%
% Statistics for decomposed arm trajectory

% Reshape data into long format
kinematicsLong = stack(kinematicsTable, ...
    {'LIFT1','LIFT2'}, ...
    'NewDataVariableName','VALUE', ...
    'IndexVariableName','PHASE');

% Convert PHASE into categorical with nice labels
kinematicsLong.PHASE = categorical(kinematicsLong.PHASE, ...
    {'LIFT1','LIFT2'}, ...
    {'Acceleration','Deceleration'});

head(kinematicsLong)

% Fit mixed model with interaction
lme_decomposed = fitlme(kinematicsLong, 'VALUE ~ GROUP*PHASE + (1+GROUP|SUBJ_NUM)');

[beta,~,stats] = fixedEffects(lme_decomposed,'DFMethod','satterthwaite');
coefNames = lme_decomposed.Coefficients.Name;  % Names of coefficients

phaseLabels = categories(kinematicsLong.PHASE);
groups = categories(kinematicsLong.GROUP);

fprintf('Group Differences within Each Phase (Satterthwaite-corrected DF)\n');
fprintf('---------------------------------------------------------------\n');

for i = 1:length(phaseLabels)
    
    % Build contrast vector L
    L = zeros(1,length(coefNames));
    
    % Always include GROUP main effect
    idx_group = strcmp(coefNames,'GROUP_NT');
    L(idx_group) = 1; % Reference is first level of muscle
    
    % Add interaction term if not reference muscle
    interactionName = ['GROUP_NT:PHASE_' char(phaseLabels{i})];
    idx_inter = strcmp(coefNames,interactionName);
    if any(idx_inter)
        L(idx_inter) = 1;
    end
    
    % Contrast estimate (linear combo of coeffs -- estimated group diff for
    % muscle)
    contrastEstimate = L * beta;
    
    % Standard error of contrast
    Covb = lme_decomposed.CoefficientCovariance;
    contrastSE = sqrt(L * Covb * L'); % To get standard error. Accounts for correlation btwn coeffs (group and group x muscle)
    
    % t-statistic
    tVal = contrastEstimate / contrastSE;
    
    % DF from Satterthwaite
    df = stats.DF;  % satterthwaite-corrected DF
    idx_nonzero = find(L~=0);
    df_contrast = min(stats.DF(idx_nonzero));
    
    % Two-tailed p-value
    pVal = 2*(1 - tcdf(abs(tVal), df_contrast));
    
    %fprintf('%-25s  t(%.1f) = %.3f, p = %.4f\n', phaseLabels{i}, df_contrast, tVal, pVal);
    fprintf('%-20s  b = %.3f, SE = %.3f, t(%.1f) = %.3f, p = %.4f\n', ...
    phaseLabels{i}, contrastEstimate, contrastSE, df_contrast, tVal, pVal);
    
end
%%
M_reduced = fitlme(kinematicsLong, ...
    'VALUE ~ PHASE + (1|SUBJ_NUM)');
compare(M_reduced, lme_decomposed)