% ============================
% Kinetics Metrics Extraction
% ============================
% This code compute entropy of CoP velocity

% ===================
% Subject Information
% ===================
orderedNT = [222 223 226 230 231 232 234 239 244 250 253 258 259 263 264 265 271 276 277 279 280 282];
orderedASD = [20 21 215 218 233 240 267 269 270 273 274 284 319 326];

 allParticipants = [222 223	226	230	234 231	239	250	244	253	232	263	259	258	266	229	264	265	271	280	279	276	282	277	218	20	21 240	233	267	269	215	270	274	284 273 319 326];

% ================================
% Paths, Constants, Initialization
% ================================
PathRawData = 'C:\Users\Action Lab\OneDrive - Northeastern University\Desktop\OneDrive - Northeastern University\Action Lab\Arm Raising\Armlifting\MATLAB\Sabrina Codes\';
%PathRawData = 'C:\Users\xvssc\OneDrive - Northeastern University\Action Lab\Arm Raising\Armlifting\MATLAB\Sabrina Codes\';
clear entropyPrepAllDataNT entropyPrepAllDataASD entropy90AllDataNT entropy90AllDataASD
clear entropy180AllDataNT entropy180AllDataASD entropyAllDataNT entropyAllDataASD

% ====================
% Entropy Computation
% ====================
for subjSeq = 1:length(allParticipants)
    subj = allParticipants(subjSeq);
    if subj == 266 || subj == 229
        continue;
    end
    if subjSeq <= 24
        subjID = find(orderedNT == subj);
        ControlLiftData = load(['UpdatedMovementOnsetInfo.mat']);
        goodTrialInfo = load(['ParsedTrialsNT.mat']);
        goodTrials = goodTrialInfo.ParsedTrialsNT{subjID};
        load('initiationNT2.mat');
        load('terminationNT4.mat');
        load('velocityMaxTimeNT3.mat')
        movementOnset = initiationNT{subjSeq};
        terminationIDX = terminationNT{subjSeq};
        peakIDX = velocityMaxTimeNT{subjSeq}*200;
    else
        subjID = find(orderedASD == subj);
        ControlLiftData = load(['UpdatedMovementOnsetInfo.mat']);
        goodTrialInfo = load(['ParsedTrialsASD.mat']);
        goodTrials = goodTrialInfo.ParsedTrialsASD{subjID};   
        load('initiationASD2.mat');
        load('terminationASD3.mat');
        load('velocityMaxTimeASD3.mat');
        movementOnset = initiationASD{subjSeq};
        terminationIDX = terminationASD{subjSeq};
        peakIDX = velocityMaxTimeASD{subjSeq}*200;
    end

    % Load in CoP data
    fileName = sprintf('%d_Control_Lift',subj);
    FILE = load([PathRawData fileName '.mat']);
    fileName2 = sprintf('qtm_%d_Control_Lift',subj);
    load([PathRawData fileName '_qtm.mat']);
    load([PathRawData fileName '.mat']);

    [rotatedCoPX,rotatedCoPY,rotatedCoPZ,cordinateSystemRot]=TransformingCoP(FILE.(fileName2),subj);    
    rotatedCoPX = rotatedCoPX'; rotatedCoPY = rotatedCoPY'; rotatedCoPZ = rotatedCoPZ';
    
    % Filter Data
    m = 4; fl = 9;
    rotatedCoPX = sgolayfilt(rotatedCoPX,m,fl);
    rotatedCoPY = sgolayfilt(rotatedCoPY,m,fl);
    rotatedCoPZ = sgolayfilt(rotatedCoPZ,m,fl);

    % Base off of peak velocity 
    peakValues = peakIDX;

    % Calculate CoP Angles relative to kinematic features
        clear entropyCOPpreparation entropyCOPlift90 entropyCOPlift180 entropyCOPall
         angvelCOPpreparation = [];
         angvelCOPacceleration = [];
         angvelCOPdeceleration = [];
        if ~isempty(terminationIDX)       
            for i = 1:length(movementOnset) 
                if ~isnan(terminationIDX(i)) && ~isnan(movementOnset(i)) && ~isnan(peakValues(i))
                    CoP90Index = movementOnset(i) + peakValues(i) + 1;
                    angleIndex = terminationIDX(i);
                    CoP180Index = angleIndex;                 
    
                    maxnStep = 10;
                    for nStep = 1 %1:maxnStep
                        APAstartWindow = movementOnset(i) - 40;
                        [pathPrep(nStep) pPrep] = calPathlength(rotatedCoPX,rotatedCoPY,APAstartWindow,movementOnset(i),nStep);
                        pathTimePrep(nStep) = length(APAstartWindow:movementOnset(i));
                        [pathLiftBegin(nStep) pLiftBegin] = calPathlength(rotatedCoPX,rotatedCoPY,movementOnset(i),CoP90Index,nStep);
                        [pathLiftEnd(nStep) pLiftEnd] = calPathlength(rotatedCoPX,rotatedCoPY,CoP90Index,CoP180Index,nStep);
                    end
    
                    % SHANNON's ENTROPY
                    % 1. CoP during prep: CoP of movementOnset-40 to movementOnset
                    histN = 15;
                    angvelCOPprep = atan2(diff(rotatedCoPY((movementOnset(i)-40):movementOnset(i))), diff(rotatedCoPX((movementOnset(i)-40):movementOnset(i))));
                    [angvelCOPprephist, binN] = hist(angvelCOPprep, histN);
                    probangvelCOPprephist = angvelCOPprephist./length(angvelCOPprep);
                    infoangvelCOPprephist = -1*probangvelCOPprephist.*log(probangvelCOPprephist);
                    infoangvelCOPprephist(isnan(infoangvelCOPprephist)) = 0;
                    entropyCOPpreparation(i) = sum(infoangvelCOPprephist);

                    % 2. CoP during first half arm raise: CoP of movementOnset to 90 degrees
                    angvelCOPlift = atan2(diff(rotatedCoPY(movementOnset(i):CoP90Index)),diff(rotatedCoPX(movementOnset(i):CoP90Index)));           
                    angvelCOPlifthist = hist(angvelCOPlift,[0:2*pi/histN:(2*pi-2*pi/histN)]+pi/histN-pi);
                    [angvelCOPlifthist, binN] = hist(angvelCOPlift,histN);            
                    probangvelCOPlifthist = angvelCOPlifthist./length(angvelCOPlift);
                    infoangvelCOPlifthist = -1* probangvelCOPlifthist.*log(probangvelCOPlifthist);
                    infoangvelCOPlifthist(isnan(infoangvelCOPlifthist)) = 0;
                    entropyCOPlift90(i) = sum(infoangvelCOPlifthist);

                    % 3. CoP during second half arm raise: CoP of movementOnset to max degrees
                    angvelCOPlift = atan2(diff(rotatedCoPY(CoP90Index:CoP180Index)),diff(rotatedCoPX(CoP90Index:CoP180Index)));           
                    angvelCOPlifthist = hist(angvelCOPlift,[0:2*pi/histN:(2*pi-2*pi/histN)]+pi/histN-pi);
                    [angvelCOPlifthist, binN] = hist(angvelCOPlift,histN);            
                    probangvelCOPlifthist = angvelCOPlifthist./length(angvelCOPlift);
                    infoangvelCOPlifthist = -1* probangvelCOPlifthist.*log(probangvelCOPlifthist);
                    infoangvelCOPlifthist(isnan(infoangvelCOPlifthist)) = 0;
                    entropyCOPlift180(i) = sum(infoangvelCOPlifthist);

                    % Angular Data 
                    angvelCOPprep = atan2(diff(rotatedCoPY((movementOnset(i)-40):movementOnset(i))), diff(rotatedCoPX((movementOnset(i)-40):movementOnset(i))));
                    angvelCOPlift1 = atan2(diff(rotatedCoPY(movementOnset(i):CoP90Index)),diff(rotatedCoPX(movementOnset(i):CoP90Index)));
                    angvelCOPlift2 = atan2(diff(rotatedCoPY(CoP90Index:CoP180Index)),diff(rotatedCoPX(CoP90Index:CoP180Index)));

                    angvelCOPpreparation = [angvelCOPpreparation angvelCOPprep];
                    angvelCOPacceleration = [angvelCOPacceleration angvelCOPlift1];
                    angvelCOPdeceleration = [angvelCOPdeceleration angvelCOPlift2];                    
                end
            end
        end

    % Collect all data
    if ~isempty(terminationIDX)  && ~all(isnan(peakValues)) 
        if ~all(isnan((terminationIDX)))

                % Processed Entropy Data
                entropyCOPprep{subjSeq} = entropyCOPpreparation;
                entropyCOPlift1{subjSeq} = entropyCOPlift90;
                entropyCOPlift2{subjSeq} = entropyCOPlift180;

                % Filtered Position and Velocity Data
                copPrepDatax{subjSeq} = rotatedCoPX(movementOnset(i)-40:movementOnset(i));
                copPrepDatay{subjSeq} = rotatedCoPY(movementOnset(i)-40:movementOnset(i));
                copLift1Datax{subjSeq} = rotatedCoPX(movementOnset(i):CoP90Index);
                copLift1Datay{subjSeq} = rotatedCoPY(movementOnset(i):CoP90Index);
                copLift2Datax{subjSeq} = rotatedCoPX(CoP90Index:CoP180Index);
                copLift2Datay{subjSeq} = rotatedCoPY(CoP90Index:CoP180Index);
                copPrepDataVelx{subjSeq} = diff(rotatedCoPX(movementOnset(i)-40:movementOnset(i)));
                copPrepDataVely{subjSeq} = diff(rotatedCoPY(movementOnset(i)-40:movementOnset(i)));
                copLift1DataVelx{subjSeq} = diff(rotatedCoPX(movementOnset(i):CoP90Index));
                copLift1DataVely{subjSeq} = diff(rotatedCoPY(movementOnset(i):CoP90Index));
                copLift2DataVelx{subjSeq} = diff(rotatedCoPX(CoP90Index:CoP180Index));
                copLift2DataVely{subjSeq} = diff(rotatedCoPY(CoP90Index:CoP180Index));

                % Angular Data
                angPrepAll{subjSeq} = angvelCOPprep;
                angLift1All{subjSeq} = angvelCOPlift1;
                angLift2All{subjSeq} = angvelCOPlift2;
        else
            disp(subjSeq)
        end
    else
        disp(subjSeq)
   end
end

%% ============
% Plotting
% =============
subjects = 38;

for n = 1:24
    prepData = entropyCOPprep{n}; prepData(prepData == 0) = NaN;
    entropyAllData(n).entropyPrep = prepData;
    Data90 = entropyCOPlift1{n}; Data90(Data90 == 0) = NaN;
    entropyAllData(n).entropy90 = Data90;
    Data180 = entropyCOPlift2{n}; Data180(Data180 == 0) = NaN;
    entropyAllData(n).entropy180 = Data180;
end

for n = 25:38
    prepDataASD = entropyCOPprep{n}; prepDataASD(prepDataASD == 0) = NaN;
    entropyAllData(n).entropyPrep = prepDataASD;
    Data90ASD = entropyCOPlift1{n}; Data90ASD(Data90ASD == 0) = NaN;
    entropyAllData(n).entropy90 = Data90ASD;
    Data180ASD = entropyCOPlift2{n}; Data180ASD(Data180ASD == 0) = NaN;
    entropyAllData(n).entropy180 = Data180ASD;
end

% Box Plot of subject averages, by group, all three phases
prepEntropyNT = []; lift1EntropyNT = []; lift2EntropyNT = [];
prepEntropyASD = []; lift1EntropyASD = []; lift2EntropyASD = [];

for n = 1:24
    prepEntropyNT = [prepEntropyNT; nanmean(entropyAllData(n).entropyPrep')];
    lift1EntropyNT = [lift1EntropyNT; nanmean(entropyAllData(n).entropy90')];
    lift2EntropyNT = [lift2EntropyNT; nanmean(entropyAllData(n).entropy180')];
end
for n = 25:38
    prepEntropyASD = [prepEntropyASD; nanmean(entropyAllData(n).entropyPrep')];
    lift1EntropyASD = [lift1EntropyASD; nanmean(entropyAllData(n).entropy90')];
    lift2EntropyASD = [lift2EntropyASD; nanmean(entropyAllData(n).entropy180')];
end

% Plot boxplots for each of the metrics overlayed with subject averages as scatter
figure; hold on; tiledlayout("horizontal");
nexttile; hold on; BoxPlotGenerator(prepEntropyNT, prepEntropyASD, 'Entropy During Preparation Phase', 'Entropy');
nexttile; hold on; BoxPlotGenerator(lift1EntropyNT, lift1EntropyASD, 'Entropy During Acceleration Phase', 'Entropy');
nexttile; hold on; BoxPlotGenerator(lift2EntropyNT, lift2EntropyASD, 'Entropy During Deceleration Phase', 'Entropy');

% Plot CoP Position and Velocity (Exemplar)
figure; hold on;
subplot(2,2,1); hold on; % NT position
for i = 1:14
    CoPx = [copLift1Datax{i} copLift2Datax{i}];
    CoPy = [copLift1Datay{i} copLift2Datay{i}];
    plot(CoPx, CoPy, 'b', 'linewidth',2);
    scatter(CoPx(1), CoPy(1),'g','filled');
    scatter(CoPx(end), CoPy(end),'k','filled');
    xlim([-0.04,0.06]); ylim([-0.08,0.06]);
end
subplot(2,2,2); hold on; % ASD position
for i = 25:38
    CoPx = [copLift1Datax{i} copLift2Datax{i}];
    CoPy = [copLift1Datay{i} copLift2Datay{i}];
    plot(CoPx, CoPy, 'r','linewidth',2);
    scatter(CoPx(1), CoPy(1),'g','filled');
    scatter(CoPx(end), CoPy(end),'k','filled');
    xlim([-0.04,0.06]); ylim([-0.08,0.06]);
end

subplot(2,2,3); hold on; % NT velocity
for i = 1:14
    CoPx = [copLift1DataVelx{i} copLift2DataVelx{i}];
    CoPy = [copLift1DataVely{i} copLift2DataVely{i}];
    plot(CoPx, CoPy, 'b', 'linewidth',2);
    scatter(CoPx(1), CoPy(1),'g','filled');
    scatter(CoPx(end), CoPy(end),'k','filled');
    xlim([-0.004,0.005]); ylim([-0.008,0.006]);
end
subplot(2,2,4); hold on; % ASD velocity
for i = 25:38
    CoPx = [copLift1DataVelx{i} copLift2DataVelx{i}];
    CoPy = [copLift1DataVely{i} copLift2DataVely{i}];
    plot(CoPx, CoPy, 'r','linewidth',2);
    scatter(CoPx(1), CoPy(1),'g','filled');
    scatter(CoPx(end), CoPy(end),'k','filled');
    xlim([-0.004,0.005]); ylim([-0.008,0.006]);
end

% Plot exemplar rose plots for each phase
NT_idx = 1:24; 
ASD_idx = 25:length(allParticipants);

NT_data = {angPrepAll(NT_idx)', angLift1All(NT_idx)', angLift2All(NT_idx)'};
ASD_data = {angPrepAll(ASD_idx)', angLift1All(ASD_idx)', angLift2All(ASD_idx)'};
explain_your_data_pattern(NT_data, ASD_data)


function explain_your_data_pattern(NT_data, ASD_data)    
    for phase_idx = 1:3
        all_angles = [];
        all_groups = [];

        for s = 1:14
            if ~isempty(NT_data{phase_idx}{s})
                this_angles = NT_data{phase_idx}{s}(:) * 180/pi;   % radians → degrees
                all_angles = [all_angles; this_angles];
                all_groups = [all_groups; repmat(s, numel(this_angles), 1)];
            end
        end

        [figure_handle,count,speeds,directions,Table,Others] = WindRose(all_angles,all_groups,'ndirs',15, 'vwinds',[1:1:14], 'labels', {'Front', 'Right', 'Back', 'Left'}, ...
            'cmap', 'jet');
        
        % Update colorbar and figure colormap
        load('NTColorMap.mat');
        colormap(figure_handle, NTColorMap);
        %colorbar;
        
        % Manually recolor the wedges
        patch_handles = findall(figure_handle, 'Type', 'Patch');
        
        nColors = size(NTColorMap, 1);
        for i = 1:length(patch_handles)
            color_idx = mod(i-1, nColors) + 1;  
            set(patch_handles(i), 'FaceColor', NTColorMap(color_idx,:), 'EdgeColor', 'none');
        end
        set(findall(gcf,'Type','Patch'),'EdgeColor','none');  % simpler vectors

    end
end
%% Statistics
phases = {'entropyPrep', 'entropy90', 'entropy180'};       
phaseNames = {'Preparation', 'Lift 1', 'Lift 2'};

% Preallocate
allSubject = [];
allGroup = [];
allPhase = [];
allEntropy = [];
groupLabel = [repmat("NT", 24, 1); repmat("ASD", 14, 1)];
for s = 1:subjects
    numTrials = length(entropyAllData(s).entropyPrep); % assuming same across muscles
    for m = 1:length(phases)
        entropys = entropyAllData(s).(phases{m});
        allSubject = [allSubject; repmat(s, numTrials, 1)];
        allGroup = [allGroup; repmat(groupLabel(s), numTrials, 1)]; % groupLabel must be pre-defined, e.g., NT=0, ASD=1
        allPhase = [allPhase; repmat(phaseNames(m), numTrials, 1)];
        allEntropy = [allEntropy; entropys(:)];
    end
end

% Create table
T = table(categorical(allSubject), categorical(allGroup), ...
          categorical(allPhase), allEntropy, ...
          'VariableNames', {'Subject', 'Group', 'Phase', 'Entropy'});

% Linear Mixed Model
M3 = fitlme(T, 'Entropy ~ Group * Phase + (Phase|Subject)');

[beta,~,stats] = fixedEffects(M3,'DFMethod','satterthwaite');
coefNames = M3.Coefficients.Name;  % Names of coefficients

phaseLabels = categories(T.Phase);
groups = categories(T.Group);

fprintf('Group Differences within Each Muscle (Satterthwaite-corrected DF)\n');
fprintf('---------------------------------------------------------------\n');

% phaseLabels = {'Phase_Preparation', 'Phase_Lift 1', 'Phase_Lift 2'};
for i = 1:length(phaseLabels)
    
    % Build contrast vector L
    L = zeros(1,length(coefNames));
    
    % Always include GROUP main effect
    idx_group = strcmp(coefNames,'Group_NT');
    L(idx_group) = 1; % Reference is first level of muscle
    
    % Add interaction term if not reference muscle
    interactionName = ['Group_NT:Phase_' char(phaseLabels{i})];
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
    df_contrast = min(stats.DF(idx_nonzero));
    
    % Two-tailed p-value
    pVal = 2*(1 - tcdf(abs(tVal), df_contrast));
    
    %fprintf('%-25s  t(%.1f) = %.3f, p = %.4f\n', phaseLabels{i}, df_contrast, tVal, pVal);
    fprintf('%-20s  b = %.3f, SE = %.3f, t(%.1f) = %.3f, p = %.4f\n', ...
    phaseLabels{i}, contrastEstimate, contrastSE, df_contrast, tVal, pVal);
    
end

M_reduced = fitlme(T, ...
    'Entropy ~ Phase + (Phase | Subject)', ...
    'FitMethod','ML');
compare(M_reduced, M3)
