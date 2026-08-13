%% ============================
% Subject Information
% ============================

allParticipants = [222 223 226 230 234 231 239 250 244 253 ...
                   232 263 259 258 266 229 264 265 271 280 ...
                   279 276 282 277 218 20 21 240 233 267 ...
                   269 215 270 274 284 273 319 326];

orderedNT = [222 223 226 230 231 232 234 239 ...
             244 250 253 258 259 263 264 265 ...
             271 276 277 279 280 282];

orderedASD = [20 21 215 218 233 240 ...
              267 269 270 273 274 284 ...
              319 326];

excluded = [266 229];

%% ============================
% Paths and Constants
% ============================

PathRawData = ...
'C:\Users\Action Lab\OneDrive - Northeastern University\Desktop\OneDrive - Northeastern University\Action Lab\Arm Raising\Armlifting\MATLAB\Sabrina Codes\';

fs = 2000;
dt = 1/fs;

nInterp = 101;

%% ============================
% Load Data
% ============================

load initiationNT2.mat
load initiationASD2.mat

load terminationNT4.mat
load terminationASD3.mat

%% ============================
% Storage
% ============================


trajNT  = {};
trajASD = {};

velNT   = {};
velASD  = {};

NTcount  = 0;
ASDcount = 0;

subjnumber = [];

% ============================
% Subject Loop
% ============================

for subjSeq = 1:length(allParticipants)

    subj = allParticipants(subjSeq);

    if ismember(subj,excluded)
        continue
    end

    % ------------------------
    % Determine Group
    % -------------------------

    if subjSeq <= 24

        subjID = find(orderedNT==subj);

        starts = initiationNT{subjSeq};
        stops  = terminationNT{subjSeq};

        isNT = true;

        NTcount = NTcount + 1;

        trajNT{NTcount} = {};
        velNT{NTcount}  = {};

    else

        subjID = find(orderedASD==subj);

        starts = initiationASD{subjSeq};
        stops  = terminationASD{subjSeq};

        isNT = false;

        ASDcount = ASDcount + 1;

        trajASD{ASDcount} = {};
        velASD{ASDcount}  = {};

    end

    % ------------------------
    % Load Subject Data
    % -------------------------

    fileName = sprintf('%d_Control_Lift',subj);
    qtmName  = sprintf('qtm_%d_Control_Lift',subj);

    FILE = load([PathRawData fileName '.mat']);

    labels = FILE.(qtmName).Trajectories.Labeled.Labels;

    % ------------------------
    % Find Markers
    % -------------------------

    barMarker = [];
    shoulderMarker = [];

    for k = 1:length(labels)

        if strcmp(labels{k},'Bar2')
            barMarker = k;
        elseif strcmp(labels{k},'Shoulder')
            shoulderMarker = k;
        end

    end

    if isempty(barMarker)

        fprintf('Bar2 not found for subject %d\n',subj);
        continue

    end

    if isempty(shoulderMarker)

        fprintf('Shoulder not found for subject %d\n',subj);
        continue

    end

    % ------------------------
    % Extract Marker Data
    % -------------------------

    barXYZ = squeeze( ...
        FILE.(qtmName).Trajectories.Labeled.Data(barMarker,:,:));

    shoulderXYZ = squeeze( ...
        FILE.(qtmName).Trajectories.Labeled.Data(shoulderMarker,:,:));

    % ------------------------
    % Trial Loop
    % -------------------------

    validCount = 0;

    for tr = 1:length(starts)
        if subjSeq == 9 && tr == 7
            continue
        end
        if subjSeq == 23 && tr == 7 
            continue
        end
        if subjSeq == 5+24 && tr == 6
            continue;
        end
        if subjSeq == 12+24
            continue;
        end


        startIdx = round(starts(tr));
        stopIdx  = round(stops(tr));

        if isnan(startIdx) || isnan(stopIdx)
            continue
        end

        if stopIdx <= startIdx
            continue
        end

        if startIdx < 1
            continue
        end

        if stopIdx > size(barXYZ,2)
            continue
        end

        % ------------------------
        % Relative Position
        % ------------------------

        barPos = barXYZ(:,startIdx:stopIdx);
        shoulderPos = shoulderXYZ(:,startIdx:stopIdx);

        % Initial arm length
        armLength = norm(barPos((1:3),1) - shoulderPos((1:3),1));

        if isnan(armLength) || armLength < eps
            continue
        end

        % Position relative to shoulder
        relPos = (barPos - shoulderPos) ./ armLength;

        x = relPos(1,:);
        y = relPos(2,:);
        z = relPos(3,:);

% ------------------------
% Relative Position
% ------------------------
% 
% barPos      = barXYZ((1:3),startIdx:stopIdx);
% shoulderPos = shoulderXYZ((1:3),startIdx:stopIdx);
% 
% % Initial arm vector
% armVec0 = barPos((1:3),1) - shoulderPos((1:3),1);
% armLength = norm(armVec0);
% 
% if isnan(armLength) || armLength < eps
%     continue
% end
% 
% % Normalize initial arm vector
% zLocal = armVec0 / armLength;
% 
% % Global vertical
% up = [0;0;1];
% 
% % Build orthogonal coordinate system
% xLocal = cross(up,zLocal);
% xLocal = xLocal / norm(xLocal);
% 
% yLocal = cross(zLocal,xLocal);
% yLocal = yLocal / norm(yLocal);
% 
% % Rotation matrix
% R = [xLocal yLocal zLocal];
% 
% % Relative coordinates
% relPos = zeros(size(barPos));
% 
% for i = 1:size(barPos,2)
% 
%     rel = barPos((1:3),i) - shoulderPos((1:3),i);
% 
%     % Rotate into shoulder coordinate system
%     relPos(:,i) = R' * rel;
% 
% end
% 
% % Normalize by arm length
% relPos = relPos / armLength;
% 
% x = relPos(1,:);
% y = relPos(2,:);
% z = relPos(3,:);
        % ------------------------
        % Velocity
        % ------------------------

        vx = gradient(x)*fs;
        vy = gradient(y)*fs;
        vz = gradient(z)*fs;

        % ------------------------
        % Normalize Trial Length
        % ------------------------

        oldT = linspace(0,1,length(x));
        newT = linspace(0,1,nInterp);

        posTrial = nan(nInterp,3);
        velTrial = nan(nInterp,3);

        posTrial(:,1) = interp1(oldT,x,newT,'linear');
        posTrial(:,2) = interp1(oldT,y,newT,'linear');
        posTrial(:,3) = interp1(oldT,z,newT,'linear');

        velTrial(:,1) = interp1(oldT,vx,newT,'linear');
        velTrial(:,2) = interp1(oldT,vy,newT,'linear');
        velTrial(:,3) = interp1(oldT,vz,newT,'linear');

        validCount = validCount + 1;

        % ------------------------
        % Store Trial
        % ------------------------

        if isNT

            trajNT{NTcount}{validCount} = posTrial;
            
            velNT{NTcount}{validCount}  = velTrial;

        else

            trajASD{ASDcount}{validCount} = posTrial;
            velASD{ASDcount}{validCount}  = velTrial;

        end

    end
    if isNT
    subjnumber = [subjnumber subjSeq];
    end
    fprintf('Subject %d : %d valid trials\n',subj,validCount)

end
%%
% ============================
% NT Position
% ============================

figure
hold on

cmap = lines(length(trajNT));
nextPlot = 1;

for s = 1:length(trajNT)
    if s == 17 || s ==19 || s == 20 || s ==22
        continue
    else
        subplot(3,6,nextPlot); hold on;
        %subtightplot(3,6,nextPlot,[0.02 0.02],[0.05 0.05],[0.03 0.03]);
        nextPlot = nextPlot+1;
        
    end

    for tr = 1:length(trajNT{s})

        X = trajNT{s}{tr};

        if isnan(X)
            continue
        end

        plot3(X(:,1),...
              X(:,2),...
              X(:,3),...
              'Color',[0.5 0.5 0.5],...
              'LineWidth',1)

        scatter3(X(1,1),...
                 X(1,2),...
                 X(1,3),...
                 40,...
                 'k',...
                 'filled')
grid on
axis equal
view(3)
    end

end

ax = findall(gcf,'Type','axes');

for k = 1:length(ax)
    pos = ax(k).Position;

    pos(3) = pos(3)*1.12;   % 12% wider
    pos(4) = pos(4)*1.12;   % 12% taller

    ax(k).Position = pos;
end

grid on
axis equal
view(3)

% xlabel('X')
% ylabel('Y')
% zlabel('Z')

sgtitle('NT Position Trajectories')

%% ============================
% NT Position - ML
% ============================

figure
hold on

cmap = lines(length(trajNT));
nextPlot = 1;

for s = 1:length(trajNT)

    if ismember(s,[17 19 20 22])
        continue
    end

    subplot(3,6,nextPlot)
    hold on
    nextPlot = nextPlot + 1;

    % ---------- Determine average forward direction for this subject ----------
    dirs = [];

    for tr = 1:length(trajNT{s})

        X = trajNT{s}{tr};

        if any(isnan(X(:))) || size(X,1) < 6
            continue
        end

        % Horizontal movement over first few samples
        d = X(6,1:2) - X(1,1:2);

        if norm(d) < 1e-6
            continue
        end

        dirs = [dirs;
                d/norm(d)];
    end

    % Average forward direction
    forward = mean(dirs,1);
    forward = forward / norm(forward);

    % Rotation angle about the vertical (Z) axis
    theta = atan2(forward(2),forward(1));

    % Rotation matrix
    R = [ cos(-theta)  -sin(-theta)   0;
          sin(-theta)   cos(-theta)   0;
               0             0        1];

    % ---------- Plot all trials ----------
    for tr = 1:length(trajNT{s})

        X = trajNT{s}{tr};

        if any(isnan(X(:)))
            continue
        end

        % Translate so every trajectory starts at the origin
        X0 = X - X(1,:);

        % Rotate into subject coordinates
        Xrot = (R * X0')';

        % Plot sagittal plane (forward vs vertical)
        plot(Xrot(:,1), Xrot(:,3), ...
            'Color', [0.5 0.5 0.5], ...
            'LineWidth', 1)

        scatter(Xrot(1,1), Xrot(1,3), ...
            40, 'k', 'filled')
    end

    grid on
    axis equal
    % xlabel('Forward (m)')
    % ylabel('Vertical (m)')
    title(sprintf('Subject %d',s))

end

sgtitle('NT Position Trajectories')


%%
figure

nextPlot = 1;

for s = 1:length(trajNT)

    if ismember(s,[17 19 20 22])
        continue
    end

    % ---------- Compute subject rotation ----------
    dirs = [];

    for tr = 1:length(trajNT{s})

        X = trajNT{s}{tr};

        if any(isnan(X(:))) || size(X,1) < 6
            continue
        end

        d = X(6,1:2) - X(1,1:2);

        if norm(d) < 1e-6
            continue
        end

        dirs = [dirs;
                d/norm(d)];
    end

    forward = mean(dirs,1);
    forward = forward/norm(forward);

    theta = atan2(forward(2),forward(1));

    R = [ cos(-theta) -sin(-theta) 0;
          sin(-theta)  cos(-theta) 0;
               0            0      1];

    % ---------- Rotate all trials ----------
    rotTraj = {};

    for tr = 1:length(trajNT{s})

        X = trajNT{s}{tr};

        if any(isnan(X(:)))
            continue
        end

        X0 = X - X(1,:);
        Xrot = (R*X0')';

        rotTraj{end+1} = Xrot;
    end

    if isempty(rotTraj)
        continue
    end

    % ---------- Stack trajectories ----------
    nTrials = numel(rotTraj);
    nPts = size(rotTraj{1},1);

    Xall = nan(nPts,nTrials);
    Yall = nan(nPts,nTrials);

    for tr = 1:nTrials
        Xall(:,tr) = rotTraj{tr}(:,1);
        Yall(:,tr) = rotTraj{tr}(:,2);
    end

    % ---------- Compute SD ----------
    sdX  = std(Xall,0,2,'omitnan');
    sdY  = std(Yall,0,2,'omitnan');
    sdXY = sqrt(sdX.^2 + sdY.^2);

    % ---------- Plot ----------
    subplot(3,6,nextPlot)
    hold on

    t = linspace(0,100,nPts);

    plot(t,sdX,'r','LineWidth',1.5)
    plot(t,sdY,'b','LineWidth',1.5)
    plot(t,sdXY,'k','LineWidth',2)

    title(sprintf('Subject %d',s))
    xlabel('% Movement')
    ylabel('SD (m)')
    grid on

    nextPlot = nextPlot + 1;

end

legend({'X','Y','X+Y'},'Location','bestoutside')
sgtitle('Subject Trajectory Variability')
%%
%% ============================
% NT Trajectory Bundle Variability
% ============================
% 
% bundleSD_NT = cell(length(trajNT),1);
% 
% for s = 1:length(trajNT)
% 
%     nTrials = length(trajNT{s});
% 
%     % Skip subjects with no trials
%     if nTrials == 0
%         continue
%     end
% 
%     % Stack trials into a 3D array
%     % Dimensions: time x coordinate x trial
%     trajArray = nan(nInterp,3,nTrials);
% 
%     for tr = 1:nTrials
%         trajArray(:,:,tr) = trajNT{s}{tr};
%     end
% 
%     % Standard deviation across trials
%     stdTrajectory = std(trajArray,0,3,'omitnan');
% 
%     % Magnitude of the SD vector at each time point
%     bundleSD_NT{s} = ...
%     sqrt( ...
%         stdTrajectory(:,1).^2 + ...
%         stdTrajectory(:,2).^2 + ...
%         stdTrajectory(:,3).^2 );
% %stdTrajectory(:,3);
% end
% 
% figure
% hold on
% 
% for s = 1:length(bundleSD_NT)
% 
%     if isempty(bundleSD_NT{s})
%         continue
%     end
% 
%     plot(linspace(0,100,nInterp),bundleSD_NT{s},'LineWidth',1.5)
% 
% end
% 
% xlabel('Normalized Movement (%)')
% ylabel('Trajectory Variability')
% title('NT Trajectory Bundle Variability')
% grid on

%% ============================
% ASD Position
% ============================

figure
hold on

cmap = lines(length(trajASD));
nextPlot = 1;

for s = 1:length(trajASD)
    if s == 12 || s ==7
        continue;
    else
        subplot(3,6,nextPlot); hold on;
        %subtightplot(3,6,nextPlot,[0.02 0.02],[0.05 0.05],[0.03 0.03]);
        nextPlot = nextPlot+1;
    end
        
    for tr = 1:length(trajASD{s})

        X = trajASD{s}{tr};

        if isnan(X)
            continue
        end

        plot3(X(:,1),...
              X(:,2),...
              X(:,3),...
              'Color',[0.5 0.5 0.5],...
              'LineWidth',1)

        scatter3(X(1,1),...
                 X(1,2),...
                 X(1,3),...
                 40,...
                 'k',...
                 'filled')
grid on
axis equal
view(3)
    end

end

ax = findall(gcf,'Type','axes');

for k = 1:length(ax)
    pos = ax(k).Position;

    pos(3) = pos(3)*1.12;   % 12% wider
    pos(4) = pos(4)*1.12;   % 12% taller

    ax(k).Position = pos;
end

grid on
axis equal
view(3)

% xlabel('X')
% ylabel('Y')
% zlabel('Z')

sgtitle('ASD Position Trajectories')

%% ============================
% ASD Position - ML
% ============================

figure
hold on

cmap = lines(length(trajASD));
nextPlot = 1;

for s = 1:length(trajASD)

    if ismember(s,[7, 12])
        continue
    end

    subplot(3,6,nextPlot)
    hold on
    nextPlot = nextPlot + 1;

    % ---------- Determine average forward direction for this subject ----------
    dirs = [];

    for tr = 1:length(trajASD{s})

        X = trajASD{s}{tr};

        if any(isnan(X(:))) || size(X,1) < 6
            continue
        end

        % Horizontal movement over first few samples
        d = X(6,1:2) - X(1,1:2);

        if norm(d) < 1e-6
            continue
        end

        dirs = [dirs;
                d/norm(d)];
    end
    
    % Average forward direction
    forward = mean(dirs,1);
    forward = forward / norm(forward);

    % Rotation angle about the vertical (Z) axis
    theta = atan2(forward(2),forward(1));

    % Rotation matrix
    R = [ cos(-theta)  -sin(-theta)   0;
          sin(-theta)   cos(-theta)   0;
               0             0        1];

    % ---------- Plot all trials ----------
    for tr = 1:length(trajASD{s})

        X = trajASD{s}{tr};

        if any(isnan(X(:)))
            continue
        end

        % Translate so every trajectory starts at the origin
        X0 = X - X(1,:);

        % Rotate into subject coordinates
        Xrot = (R * X0')';

        % Plot sagittal plane (forward vs vertical)
        plot(Xrot(:,2), Xrot(:,3), ...
            'Color', [0.5 0.5 0.5], ...
            'LineWidth', 1)

        scatter(Xrot(1,2), Xrot(1,3), ...
            40, 'k', 'filled')
    end

    grid on
    axis equal
    % xlabel('Forward (m)')
    % ylabel('Vertical (m)')
    title(sprintf('Subject %d',s))

end

sgtitle('ASD Position Trajectories')

% %% ============================
% % ASD Trajectory Bundle Variability
% % ============================
% 
% bundleSD_ASD = cell(length(trajASD),1);
% 
% for s = 1:length(trajASD)
% 
%     nTrials = length(trajASD{s});
% 
%     % Skip subjects with no trials
%     if nTrials == 0
%         continue
%     end
% 
%     % Stack trials into a 3D array
%     % Dimensions: time x coordinate x trial
%     trajArray = nan(nInterp,3,nTrials);
% 
%     for tr = 1:nTrials
%         trajArray(:,:,tr) = trajASD{s}{tr};
%     end
% 
%     % Standard deviation across trials
%     stdTrajectory = std(trajArray,0,3,'omitnan');
% 
%     % Magnitude of the SD vector at each time point
%     bundleSD_ASD{s} = ...
%     sqrt( ...
%         stdTrajectory(:,1).^2 + ...
%         stdTrajectory(:,2).^2 + ...
%         stdTrajectory(:,3).^2 );
% 
% end
% 
% figure
% hold on
% 
% for s = 1:length(bundleSD_ASD)
% 
%     if isempty(bundleSD_ASD{s})
%         continue
%     end
% 
%     plot(linspace(0,100,nInterp),bundleSD_ASD{s},'LineWidth',1.5)
% 
% end
% 
% xlabel('Normalized Movement (%)')
% ylabel('Trajectory Variability')
% title('ASD Trajectory Bundle Variability')
% grid on

%%
figure

nextPlot = 1;

for s = 1:length(trajASD)

    if ismember(s,[7 12])
        continue
    end

    % ---------- Compute subject rotation ----------
    dirs = [];

    for tr = 1:length(trajASD{s})

        X = trajASD{s}{tr};

        if any(isnan(X(:))) || size(X,1) < 6
            continue
        end

        d = X(6,1:2) - X(1,1:2);

        if norm(d) < 1e-6
            continue
        end

        dirs = [dirs;
                d/norm(d)];
    end

    forward = mean(dirs,1);
    forward = forward/norm(forward);

    theta = atan2(forward(2),forward(1));

    R = [ cos(-theta) -sin(-theta) 0;
          sin(-theta)  cos(-theta) 0;
               0            0      1];

    % ---------- Rotate all trials ----------
    rotTraj = {};

    for tr = 1:length(trajASD{s})

        X = trajASD{s}{tr};

        if any(isnan(X(:)))
            continue
        end

        X0 = X - X(1,:);
        Xrot = (R*X0')';

        rotTraj{end+1} = Xrot;
    end

    if isempty(rotTraj)
        continue
    end

    % ---------- Stack trajectories ----------
    nTrials = numel(rotTraj);
    nPts = size(rotTraj{1},1);

    Xall = nan(nPts,nTrials);
    Yall = nan(nPts,nTrials);

    for tr = 1:nTrials
        Xall(:,tr) = rotTraj{tr}(:,1);
        Yall(:,tr) = rotTraj{tr}(:,2);
    end

    % ---------- Compute SD ----------
    sdX  = std(Xall,0,2,'omitnan');
    sdY  = std(Yall,0,2,'omitnan');
    sdXY = sqrt(sdX.^2 + sdY.^2);

    % ---------- Plot ----------
    subplot(3,6,nextPlot)
    hold on

    t = linspace(0,100,nPts);

    plot(t,sdX,'r','LineWidth',1.5)
    plot(t,sdY,'b','LineWidth',1.5)
    plot(t,sdXY,'k','LineWidth',2)

    title(sprintf('Subject %d',s))
    xlabel('% Movement')
    ylabel('SD (m)')
    grid on

    nextPlot = nextPlot + 1;

end

legend({'X','Y','X+Y'},'Location','bestoutside')
sgtitle('Subject Trajectory Variability')
%% ============================
%% ============================================
% Group Mean + Covariance Ellipsoids
% ============================================

figure
hold on

plotCovarianceTrajectory(trajNT,[0 0.45 0.74])
plotCovarianceTrajectory(trajASD,[0.85 0.2 0.2])

grid on
axis equal
view(3)

xlabel('X')
ylabel('Y')
zlabel('Z')

title('Mean Trajectory with Covariance Ellipsoids')

legend({'NT','ASD'})

%%
%%
%% ============================================
% Trajectory Density Maps (XZ Plane)
% ============================================

figure

nBins = 100;      % number of histogram bins

% -----------------------------
% NT
% -----------------------------

allX = [];
allZ = [];

for s = 1:length(trajNT)

    if isempty(trajNT{s}) 
        continue
    end

    for tr = 1:length(trajNT{s})

        X = trajNT{s}{tr};
        if isnan(X)
            continue
        end

        allX = [allX; X(:,1)];
        allZ = [allZ; X(:,3)];

    end

end

subplot(1,2,1)

[N,Xedges,Zedges] = histcounts2(allX,allZ,nBins);

imagesc(Xedges,Zedges,N')

axis xy
axis equal
box on

xlabel('Anterior-Posterior (mm)')
ylabel('Vertical (mm)')

title('NT Trajectory Density')

colormap(gca,parula)
colorbar


% -----------------------------
% ASD
% -----------------------------

allX = [];
allZ = [];

for s = 1:length(trajASD)

    if isempty(trajASD{s})
        continue
    end

    for tr = 1:length(trajASD{s})

        X = trajASD{s}{tr};
        if isnan(X)
            continue
        end        

        allX = [allX; X(:,1)];
        allZ = [allZ; X(:,3)];

    end

end

subplot(1,2,2)

[N,Xedges,Zedges] = histcounts2(allX,allZ,nBins);

imagesc(Xedges,Zedges,N')

axis xy
axis equal
box on

xlabel('Anterior-Posterior (mm)')
ylabel('Vertical (mm)')

title('ASD Trajectory Density')

colormap(gca,parula)
colorbar

%%

function drawCovarianceEllipsoid(mu,C,color)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eigen decomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[V,D] = eig(C);

radii = sqrt(diag(D));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unit sphere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x,y,z] = sphere(20);

xyz = [x(:)'; y(:)'; z(:)'];

xyz = V*diag(radii)*xyz;

X = reshape(xyz(1,:),size(x)) + mu(1);
Y = reshape(xyz(2,:),size(y)) + mu(2);
Z = reshape(xyz(3,:),size(z)) + mu(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

surf(X,...
     Y,...
     Z,...
     'FaceColor',color,...
     'EdgeColor','none',...
     'FaceAlpha',0.25)

end

function plotCovarianceTrajectory(trajGroup,color)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute subject mean trajectories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nSub = length(trajGroup);

% first non-empty subject
idx = find(~cellfun(@isempty,trajGroup),1);

nInterp = size(trajGroup{idx}{1},1);

subjectMean = [];

for s = 1:nSub

    if isempty(trajGroup{s})
        continue
    end

    nTrials = length(trajGroup{s});

    temp = nan(nInterp,3,nTrials);

    for tr = 1:nTrials
        temp(:,:,tr) = trajGroup{s}{tr};
    end

    subjectMean(:,:,end+1) = mean(temp,3,'omitnan');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Group mean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

groupMean = mean(subjectMean,3,'omitnan');

plot3(groupMean(:,1),...
      groupMean(:,2),...
      groupMean(:,3),...
      'Color',color,...
      'LineWidth',4)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ellipsoids every 10%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pts = round(linspace(1,nInterp,11));

for p = pts

    P = squeeze(subjectMean(p,:,:))';

    % remove NaN subjects
    P = P(~any(isnan(P),2),:);

    if size(P,1)<3
        continue
    end

    C = cov(P);

    drawCovarianceEllipsoid(mean(P),C,color)

end

end