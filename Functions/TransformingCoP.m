function [newCoPTransformedX,newCoPTransformedY,newCoPTransformedZ,cordinateSystemRot] = TransformingCoP(qtm_in, subject)

% SabrinaASD_FPdiagnostics_Rocking(Sabrina_Rocking_Diagonal0004)
%SabrinaASD_FPdiagnostics_Rocking(Bottle_4_corners0001)

SampleRate = qtm_in.FrameRate;
MarkerData = (qtm_in.Trajectories.Labeled.Data(:,:,:)) ./ 1000; % meters
MarkerData = permute(MarkerData,[3 2 1]); % Frame x Dimension x Marker
FPSampMod = qtm_in.Force(1).SamplingFactor;
FPSampleRate = SampleRate * FPSampMod;
FPForce = qtm_in.Force(1).Force';
FPCOP = qtm_in.Force(1).COP' ./ 1000;
fpLoc = qtm_in.Force(1).ForcePlateLocation ./ 1000; %meters;

if subject == 253
    FPForce = qtm_in.Force(2).Force';
    FPCOP = qtm_in.Force(2).COP'./ 1000;
    fpLoc = qtm_in.Force(2).ForcePlateLocation ./ 1000;
end



    labels = qtm_in.Trajectories.Labeled.Labels;
    
for x = 1:length(labels)
currLabel = qtm_in.Trajectories.Labeled.Labels(x);
        
        if any(strcmp(currLabel,'ToeR'))
            toeRMarker = x;
        end
        if any(strcmp(currLabel,'ToeL'))
            toeLMarker = x;
        end
end

% downsample fp data for faster plotting and same dimension as marker data, with a touch of smoothing.
FPForce = interp1(1:length(FPCOP),FPForce,linspace(1,length(FPCOP),length(MarkerData)),'makima');
FPCOP = interp1(1:length(FPCOP),FPCOP,linspace(1,length(FPCOP),length(MarkerData)),'makima');

% Add some Z value to FPCOP original, to have 'positive z' - to check whether Z is
% inverted in force plate r.f.
FPCOP(:,3) = FPCOP(:,3) + 0.1;

% Plot force plate in QTM coords
curFig = Create_Reuse_Figure([],'FP frame diagnostics',[2200 200 800 600]);
sp1 = subplot(2,2,[1 2]); hold on; % QTM REFERENCE FRAME FOR NOW
xlim([-0.6 0.6]); xlabel('QTM X (m)','fontweight','bold');
ylim([-0.6 0.8]); ylabel('QTM Y (m)','fontweight','bold');
zlim([-0.1 0.3]); zlabel('QTM Z (m)','fontweight','bold');
axis equal
view(0,90);

FPcoordsX = [0 0.4 0.4 0];
FPcoordsY = [0 0 0.6 0.6];
FPcoordsZ = [0 0 0 0];
% Normally I'd avoid using patch because it automatically deteriorates quality of
% graphics. 'Rectangle' function had a simpler input and no quality effect, but it only draws
% rectangles with sides parallel to X and Y. We will be rotating a rectangle, so that wouldn't 
% work and we need vertice coordinates instead. They could be rotated similarly to other data.
patch(FPcoordsX,FPcoordsY,FPcoordsZ,'FaceColor',[0.4 0.4 0.4],'EdgeColor','k','LineWidth',2);

%rectangle('Position',[0 0 0.4 0.6],'FaceColor',[0.4 0.4 0.4],'LineWidth',2);



ForeFeetRight = MarkerData(:,:,toeRMarker); % In case you have more markers
% Rocking. Plot forefoot markers, they're static, so use mean.
[XplotRight,YplotRight,ZplotRight] = FindMeans(ForeFeetRight);
scatter3(XplotRight,YplotRight,ZplotRight,40,'MarkerFaceColor','r','MarkerEdgeColor','none');

ForeFeetLeft = MarkerData(:,:,toeLMarker); % In case you have more markers
% Rocking. Plot forefoot markers, they're static, so use mean.
[XplotLeft,YplotLeft,ZplotLeft] = FindMeans(ForeFeetLeft);
scatter3(XplotLeft,YplotLeft,ZplotLeft,40,'MarkerFaceColor','b','MarkerEdgeColor','none');


% Create transformation (R and t) from FP to QTM
    % express r.f. basis of FP w.r.t. the r.f. of QTM
R_q_fp = [-1 0 0;...
    0 1 0;...
    0 0 1];  % X is inverted. The rest appears to be the same.
trans = [0.2 0.3 0]; % This is again FP origin w.r.t. QTM origin
FPCOP_q = (R_q_fp * FPCOP')' + trans;
irange = 1:length(FPCOP);
plot3(FPCOP_q(irange,1),FPCOP_q(irange,2),FPCOP_q(irange,3),'g','LineWidth',1);
% The plot makes sense!!!!!


%% Now try the new r.f.


sp2 = subplot(2,2,[3 4]); hold on; % This will be the new r.f.
xlim([-0.6 0.6]); xlabel('SA X (m)','fontweight','bold');
ylim([-0.6 0.8]); ylabel('SA Y (m)','fontweight','bold');
zlim([-0.1 0.3]); zlabel('SA Z (m)','fontweight','bold');
axis equal
view(0,90);

% Initialize the rotated r.f. Compared to QTM, alpha counterclockwise would be about 90 + 57 deg
alRot_q_sa = 130.53866648;

R_q_sa = axang2rotm([0 0 1 alRot_q_sa/180*pi]); % Rotation matrix, about Z axis, counterclockwise, by a given angle (deg2rad).
% This Matlab function uses clockwise angle! So added minus.
%R_q_sa = eye(3); % no rotation just for checking translation only.
trans_q_sa = -[0.2; 0.3; 0]; % It's easier for me to write "How I see new r.f. from Qualisys". 
%%%%%But I actually need% the opposite. So I added a minus. %%%% Do I?.....
%trans_q_sa = [0; 0; 0]; % no translation, just for cheking.
FPcoords = [FPcoordsX; FPcoordsY; FPcoordsZ]; %Stack them for shorter matrix notation.

FPcoords_sa = inv(R_q_sa) * (FPcoords + trans_q_sa);
patch(FPcoords_sa(1,:),FPcoords_sa(2,:),FPcoords_sa(3,:),'FaceColor',[0.4 0.4 0.4],'EdgeColor','k','LineWidth',2);

MarkVecsR_q = [XplotRight'; YplotRight'; ZplotRight'];
MarkVecsR_sa = inv(R_q_sa) * (MarkVecsR_q + trans_q_sa);
scatter3(MarkVecsR_sa(1,:),MarkVecsR_sa(2,:),MarkVecsR_sa(3,:),40,'MarkerFaceColor','r','MarkerEdgeColor','none');

MarkVecsL_q = [XplotLeft'; YplotLeft'; ZplotLeft'];
MarkVecsL_sa = inv(R_q_sa) * (MarkVecsL_q + trans_q_sa);
scatter3(MarkVecsL_sa(1,:),MarkVecsL_sa(2,:),MarkVecsL_sa(3,:),40,'MarkerFaceColor','b','MarkerEdgeColor','none');


FPCOP_sa = (inv(R_q_sa) * (FPCOP_q' + trans_q_sa))';
plot3(FPCOP_sa(irange,1),FPCOP_sa(irange,2),FPCOP_sa(irange,3),'g','LineWidth',1);

newCoPTransformedX = FPCOP_sa(irange,1);
newCoPTransformedY= FPCOP_sa(irange,2);
newCoPTransformedZ = FPCOP_sa(irange,3);
cordinateSystemRot = FPcoords_sa;

end




function [Xplot,Yplot,Zplot] = FindMeans(source,varargin)

if any(strcmpi(varargin,'std'))
    flag_std = 1;
else
    flag_std = 0;
end

if flag_std
    Xplot = squeeze(std(source(:,1,:),[],1,'omitnan'));
    Yplot = squeeze(std(source(:,2,:),[],1,'omitnan'));
    Zplot = squeeze(std(source(:,3,:),[],1,'omitnan'));
else
    Xplot = squeeze(mean(source(:,1,:),1,'omitnan'));
    Yplot = squeeze(mean(source(:,2,:),1,'omitnan'));
    Zplot = squeeze(mean(source(:,3,:),1,'omitnan'));
end

end
