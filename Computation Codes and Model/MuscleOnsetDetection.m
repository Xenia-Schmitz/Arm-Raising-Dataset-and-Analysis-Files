% =================
% EMG APA Detection
% =================
% This code extracts muscle onset for all trials

% =================================
% Subject Information and Constants
% =================================
NT_subj = [222 223 226 230 231 232 234 239 244 250 253 258 259 263 264 265 271 276 277 279 280 282];
ASD_subj = [20 21 215 218 233 240 267 269 270 273 274 284 319 326];

fs = 1000; % Sampling frequency
min_duration = round(0.08*fs);


for n = 1:length(ASD_subj) % Change to NT_subj for NT analysis

    subjID = ASD_subj(n);  % Change to NT_subj for NT analysis
    disp('Loading raw EMG data (reaching to 8 target in sagittal plane)')

    % Load Data
    file = load(['subj' num2str(subjID) '_data.mat']);
    fileName = fieldnames(file);
    data = file.subj_data;
    emgchannels = file.emgchannels;
    kinchannels = file.kinchannels;

    % Realign EMG data (22ms delay relative to kinematic data)
    for i = 1:length(data)            
        data(i).emgtime = data(i).emgtime - 0.022;
    end

    % Identify which elements have non-empty 'info'
    hasInfo = ~arrayfun(@(x) isempty(x.info), data);
    data = data(hasInfo);

    % Process EMG
    info = getInfo(data);

    par.type = 'EmgData';
    par.chlabels = emgchannels;
    
    sa = SynergyAnalyzer(data,info,par);
    sa.opt.verbose = 1;

    for i=1:length(data)
      itrial = sa.findTrials('type',[2 1 i]); 
      if ~isempty(itrial)
            ind(i) = itrial(1);
      end
    end
    
    % FIR bandpass 30-300Hz
    sa.opt.emgFilter.type = 'fir1 band';
    sa.opt.emgFilter.resample = 0;
    sa.opt.emgFilter.par = [60 30/(2000/2) 300/(2000/2)];
    sah = sa.dataFilter;
    
    % Notch 60Hz
    sah.opt.emgFilter.type = 'notch';
    sah.opt.emgFilter.resample = 0;
    sah.opt.emgFilter.par = [60/(2000/2) 30]; 
    san = sah.dataFilter;
    
    % Rectify
    san.opt.emgFilter.type = 'rectify';
    san.opt.emgFilter.resample = 0;
    san.opt.emgFilter.par = 1; 
    sar = san.dataFilter;
    
    % Butterworth low-pass
    sar.opt.emgFilter.type = 'butter low';
    sar.opt.emgFilter.resample = 0;
    sar.opt.emgFilter.par = [4 10/(2000/2)]; 
    sal = sar.dataFilter;
    
    % Resample
    sal.opt.emgFilter.resample = 1;
    sal.opt.emgFilter.resample_period = .01;
    saf = sal.dataFilter;
    
    % Average trials
    saf.opt.average.gr = saf.groupTrials('type3',[1:length(data)]');
    saf.opt.average.trange = [-.622 .978];
    sav = saf.average;
    
    % Normalize
    sav.opt.normalize.type = 2;
    sav = sav.normalize;

    % APA Onset Detection

    % Muscle assignment
    bicepsEMG = []; tricepsEMG = []; 
    rectusAbdEMG = []; LatissimusEMG = []; ObliqueEMG = []; ErectorEMG = [];
    for i = 1:length(data)
        bicepsEMG = [bicepsEMG sav.data(1,i).data(1,:)'];
        tricepsEMG = [tricepsEMG sav.data(1,i).data(2,:)'];
        rectusAbdEMG = [rectusAbdEMG sav.data(1,i).data(4,:)'];
        LatissimusEMG = [LatissimusEMG sav.data(1,i).data(5,:)'];
        ObliqueEMG = [ObliqueEMG sav.data(1,i).data(6,:)'];
        ErectorEMG = [ErectorEMG sav.data(1,i).data(7,:)'];
    end
    

    onset_time = [];

    for i = 1:length(data)
        for j = 1:9
            emg_signal = sav.data(1,i).data(j,:);
            time_signal = sav.data(1,i).time;
            
            % Spline interpolation
            splineFit = spline(time_signal, emg_signal); 
            fine_t = linspace(min(time_signal), max(emg_signal), 1625);
            fine_emg = ppval(splineFit, fine_t);

            % Baseline and threshold
            baseline = mean(fine_emg(248:259+round(0.05*fs)));
            threshold = baseline + 2.5 * std(fine_emg(248:259+round(0.05*fs)));
            above_threshold = fine_emg > threshold;

            % Find onset (first time signal crosses threshold)
            onset_idx_candidates = find(fine_emg > threshold);
            onset_idx = [];

            for k = 1:length(onset_idx_candidates)
                candidate = onset_idx_candidates(k);
                if candidate + min_duration <= length(above_threshold)
                    if above_threshold(candidate:candidate+min_duration)
                        onset_idx = [onset_idx candidate];
                    end
                end
            end

            thresholdValues(i,j) = threshold;
            
            % Apply timing constraints
            if isempty(onset_idx)
                onset_time(j,i) = NaN; % No onset detected
            else
                if fine_t(onset_idx(1)) < -0.25
                    newFirst = find(fine_t(onset_idx) >= -0.25, 1, 'first');
                    newFirst = onset_idx(newFirst);
                    if ~isempty(newFirst)
                        onset_time(j,i) = fine_t(newFirst(1));
                        onset_idx(1) = newFirst(1);
                    else
                        onset_time(j,i) = NaN;
                    end
                elseif fine_t(onset_idx(1)) > 0.2
                    onset_time(j,i) = NaN;
                else
                    onset_time(j,i) = fine_t(onset_idx(1));
                end
            end
            if all((fine_emg(453:1075))-baseline < 0.2) %453:1075
                onset_time(j,i) = NaN;
            end
            if any(fine_emg(1:275) >= 0.6)
                onset_time(j,i) = NaN;
            end
        end
    end

   % Plot all trials with APA threshold
    yaxis = [0:3:27]; 
    numLines = length(data);

    cmap = [linspace(1, 0.4, numLines)', ...
            linspace(0.6, 0, numLines)', ...
            linspace(0.6, 0, numLines)'];
    for i = 1:length(data)
        subplot(1,length(data),i); hold on;
        time_emg = sav.data(i).time(1,:);
        time_kin = data(i).postime(1,:);
        if i == 1
            yticks(yaxis);
            yticklabels(flip(emgchannels(1:9)'));
        else
            set(gca, 'YTickLabel', []);
            set(gca, 'YTick', []);
        end
        x_coords = [sav.data(i).time(1,19) sav.data(i).time(1,41) sav.data(i).time(1,41) sav.data(i).time(1,19)];
        y_coords = [-0.1 -0.1 27 27];
        fill(x_coords, y_coords, 'k','FaceAlpha', 0.25, 'EdgeColor', 'none');
        for k = 1:(length(yaxis)-2)
            index = 9-k;
            APA_index = 1;
            if all(isnan(sav.data(1,i).data(index,:)))
                yline(0);
            else
            yline(yaxis(k+1)+thresholdValues(i,index),'r');

            plot(time_emg, yaxis(k+1)+1*sav.data(1,i).data(index,:),'Color', cmap(3, :),'linewidth',1.5);
            fill([time_emg fliplr(time_emg)], [(yaxis(k+1)+1*sav.data(i).data(index,:)) yaxis(k+1)*ones(size(time_emg))], cmap(3, :), 'FaceAlpha', 0.3); 
            plot([onset_time(index,i) onset_time(index,i)], [yaxis(k+1) yaxis(k+1)+1], '-.k', 'LineWidth',2);
            title(string(i)); xlabel('Time [s]');
            xlim([-0.6, 1]);
            APA_index = APA_index+1;
            end
        end
        ax = gca; ax.FontSize = 14;
    end

    % Save APA timing
    APAtimeAllParticipants(n).subjID = subjID;
    APAtimeAllParticipants(n).timing = onset_time;
end