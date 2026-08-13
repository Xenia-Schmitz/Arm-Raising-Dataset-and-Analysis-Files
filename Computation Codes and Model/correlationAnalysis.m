% ==================================================================
% Correlations Between Arm-Raising Performance and Clinical Measures
% ==================================================================

clear; clc;

% ==========
% Load Data
% ==========
movement   = readtable('AllExtractedMetrics2.xlsx');
behavioral = readtable('AllSubjectMetrics.xlsx');
movement.TIBIALIS_ANTERIOR = str2double(movement.TIBIALIS_ANTERIOR);

% =============================================================
% Aggregate movement metrics to subject-level means (omit NaNs)
% =============================================================

idVars   = {'SUBJ_NUM','GROUP'};
moveVars = movement.Properties.VariableNames([2,4:end]);

movement = groupsummary( ...
    movement, ...
    idVars, ...
    @(x) mean(x,'omitnan'), ...
    moveVars);

% Clean variable names
movement.Properties.VariableNames = ...
    strrep(movement.Properties.VariableNames,'','');

% =================
% Define variables
% =================
behavior_vars = behavioral.Properties.VariableNames(6:end);
movement_vars = movement.Properties.VariableNames(6:end);

n_behav = numel(behavior_vars);
n_move  = numel(movement_vars);

% ======================
% NT group correlations
% ======================
behNT_all = behavioral(strcmp(behavioral.GROUP,'NT'), :);
movNT_all = movement(strcmp(movement.GROUP,'NT'), :);

NT = innerjoin(behNT_all, movNT_all, 'Keys','SUBJ_NUM');

R_NT = NaN(n_behav,n_move);
P_NT = NaN(n_behav,n_move);

for i = 1:n_behav
    b = NT{:,behavior_vars{i}};
    for j = 1:n_move
        m = NT{:,movement_vars{j}};
        v = ~isnan(b) & ~isnan(m);

        if sum(v) >= 3 && std(b(v)) > 0 && std(m(v)) > 0
            [R_NT(i,j),P_NT(i,j)] = corr(b(v),m(v),'Type','Spearman');
        end
    end
end

% =======================
% ASD group correlations
% =======================
behASD_all = behavioral(strcmp(behavioral.GROUP,'ASD'), :);
movASD_all = movement(strcmp(movement.GROUP,'ASD'), :);

ASD = innerjoin(behASD_all, movASD_all, 'Keys','SUBJ_NUM');

R_ASD = NaN(n_behav,n_move);
P_ASD = NaN(n_behav,n_move);

for i = 1:n_behav
    b = ASD{:,behavior_vars{i}};
    for j = 1:n_move
        m = ASD{:,movement_vars{j}};
        v = ~isnan(b) & ~isnan(m);

        if sum(v) >= 3 && std(b(v)) > 0 && std(m(v)) > 0
            [R_ASD(i,j),P_ASD(i,j)] = corr(b(v),m(v),'Type','Spearman');
        end
    end
end

% ==========================
% Fisher z-tests (NT vs ASD)
% ==========================

Z_diff = NaN(n_behav,n_move);
P_z    = NaN(n_behav,n_move);

for i = 1:n_behav
    for j = 1:n_move
        r1 = R_NT(i,j);
        r2 = R_ASD(i,j);

        n1 = sum(~isnan(NT{:,behavior_vars{i}}) & ...
                 ~isnan(NT{:,movement_vars{j}}));

        n2 = sum(~isnan(ASD{:,behavior_vars{i}}) & ...
                 ~isnan(ASD{:,movement_vars{j}}));

        if ~isnan(r1) && ~isnan(r2) && n1 > 3 && n2 > 3
            z = (atanh(r1) - atanh(r2)) / sqrt(1/(n1-3) + 1/(n2-3));
            P_z(i,j)    = 2*(1-normcdf(abs(z)));
            Z_diff(i,j) = z;
        end
    end
end

% =========================================================
% ALL participants correlations (ADOS restricted to ASD)
% =========================================================

ALL = innerjoin(behavioral, movement, 'Keys','SUBJ_NUM');

R_ALL = NaN(n_behav,n_move);
P_ALL = NaN(n_behav,n_move);

for i = 1:n_behav
    b = ALL{:,behavior_vars{i}};
    for j = 1:n_move
        m = ALL{:,movement_vars{j}};

        % Restrict ADOS measures to ASD only
        if contains(behavior_vars{i},'ADOS','IgnoreCase',true)
            idx = strcmp(ALL.GROUP_behavioral,'ASD');
            b2 = b(idx);
            m2 = m(idx);
        else
            b2 = b;
            m2 = m;
        end

        v = ~isnan(b2) & ~isnan(m2);

        if sum(v) >= 3
            [R_ALL(i,j),P_ALL(i,j)] = corr(b2(v),m2(v),'Type','Spearman', 'rows','complete');
            if isnan(R_ALL(i,j))
                disp('pause')
            end
        else
            disp('pause')
        end
    end
end

% ===============
% Create Table 2
% ===============
cellOut = cell(n_behav,n_move);

for i = 1:n_behav
    for j = 1:n_move
        if ~isnan(R_ALL(i,j))
            cellOut{i,j} = sprintf('%.2f (%.3f)',R_ALL(i,j),P_ALL(i,j));
        else
            disp('pause')
            cellOut{i,j} = '';
        end
    end
end

Table2 = cell2table( ...
    cellOut, ...
    'RowNames',behavior_vars, ...
    'VariableNames',movement_vars);

% ===============
% Summary output
% ===============
fprintf('\nSummary:\n');
fprintf('NT significant correlations:  %d\n',sum(P_NT(:)<alpha));
fprintf('ASD significant correlations: %d\n',sum(P_ASD(:)<alpha));
fprintf('ALL significant correlations: %d\n',sum(P_ALL(:)<alpha));
fprintf('Fisher z group differences:   %d\n\n',sum(P_z(:)<alpha));

disp('Analysis complete. Tables written to disk.');

fprintf('\n=============================================\n');
fprintf('Significant correlations across ALL participants (Spearman, p < %.2f)\n', alpha);
fprintf('=============================================\n\n');

[sig_row_all, sig_col_all] = find(P_ALL < alpha);

if isempty(sig_row_all)
    fprintf('None.\n');
else
    for k = 1:length(sig_row_all)
        i = sig_row_all(k);
        j = sig_col_all(k);

        fprintf('Behavior: %-20s | Movement: %-30s | rho = %.2f | p = %.3f\n', ...
            behavior_vars{i}, movement_vars{j}, R_ALL(i,j), P_ALL(i,j));
    end
end

fprintf('\n=============================================\n');
fprintf('Significant Fisher''s z tests (NT vs ASD)\n');
fprintf('Only for correlations significant across ALL participants\n');
fprintf('=============================================\n\n');

found = false;

for k = 1:length(sig_row_all)
    i = sig_row_all(k);
    j = sig_col_all(k);

    if P_z(i,j) < alpha
        found = true;

        fprintf(['Behavior: %-20s | Movement: %-30s | ' ...
                 'z = %.2f | p = %.3f | ' ...
                 'NT rho = %.2f | ASD rho = %.2f\n'], ...
            behavior_vars{i}, ...
            movement_vars{j}, ...
            Z_diff(i,j), ...
            P_z(i,j), ...
            R_NT(i,j), ...
            R_ASD(i,j));
    end
end

if ~found
    fprintf('None.\n');
end



%% %%%%%%%
% ==================================================================
% Correlations Between Arm-Raising Performance and Clinical Measures
% WITH FDR CORRECTION (BH)
% ==================================================================

clear; clc;

alpha = 0.05;

% ==========
% Load Data
% ==========
movement   = readtable('AllExtractedMetrics2.xlsx');
behavioral = readtable('AllSubjectMetrics.xlsx');

movement.TIBIALIS_ANTERIOR = str2double(movement.TIBIALIS_ANTERIOR);

% =============================================================
% Aggregate movement metrics to subject-level means (omit NaNs)
% =============================================================
idVars   = {'SUBJ_NUM','GROUP'};
moveVars = movement.Properties.VariableNames([2,4:end]);

movement = groupsummary( ...
    movement, ...
    idVars, ...
    @(x) mean(x,'omitnan'), ...
    moveVars);

% =================
% Define variables
% =================
behavior_vars = behavioral.Properties.VariableNames(6:end);
movement_vars = movement.Properties.VariableNames(6:end);

n_behav = numel(behavior_vars);
n_move  = numel(movement_vars);

% ======================
% NT correlations
% ======================
behNT_all = behavioral(strcmp(behavioral.GROUP,'NT'), :);
movNT_all = movement(strcmp(movement.GROUP,'NT'), :);
NT = innerjoin(behNT_all, movNT_all, 'Keys','SUBJ_NUM');

R_NT = NaN(n_behav,n_move);
P_NT = NaN(n_behav,n_move);

for i = 1:n_behav
    b = NT{:,behavior_vars{i}};
    for j = 1:n_move
        m = NT{:,movement_vars{j}};
        v = ~isnan(b) & ~isnan(m);

        if sum(v) >= 3 && std(b(v)) > 0 && std(m(v)) > 0
            [R_NT(i,j),P_NT(i,j)] = corr(b(v),m(v),'Type','Spearman');
        end
    end
end

% ======================
% ASD correlations
% ======================
behASD_all = behavioral(strcmp(behavioral.GROUP,'ASD'), :);
movASD_all = movement(strcmp(movement.GROUP,'ASD'), :);
ASD = innerjoin(behASD_all, movASD_all, 'Keys','SUBJ_NUM');

R_ASD = NaN(n_behav,n_move);
P_ASD = NaN(n_behav,n_move);

for i = 1:n_behav
    b = ASD{:,behavior_vars{i}};
    for j = 1:n_move
        m = ASD{:,movement_vars{j}};
        v = ~isnan(b) & ~isnan(m);

        if sum(v) >= 3 && std(b(v)) > 0 && std(m(v)) > 0
            [R_ASD(i,j),P_ASD(i,j)] = corr(b(v),m(v),'Type','Spearman');
        end
    end
end

% ==========================
% Fisher z-tests (NT vs ASD)
% ==========================
Z_diff = NaN(n_behav,n_move);
P_z    = NaN(n_behav,n_move);

for i = 1:n_behav
    for j = 1:n_move

        r1 = R_NT(i,j);
        r2 = R_ASD(i,j);

        n1 = sum(~isnan(NT{:,behavior_vars{i}}) & ~isnan(NT{:,movement_vars{j}}));
        n2 = sum(~isnan(ASD{:,behavior_vars{i}}) & ~isnan(ASD{:,movement_vars{j}}));

        if ~isnan(r1) && ~isnan(r2) && n1 > 3 && n2 > 3
            z = (atanh(r1) - atanh(r2)) / sqrt(1/(n1-3) + 1/(n2-3));
            P_z(i,j)    = 2*(1-normcdf(abs(z)));
            Z_diff(i,j) = z;
        end
    end
end
%%
% =========================================================
% ALL participants correlations
% =========================================================
ALL = innerjoin(behavioral, movement, 'Keys','SUBJ_NUM');

R_ALL = NaN(n_behav,n_move);
P_ALL = NaN(n_behav,n_move);

for i = 1:n_behav
    b = ALL{:,behavior_vars{i}};
    for j = 1:n_move
        m = ALL{:,movement_vars{j}};

        if contains(behavior_vars{i},'ADOS','IgnoreCase',true)
            idx = strcmp(ALL.GROUP_behavioral, 'ASD');
            b2 = b(idx);
            m2 = m(idx);
        else
            b2 = b;
            m2 = m;
        end

        v = ~isnan(b2) & ~isnan(m2);

        if sum(v) >= 3
            [R_ALL(i,j),P_ALL(i,j)] = corr(b2(v),m2(v),'Type','Spearman');
        end
    end
end
%%
% ======================
% FDR CORRECTION SECTION
% ======================

% ---- ALL ----
P_ALL_FDR = NaN(size(P_ALL));
vec = P_ALL(:);
v = ~isnan(vec);
vec_fdr = NaN(size(vec));
vec_fdr(v) = mafdr(vec(v),'BHFDR',true);
P_ALL_FDR = reshape(vec_fdr,size(P_ALL));

% ---- NT ----
P_NT_FDR = NaN(size(P_NT));
vec = P_NT(:);
v = ~isnan(vec);
vec_fdr = NaN(size(vec));
vec_fdr(v) = mafdr(vec(v),'BHFDR',true);
P_NT_FDR = reshape(vec_fdr,size(P_NT));

% ---- ASD ----
P_ASD_FDR = NaN(size(P_ASD));
vec = P_ASD(:);
v = ~isnan(vec);
vec_fdr = NaN(size(vec));
vec_fdr(v) = mafdr(vec(v),'BHFDR',true);
P_ASD_FDR = reshape(vec_fdr,size(P_ASD));

% ---- Fisher Z ----
P_z_FDR = NaN(size(P_z));
vec = P_z(:);
v = ~isnan(vec);
vec_fdr = NaN(size(vec));
vec_fdr(v) = mafdr(vec(v),'BHFDR',true);
P_z_FDR = reshape(vec_fdr,size(P_z));

% ======================
% Summary output (FDR)
% ======================
fprintf('\nSummary (FDR corrected):\n');
fprintf('NT significant correlations:  %d\n',sum(P_NT_FDR(:) < alpha));
fprintf('ASD significant correlations: %d\n',sum(P_ASD_FDR(:) < alpha));
fprintf('ALL significant correlations: %d\n',sum(P_ALL_FDR(:) < alpha));
fprintf('Fisher z group differences:   %d\n\n',sum(P_z_FDR(:) < alpha));

% ======================
% ALL significant report
% ======================
fprintf('\n=============================================\n');
fprintf('Significant correlations (ALL, FDR corrected)\n');
fprintf('=============================================\n\n');

[sig_row_all, sig_col_all] = find(P_ALL_FDR < alpha);

if isempty(sig_row_all)
    fprintf('None.\n');
else
    for k = 1:length(sig_row_all)
        i = sig_row_all(k);
        j = sig_col_all(k);

        fprintf('Behavior: %-20s | Movement: %-30s | rho = %.2f | p_FDR = %.3f\n', ...
            behavior_vars{i}, movement_vars{j}, R_ALL(i,j), P_ALL_FDR(i,j));
    end
end

% ======================
% Fisher z (FDR)
% ======================
fprintf('\n=============================================\n');
fprintf('Significant Fisher z tests (FDR corrected)\n');
fprintf('=============================================\n\n');

found = false;

for i = 1:n_behav
    for j = 1:n_move

        if P_z_FDR(i,j) < alpha
            found = true;

            fprintf(['Behavior: %-20s | Movement: %-30s | ' ...
                     'z = %.2f | p_FDR = %.3f | NT rho = %.2f | ASD rho = %.2f\n'], ...
                behavior_vars{i}, movement_vars{j}, ...
                Z_diff(i,j), P_z_FDR(i,j), R_NT(i,j), R_ASD(i,j));
        end
    end
end

if ~found
    fprintf('None.\n');
end

disp('Analysis complete with FDR correction.');