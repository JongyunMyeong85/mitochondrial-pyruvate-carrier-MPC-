%% call temp_exp again
% spatio analysis

%% asking what tyype of analysis
close all; 
clear
set(0,'defaultAxesFontName', 'arial');
global set_fg set_ax 
set_fg = {'units', 'inches', 'color', 'w'};
set_ax = {'units', 'inches', 'fontsize', 10, 'box', 'off', 'color', 'none'};

answer_analize = questdlg('What type of summarized data do you want to load?', 'Load data', 'Load from Temporal Data', 'Load from Spatio Data', 'Load from Temporal Data');
% close all;, clear;
ComputerName = getComputerName;

if strcmp(answer_analize, 'Load from Spatio Data')
    if strcmp(ComputerName, 'asuslaptop')
        cd('C:\Users\jongy\My Drive\making paper\Synapse\spatial_analysis');                       % Laptop change CD
    elseif strcmp(ComputerName, 'desktop-d5sp5os')
        cd('H:\My Drive\making paper\Synapse\spatial_analysis');                                                 % change CD
    end
    [Name_Exp_file_mat, path] = uigetfile;
elseif strcmp(answer_analize, 'Load from Temporal Data')
    if strcmp(ComputerName, 'asuslaptop')
        cd('C:\Users\jongy\My Drive\making paper\Synapse\temporal_analysis');                       % Laptop change CD
    elseif strcmp(ComputerName, 'desktop-d5sp5os')
        cd 'H:\My Drive\making paper\Synapse\temporal_analysis'
    end
    [Name_Exp_file_mat, path] = uigetfile;
end

load(Name_Exp_file_mat);
Name_Exp_file = erase(Name_Exp_file_mat, '.mat');
temp_Exp = eval(Name_Exp_file);
denoised = temp_Exp.info.denoised;
if isfield(temp_Exp.info, 'threshEv')
    threshEv = temp_Exp.info.threshEv;
else
    threshEv = str2double(inputdlg("There is no threEv. Enter the 'threshEv'", "threshEv", [1, 45], {'3.5'}));
end

%% load Stuff data
for i = 1:length(temp_Exp.info.date)
    temp_Load_data_info = temp_Exp.info.date{i};
    temp_Load_data_info = strsplit(temp_Load_data_info, '_');
    temp_Load_data_date = strjoin(temp_Load_data_info(2:4), '.');                                           % '01.31.2021'
    temp_Load_name = strjoin({'StuffLoc', temp_Load_data_date, temp_Load_data_info{5}, temp_Load_data_info{6}}, '_');
    cd(['F:\DATA\Matlab\analysis\', temp_Load_data_info{4}, '\',  temp_Load_data_date]);
    % change denoised folder
    if denoised
        cd denoised
    end

    if isfile([temp_Load_name, '.mat'])
        load([temp_Load_name, '.mat']);
    elseif isfile([replace(temp_Load_name, 'StuffLoc', 'ROI'), '.mat'])
        warndlg(['There is no ', temp_Load_name, ', ', newline, 'but there is ', ...
            replace(temp_Load_name, 'StuffLoc', 'ROI'), newline, 'There is no problem'], 'No StuffLoc file');
    else
        errordlg('There is no corresponding files ("StuffLoc" and "ROI")','No files');
        return
    end

end
% end
variables = who;
variables_StuffLoc = variables(contains(variables, 'StuffLoc'));
variables_StuffLoc = variables_StuffLoc(contains(variables_StuffLoc, 'Data'));
Num_variables_StuffLoc = length(variables_StuffLoc);

%% SETTING THRESHOLD
diameter = [];                                                                                                                                  % SET Divide diameter for relationship between endocytosis and distance from center
ColorLinear = {'#D95319', '#0072BD', '#77AC30', '#EDB120', '#7E2F8E', 'b','#A2142F', '#4DBEEE', 'k','m', 'r','g','y','c'};                            % SET Color
nbins = (0:0.2:2);
temp_Prob = [];
temp_std = [];
temp_area = [];
n = 1; m = 1;
% Synapse = struct('Cent2distal', [], 'NumClusters', []);
Synapse(200).Cent2distal= [];
Synapse(200).NumClusters = [];
NumClusters = [];
NumClusters_less22 = [];
NumClusters_than2 = [];
Num_Cent2distal = 0;                                                                        %
ReusedTime = [];
Cluster2Cluster = [];
Center2Cluster = [];
Mean_Cent2distal = [];
Dist_Event2Event = [];                                                              % distance between Event2Event
Num_event_in_cluster = [];                                                          % events number in single cluster
ratio_more_2event = [];                                                             % ratio of clusters having multiple events
dist_AZ2cluster = [];                                                               
if isfield(temp_Exp.Period, 'Dist2Marker') ||isfield(temp_Exp.Period, 'Dist2Marker_before') || isfield(temp_Exp.Period, 'Dist2Marker_after')
    Dist_Event2Event_marker = [];                                                              % distance between Event2Event
    Dist_Event2Event_before_near = [];                                                              % distance between Event2Event
    Dist_Event2Event_before_far = [];                                                              % distance between Event2Event
    NumClusters_before_near = [];                                                              % distance between Event2Event
    NumClusters_before_far = [];                                                              % distance between Event2Event
    Dist_Event2Event_after_near = [];                                                              % distance between Event2Event
    Dist_Event2Event_after_far = [];                                                              % distance between Event2Event
    Marker_threshold = inputdlg('Enter Marker distance threshold', 'Marker Threshold', [1, 35], {'0.5'});
    Marker_threshold = str2double(Marker_threshold{1});
end

ii = 1;
%% EXPORT DATA
for j = 1:Num_variables_StuffLoc
    temp_StuffLoc = eval(variables_StuffLoc{j});                                    % assign one StuffLoc file to temp_StuffLo
    % check denoised and threshold
    if  temp_StuffLoc.InfoExp.threshEv ~= threshEv || temp_StuffLoc.InfoExp.denoised ~= denoised
        fprintf('Un-maching denoised or threEV\n(Data name: %s)\n', variables_StuffLoc{j});
        return
    end
    % define pixel
    if j == 1
        PixelSize = temp_StuffLoc.InfoExp.PixelSize;
    end

    %% GETTING DATA
    % getting probability
    temp_Prob = [temp_Prob, temp_StuffLoc.perSyn.p];                                % 'p' is probability
    % getting PSD
    temp_std = [temp_std, [temp_StuffLoc.allEv.loc_std]*PixelSize];
    % getting area
    temp_area = [temp_area, [temp_StuffLoc.perSyn.area]];
    for k = 1:length(temp_StuffLoc.allEv)
        % getting success distance
        if  k ~= length(temp_StuffLoc.allEv) && temp_StuffLoc.allEv(k).synNum == temp_StuffLoc.allEv(k+1).synNum
            temp_dist = [temp_StuffLoc.allEv(k).x,temp_StuffLoc.allEv(k).y; temp_StuffLoc.allEv(k+1).x,temp_StuffLoc.allEv(k+1).y];
            dist(n,m) = pdist(temp_dist) * PixelSize;
            n = n+1;
        else
            m = m+1;
            n = 1;
        end
        temp_synNum = temp_StuffLoc.allEv(k).synNum;
        temp_eventN = temp_StuffLoc.allEv(k).EventN;
        MVR = F_idx_MVR(variables_StuffLoc{j}, temp_StuffLoc.allEv(k).synNum,...
            temp_StuffLoc.allEv(k).EventN, temp_Exp);
        if isnan(MVR)
            continue
        end
        temp_Cent2Event(ii).MVR = MVR;
        temp_Cent2Event(ii).centPhy_dist = temp_StuffLoc.allEv(k).dist2centPhy;          %#ok<*SAGROW> % distance from center to event
        temp_Cent2Event(ii).centF_dist = temp_StuffLoc.allEv(k).dist2centF;          % distance from center to event
        ii = ii + 1;
    end

    for o = 1: length(temp_StuffLoc.perSyn)
        % getting distance between event to event
        temp_Event2Event = temp_StuffLoc.perSyn(o).XY;
        temp_Event2Event = pdist(temp_Event2Event);
        Dist_Event2Event = [Dist_Event2Event, temp_Event2Event * PixelSize;];

        % gettind distance between center to distal ,
        temp_Cent2distal = pdist([[temp_StuffLoc.perSyn(o).centerPhyX, temp_StuffLoc.perSyn(o).centerPhyY]; temp_StuffLoc.perSyn(o).area_point])...
            *PixelSize;
        temp_Cent2distal = temp_Cent2distal(1:length(temp_StuffLoc.perSyn(o).area_point));                                            % last point is same with first point
        Synapse(temp_StuffLoc.perSyn(o).numEv).Cent2distal = [Synapse(temp_StuffLoc.perSyn(o).numEv).Cent2distal, temp_Cent2distal];
        Num_Cent2distal = Num_Cent2distal + 1;

        % getting number of cluster ((x: Event Number, y: Cluster number)
        temp_NumClusters = temp_StuffLoc.perSyn(o).numClusters;
        Synapse(temp_StuffLoc.perSyn(o).numEv).NumClusters = [Synapse(temp_StuffLoc.perSyn(o).numEv).NumClusters, temp_NumClusters];
        NumClusters = [NumClusters, temp_NumClusters];
        if temp_StuffLoc.perSyn(o).numEv <= 22
        NumClusters_less22 = [NumClusters_less22, temp_NumClusters];
        end
        % more than 2 event cluster number
        temp_NumCluster_than2 = nnz([temp_StuffLoc.perSyn(o).cluster.maxRad] ~= 0); 
        NumClusters_than2 = [NumClusters_than2, temp_NumCluster_than2];
        % getting distance of Cluster to Cluster & Center to Cluster
        temp_temp_Cluster2Cluster = [];
        temp_temp_Center2Cluster = [temp_StuffLoc.perSyn(o).centerPhyX, temp_StuffLoc.perSyn(o).centerPhyY];

        for q = 1:length(temp_StuffLoc.perSyn(o).cluster)
            temp_temp_Cluster2Cluster = [temp_temp_Cluster2Cluster; temp_StuffLoc.perSyn(o).cluster(q).centroid];
            temp_temp_Center2Cluster = [temp_temp_Center2Cluster;  temp_StuffLoc.perSyn(o).cluster(q).centroid];

            % finding reuse time
            temp_T = find(temp_StuffLoc.perSyn(o).T == q);
            temp_temp_syn = find([temp_StuffLoc.allEv.synNum] == temp_StuffLoc.perSyn(o).synNum);
            temp_synNum = [temp_StuffLoc.allEv(temp_temp_syn(temp_T)).EventN];
            if length(temp_synNum) ~= 1
                ReusedTime = [ReusedTime, [temp_synNum(2:end)- temp_synNum(1:end-1)]];
            end
            temp_Num_event_in_cluster = size(temp_StuffLoc.perSyn(o).cluster(q).loc, 1);      % number of events per cluster

            Num_event_in_cluster = [Num_event_in_cluster, temp_Num_event_in_cluster];
            AZ = [temp_StuffLoc.perSyn(o).centerPhyX, temp_StuffLoc.perSyn(o).centerPhyY];
            dist_AZ2cluster = [dist_AZ2cluster, pdist([AZ; temp_StuffLoc.perSyn(o).cluster(q).centroid]) * 86.66];            % pixel size is 86.666 nm
        end
        temp_Cluster2Cluster = pdist(temp_temp_Cluster2Cluster);
        temp_Center2Cluster = pdist(temp_temp_Center2Cluster);
        temp_Center2Cluster = temp_Center2Cluster(1:length(temp_StuffLoc.perSyn(o).cluster));
        Cluster2Cluster = [Cluster2Cluster, temp_Cluster2Cluster];
        Center2Cluster = [Center2Cluster, temp_Center2Cluster];
        temp_ratio_more_2event = nnz([temp_StuffLoc.perSyn(o).cluster.maxRad] ~= 0)/temp_StuffLoc.perSyn(o).numClusters;          % ratio of cluster having multiple events
        ratio_more_2event = [ratio_more_2event, temp_ratio_more_2event];
        % check marker
        if isfield(temp_Exp.Period, 'Dist2Marker') ||isfield(temp_Exp.Period, 'Dist2Marker_before') || isfield(temp_Exp.Period, 'Dist2Marker_after')
            File_name = strrep(variables_StuffLoc{j}, 'StuffLoc', 'ROI');
            IndexData = temp_StuffLoc.perSyn(o).synNum;
            idx = and(strcmp(File_name, {temp_Exp.Period.File_name}), IndexData == [temp_Exp.Period.IndexData]);
            if sum(idx) == 1
                if isfield(temp_Exp.Period, 'Dist2Marker_before') && temp_Exp.Period(idx).Dist2Marker_before <= Marker_threshold
                    Dist_Event2Event_before_near = [Dist_Event2Event_before_near, temp_Event2Event*temp_StuffLoc.InfoExp.PixelSize;];
                    NumClusters_before_near = [NumClusters_before_near, temp_NumClusters];
                else
                    Dist_Event2Event_before_far = [Dist_Event2Event_before_far, temp_Event2Event*temp_StuffLoc.InfoExp.PixelSize;];
                    NumClusters_before_far = [NumClusters_before_far, temp_NumClusters];
                end

                if isfield(temp_Exp.Period, 'Dist2Marker_after') && temp_Exp.Period(idx).Dist2Marker_after <= Marker_threshold
                    Dist_Event2Event_after_near = [Dist_Event2Event_after_near, temp_Event2Event*temp_StuffLoc.InfoExp.PixelSize;];
                else
                    Dist_Event2Event_after_far = [Dist_Event2Event_after_far, temp_Event2Event*temp_StuffLoc.InfoExp.PixelSize;];
                end
            end
        end

    end

end

%%

% probability
Mean_Prob = mean(temp_Prob);                                                    % mean of probability
SEM_Prob = std(temp_Prob)/sqrt(length(temp_Prob));                                                      % standard deviation of probability
Num_Prob = length(temp_Prob);
%
% PSD size
Mean_std = mean(temp_std);                                                      % mean of std
SEM_std = std(temp_std)/sqrt(length(temp_std));                                                        % standard of std
Num_std = length(temp_std);

% area size
Mean_area = mean(temp_area);
SEM_area = std(temp_area)/sqrt(length(temp_area));
Num_area = length(temp_area);

% successive distance
dist_zero = dist;
dist(dist == 0) = NaN;
Mean_Succdist = nanmean(dist, 'all');
SEM_Succdist = std(dist, 0, 'all', 'omitnan')/sqrt(nnz(dist_zero));
Mean_succdist_seq = nanmean(dist,2);
Std_succdist_seq = std(dist,0,2,'omitnan');
Num_Succdist = size(dist,2);

% mass center to distal
for p = 1:200
    %Mean_Cent2distal_seq(p) = mean(Synapse(p).Cent2distal);
    Synapse(p).Mean_Cent2Distal = mean(Synapse(p).Cent2distal);
    %Std_Cent2distal_seq(p) = std(Synapse(p).Cent2distal);
    Synapse(p).Std_Cent2Distal = std(Synapse(p).Cent2distal);
    Mean_Cent2distal = [Mean_Cent2distal, Synapse(p).Cent2distal];

    % num of cluster
    Synapse(p).Mean_NumClusters = mean(Synapse(p).NumClusters);
    Synapse(p).Std_NumClusters = std(Synapse(p).NumClusters);
    Synapse(p).Num_Numclusters = length(Synapse(p).NumClusters);
end
Mean_Cent2distal = mean(Mean_Cent2distal);
Std_Cent2distal = std(Mean_Cent2distal);

% Cluster2Cluster
Cluster2Cluster = Cluster2Cluster*temp_StuffLoc.InfoExp.PixelSize;
Num_Cluster2Cluster = length(Cluster2Cluster);
Mean_Cluster2Cluster = mean(Cluster2Cluster);
Std_Cluster2Cluster = std(Cluster2Cluster);

% Center2Cluster
Center2Cluster = Center2Cluster*temp_StuffLoc.InfoExp.PixelSize;
Num_Center2Cluster = length(Center2Cluster);
Mean_Center2Cluster = mean(Center2Cluster);
Std_Center2Cluster = std(Center2Cluster);

% Reused Time
Num_ReusedTime = length(ReusedTime);
Mean_ReusedTime = mean(ReusedTime);
Std_ReusedTime = std(ReusedTime);

% exponential fit" vs "distance between Cent2Event
Cent2Event = temp_Cent2Event;
% IgnoreRow = find(([Cent2Event.roi_traceF_fitB] <= min_tau) | (max_tau < [Cent2Event.roi_traceF_fitB]));                                      % find positive B value
% Cent2Event(IgnoreRow) =  [];
Median_dist_centPhy = median([Cent2Event.centPhy_dist]);
Max_dist_centPhy = max([Cent2Event.centPhy_dist]);
% diameter = ([0,Median_dist,Max_dist]);
if isempty(diameter) == 0
    for s = 1:length(diameter)
        if s == length(diameter)
            %        sorted_temp_roi_trace_fitA = [Cent2Event(find(diameter(s) < [Cent2Event.dist])).roi_trace_fitA];         %ss  = 1:length(Cent2Event)
            sorted_temp_roi_traceF_fitA = [Cent2Event(find(diameter(s) < [Cent2Event.dist])).roi_traceF_fitA];         %#ok<*FNDSB> %ss  = 1:length(Cent2Event)
            %        sorted_temp_roi_trace_fitB = [Cent2Event(find(diameter(s) < [Cent2Event.dist])).roi_trace_fitB];         %ss  = 1:length(Cent2Event)
            sorted_temp_roi_traceF_fitB = [Cent2Event(find(diameter(s) < [Cent2Event.dist])).roi_traceF_fitB];
            sorted_temp_N_roi_trace = [Cent2Event(find(diameter(s) < [Cent2Event.dist])).roi_N_trace];
            sorted_temp_N_roi_traceF = [Cent2Event(find(diameter(s) < [Cent2Event.dist])).roi_N_traceF];
            sorted_temp_EventNum = length(sorted_temp_roi_traceF_fitB);
            temp_diameter = sprintf([' ', num2str(diameter(s)), ' <  Diameter']);
        else
            %        sorted_temp_roi_trace_fitA = [Cent2Event(find(diameter(s) < [Cent2Event.dist] & [Cent2Event.dist] <= diameter(s+1))).roi_trace_fitA];         %ss  = 1:length(Cent2Event)
            sorted_temp_roi_traceF_fitA = [Cent2Event(find(diameter(s)< [Cent2Event.dist] & [Cent2Event.dist] <= diameter(s+1))).roi_traceF_fitA];         %ss  = 1:length(Cent2Event)
            %        sorted_temp_roi_trace_fitB = [Cent2Event(find(diameter(s)< [Cent2Event.dist] & [Cent2Event.dist] <= diameter(s+1))).roi_trace_fitB];         %ss  = 1:length(Cent2Event)
            sorted_temp_roi_traceF_fitB = [Cent2Event(find(diameter(s)< [Cent2Event.dist] & [Cent2Event.dist] <= diameter(s+1))).roi_traceF_fitB];
            sorted_temp_N_roi_trace = [Cent2Event(find(diameter(s)< [Cent2Event.dist] & [Cent2Event.dist] <= diameter(s+1))).roi_N_trace];
            sorted_temp_N_roi_traceF = [Cent2Event(find(diameter(s)< [Cent2Event.dist] & [Cent2Event.dist] <= diameter(s+1))).roi_N_traceF];
            sorted_temp_EventNum = length(sorted_temp_roi_traceF_fitB);
            temp_diameter = sprintf([' ', num2str(diameter(s)), ' < Diameter ≤ ' , num2str(diameter(s+1))]);
        end
        Diameter(s).criteria = temp_diameter;
        %Diameter(s).sorted_temp_roi_trace_fitA = sorted_temp_roi_trace_fitA;
        Diameter(s).sorted_temp_roi_traceF_fitA = sorted_temp_roi_traceF_fitA;
        %Diameter(s).sorted_temp_roi_trace_fitB = sorted_temp_roi_trace_fitB;
        Diameter(s).sorted_temp_roi_traceF_fitB = sorted_temp_roi_traceF_fitB;
        Diameter(s).sorted_temp_N_roi_trace = mean(sorted_temp_N_roi_trace,2);
        Diameter(s).sorted_temp_N_roi_traceF = mean(sorted_temp_N_roi_traceF,2);
        Diameter(s).EventNum = sorted_temp_EventNum;
        Diameter(s).N_sorted_roi_trace_fitB = mean(Diameter(s).sorted_temp_roi_traceF_fitB);

    end
end

%% Load Control
cd('H:\My Drive\making paper\Synapse\spatial_analysis');                                                 % change CD
Control_name = split(Name_Exp_file, '_');
Control_name = ['Control_', Control_name{end}];
load(Control_name);
Control = eval(Control_name);


%% Draw graphs
%(Probability)
fg_1 = figure(1);
set(fg_1, 'units', 'inches', 'position', [0.5, 15.5, 2.5, 2.5], 'color', 'w');
ax_1 = axes(fg_1);
his_Pro = histogram(ax_1, temp_Prob, [0.025:0.025:0.3], 'Normalization', 'probability', 'facecolor', [0.5, 0.5, 0.5], 'edgecolor', 'w');                                                 % draw graph
set(ax_1, 'units', 'inches', 'Xlim', [0,  0.3], 'ylim', [0, 0.4], 'fontsize', 10, 'box', 'off', 'color', 'none');
xticks(ax_1, [0:0.1:1]);
txt_Prob = strcat([num2str(Mean_Prob,3), ' ± ', num2str(SEM_Prob,3), ', n = ', num2str(Num_Prob)]);
disp(txt_Prob);
xlabel('Probability', 'fontsize', 12);
ylabel('norm.count', 'fontsize', 12);
fg_1_title = 'Probability';

% Draw graphs (size single-particle)
fg_2 = figure(2);
set(fg_2, 'units', 'inches', 'position', [3.5, 15.5, 2.5, 2.5], 'color', 'w');
ax_2 = axes(fg_2);
his_PSD = histogram(ax_2, temp_std, [5:2:50], 'Normalization', 'probability', 'facecolor', [0.5, 0.5, 0.5], 'edgecolor', 'w');
set(ax_2, 'units', 'inches', 'xlim', [0, 40], 'ylim', [0, 0.15], 'fontsize', 10, 'box', 'off', 'color', 'none');
xticks(ax_2, [0:10:50]);
yticks(ax_2, [0:0.1:1]);
txt_std = strcat([num2str(Mean_std,3), ' ± ', num2str(SEM_std,3), 'nm, n = ', num2str(Num_std)]);
disp(txt_std);
xlabel('Δx (nm)', 'fontsize', 12);
ylabel('norm.ount', 'fontsize', 12);
fg_2_title = 'PSD';

% Draw histagram graph (area)
fg_3 = figure(3);
set(fg_3, 'units', 'inches', 'position', [6.5, 15.5, 2.5, 2.5], 'color', 'w');
ax_3 = axes(fg_3);
his_Area = histogram(ax_3, temp_area, [0:0.02:0.3], 'Normalization', 'Probability', 'facecolor', [0.5, 0.5, 0.5], 'edgecolor', 'w');
set(ax_3, 'units', 'inches', 'xlim', [0 0.3], 'ylim', [0, 0.2], 'fontsize', 10, 'box', 'off', 'color', 'none');
xticks(ax_3, [0:0.1:0.6]);
yticks(ax_3, [0:0.1:1]);
xlabel(ax_3, 'area size (µm^{2})', 'fontsize', 12);
ylabel(ax_3, 'norm. count', 'fontsize', 12);
fg_3_title = 'his_area size';
txt_Area = strcat([num2str(Mean_area,2), ' ± ', num2str(SEM_area,3), ', n = ', num2str(Num_area)]);
disp(txt_Area);

fg_301 = figure(301);
set(fg_301, 'units', 'inches', 'position', [6.5, 19, 2.5, 2.5], 'color', 'w');
ax_301 = axes(fg_301);
bar(ax_301, 1, mean(temp_area), 'facecolor', 'k', 'edgecolor', 'none');
hold(ax_301, 'on');
errorbar(ax_301, 1, mean(temp_area), SEM(temp_area), 'k.', 'markersize', 1, 'capsize', 8);
set(ax_301, 'units', 'inches', 'ylim', [0, 0.1], 'tickdir', 'both', 'box', 'off', 'color', 'none');
yticks(ax_301, [0:0.05:0.1]);
yticklabels(ax_301, {'0.00', '0.05', '0.10'});
ylabel(ax_301, 'area size (µm^{2})', 'fontsize', 12);
fg_301_title = 'area size';

% Draw graphs (success distance)
fg_4 = figure(4);
set(fg_4, 'units', 'inches', 'position', [9.5, 15, 2.5, 2.5], 'color', 'w');
ax_4 = axes(fg_4);
errorbar(ax_4, Mean_succdist_seq, Std_succdist_seq);
set(ax_4, 'units', 'inches', 'xlim', [0 40], 'ylim', [0, 400], 'fontsize', 10, 'box', 'off', 'color', 'none');
xlabel(ax_4, 'times (s)');
ylabel(ax_4, 'distance (nm)');
fg_4_title = 'Succeed distance';
txt_Succdist = strcat([num2str(Mean_Succdist,'%.1f'), ' ± ', num2str(SEM_Succdist,'%.1f'), 'nm, n = ', num2str(Num_Cent2distal)]);

% Draw graphs (center to distal)
fg_5 = figure(5);
set(fg_5, 'units', 'inches', 'position', [12.5, 15.5, 2.5, 2.5], 'color', 'w');
ax_5 = axes(fg_5);
errorbar(ax_5, [10:30], [Synapse(10:30).Mean_Cent2Distal], [Synapse(10:30).Std_Cent2Distal]);
set(ax_5, 'units', 'inches', 'xlim', [9, 31], 'ylim', [0, 400], 'fontsize', 10, 'box', 'off', 'color', 'none');
xlabel(ax_5, 'events number', 'fontsize', 12);
ylabel(ax_5, 'distance (nm)');
text(ax_5, 5, 1.1 * max([Synapse.Mean_Cent2Distal] + [Synapse.Std_Cent2Distal]), txt_std, 'fontsize', 15);
fg_5_title ='center2distal'; 
txt_std = strcat([num2str(Mean_Cent2distal,'%.1f'), ' ± ', num2str(Std_Cent2distal,'%.1f'), 'nm, n = ', num2str(Num_Cent2distal)]);

% events vs numbcluster
fg_6 = figure(6);
set(fg_6, 'units', 'inches', 'position',  [15.5, 15.5, 2.5, 2.5], 'color', 'w');
ax_6 = axes(fg_6);
errorbar(ax_6, [5:30],[Synapse([5:30]).Mean_NumClusters], [Synapse([5:30]).Std_NumClusters]);
set(ax_6, 'units', 'inches', 'xlim', [9, 31], 'ylim', [0, 20], 'fontsize', 10, 'box', 'off', 'color', 'none');
xticks(ax_6, [10:10:30]);
yticks(ax_6, [0:10:20]);
xlabel(ax_6, 'events number', 'fontsize', 12);
ylabel(ax_6, 'release site number')
fg_6_title = 'errorbar_NumClusters';

% mean NumClusters
fg_601 = figure(601);
set(fg_601, 'units', 'inches', 'position', [18.5, 15.5, 2.5, 2.5], 'color', 'w');
ax_601 = axes(fg_601);
bar(ax_601, 1, mean(NumClusters));
hold(ax_601, 'on');
errorbar(ax_601, 1, mean(NumClusters), SEM(NumClusters), 'k.', 'markersize', 1, 'capsize', 8);
set(ax_601, 'units', 'inches', 'ylim', [0, 15], 'TickDir', 'both', 'fontsize', 10, 'box', 'off', 'color', 'none');
yticks(ax_601, [0:5:20]);
ylabel(ax_601, 'release site number', 'fontsize', 12);
fg_601_title = 'release site number';

% mean NumClusters <22 events
fg_602 = figure(602);
set(fg_602, 'units', 'inches', 'position', [18.5, 18.5, 2.5, 2.5], 'color', 'w');
ax_602 = axes(fg_602);
bar(ax_602, 1, mean(NumClusters_less22));
hold(ax_602, 'on');
errorbar(ax_602, 1, mean(NumClusters_less22), SEM(NumClusters_less22), 'k.', 'markersize', 1, 'capsize', 8);
set(ax_602, 'units', 'inches', 'ylim', [0, 15], 'TickDir', 'both', 'fontsize', 10, 'box', 'off', 'color', 'none');
yticks(ax_602, [0:5:20]);
ylabel(ax_602, 'release site number', 'fontsize', 12);
fg_602_title = 'release site number_less_22';

% mean number of cluster having more than 2 events;
fg_603 = figure(603);
set(fg_603, 'units', 'inches', 'position', [21.5, 18.5, 2.5, 2.5], 'color', 'w');
ax_603 = axes(fg_603);
bar(ax_603, 1, mean(NumClusters_than2));
hold(ax_603, 'on');
errorbar(ax_603, 1, mean(NumClusters_than2), SEM(NumClusters_than2), 'k.', 'markersize', 1, 'capsize', 8);
set(ax_603, 'units', 'inches', 'ylim', [0, 15], 'TickDir', 'both', 'fontsize', 10, 'box', 'off', 'color', 'none');
yticks(ax_603, [0:5:20]);
ylabel(ax_603, 'release site number than 2 events', 'fontsize', 12);
fg_603_title = 'release site number than 2 events';

% Draw graphs (Cluster2Cluster)
fg_7 = figure(7);
set(fg_7, 'position', [510*6 + 10, 490*1, 500, 500]);
ax_7 = axes(fg_7);
his_Clu2Clu = histogram(ax_7, Cluster2Cluster, [0:50:500]);
set(ax_7, 'fontsize', 20, 'xlim', [0 510], 'ylim', [0, round(max(his_Clu2Clu.Values)/2, 1, 'significant')*3]);
txt_Cluster2Cluster = strcat([num2str(Mean_Cluster2Cluster,'%.1f'), ' ± ', num2str(Std_Cluster2Cluster,'%.1f'),...
    'nm, n = ', num2str(Num_Cluster2Cluster)]);
text(ax_7, 20, 1.1 * max(his_Clu2Clu.Values), txt_Cluster2Cluster);
disp(['Cluster2Cluster2: ',txt_Cluster2Cluster]);
title('Cluster2Cluster2');

% ceter to clusters
fg_8 = figure(8);
set(fg_8, 'units', 'inches', 'position', [21.5, 15.5, 2.5, 2.5], 'color', 'w');
ax_8 = axes(fg_8);
his_Cen2Clu = histogram(ax_8, Center2Cluster, [0:50:500], 'normalization', 'probability', 'facecolor', [0.5, 0.5, 0.5], 'edgecolor', 'w');
set(ax_8, 'units', 'inches', 'xlim', [0 500], 'ylim', [0, 0.3001], 'fontsize', 10, 'box', 'off', 'color', 'none');
yticks(ax_8, [0:0.1:1]);
yticklabels(ax_8, {'0.0', '0.1', '0.2', '0.3'});
xlabel(ax_8, 'dist cluster to cluster', 'fontsize', 12);
ylabel(ax_8, 'norm.count', 'fontsize', 12);
txt_Center2Cluster = strcat(['Center2Cluster2: ', num2str(Mean_Center2Cluster,'%.1f'), ' ± ', num2str(Std_Center2Cluster,'%.1f'),...
    'nm, n = ', num2str(Num_Center2Cluster)]);
disp(txt_Center2Cluster);
fg_8_title = 'dis_cluster2cluster';

% Draw graphs (Reused time)
fg_9 = figure(9);
set(fg_9, 'units', 'inches', 'position', [24.5, 15.5, 2.5, 2.5], 'color', 'w');
ax_9 = axes(fg_9);
RUT = histogram(ax_9, ReusedTime, [0:15:150], 'Normalization', 'Probability');
set(ax_9, 'fontsize', 20, 'xlim', [0 150], 'ylim', [0, 0.6]);
grid on;
txt_ReusedTime = strcat([num2str(Mean_ReusedTime,'%.1f'), ' ± ', num2str(Std_ReusedTime,'%.1f'),...
    'nm, n = ', num2str(Num_ReusedTime)]);
disp(['ReusedTime: ', txt_ReusedTime]);
text(ax_9, 10, 1.1 * max(RUT.Values), 5, txt_ReusedTime, 'fontsize', 15);
title(ax_9, 'Reused Time');

fg_101 = figure(101);
set(fg_101, 'unit', 'inch', 'position', [0.5, 19, 2.5, 2.5], 'color', 'w');
ax_101 = axes(fg_101);
histogram(ax_101, [Cent2Event.centPhy_dist], [0:40:400], 'normalization', 'probability', 'facecolor', [0.5, 0.5, 0.5], 'edgecolor', 'w');
set(ax_101, 'units', 'inches', 'xlim', [-10, 600], 'ylim', [0, 0.30001], 'fontsize', 10, 'box', 'off');
xticks(ax_101, [0:200:1000]);
yticks(ax_101, [0:0.1:1]);
yticklabels(ax_101, {'0.0', '0.1', '0.2', '0.3'});
xlabel(ax_101, 'event to AZ center (nm)', 'fontsize', 12);
ylabel(ax_101, 'norm.count', 'fontsize', 12);
fg_101_title = 'event2centPhy';

% center to event cdf graph
fg_1011 = figure(1011);
set(fg_1011, 'units', 'inches', 'position', [3.5, 19, 2.5, 2.5], 'color', 'w');
ax_1011 = axes(fg_1011);
cdfplot([Cent2Event.centPhy_dist]);
set(ax_1011, 'units', 'inches', 'xlim', [0, 400], 'xgrid', 'off', 'ygrid', 'off', 'fontsize', 10, 'box', 'off');
xticks(ax_1011, [0:200:1000]);
yticks(ax_1011, [0:0.5:1]);
yticklabels(ax_1011, {'0.0', '0.5', '1.0'});
xlabel(ax_1011, 'event to AZ center (nm)', 'fontsize', 12);
ylabel(ax_1011, 'norm.count', 'fontsize', 12);
fg_1011_title = 'event2centPhy cdf';
ax_1011.Title.Visible = 'off';

% center to event mean + SEM
fg_1012 = figure(1012);
set(fg_1012, 'units', 'inches', 'position', [6.5, 19, 2.5, 2.5], 'color', 'w');
ax_1012 = axes(fg_1012);
errorbar(ax_1012, 1, mean([Cent2Event.centPhy_dist]), SEM([Cent2Event.centPhy_dist]))
set(ax_1012, 'units', 'inches', 'fontsize', 10, 'box', 'off');
ylabel(ax_1012, 'distance (nm)', 'fontsize', 12);
title(ax_1012, 'event2centPhy average');
fg_1012_title = 'event2centPhy average';

% ceter to MVR
idx_MVR = [Cent2Event.MVR];
fg_1013 = figure(1013);
set(fg_1013, set_fg{:}, 'position', [9.5, 19, 2.5, 2.5]);
ax_1013 = axes(fg_1013);
errorbar(ax_1013, 1, mean([Cent2Event(idx_MVR).centPhy_dist]), SEM([Cent2Event(idx_MVR).centPhy_dist]))
set(ax_1013, set_ax{:});
title(ax_1013, 'event2centPhy MVR average');
fg_1013_title = 'event2centPhy MVR average';

% physiological center
fg_102 = figure(102);
set(fg_102, 'unit', 'inch', 'position', [3.5, 19, 2.5, 2.5], 'color', 'w');
ax_102 = axes(fg_102);
histogram(ax_102, [Cent2Event.centF_dist], [0:40:400], 'normalization', 'probability', 'facecolor', [0.5, 0.5, 0.5], 'edgecolor', 'w');
set(ax_102, 'units', 'inches', 'xlim', [-10, 600], 'ylim', [0, 0.30001], 'fontsize', 10, 'box', 'off');
xticks(ax_102, [0:200:1000]);
yticks(ax_102, [0:0.1:1]);
yticklabels(ax_102, {'0.0', '0.1', '0.2', '0.3'});
xlabel(ax_102, 'event to AZ center (nm)', 'fontsize', 12);
ylabel(ax_101, 'norm.count', 'fontsize', 12);
fg_102_title = 'event2centF';

% fg_102 = figure(102);
% his_Cent2Event = histogram([Cent2Event(:).centPhy_dist], [0:40:500], 'Normalization', 'probability');

% center to event cdf graph
fg_1021 = figure(1021);
set(fg_1021, 'units', 'inches', 'position', [12.5, 19, 2.5, 2.5], 'color', 'w');
ax_1021 = axes(fg_1021);
cdfplot([Cent2Event.centF_dist]);
set(ax_1021, 'units', 'inches', 'xlim', [0, 400], 'xgrid', 'off', 'ygrid', 'off', 'fontsize', 10, 'box', 'off');
xticks(ax_1021, 0:200:1000);
yticks(ax_1021, 0:0.5:1);
yticklabels(ax_1021, {'0.0', '0.5', '1.0'});
xlabel(ax_1021, 'event to AZ center (nm)', 'fontsize', 12);
ylabel(ax_1021, 'norm.count', 'fontsize', 12);
fg_1021_title = 'event2centF cdf';
ax_1021.Title.Visible = 'off';

% centerF to event mean + SEM
fg_1022 = figure(1022);
set(fg_1022, 'units', 'inches', 'position', [15, 19, 2.5, 2.5], 'color', 'w');
ax_1022 = axes(fg_1022);
errorbar(ax_1022, 1, mean([Cent2Event.centF_dist]), SEM([Cent2Event.centF_dist]))
set(ax_1022, 'units', 'inches', 'fontsize', 10, 'box', 'off');
ylabel(ax_1022, 'distance (nm)', 'fontsize', 12);
title(ax_1022, 'event2centF average');
fg_1022_title = 'event2centF average';


% Event to Event
% figure(111);
% his_Event2Event_Cont = histogram(Control.Dist_Event2Event, [0:40:600], 'Normalization', 'probability');
fg_112 = figure(112);
set(fg_112, 'units', 'inches', 'position', [0.5, 12, 2.5, 2.5], 'color', 'w');
ax_112 = axes(fg_112);
his_Event2Event = histogram(ax_112, Dist_Event2Event, [0:40:600], 'Normalization', 'probability');
his_Event2Event.FaceColor = [0.5, 0.5, 0.5];
his_Event2Event.EdgeColor = 'w';
set(ax_112, 'units', 'inches', 'xlim', [0, 400], 'ylim', [0, 0.3001], 'fontsize', 10, 'box', 'off', 'color', 'none');
xticks(ax_112, [0:200:600]);
yticks(ax_112, [0:0.1:1]);
yticklabels(ax_112, {'0.0', '0.1', '0.2', '0.3', '0.4'});
xlabel(ax_112, 'distance between events (nm)', 'fontsize', 12);
ylabel(ax_112, 'norm.count', 'fontsize', 12);
fg_112_title = 'dis_event2event';

% distance event2event cdf
fg_1121 = figure(1121);
set(fg_1121, 'units', 'inches', 'position', [9.5, 12, 2.5, 2.5], 'color', 'w');
ax_1121 = axes(fg_1121);
cdfplot(Dist_Event2Event);
set(ax_1121, 'units', 'inches', 'xlim', [0, 400], 'xgrid', 'off', 'ygrid', 'off', 'fontsize', 10, 'box', 'off');
xticks(ax_1121, [0:200:600]);
yticks(ax_1121, [0:0.5:1]);
yticklabels(ax_1121, {'0.0', '0.5', '1.0'});
xlabel(ax_1121, 'distance between events (nm)', 'fontsize', 12);
ylabel(ax_1121, 'norm.count', 'fontsize', 12);
fg_1121_title = 'dis_event2event cdf';
ax_1121.Title.Visible = 'off';

% distance event2event mean + SEM
fg_1122 = figure(1122);
set(fg_1122, 'units', 'inches', 'position', [12.5, 12, 2.5, 2.5], 'color', 'w');
ax_1122 = axes(fg_1122);
errorbar(ax_1122, 1, mean(Dist_Event2Event), SEM(Dist_Event2Event), 'k.', 'markersize', 1, 'capsize', 8);
set(ax_1122, 'units', 'inches', 'fontsize', 10, 'box', 'off');
ylabel(ax_1122, 'distance (nm)', 'fontsize', 12);
title(ax_1122, 'dist event2evet');
fg_1122_title = 'dist_event2event_mean';

% number of events in single cluster
fg_113 = figure(113);
set(fg_113, 'units', 'inches', 'position', [3.5, 12, 2.5, 2.5], 'color', 'w');
ax_113 = axes(fg_113);
histogram(ax_113, Num_event_in_cluster, [0:8], 'normalization', 'probability', 'facecolor', [0.5, 0.5, 0.5], 'edgecolor', 'w');
set(ax_113, 'units', 'inches', 'xlim', [0, 8], 'ylim', [0, 0.8], 'fontsize', 10, 'box', 'off', 'color', 'none');
xticks(ax_113, [0:3:6]);
yticks(ax_113, [0:0.4:0.8]);
yticklabels(ax_113, {'0.0', '0.4', '0.8'});
xlabel(ax_113, 'event in release site', 'fontsize', 12);
ylabel(ax_113, 'norm.count', 'fontsize', 12);
fg_113_title = 'his_event in release site';

fg_1131 = figure(1131);
set(fg_1131, 'units', 'inches', 'position', [3.5, 9, 2.5, 2.5], 'color', 'w');
ax_1131 = axes(fg_1131);
x = [1:4];
y = nan(size(x));
for i = x
    if i == x(end)
        y(i) = mean(dist_AZ2cluster(Num_event_in_cluster >= x(i)));
        Exp.AZ2cluster_Nevent(i).dist = dist_AZ2cluster(Num_event_in_cluster >= x(i));
    else
        y(i) = mean(dist_AZ2cluster(Num_event_in_cluster == x(i)));
        Exp.AZ2cluster_Nevent(i).dist = dist_AZ2cluster(Num_event_in_cluster == x(i));
    end
end
bar(ax_1131, x, y, 'facecolor', 'k', 'edgecolor', 'none');
set(ax_1131, 'units', 'inches', 'xlim', [0, 5], 'ylim', [0, 200], 'fontsize', 10, 'box', 'off', 'color', 'none');
xticklabels(ax_1131, {'0.008', '0.017', '0.025', '>0.033'});
xlabel(ax_1131, 'release probability (Pr-site)', 'fontsize', 12);
ylabel(ax_1131, 'distance to AZ center (nm)', 'fontsize', 12);
fg_1131_title = 'his_dist_RS2AZ';

% Pr in RS depend on the dist AZ2RS
fg_1132 = figure(1132);
set(fg_1132, 'units', 'inches', 'position', [3.5, 9, 2.5, 2.5], 'color', 'w');
ax_1132 = axes(fg_1132);
x = [1:4];
y = nan(size(x));
for i = x
    idx_1132 = and((i-1) * 100 <=  dist_AZ2cluster , dist_AZ2cluster < i*100);
    y(i) = mean(Num_event_in_cluster(idx_1132)/120);
    Exp.AZ2cluster_Nevent(i).Pr = Num_event_in_cluster(idx_1132)/120;
end
pl_1132 = plot(ax_1132, x-0.5, y, '.:', 'markersize', 20);
set(ax_1132, 'units', 'inches', 'xlim', [0, 4], 'ylim', [0.008, 0.02], 'fontsize', 10, 'box', 'off', 'color', 'none');
yticks(ax_1132, [0.01:0.005:0.02]);
xticklabels(ax_1132, {'0', '100', '200', '300', '400'});
yticklabels(ax_1132, {'0.010', '0.015', '0.020'});
xlabel(ax_1132, 'distance to AZ center (nm)', 'fontsize', 12);
ylabel(ax_1132, 'probability', 'fontsize', 12);
fg_1132_title = 'Pr_RS2AZ_distance';

fg_114 = figure(114);
set(fg_114, 'units', 'inches', 'position', [6.5, 12, 2.5, 2.5], 'color', 'w');
ax_114 = axes(fg_114);
bar(ax_114, 1, mean(Num_event_in_cluster), 'facecolor', 'k', 'edgecolor', 'none');
hold(ax_114, 'on');
errorbar(ax_114, 1, mean(Num_event_in_cluster), SEM(Num_event_in_cluster), 'k.', 'markersize', 1, 'capsize', 8);
set(ax_114, 'units', 'inches', 'ylim', [0, 2], 'tickdir', 'both', 'box', 'off', 'color','none');
yticks(ax_114, 0:1:2);
ylabel(ax_114, ['events number', newline, 'in a release site'], 'fontsize', 12);
fg_114_title = 'bar_event in release site';

fg_115 = figure(115);
set(fg_115, 'unit', 'inches', 'position', [6.5, 12, 2.5, 2.5], 'color', 'w');
ax_115 = axes(fg_115);
cdfplot(Num_event_in_cluster);
set(ax_115, 'units', 'inches', 'xgrid', 'off', 'ygrid', 'off', 'tickdir', 'both', 'box', 'off', 'color',' none');
xlabel(ax_115, '#event in release site', 'fontsize', 12);
ylabel(ax_115, 'norm.count', 'fontsize', 12);
fg_115_title = 'cdf_event in release site';
ax_115.Title.Visible = 'off';

% #event in release site more thatn two events
fg_116 = figure(116);
set(fg_116, 'units', 'inches', 'position', [6.5, 9, 1.5, 2.5], 'color', 'w');
ax_116 = axes(fg_116);
bar(ax_116, 1, nnz(Num_event_in_cluster >=2)/numel(Num_event_in_cluster), 'facecolor', 'k', 'edgecolor', 'none');
ylabel(ax_116, 'ratio (release site events more than 2/ release site)', 'fontsize', 12);
fg_116_title = '#event in release site more than 2';

fg_1161 = figure(1161);
set(fg_1161, set_fg{:}, 'position', [8.5, 9, 1.5, 2.5]);
ax_1161 = axes(fg_1161);
bar(ax_1161, 1, mean(ratio_more_2event), 'facecolor', 'k', 'edgecolor', 'none');
hold(ax_1161, 'on');
errorbar(ax_1161, 1, mean(ratio_more_2event), SEM(ratio_more_2event), '.k');
title(ax_1161, 'ratio cluster more than 2events');

if isfield(temp_Exp.Period, 'Dist2Marker_before')
    fg_201 = figure(201);
    set(fg_201, 'units', 'inches', 'position', [15.5, 12, 2.5, 2.5], 'color', 'w');
    ax_201 = axes(fg_201);
    hold(ax_201, 'on')
    cdf_event2event_before_near = cdfplot(Dist_Event2Event_before_near);
    cdf_event2event_before_near.Color = [1 0.7529 0];
    cdf_event2event_before_near.LineWidth = 1;
    cdf_event2event_before_far = cdfplot(Dist_Event2Event_before_far);
    cdf_event2event_before_far.Color = [0.4392, 0.6784, 0.2784];
    cdf_event2event_before_far.LineWidth = 1;
    set(ax_201, 'units', 'inches', 'xlim', [0, 600], 'xgrid', 'off', 'ygrid', 'off', 'fontsize', 10, 'box', 'off');
    xticks(ax_201, 0:300:600);
    yticks(ax_201, [0:0.5:1]);
    yticklabels(ax_201, {'0.0', '0.5', '1.0'});
    xlabel(ax_201, 'distance between events (nm)', 'fontsize', 12);
    ylabel(ax_201, 'norm.count', 'fontsize', 12);
    fg_201_title = 'event2event near marker';
    title(ax_201, fg_201_title);
    [a201, b201] = kstest2(Dist_Event2Event_before_near, Dist_Event2Event_before_far);
    display(['Event2Event by marker: ', num2str(b201)]);

    fg_202 = figure(202);
    set(fg_202, 'units', 'inches', 'position', [18.5, 12, 2.5, 2.2], 'color', 'w');
    ax_202 = axes(fg_202);
    b_202 = bar(ax_202, [1, 2], [mean(NumClusters_before_near), mean(NumClusters_before_far)], 'facecolor', 'flat', 'edgecolor', 'none', 'barwidth', 0.8);
    b_202.CData = [255, 192, 0;112, 173, 71] / 255;
    hold(ax_202, 'on');
    errorbar(ax_202, [1, 2], [mean(NumClusters_before_near), mean(NumClusters_before_far)], [SEM(NumClusters_before_near), SEM(NumClusters_before_far)], '.k', 'markersize', 1, 'capsize', 20);
    set(ax_202,'units', 'inches', 'xlim', [0.5, 2.5], 'box', 'off', 'color', 'none');
    yticks(ax_202, [0:5:10]);
    ylabel(ax_202, 'release site number', 'fontsize', 12);
    ax_202.XTick = [];
    ax_202.Position(4) = 1.6;
    [~, b_202] = ttest2(NumClusters_before_near, NumClusters_before_far);
    display(['release site number by marker: ', num2str(b_202)]);

end
%% saving figures
cd('H:\My Drive\making paper\figure');
split_Name = strsplit(Name_Exp_file, '_');
save_fig_name = strjoin(split_Name(1:end-1), '_');  
save_num = [2, 3, 301, 4, 5, 6, 601, 603, 8, 101, 1011, 1012, 1013, 102, 1021, 1022, 112, 113, 114, 115, 116, 1121, 1122, 113, 1131, 1132, 114];
save_fg = cell(length(save_num), 1);
for i = 1:length(save_num)
    n = save_num(i);
    save_fg{i} = eval(sprintf('erase(strjoin({save_fig_name, fg_%d_title}, "_"), "\");', n));
    eval(['savefig(fg_', num2str(n), ', save_fg{i});']);
end
%% saving the data
answer = questdlg('Would you like to save/overlap the data?');                         % ask whether save this data
if strcmp(answer,'Yes')
    cd('H:\My Drive\making paper\Synapse\spatial_analysis');   
    Exp.info = temp_Exp.info;
    Exp.info.date = variables_StuffLoc;
    Exp.info.name = 'stuffFrlomLoc';
    Exp.info.date = temp_Exp.info.date;
    Exp.info.denoised = denoised;
    Exp.info.threshEv = threshEv;
    Exp.info.PixelSize = PixelSize;
    Exp.Pro = temp_Prob;
    Exp.Mean = mean(Exp.Pro);
    Exp.SEM = std(Exp.Pro)/sqrt(length(Exp.Pro));
    Exp.std = temp_std;
    Exp.Area = temp_area;
    Exp.dist = dist;
    Exp.Dist_Event2Event = Dist_Event2Event;
    Exp.Cent2Event = Cent2Event;
    Exp.Mean_succdist_seq = Mean_succdist_seq;
    Exp.Std_succdist_seq = Std_succdist_seq;
    Exp.Mean_Cent2Distal = [Synapse.Mean_Cent2Distal];
    Exp.Std_Cent2Distal = [Synapse.Std_Cent2Distal];
    Exp.Mean_NumClusters = [Synapse.Mean_NumClusters];
    Exp.Std_NumClusters = [Synapse.Std_NumClusters];
    Exp.Cluster2Cluster = Cluster2Cluster;
    Exp.Center2Cluster = Center2Cluster;
    Exp.ReusedTime = ReusedTime;
    Exp.Num_event_in_cluster = Num_event_in_cluster;
    Exp.NumClusters = NumClusters;
    Exp.NumClusters_less22 = NumClusters_less22;
    Exp.ratio_more_2event = ratio_more_2event;
    eval([Name_Exp_file, ' = Exp;']);
    eval(['save(''', Name_Exp_file_mat, ''', ''', Name_Exp_file, ''');']);
end

%% finding MVR index
function idx_MVR = F_idx_MVR(File_name, synNum, EventN, temp_Exp)
    File_name = replace(File_name, 'StuffLoc', 'ROI');
    syn_idx = strcmp(File_name, {temp_Exp.Period.File_name}) & synNum == [temp_Exp.Period.IndexData];
    if sum(syn_idx) == 0
        disp(['File name is ', File_name]);
        disp(['syn num is ', num2str(synNum)]);
        disp(['eventNum is ', num2str(EventN)]);
        idx_MVR = nan;
        return
    end
    if ismember(EventN, temp_Exp.Period(syn_idx).EventN_MVR)
        idx_MVR = true;
    elseif ~ismember(EventN, temp_Exp.Period(syn_idx).EventN_MVR)
        idx_MVR = false;
    else
        fprintf('Something wrong');
    end
end