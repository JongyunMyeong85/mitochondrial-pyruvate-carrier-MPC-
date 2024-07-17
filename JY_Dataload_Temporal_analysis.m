%% data load from summarized data
% It load ROI data from name of origianl data
% You can make graphs through this code

%% Ask what type of analysis
list = {'Load & Analize', 'Analize', 'Add & Load & Analize', 'Add & new'};
[indx, ~] = listdlg('ListString', list, 'selectionMode', 'singe');
if isempty(indx)
    disp('Don''t you want to analzie?')
    return
end
answer_analize = list{indx};
% clean var
if strcmp(answer_analize, 'Analize') || strcmp(answer_analize, 'Load & Analize')
    clearvars -except answer_analize UK5099
end
close all; set(0,'defaultAxesFontName', 'arial');
ColorLinear = [[0.4, 0.4, 0.4]; [0.851, 0.3255, 0.098]; [0, 0.4471, 0.7412]; [0.4667, 0.6745, 0.1882]; [0.9294, 0.6941, 0.1225]; [0.4940 0.1840 0.5560]; [0 0 1]; [0.6350 0.0780 0.1840]; [0.3010 0.7450 0.9330]; [0 0 0]; [1 0 1]; [1 0 0]; [0 1 0]; [1 1 0]; [0 1 1]];                            % SET Color
set(0, 'defaultaxesFontName', 'arial');
set_fg = {'units', 'inches', 'color', 'w'};
set_ax = {'units', 'inches', 'fontsize', 10, 'box', 'off', 'color', 'none'};
threhold_Pro_low = 10;                                                                          % threshold for low and high probability synapse

% load summary file
if strcmp(getComputerName, 'desktop-d5sp5os')
    cd('H:\My Drive\making paper\Synapse\temporal_analysis');                                       % change CD
elseif strcmp(getComputerName, 'jongyunmyeonglabtop')
    cd('C:\Users\jongy\My Drive\making paper\Synapse\temporal_analysis');                       % Laptop change CD
end
%
if ~strcmp(answer_analize, 'Add & new')
    [Name_Exp_file_mat, path] = uigetfile;
    cd(path);
    load(Name_Exp_file_mat);
    Name_Exp_file = erase(Name_Exp_file_mat, '.mat');
    temp_Exp = eval(Name_Exp_file);                                                 % assign file contens to temp_Exp_file
    threshold = split(Name_Exp_file, '_');
    threshold = threshold{end};
    threshold = extractBefore(threshold, 'thresh');
    NumberAP = temp_Exp.info.NumberAP;
    denoised = temp_Exp.info.denoised;

    if strcmp(answer_analize, 'Load & Analize') || strcmp(answer_analize, 'Add & Load & Analize')
        for i = 1:length(temp_Exp.info.date)
            temp_Load_data_info = temp_Exp.info.date{i};
            temp_Load_data_info = strsplit(temp_Load_data_info, '_');
            temp_Load_data_date = strjoin(temp_Load_data_info(2:4), '.');                                           % '01.31.2021'
            temp_Load_name = strjoin({temp_Load_data_info{1}, temp_Load_data_date, temp_Load_data_info{5:6}}, '_');
            if strcmp(getComputerName, 'asuslaptop')
                cd(['C:\DATA\Matlab\analysis\', temp_Load_data_info{4}, '\',  temp_Load_data_date]);
            elseif strcmp(getComputerName, 'desktop-d5sp5os')
                % check denoised
                if denoised
                    cd(fullfile('F:\DATA\Matlab\analysis\', temp_Load_data_info{4}, temp_Load_data_date, 'denoised'));
                else
                    cd(fullfile('F:\DATA\Matlab\analysis\', temp_Load_data_info{4}, temp_Load_data_date));
                end
            end
            load([temp_Load_name, '.mat']);                                                 % load ROI files
        end
    end
%     answer_reload = questdlg('Would you like to reload data?');                         % ask whether save this data
else                                                                                                % making "New" file
    answer_file_name = inputdlg(["flie name", "Number AP", "threshold", "denoised(true/false)"], "Input", [1, 35], ["Control", "200", "10", "false"]);
    path = 'H:\My Drive\making paper\Synapse\temporal_analysis';
    threshold = answer_file_name{3};
    if strcmp(answer_file_name{4}, 'true')
        denoised = true;
    elseif strcmp(answer_file_name{4}, 'false')
        denoised = false;
    else
        disp("Please define the dnoised status");
        return
    end    
    
    Name_Exp_file = [answer_file_name{1}, '_', threshold, 'thresh'];
    Name_Exp_file_mat = [Name_Exp_file, '.mat'];
    NumberAP = str2num(answer_file_name{2});
    temp_Exp = struct;
end

variables = who;
Name_variables_ROI = char;                                                                      % export ROI files
variables_ROI = variables(contains(variables, 'ROI'));
variables_ROI = variables_ROI(contains(variables_ROI, 'Data'));
Num_variables_ROI = length(variables_ROI);

cd(path);
% Control_name = eval(['''Control_', threshold, 'thresh'';']);
Control_name = 'Control_10thresh';
temp_Control = load(Control_name);
temp_Control = temp_Control.(Control_name);
% Ask MVR
answer_MVR = questdlg('Do you want analyze MVR? as well?');                                     % ask MVR analysis

% threshold
threshold_5th = questdlg('What threshold for 5th fluorsence', 'threshold_5th_F', '0.9', '1', 'The other', '1');
if strcmp(threshold_5th, 'The other')
    threshold_5th = inputdlg({'What threshold for 5th fluorsence'}, 'threshold_5th_F', [1 35], {'1.01'});
    threshold_5th = threshold_5th{1};
end
threshold_5th = str2num(threshold_5th);

%%
if ~strcmp(answer_analize, 'Analize')
    %% export needed data
    for i=1:Num_variables_ROI
        Name_variables_ROI = sprintf([Name_variables_ROI,' ',variables_ROI{i}]);
    end
      
    %% Thresh holds
    max_thresh_event = NumberAP;
    min_thresh_event = str2num(extractAfter(variables_ROI{1}, 'thresh')); %#ok<*ST2NM>              % general 10 events or 6 , 5
    %% preparing variables
    p = 0;                                                                                          % period index
    All_Both = [];                                                                                  % Both
    All_Between = [];                                                                                 % Event
    All_Fail = [];                                                                                  % Fail
    All_Event = [];
    Total.Event_Number = [];
    Total.Lastframe = [];
    temp_Prob = [];                                                                                 % probabiltiy
    All_eleN_T = [];
    Sim_All_eleN_T = [];
    SynN = [];
    for i = 1:40                                                                                    % No more 40 event cluster
        Term_pre_T.("cluster" + num2str(i) + "N") = [];
        Term_post_T.("cluster" + num2str(i) + "N") = [];
    end
    %% Export data
    for i = 1:Num_variables_ROI
        temp_ROI = eval(variables_ROI{i});                                                                  % defind temporal ROI file
        %check denoised data
        if isfield(temp_ROI.InfoExp, 'denoised') && temp_ROI.InfoExp.denoised ~= denoised
            fprintf("Not match denoised data. Please check file: %s\n", variables_ROI{i});
            return
        end

        % camera exposer time
        if i == 1
            threshEv = temp_ROI.InfoExp.threshEv;
            timeframe = temp_ROI.InfoExp.timeframe;
            Frequency = temp_ROI.InfoExp.Frequency;
            interEstimulation = temp_ROI.InfoExp.interEstimulation;
            timeStartEstimulation = temp_ROI.InfoExp.timeStartEstimulation;                                 % ex 0.02 
        end
        
        for ii = length(temp_ROI.all):-1:1                                                                  % export data to "Period data"
            for iii = length(temp_ROI.all(ii).trace):-1:1                                                 
                % 5th intensity
                temp_ROI.all(ii).trace(iii).Int5th_roi_traceF = temp_ROI.all(ii).trace(iii).roi_traceF(interEstimulation)...
                    - mean(temp_ROI.all(ii).trace(iii).roi_traceF(1:interEstimulation - 1));
                temp_ROI.all(ii).trace(iii).Int5th_N_roi_traceF = temp_ROI.all(ii).trace(iii).N_roi_traceF(interEstimulation);
                temp_ROI.all(ii).trace(iii).Int5th_roi_trace = temp_ROI.all(ii).trace(iii).roi_trace(interEstimulation)...
                    - mean(temp_ROI.all(ii).trace(iii).roi_trace(1:interEstimulation-1));
                temp_ROI.all(ii).trace(iii).Int5th_N_roi_trace = temp_ROI.all(ii).trace(iii).N_roi_trace(interEstimulation);
                % apply 5th threshold.
                if  temp_ROI.all(ii).trace(iii).N_roi_traceF(timeStartEstimulation/timeframe + 1) < threshold_5th                       
                    temp_ROI.all(ii).trace(iii) = [];
                end
            end
            
            if (length(temp_ROI.all(ii).trace) > max_thresh_event) || (length(temp_ROI.all(ii).trace) < min_thresh_event)
                temp_ROI.all(ii) = [];
            end
        end
        for iiii = 1:length(temp_ROI.all)                                                                         % export data to "Period data"
            p = p+1;            
            Period(p).File_name = variables_ROI{i};
            Period(p).IndexData = temp_ROI.all(iiii).indexData;
            Period(p).EventN = [temp_ROI.all(iiii).trace.FrameN];
            Period(p).NEvent = length(temp_ROI.all(iiii).trace); %#ok<*SAGROW>
            Period(p).Between = Period(p).EventN(2:end) - Period(p).EventN(1:end-1);                                 % between
            if Period(p).EventN(end) ~= NumberAP
                Period(p).Between = [Period(p).Between, nan];
            end
            Period(p).Fail = [];
            Period(p).Prob = length(temp_ROI.all(iiii).trace)/NumberAP;                             % Probabiltiy
            temp_Prob = [temp_Prob, Period(p).Prob];
            for iiiii = 0:length(Period(p).EventN)                                                                  % Fail defind
                if (iiiii ~= 0) && (Period(p).EventN(iiiii) == NumberAP)
                else
                    if (iiiii == 0) && (Period(p).EventN(1) ~= 1)                                                          % initial fail Between
                        Period(p).Fail = Period(p).EventN(1)-1 : -1 : 1;
                    elseif (iiiii == length(Period(p).Between))
                        Period(p).Fail = [Period(p).Fail, NaN(1, NumberAP - 1 - Period(p).EventN(end))];
                    elseif (1 <= iiiii) && (Period(p).Between(iiiii) ~= 1)
                        Period(p).Fail = [Period(p).Fail, Period(p).Between(iiiii)-1 : -1 : 1];
                    end

                    if (iiiii == 0) && (Period(p).EventN(1) ~= 1)                                                           % Both
                        Period(p).Both = Period(p).EventN(1)-1:-1:1;
                    elseif (iiiii == length(Period(p).Between))
                        Period(p).Both = [Period(p).Both, NaN(1, NumberAP - Period(p).EventN(end))];
                    elseif 1 <= iiiii
                        if nnz(strcmp(fieldnames(Period), 'Both')) == 1
                            Period(p).Both = [Period(p).Both, Period(p).Between(iiiii):-1:1];
                        else
                            Period(p).Both = Period(p).Between(iiiii):-1:1;
                        end
                    end
                end
            end
            %% tempeoral cluster
            temp_clu_values = Period(p).EventN';
            Z = linkage(temp_clu_values(:,1), 'single');
            T = cluster(Z,'cutoff', 2,'criterion','distance');                     %
            jj = 1;
            % reordering
            temp_T = T;
            for j = 1:length(T)
                if j ~= length(T) && T(j) == T(j+1)
                    temp_T(j) = jj;
                else
                    temp_T(j) = jj;
                    jj = jj + 1;
                end
            end
            T = temp_T;
            Period(p).T = T;
            Period(p).Num_T = max(T);
            Period(p).eleN_T = zeros(1, Period(p).Num_T);
            for v = 1:Period(p).Num_T
                Period(p).eleN_T(v) = nnz(Period(p).T == v);
            end
            % count number in the cluster
            Period(p).T_temporal_order = [];
            for k = 1:max(T)
                loc_TT = find(T == k);
                Period(p).T_temporal_order(loc_TT) = (1:nnz(T == k));
            end
            
            % term between temporal cluster
            for jjj = 1:max(T)
                if jjj == 1
                    EventN_C = min(nnz(T == jjj), 40);                       % Event numbers in cluster
%                     eval(['Term_post_T.cluster', num2str(EventN_C),'N = [Term_post_T.cluster', num2str(EventN_C), 'N Period(p).Between(find(T == jjj, 1, ''last''))];']);
                    Term_post_T.("cluster" + num2str(EventN_C) + 'N') = [Term_post_T.("cluster" + num2str(EventN_C) + "N"), Period(p).Between(find(T == jjj, 1, 'last'))];
                elseif jjj == max(T)
                    EventN_C = min(nnz(T == jjj), 40);                       % Event numbers in cluster
                    eval(['Term_pre_T.cluster', num2str(EventN_C),'N = [Term_pre_T.cluster', num2str(EventN_C), 'N Period(p).Between(find(T == jjj, 1, ''first'')-1)];']);
                    
                else
                    EventN_C = min(nnz(T == jjj), 40);                       % Event numbers in cluster
                    eval(['Term_pre_T.cluster', num2str(EventN_C),'N = [Term_pre_T.cluster', num2str(EventN_C), 'N Period(p).Between(find(T == jjj, 1, ''first'')-1)];']);
                    eval(['Term_post_T.cluster', num2str(EventN_C),'N = [Term_post_T.cluster', num2str(EventN_C), 'N Period(p).Between(find(T == jjj, 1, ''last''))];']);
                end
            end
            Period(p).max_eleN_T = max(Period(p).eleN_T);
            All_eleN_T = [All_eleN_T, Period(p).eleN_T];
            % traces
            Period(p).trace = temp_ROI.all(iiii).trace;
            Period(p).mean_roi_trace = temp_ROI.mean_roi_trace;
            Period(p).roi_traceF = temp_ROI.all(iiii).roi_traceF;
            Period(p).roi_trace = temp_ROI.all(iiii).roi_trace;
            Period(p).Background = temp_ROI.all(iiii).Background;
            % store
            All_Both = [All_Both, Period(p).Both];
            All_Between = [All_Between, Period(p).Between];
            All_Fail = [All_Fail, Period(p).Fail];
            All_Event = [All_Event, Period(p).EventN];
            Total.Event_Number = [Total.Event_Number, length(temp_ROI.all(iiii).trace)];
            Total.Lastframe = [Total.Lastframe, temp_ROI.all(iiii).trace(end).FrameN];
            % MARKER
            if isfield(temp_ROI.all, 'Dist2Marker')
                Period(p).Dist2Marker = temp_ROI.all(iiii).Dist2Marker;
            end
            if isfield(temp_ROI.all, 'Dist2Marker_before')
                Period(p).Dist2Marker_before = temp_ROI.all(iiii).Dist2Marker_before;
            end
            if isfield(temp_ROI.all, 'Dist2Marker_after')
                Period(p).Dist2Marker_after = temp_ROI.all(iiii).Dist2Marker_after;
            end
        end
        temp_Exp.Dish(i).File_name = variables_ROI{i};
        temp_Exp.Dish(i).NEvent = arrayfun(@(x) size(x.trace, 2), temp_ROI.all);
        temp_Exp.Dish(i).Pro_ratio_low = nnz(temp_Exp.Dish(i).NEvent == threhold_Pro_low)/numel(temp_Exp.Dish(i).NEvent);
        %% number test %%
        if exist('Period', 'var')
            if ((length(Period(p).Fail) + length(Period(p).Between) ~= length(Period(p).Both) + 1) && length(Period(p).Both) ~= NumberAP -1)
                display(variables_ROI{i});
                display(p);
                break
            end
        end
    end
    
    %% After continues Event
    Continue(50).Fail = [];
    Continue(5).Success = [];
    Fewsecond_Event = cell(6,1);
    Fewsecond_Event{1} = [All_Between];                                                                    % single continue
    for k = 1:length(Period)
        for kk = 1:length(Period(k).Between)-1
            if Period(k).Between(kk) == 1
                Fewsecond_Event{2} = [Fewsecond_Event{2}, Period(k).Between(kk+1)];                             % double continue
                if (Period(k).Between(kk+1) == 1) && (length(Period(k).Between)-1 ~= kk)
                    Fewsecond_Event{3} = [Fewsecond_Event{3}, Period(k).Between(kk+2)];                             % double continue
                    if (Period(k).Between(kk+2) == 1) && (length(Period(k).Between)-2 ~= kk)
                        Fewsecond_Event{4} = [Fewsecond_Event{4}, Period(k).Between(kk+3)];
                        if (Period(k).Between(kk+3) == 1) && (length(Period(k).Between)-3 ~= kk)
                            Fewsecond_Event{5} = [Fewsecond_Event{5}, Period(k).Between(kk+3)];
                            if (Period(k).Between(kk+4) == 1) && (length(Period(k).Between)-4 ~= kk)
                                Fewsecond_Event{6} = [Fewsecond_Event{6}, Period(k).Between(kk+4)];
                            end
                        end
                    end
                end
            end
        end

        %% consicuvie Fail
        for kkk = 1:50
            temp_Continue_Fail = Period(k).Between;
            temp_Continue_Fail = temp_Continue_Fail(temp_Continue_Fail >= kkk);
            if NumberAP - Period(k).EventN(end)-kkk+1 > 0
                temp_Continue_Fail = [temp_Continue_Fail, nan] - kkk +1;
            else
                temp_Continue_Fail = temp_Continue_Fail - kkk +1;
            end
            Continue(kkk).Fail = [Continue(kkk).Fail, nnz(temp_Continue_Fail == 1)/length(temp_Continue_Fail)];
        end
%% consicuvie Success
        temp_Between = Period(k).Between;
        for kkkk = 1:10
            Continue(kkkk).Success = [Continue(kkkk).Success, nnz(temp_Between == 1)/length(temp_Between)];
            temp_Between = temp_Between([false, temp_Between(1:end-1) == 1]);
            if isempty(temp_Between)
                break
            end
        end
    end
    
    for kkk = 1:50
        Continue(kkk).Fail = rmmissing(Continue(kkk).Fail);
        Continue(kkk).Pr_Fail = mean(Continue(kkk).Fail);
        Continue(kkk).SEM_Fail = std(Continue(kkk).Fail)/sqrt(length(Continue(kkk).Fail));
        Continue(kkk).Success = rmmissing(Continue(kkk).Success);
        Continue(kkk).Pr_Success = mean(Continue(kkk).Success);
        Continue(kkk).SEM_Success = std(Continue(kkk).Success)/sqrt(length(Continue(kkk).Success));
    end
    Fewsecond_Event_Value(1) = nnz(Fewsecond_Event{1} == 1)/numel(Fewsecond_Event{1});
    Fewsecond_Event_Value(2) = nnz(Fewsecond_Event{2} == 1)/numel(Fewsecond_Event{2});
    Fewsecond_Event_Value(3) = nnz(Fewsecond_Event{3} == 1)/numel(Fewsecond_Event{3});
    Fewsecond_Event_Value(4) = nnz(Fewsecond_Event{4} == 1)/numel(Fewsecond_Event{4});
    Fewsecond_Event_Value(5) = nnz(Fewsecond_Event{5} == 1)/numel(Fewsecond_Event{5});
    Fewsecond_Event_Value(6) = nnz(Fewsecond_Event{6} == 1)/numel(Fewsecond_Event{6});
    
    %% After continues Fail
    for d = 0 : 99
        temp_Fewsecond_Fail = All_Between - d;
        Fewsecond_Fail(d+1) = nnz(temp_Fewsecond_Fail == 1)/nnz(temp_Fewsecond_Fail >= 1);
    end
    
    %% cluster vs Nevent and Number of event of burst
    TvsEventN.T = struct;                     % temporal cluster vs 
    TvsEventN.Total_eleN_T = [];                      % number of event during burst
    
    for d = 1:length(Period)
       X_T = Period(d).Num_T;
       if eval(['sum(strcmp(fieldnames(TvsEventN.T) ,''N', num2str(X_T), ''')) == 1'])
%        eval(['TvsEventN.T.N', num2str(X_T), ' = [TvsEventN.T.N', num2str(X_T), ', Period(d).NEvent];']);
       TvsEventN.T.("N" + num2str(X_T)) = [TvsEventN.T.("N" + num2str(X_T)), Period(d).NEvent];
       else
%        eval(['TvsEventN.T.N', num2str(X_T), ' = Period(d).NEvent;']);
       TvsEventN.T.("N" + num2str(X_T)) = Period(d).NEvent;

       end
       TvsEventN.Total_eleN_T = [TvsEventN.Total_eleN_T, Period(d).eleN_T];
    end
    % mean and SEM
    for dd = 1:50
        if eval(['sum(strcmp(fieldnames(TvsEventN.T), "N', num2str(dd), '")) == 1'])
            if dd == 1
                TvsEventN.mean_TvsEventN = mean(TvsEventN.T.N1);
                TvsEventN.SEM_TvsEventN = SEM(TvsEventN.T.N1);
            else
                eval(['TvsEventN.mean_TvsEventN = [TvsEventN.mean_TvsEventN, mean(TvsEventN.T.N', num2str(dd), ')];']);
                eval(['TvsEventN.SEM_TvsEventN = [TvsEventN.SEM_TvsEventN, SEM(TvsEventN.T.N', num2str(dd), ')];']);
            end
        elseif dd == 1
           TvsEventN.mean_TvsEventN = nan;
           TvsEventN.SEM_TvsEventN = nan;
        else
           TvsEventN.mean_TvsEventN = [TvsEventN.mean_TvsEventN , nan];
           TvsEventN.SEM_TvsEventN = [TvsEventN.SEM_TvsEventN , nan];
        end
    end
    
    %% Another analysis
    for l = 1:length(Period)
        Period(l).Pr_Both = length(Period(l).EventN)/NumberAP;
        Period(l).Pr_Between = nnz(Period(l).Between == 1)/length(Period(l).EventN);
        Period(l).Pr_Fail = nnz(Period(l).Fail == 1)/(NumberAP - 1 - length(Period(l).EventN));
    end
    Mean_Pr_Both = mean([Period.Pr_Both])
    SEM_Pr_Both = SEM([Period.Pr_Both]);
    Mean_Pr_Between = mean([Period.Pr_Between])
    SEM_Pr_Between = SEM([Period.Pr_Between]);
    Mean_Pr_Fail = mean([Period.Pr_Fail])
    SEM_Pr_Fail = SEM([Period.Pr_Fail]);
    SEM_Pr_Both = std([Period.Pr_Both])/sqrt(length([Period.Pr_Both]));
    SEM_Pr_Between = std([Period.Pr_Between])/sqrt(length([Period.Pr_Between]));
    
    %% simulation
    Pr_group = Total.Event_Number;
    temp_Simulation_Period_Between = []; temp_Simulation_Fail_Between = [];
    All_Simulation_Both = [];
    All_Simulation_Between = [];
    All_Simulation_Fail = [];
    All_Simulation_Event = [];
    Simulation_Period = struct('Fail', [], 'Both', []);
    for i = 1:length(Pr_group)
        temp_group = randperm(Total.Lastframe(i)-1);
        temp_group = find(temp_group <= Pr_group(i)-1);
        temp_group = [temp_group, Total.Lastframe(i)];
        Simulation_Period(i).EventN = temp_group;
        Simulation_Period(i).Between = temp_group(2:end)-temp_group(1:end-1);
        if temp_group(end) ~= NumberAP
            Simulation_Period(i).Between = [Simulation_Period(i).Between, nan];
        end
        for ii = 0:length(Simulation_Period(i).EventN)
            if (ii ~= 0) && (Simulation_Period(i).EventN(ii) == NumberAP)
            else
                if (ii == 0) && (Simulation_Period(i).EventN(1) ~= 1)                                                     % fail
                    Simulation_Period(i).Fail = [Simulation_Period(i).EventN(1)-1: -1 : 1];
                elseif ii == length(Simulation_Period(i).Between)
                    Simulation_Period(i).Fail = [Simulation_Period(i).Fail, NaN(1, NumberAP - 1 - Period(i).EventN(end))];
                elseif (1 <= ii) && (Simulation_Period(i).Between(ii) ~= 1)
                    Simulation_Period(i).Fail = [Simulation_Period(i).Fail, Simulation_Period(i).Between(ii)-1:-1:1];
                end

                if ii == 0 && (Simulation_Period(i).EventN(1) ~= 1)
                    Simulation_Period(i).Both = Simulation_Period(i).EventN(1)-1:-1:1;
                elseif ii == length(Period(i).Between)
                    Simulation_Period(i).Both = [Simulation_Period(i).Both, NaN(1, NumberAP - Simulation_Period(i).EventN(end))];
                elseif 1 <= ii
                    if nnz(strcmp(fieldnames(Simulation_Period), 'Both')) ==1
                        Simulation_Period(i).Both = [Simulation_Period(i).Both, Simulation_Period(i).Between(ii):-1:1];
                    else
                        Simulation_Period(i).Bot = Simulation(i).Between(ii):-1:1;
                    end
                end
            end

        end
        All_Simulation_Both = [All_Simulation_Both, Simulation_Period(i).Both];                                  % All Both
        All_Simulation_Between = [All_Simulation_Between, Simulation_Period(i).Between];                                  % All Both
        All_Simulation_Fail = [All_Simulation_Fail, Simulation_Period(i).Fail];                                  % All Both
        All_Simulation_Event = [All_Simulation_Event, Simulation_Period(i).EventN];                                  % All Both
        
        %% tempeoral cluster
        Sim_temp_clu_values = [];
        Sim_temp_clu_values = Simulation_Period(i).EventN';
        Sim_Z = linkage(Sim_temp_clu_values, 'single');
        Sim_T = cluster(Sim_Z,'cutoff', 2,'criterion','distance');                     %
        Simulation_Period(i).T = Sim_T;
        Simulation_Period(i).Num_T = max(Sim_T);
        Simulation_Period(i).eleN_T = zeros(1, Simulation_Period(i).Num_T);
        for v = 1:Simulation_Period(i).Num_T
            Simulation_Period(i).eleN_T(v) = nnz(Simulation_Period(i).T == v);
        end
        Simulation_Period(i).max_eleN_T = max(Simulation_Period(i).eleN_T);
        Sim_All_eleN_T = [Sim_All_eleN_T, Simulation_Period(i).eleN_T];
    end
    
    %% After continues Event
    Fewsecond_Simulation_Event = cell(6,1);
    Fewsecond_Simulation_Event{1} = [All_Simulation_Between];
    for k = 1:length(Simulation_Period)
        for kk = 1:length(Simulation_Period(k).Between)-1
            if Simulation_Period(k).Between(kk) == 1
                Fewsecond_Simulation_Event{2} = [Fewsecond_Simulation_Event{2}, Simulation_Period(k).Between(kk+1)];
                if (Simulation_Period(k).Between(kk+1) == 1) && (length(Simulation_Period(k).Between)-1 ~= kk)
                    Fewsecond_Simulation_Event{3} = [Fewsecond_Simulation_Event{3}, Simulation_Period(k).Between(kk+2)];
                    if (Simulation_Period(k).Between(kk+2) == 1) && (length(Simulation_Period(k).Between)-2 ~=kk)
                        Fewsecond_Simulation_Event{4} = [Fewsecond_Simulation_Event{4}, Simulation_Period(k).Between(kk+3)];
                        if (Simulation_Period(k).Between(kk+2) == 1) && (length(Simulation_Period(k).Between)-2 ~=kk)
                            Fewsecond_Simulation_Event{4} = [Fewsecond_Simulation_Event{4}, Simulation_Period(k).Between(kk+3)];
                            if (Simulation_Period(k).Between(kk+3) == 1) && (length(Simulation_Period(k).Between)-3 ~=kk)
                                Fewsecond_Simulation_Event{5} = [Fewsecond_Simulation_Event{5}, Simulation_Period(k).Between(kk+4)];
                                if (Simulation_Period(k).Between(kk+4) == 1) && (length(Simulation_Period(k).Between)-4 ~=kk)
                                    Fewsecond_Simulation_Event{6} = [Fewsecond_Simulation_Event{6}, Simulation_Period(k).Between(kk+5)];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    Fewsecond_Simulation_Event_Value(1) = nnz(Fewsecond_Simulation_Event{1} == 1)/numel(Fewsecond_Simulation_Event{1});
    Fewsecond_Simulation_Event_Value(2) = nnz(Fewsecond_Simulation_Event{2} == 1)/numel(Fewsecond_Simulation_Event{2});
    Fewsecond_Simulation_Event_Value(3) = nnz(Fewsecond_Simulation_Event{3} == 1)/numel(Fewsecond_Simulation_Event{3});
    Fewsecond_Simulation_Event_Value(4) = nnz(Fewsecond_Simulation_Event{4} == 1)/numel(Fewsecond_Simulation_Event{4});
    Fewsecond_Simulation_Event_Value(5) = nnz(Fewsecond_Simulation_Event{5} == 1)/numel(Fewsecond_Simulation_Event{5});
    Fewsecond_Simulation_Event_Value(6) = nnz(Fewsecond_Simulation_Event{6} == 1)/numel(Fewsecond_Simulation_Event{6});
    
    %% Continue fail
    for d = 0 : 99
        temp_Fewsecond_Simulation_Fail = All_Simulation_Between - d;
        Fewsecond_Simulation_Fail(d+1) = nnz(temp_Fewsecond_Simulation_Fail == 1)/nnz(temp_Fewsecond_Simulation_Fail >= 1);
    end
    
    %% Another analysis
for l = 1:length(Simulation_Period)
    Simulation_Period(l).Pr_Both = length(Simulation_Period(l).EventN)/NumberAP;
    if Simulation_Period(l).EventN(end) == NumberAP
        Simulation_Period(l).Pr_Between = nnz(Simulation_Period(l).Between == 1)/(length(Simulation_Period(l).EventN) - 1);
    else
        Simulation_Period(l).Pr_Between = nnz(Simulation_Period(l).Between == 1)/length(Simulation_Period(l).EventN);
    end
        Simulation_Period(l).Pr_Between = nnz(Simulation_Period(l).Between == 1)/length(Simulation_Period(l).Between);
    Simulation_Period(l).Pr_Fail = nnz(Simulation_Period(l).Fail == 1)/(NumberAP - 1 -length(Simulation_Period(l).EventN));
    
end
    Mean_Simulation_Pr_Both = mean([Simulation_Period.Pr_Both]);
    SEM_Simulation_Pr_Both = SEM([Simulation_Period.Pr_Both]);
    Mean_Simulation_Pr_Between = mean([Simulation_Period.Pr_Between]);
    SEM_Simulation_Pr_Between = SEM([Simulation_Period.Pr_Between]);
    Mean_Simulation_Pr_Fail = mean([Simulation_Period.Pr_Fail]);
    SEM_Simulation_Pr_Fail = SEM([Simulation_Period.Pr_Fail]);
    SEM_Simulation_Pr_Both = std([Simulation_Period.Pr_Both])/sqrt(length([Simulation_Period.Pr_Both]));
    SEM_Simulation_Pr_Between = std([Simulation_Period.Pr_Between])/sqrt(length([Simulation_Period.Pr_Between]));
    SEM_Simulation_Pr_Fail = std([Simulation_Period.Pr_Fail])/sqrt(length([Simulation_Period.Pr_Fail]));
    %% Both probabiltiy
    Both_Pro(1) = length(find(All_Both == 1))/length(All_Both);                                                    % 1st probabiltiy
    for b = 1:270
        temp_All_Both = All_Both;
        temp_event_loc = ceil(rand(length(Period)*b, 1) * length(All_Both));
        temp_All_Both(temp_event_loc) = [];
        Both_Pro(b+1) = length(find(temp_All_Both == 1))/length(temp_All_Both);
    end
   
    %% copy yo temp_Exp
    temp_Exp.info.date = variables_ROI;
    temp_Exp.info.Threshold = num2str(temp_ROI.threshEvNum);
    temp_Exp.info.NumberAP = NumberAP;
    temp_Exp.info.timeframe = timeframe;
    temp_Exp.info.interEstimulation = interEstimulation;
    temp_Exp.info.Frequency = temp_ROI.InfoExp.Frequency;
    temp_Exp.info.descartFrame = temp_ROI.InfoExp.descartFrame;
    temp_Exp.info.threshEv = threshEv;
    temp_Exp.info.denoised = denoised;
    temp_Exp.info.timeStartEstimulation = timeStartEstimulation;
    temp_Exp.Pro = temp_Prob;
    temp_Exp.Period = Period;
    temp_Exp.Mean.All = mean([temp_Exp.Period.Pr_Both]);
    temp_Exp.SEM.All = SEM([temp_Exp.Period.Pr_Both]);
    temp_Exp.Mean.Between = mean([temp_Exp.Period.Pr_Between]);
    temp_Exp.SEM.Between = SEM([temp_Exp.Period.Pr_Between]);
    temp_Exp.Mean.Fail = mean([temp_Exp.Period.Pr_Fail]);
    temp_Exp.SEM.Fail = SEM([temp_Exp.Period.Pr_Fail]);
    temp_Exp.Both_Pro = Both_Pro;
    temp_Exp.Simulation_Period = Simulation_Period;
    temp_Exp.All_Both = All_Both;
    temp_Exp.All_Between = All_Between;
    temp_Exp.All_Fail = All_Fail;
    temp_Exp.All_Event = All_Event;
    temp_Exp.All_Simulation_Between = All_Simulation_Between;
    temp_Exp.All_Simulation_Both = All_Simulation_Both;
    temp_Exp.All_Simulation_Event = All_Simulation_Event;
    temp_Exp.All_Simulation_Fail = All_Simulation_Fail;
    temp_Exp.Continue = Continue;
    temp_Exp.All_eleN_T = All_eleN_T;
    temp_Exp.Sim_All_eleN_T = Sim_All_eleN_T;
    temp_Exp.Term_pre_T = Term_pre_T;
    temp_Exp.Term_post_T = Term_post_T;
    temp_Exp = JY_cluster(temp_Exp);
else
    All_Both = temp_Exp.All_Both;
    All_Between = temp_Exp.All_Between;
    All_Fail = temp_Exp.All_Fail;
    All_Event = temp_Exp.All_Event;
    All_Simulation_Between = temp_Exp.All_Simulation_Between;
    All_Simulation_Both = temp_Exp.All_Simulation_Both;
    All_Simulation_Event = temp_Exp.All_Simulation_Event;
    All_Simulation_Fail = temp_Exp.All_Simulation_Fail;
    Period = temp_Exp.Period;
    Simulation_Period = temp_Exp.Simulation_Period;
    Mean_Pr_Both = temp_Exp.Mean.All;
    Mean_Pr_Between = temp_Exp.Mean.Between;
    Mean_Pr_Fail = temp_Exp.Mean.Fail;
    SEM_Pr_Both = temp_Exp.SEM.All;
    SEM_Pr_Between = temp_Exp.SEM.Between;
    SEM_Pr_Fail = temp_Exp.SEM.Fail;
    Mean_Simulation_Pr_Both = mean([temp_Exp.Simulation_Period.Pr_Both]);
    Mean_Simulation_Pr_Between = mean([temp_Exp.Simulation_Period.Pr_Between]);
    Mean_Simulation_Pr_Fail = mean([temp_Exp.Simulation_Period.Pr_Fail]);
    SEM_Simulation_Pr_Both = SEM([temp_Exp.Simulation_Period.Pr_Both]);
    SEM_Simulation_Pr_Between = SEM([temp_Exp.Simulation_Period.Pr_Between]);
    SEM_Simulation_Pr_Fail = SEM([temp_Exp.Simulation_Period.Pr_Fail]);
    timeframe = temp_Exp.info.timeframe;
    interEstimulation = temp_Exp.info.interEstimulation;
    Frequency = temp_Exp.info.Frequency;
    timeStartEstimulation = temp_Exp.info.timeStartEstimulation;
end

%% All_eleN_T
roi_traceF = [];
N_roi_traceF = [];
signal = [];                                                                                                            % signal at 5th F
noise = [];                                                                                                             % noise at 5th F
% Peak_val = nan(1, length(temp_Exp.Period) * )
for i = 1:length(temp_Exp.Period)
    % peak values at 5th fluoresence
    idx = (1:1/timeframe:NumberAP/timeframe) - 1;                                                                       % [0, 20, 40, ....]
    idx_peak = idx + interEstimulation;                                                                                 % [5, 25, 45, ....]
    idx_basal = idx + (1:interEstimulation-1)';                                                                         % [1, 2, 3, 4; 21, 22, 23, 24;...]
    temp_Peak_Val = temp_Exp.Period(i).roi_traceF(idx_peak)';                                                            %
    temp_Basal_Val = temp_Exp.Period(i).roi_traceF(idx_basal);
    temp_Basal_Val = mean(temp_Basal_Val, 1);
    temp_Peak_Val = temp_Peak_Val - temp_Basal_Val;
    idx_Event = false(1, NumberAP);
    idx_Event(temp_Exp.Period(i).EventN) = true;
    signal = [signal, temp_Peak_Val(idx_Event)]; %#ok<*AGROW>
    noise = [noise, temp_Peak_Val(~idx_Event)];
    
    temp_Exp.Period(i).Pr_Single = sum(temp_Exp.Period(i).eleN_T == 1)/temp_Exp.Period(i).NEvent;
    temp_Exp.Period(i).Pr_NEventBurst = sum(temp_Exp.Period(i).eleN_T(temp_Exp.Period(i).eleN_T ~= 1))/temp_Exp.Period(i).NEvent;
    if temp_Exp.Period(i).Pr_Single + temp_Exp.Period(i).Pr_NEventBurst ~= 1
        sprintf('somthing wrong in synapse number: %d', i);
    end
    temp_roi_traceF = [temp_Exp.Period(i).trace.roi_traceF];                   % single AP trace
    roi_traceF = [roi_traceF, temp_roi_traceF - mean(temp_roi_traceF(1:interEstimulation-1, :), 1)];
    N_roi_traceF = [N_roi_traceF, temp_roi_traceF./mean(temp_roi_traceF(1:interEstimulation-1, :), 1)];
end



%simulation
for i = 1:length(temp_Exp.Simulation_Period)
    temp_Exp.Simulation_Period(i).Pr_Single = sum(temp_Exp.Simulation_Period(i).eleN_T == 1)/length(temp_Exp.Simulation_Period(i).EventN);
    temp_Exp.Simulation_Period(i).Pr_NEventBurst = sum(temp_Exp.Simulation_Period(i).eleN_T(temp_Exp.Simulation_Period(i).eleN_T ~= 1))/length(temp_Exp.Simulation_Period(i).EventN);
    if temp_Exp.Simulation_Period(i).Pr_Single + temp_Exp.Simulation_Period(i).Pr_NEventBurst ~= 1
        fprintf('somthing wrong in synapse number: %d\n', i);
    end
end
% average and SEM
temp_Exp.Mean.Pr_Single = mean([temp_Exp.Period.Pr_Single]);
temp_Exp.Mean.Pr_NEventBurst = mean([temp_Exp.Period.Pr_NEventBurst]);
temp_Exp.Simulation_Mean.Pr_Single = mean([temp_Exp.Simulation_Period.Pr_Single]);
temp_Exp.Simulation_Mean.Pr_NEventBurst = mean([temp_Exp.Simulation_Period.Pr_NEventBurst]);
temp_Exp.SEM.Pr_Single = SEM([temp_Exp.Period.Pr_Single]);
temp_Exp.SEM.Pr_NEventBurst = SEM([temp_Exp.Period.Pr_NEventBurst]);
temp_Exp.Simulation_SEM.Pr_Single = SEM([temp_Exp.Simulation_Period.Pr_Single]);
temp_Exp.Simulation_SEM.Pr_NEventBurst = SEM([temp_Exp.Simulation_Period.Pr_NEventBurst]);

%% load control TvsEleveN
path_TvsEventN = replace(path, 'temporal_analysis', 'TvsEventN');
cd(path_TvsEventN);

% save file
if threshold_5th >= 1
    answer_3 = questdlg('Would you like to overlap ''T vs Elevent Number'' to original data?');                         % ask whether save this data
    cd(path_TvsEventN);
    if strcmp(answer_3, 'Yes')
        if exist('TvsEventN', 'var')
            eval([Name_Exp_file, ' = TvsEventN;']);
            eval(['save(''', Name_Exp_file_mat, ''', ''', Name_Exp_file, ''');']);
        else
            load(Name_Exp_file);
            TvsEventN = eval(Name_Exp_file);
        end
    else
        load(Name_Exp_file);
        TvsEventN = eval(Name_Exp_file);
    end
end

load(Control_name);
eval(['load(''Control_',threshold, 'thresh.mat'');']);
control = eval([Control_name]); %#ok<*EVLCS> 
Pro = temp_Exp.Pro;
%% draw figure
% define variables

fg_100 = figure(100);
set(fg_100, set_fg{:}, 'position', [12, 1, 1, 1.8], 'color', 'w');
ax_100 = axes(fg_100);
bar(ax_100, 1, temp_Exp.Mean.Between, 'FaceColor', 'none', 'edgecolor', 'k', 'BarWidth', 0.6);
hold(ax_100, 'on');
bar(ax_100, 2, temp_Exp.Mean.Fail, 'FaceColor', 'none', 'edgecolor', 'r', 'BarWidth', 0.6);
set(ax_100, set_ax{:}, 'xlim', [0.5, 2.5], 'Ylim', [0, 0.12], 'YTick', [0: 0.06: 0.12]);    
ylabel(ax_100, 'Probability', 'fontsize', 12);
errorbar(ax_100, 1, temp_Exp.Mean.Between, temp_Exp.SEM.Between, 'k', 'LineStyle', 'none', 'CapSize', 5);
errorbar(ax_100, 2, temp_Exp.Mean.Fail, temp_Exp.SEM.Fail, 'r', 'LineStyle', 'none', 'CapSize', 5);
line_ax_100 = refline(ax_100, 0, Mean_Pr_Both);
set(line_ax_100, 'linestyle', ':', 'color', 'k');
ax_100.XTick = [1, 2];
xticklabels(ax_100, {'event', 'failure'});
xtickangle(ax_100, 30);
yticklabels(ax_100, ["0.00", "0.06", "0.12"]);
% hold(ax_100, 'off');
ax_100.Position([3, 4]) = [0.5, 1.2];
ps_ax_100 = get(ax_100, 'position');
fg_100_title = 'Pr_event vs failure';
[~, p_val_eventfailure] = ttest2([temp_Exp.Period.Pr_Between], [temp_Exp.Period.Pr_Fail]);
fprintf('t-test between post-event and post-failure: %s\n', num2str(p_val_eventfailure));

%% real experiment
fg_1 = figure(1);
set(fg_1, set_fg{:}, 'position', [10, 12, 8, 2])
tl_1 = tiledlayout(fg_1, 1, 3);
ax_1_1 = nexttile(tl_1, 1);
his_Both = histogram(ax_1_1, temp_Exp.All_Both, [0.5:30.5],'Normalization', 'probability');
set(ax_1_1, set_ax{:}, 'xlim', [0 31], 'ylim', [0 0.25]);
xlabel(ax_1_1, 'Time interval (s)', 'FontSize', 12);
ylabel(ax_1_1, 'Norm. Count', 'FontSize', 12);
grid(ax_1_1, 'on');
text(ax_1_1, 1, his_Both.Values(1)+0.01, num2str(his_Both.Values(1),'%.3f'), 'FontSize', 12)

ax_1_2 = nexttile(tl_1, 2);
his_event_Between = histogram(ax_1_2, temp_Exp.All_Between,[0.5:30.5],'Normalization', 'probability');
set(ax_1_2, set_ax{:}, 'xlim', [0 31], 'ylim', [0 0.25]);
xlabel(ax_1_2, 'reused time (s)', 'fontsize', 12);
ylabel(ax_1_2, 'Norm. Count', 'fontsize', 12);
grid(ax_1_2, 'on');
text(ax_1_2, 1, his_event_Between.Values(1)+0.01, num2str(his_event_Between.Values(1),'%.3f'), 'FontSize', 12)

% interval after Fail Event%
ax_1_3 = nexttile(tl_1, 3);
his_fail = histogram(ax_1_3, temp_Exp.All_Fail, [0.5:30.5], 'Normalization', 'probability');
set(ax_1_3, set_ax{:}, 'xlim', [0 31], 'ylim', [0 0.25]);
xlabel(ax_1_3, 'Reuse Time (s)', 'FontSize', 12);
ylabel(ax_1_3,  'Norm. Count', 'FontSize', 12);
grid(ax_1_3, 'on');
text(ax_1_3, 1, his_fail.Values(1)+0.01, num2str(his_fail.Values(1),'%.3f'), 'FontSize', 12);

% fg_2 = figure(2);
% ax_2 = axes(fg_2);
% plot(ax_2, [1:19], temp_Exp.Both_Pro(1:19), 'ko-');
% set(ax_2, 'FontSize', 20, 'xlim', [0, 20], 'ylim', [0, 0.2]);
% grid(ax_2, 'on');
% xlabel(ax_2, 'Time (s)', 'FontSize', 20);
% ylabel(ax_2,  'Probability', 'FontSize', 20);
% xticks(ax_2, [0, 10, 20]);
% yticks(ax_2, [0, 0.1, 0.2]);

% probability all seconds
fg_4 = figure(4);
set(fg_4, set_fg{:}, 'position', [28, 12, 2, 2]);
ax_4 = axes(fg_4);
his_all = histogram(ax_4, temp_Exp.All_Event, [0.5:1:200.5]);
bar_4 = bar(ax_4, his_all.BinEdges(1:end-1) + his_all.BinWidth/2, his_all.Values/length(temp_Exp.Period), 'facecolor', 'k', 'edgecolor', 'none', 'FaceAlpha', 0.5);
set(ax_4, set_ax{:}, 'xlim', [0 201], 'ylim', [0, 0.15]);
xlabel(ax_4, 'Reuse Time (s)', 'FontSize', 12);
ylabel(ax_4,  'Norm. Count', 'FontSize', 12);
text(ax_4, 1, his_fail.Values(1)+0.01, num2str(his_fail.Values(1),'%.3f'), 'FontSize', 12);
title(ax_4, 'Pr at time')

% probability at each period
interval = 20;                                                              % interval to making average Pr
x = (0:interval:NumberAP-1) + interval/2;
y = reshape(bar_4.YData, interval, []);
y_mean = mean(y, 1);
y_err = SEM(y, 1);

fg_41 = figure(41);
set(fg_41, set_fg{:}, 'position', [31, 12, 2, 2]);
ax_41 = axes(fg_41);
errorbar(ax_41, x, y_mean, y_err, 'k');
set(ax_41, set_ax{:}, 'xlim', [0, NumberAP], 'ylim', [0, 0.15]);
xlabel(ax_41, 'time (s)', 'fontsize', 12);
ylabel(ax_41, 'probability', 'fontsize', 12);
fg_41_title = 'Pr all time';

fg_5 = figure(5);
set(fg_5, set_fg{:}, 'position', [25, 12, 2, 2]);
ax_5 = axes(fg_5);
er_success = errorbar(ax_5, [1:6], [temp_Exp.Continue(1:6).Pr_Success], [temp_Exp.Continue(1:6).SEM_Success], 'k'); %#ok<NASGU>
set(ax_5, 'FontSize', 20, 'xlim', [0.5, 6.5], 'ylim', [0, 0.4], 'box', 'off', 'color', 'none');
xticks(ax_5, [1:6]);
yticks(ax_5, [0:0.2:0.4]);
xlabel(ax_5, 'Consecutive Event', 'FontSize', 20);
ylabel(ax_5, 'Probability', 'FontSize', 20);
grid(ax_5, 'on');

% coninuous Fail
fg_6 = figure(6);
set(fg_6, 'units', 'inch', 'position', [5, 1, 1.8, 1.6], 'color', 'w');
ax_6 = axes(fg_6);
errorbar(ax_6, [1:10], [temp_Control.Continue(1:10).Pr_Fail], [temp_Control.Continue(1:10).SEM_Fail], '.-', 'Color', [1, 0.6, 0.6], 'markersize', 8);
hold(ax_6, 'on');
% er_success = errorbar(ax_6, [1:25], [temp_Exp.Continue(1:25).Pr_Fail], [temp_Exp.Continue(1:25).SEM_Fail], 'r');
set(ax_6, 'unit', 'inch', 'xlim', [0, 11], 'ylim', [0, 0.2], ...
    'fontsize', 10, 'box', 'off', 'color', 'none'); 
errorbar(ax_6, [1:10], [temp_Exp.Continue(1:10).Pr_Fail], [temp_Exp.Continue(1:10).SEM_Fail], '.-', 'Color', [1, 0, 0]);
xticks(ax_6, [0:5:25]);
xtickangle(ax_6, 0);
xlabel(ax_6, 'Consecutive Failure', 'FontSize', 12);
yticks(ax_6, [0:0.1:0.2]);
yticklabels(ax_6, {'0.0', '0.1', '0.2'});
ylabel(ax_6, 'Probability', 'FontSize', 12);
ax_6.Position(3) = 1.1;
ax_6.Position(4) = 0.7;
legend(ax_6, 'visible', 'off')

% coninuous Fail
YValues = [temp_Exp.Continue.Pr_Fail];
options = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0, 0, 0.025], 'Upper', [inf, inf, inf], 'StartPoint', [0.2, 0.2, 0.1]);
fit_type = fittype('a1*exp(-(x/b1)) + c1', ...
    'dependent', {'YValues'}, 'independent', {'x'}, ...
    'coefficients', {'a1', 'b1', 'c1'}, 'options', options);
myfit = fit((1:15)', YValues(1:15)', fit_type);            %

fg_6_1 = figure(61);
set(fg_6_1, set_fg{:}, 'position', [5, 5, 1.8, 1.6]);
ax_6_1 = axes(fg_6_1);
hold(ax_6_1, 'on');
errorbar(ax_6_1, [1:15], [temp_Exp.Continue(1:15).Pr_Fail], [temp_Exp.Continue(1:15).SEM_Fail], '.-', 'Color', [1, 0.6, 0.6]);
plot(myfit);
set(ax_6_1, set_ax{:}, 'xlim', [0, 15], 'ylim', [0, 0.2]); 
xticks(ax_6_1, [0:5:15]);
xtickangle(ax_6_1, 0);
xlabel(ax_6_1, 'Consecutive Failure/Time (s)', 'FontSize', 12);
yticks(ax_6_1, [0:0.1:0.2]);
yticklabels(ax_6_1, {'0.0', '0.1', '0.2'});
ylabel(ax_6_1, 'Probability', 'FontSize', 12);
ax_6_1.Position(3) = 1.02;
ax_6_1.Position(4) = 0.95;
legend(ax_6_1, 'visible', 'off')


ax_6_2 = axes(fg_6_1);
set(ax_6_2, 'unit', 'inch', 'XAxisLocation', 'top', 'xlim', [0, 15], 'ylim', [0, 0.2], 'ytick', [],...
    'fontsize', 10, 'box', 'off', 'color', 'none', 'position', ax_6_1.Position, 'TickDir', 'in'); 
xlabel(ax_6_2, 'Time (s)', 'fontsize', 12);


Absolute_Pr_Fail = [];
for i = 21:25
    Absolute_Pr_Fail =  [Absolute_Pr_Fail, temp_Exp.Continue(i).Fail];
end
mean(Absolute_Pr_Fail)

% coninuous Event
fg_7 = figure(7);
set(fg_7, 'units', 'inch', 'position', [15, 5, 3.7, 1.6], 'color', 'w');
ax_7 = axes(fg_7);
set(ax_7, 'unit', 'inch', 'xlim', [0.5, 7.5], 'ylim', [0, 0.3], 'fontsize', 10, 'box', 'off', 'color', 'none'); 
hold(ax_7, 'on')
er_success = errorbar(ax_7, [1:7], [temp_Exp.Continue(1:7).Pr_Success], [temp_Exp.Continue(1:7).SEM_Success], 'k.-', 'markersize', 8); %#ok<*NASGU>
xlabel(ax_7, 'Consecutive Events/Time (s)', 'FontSize', 12);
yticks(ax_7, [0:0.1:0.3]);
yticklabels(ax_7, ["0.0", "0.1", "0.2", "0.3"]);
ylabel(ax_7, 'Probability', 'FontSize', 12);
ax_7.Position(3) = 2.83;
ax_7.Position(4) = 0.95;
hold(ax_7, 'off');
fg_7_title = 'Pr Consecutive Event';

fg_71 = figure(71);
set(fg_71, set_fg{:}, 'position',  [17, 5, 3, 1.4]);
ax_71_1 = axes(fg_71);
errorbar(ax_71_1, [1:6], [temp_Exp.Continue(1:6).Pr_Success], [temp_Exp.Continue(1:6).SEM_Success], 'color', 'k');
set(ax_71_1, 'FontSize', 10, 'xlim', [0.5, 6.5], 'ylim', [0, 0.4], 'box', 'off');
ax_71_pos = ax_71_1.Position;
ax_71_2 = axes('Position', ax_71_pos, 'YAxisLocation', 'right', 'XAxisLocation', 'top',  'Color', 'none');
hold on;
errorbar([1:25], [temp_Exp.Continue(1:25).Pr_Fail], [temp_Exp.Continue(1:25).SEM_Fail], 'Parent', ax_71_2);
set(ax_71_2, 'FontSize', 10, 'xlim', [0, 26], 'ylim', [0, 0.4], 'box', 'off')

% fg_72 = figure(72);
% set(fg_72, set_fg{:}, 'position', [15, 7, 1.75, 1.4]);
% ax_7_2 = axes(fg_72);
% er_success = errorbar(ax_7_2, [2:7], [temp_Exp.Continue(2:7).Pr_Success], [temp_Exp.Continue(2:7).SEM_Success], 'k.-', 'markersize', 8);
% hold(ax_7_2, 'on');
% errorbar(ax_7_2, [2:7], [temp_Control.Continue(2:7).Pr_Success], [temp_Control.Continue(2:7).SEM_Success], '.-', 'Color', [0.6, 0.6, 0.6], 'markersize', 8);
% set(ax_7_2, 'unit', 'inch', 'xlim', [1.5, 7.5], 'ylim', [0, 0.4], 'fontsize', 10, 'box', 'off', 'color', 'none'); 
% xticks(ax_7_2, [2:7]);
% xlabel(ax_7_2, 'Time (s)', 'FontSize', 12);
% yticks(ax_7_2, [0:0.2:0.4]);
% yticklabels(ax_7_2, {'0.0', '0.2', '04'});
% ylabel(ax_7_2, 'Probability', 'FontSize', 12);
% ax_7_2.Position(3) = 1.1;
% ax_7_2.Position(4) = 0.9;
% 
% fg_73 = figure(73);
% set(fg_73, 'units', 'inch', 'position', [15, 9, 1.8, 1.4], 'color', 'w');
% ax_7_3 = axes(fg_73);
% er_success = errorbar(ax_7_3, [2:6], [temp_Exp.Continue(1:5).Pr_Success], [temp_Exp.Continue(1:5).SEM_Success], 'k.-', 'markersize', 8);
% % er_success = errorbar(ax_7_3, [2:6], [0.0887, 0.115, 0.179, 0.304, 0.374], [0.0119, 0.0178, 0.0498, 0.1054, 0.1578], 'k.-', 'markersize', 8); % DPCPX
% % er_success = errorbar(ax_7_3, [2:6], [0.0887, 0.0791, 0.1631, 0, nan], [0.0121, 0.0178, 0.07754, 0, nan], 'k.-', 'markersize', 8);   % SCH
% % er_success = errorbar(ax_7_3, [2:6], [0.0569, 0.2312, 0, nan, nan], [0.0121, 0.2312, 0, nan, nan], 'k.-', 'markersize', 8);   % SCH+MPEP
% hold(ax_7_3, 'on');
% errorbar(ax_7_3, [2:6], [temp_Control.Continue(2:6).Pr_Success], [temp_Control.Continue(2:6).SEM_Success], '.-', 'Color', [0.6, 0.6, 0.6], 'markersize', 8);
% set(ax_7_3, 'unit', 'inch', 'xlim', [1.5, 6.5], 'ylim', [0, 0.4], 'fontsize', 10, 'box', 'off', 'color', 'none'); 
% xticks(ax_7_3, [2:2:6]);
% % xlabel(ax_7_3, 'Time (s)', 'FontSize', 12);
% xlabel(ax_7_3, 'Consecutive Events', 'FontSize', 12);
% yticks(ax_7_3, [0:0.2:0.4]);
% yticklabels(ax_7_3, {'0.0', '0.2', '0.4'});
% ylabel(ax_7_3, 'Probability', 'FontSize', 12);
% ax_7_3.Position(3) = 1.1;
% ax_7_3.Position(4) = 0.9;
% line_ax_7_3 = refline(ax_7_3, 0, temp_Exp.Mean.Between);
% set(line_ax_7_3, 'linestyle', ':', 'color', 'k');

% mean culums
fg_8 = figure(8); 
set(fg_8, 'units', 'inches', 'position', [10, 5, 4, 2]);
tl_8 = tiledlayout(fg_8, 1, 2, 'TileSpacing','compact');
Name_file = replace(Name_Exp_file, '_', ' ');                           % name without '_'
% title(tl_8, Name_file, 'fontsize', 20, 'fontweight', 'bold');
ax_8_1 = nexttile(tl_8, 1);
X = categorical({'All', 'Event', 'Failure'});
b_8_1 = bar(ax_8_1, X, [Mean_Pr_Both, Mean_Pr_Between, Mean_Pr_Fail], 'BarWidth', 0.6);
title(ax_8_1, 'Experiment', 'FontWeight', 'normal');
hold(ax_8_1, 'on');
b_8_1.FaceColor = 'flat';
b_8_1.CData(1,:) = [255 255 255]/255;
b_8_1.CData(2,:) = [0, 0, 0]/255;
b_8_1.CData(3,:) = [255, 0, 0]/255;
errorbar(ax_8_1, X, [Mean_Pr_Both, Mean_Pr_Between, Mean_Pr_Fail], [SEM_Pr_Both, SEM_Pr_Between, SEM_Pr_Fail], 'k', 'LineStyle', 'none');
set(ax_8_1, 'Units', 'inches', 'Ylim', [0, 0.15], 'YTick', [0: 0.05: 0.15], ...
    'YGrid', 'off', 'fontsize', 10, 'box', 'off', 'XTickLabelRotation', 0);
ylabel(ax_8_1, 'Probability', 'fontsize', 12);
% asterisk
[~, ~, ~, c] = ANOVA1([Period.Pr_Both], [Period.Pr_Between], [Period.Pr_Fail]);
asterisks = test_asterisk(c(1, end));
text(ax_8_1, X(2), Mean_Pr_Between + SEM_Pr_Between, asterisks, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10);
asterisks = test_asterisk(c(2, end));
asterisks = compose([num2str(asterisks), '\n', num2str(test_asterisk(c(3, end)))]);
text(ax_8_1, X(3), Mean_Pr_Fail + SEM_Pr_Fail, asterisks, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10);
hold(ax_8_1, 'off');
% simulation
ax_8_2 = nexttile(tl_8, 2);
b_8_2 = bar(ax_8_2, X, [Mean_Simulation_Pr_Both, Mean_Simulation_Pr_Between, Mean_Simulation_Pr_Fail], 'BarWidth', 0.6);
title(ax_8_2, 'Virtual Synapse', 'FontWeight', 'normal');
hold(ax_8_2, 'on');
b_8_2.FaceColor = 'flat';
b_8_2.CData(1,:) = [255 255 255]/255;
b_8_2.CData(2,:) = [0, 0, 0]/255;
b_8_2.CData(3,:) = [255, 0, 0]/255;
errorbar(ax_8_2, X, [Mean_Simulation_Pr_Both, Mean_Simulation_Pr_Between, Mean_Simulation_Pr_Fail], [SEM_Simulation_Pr_Both, SEM_Simulation_Pr_Between, SEM_Simulation_Pr_Fail], 'k', 'LineStyle', 'none');
set(ax_8_2, 'Ylim', [0, 0.15], 'YTick', [0:0.05:0.15], 'YGrid', 'on', 'fontsize', 10, 'box', 'off', 'XTickLabelRotation', 0);
ylabel(ax_8_2, 'Probability', 'fontsize', 12);
% asterisk
[p, tbl, stats, c] = ANOVA1([Simulation_Period.Pr_Both], [Simulation_Period.Pr_Between], [Simulation_Period.Pr_Fail]);
asterisks = test_asterisk(c(1, end));
text(ax_8_2, X(2), Mean_Simulation_Pr_Between + SEM_Simulation_Pr_Between, asterisks, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14);
asterisks = test_asterisk(c(2, end));
asterisks = compose([num2str(asterisks), '\n', num2str(test_asterisk(c(3, end)))]);
text(ax_8_2, X(3), Mean_Simulation_Pr_Fail + SEM_Simulation_Pr_Fail, asterisks, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14);
hold(ax_8_2, 'off');

% general probability
fg_83 = figure(83);
set(fg_83, set_fg{:}, 'position', [0.5, 18, 1.5, 2]);
ax_83 = axes(fg_83);
bar(ax_83, 1, Mean_Pr_Both, 'facecolor', 'k', 'edgecolor', 'none');
hold(ax_83, 'on');
errorbar(ax_83, 1, Mean_Pr_Both, SEM_Pr_Both, '.k', 'markersize', 0.1, 'capsize', 10);
set(ax_83, 'units', 'inches', 'fontsize', 10, 'box', 'off', 'color', 'none');
ylabel(ax_83, 'probability', 'fontsize', 12);
fg_83_title = 'Probability';
title(ax_83, fg_83_title);
ax_83.Title.Visible = 'off';

% amplitudes of noise and signal
fg_84 = figure(84);
set(fg_84, 'units', 'inches', 'position', [2.5, 18, 1.5, 2], 'color', 'w');
ax_84 = axes(fg_84);
his_noise = histogram(ax_84, noise, [-20:2:20]);
set(ax_84, 'units', 'inches', 'box', 'off', 'color', 'none')

fg_85 = figure(85);
set(fg_85, 'units', 'inches', 'position', [5.5, 18, 3, 3], 'color', 'w');
ax_85 = axes(fg_85);
his_signal = histogram(ax_85, signal, [-5:0.5:20], 'normalization', 'probability');
YValues = his_signal.Values;
x = his_signal.BinEdges(1:end-1) + his_signal.BinWidth/2;
bar(ax_85, x, YValues, 'facecolor', 'k', 'edgecolor', 'none', 'FaceAlpha', 0.2, 'barwidth', 0.8);
[max_peak, Idx_peak] = max(YValues);
options = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [max_peak*0.8, 0, 0, 0, x(Idx_peak), 0], ...
    'Upper', [max_peak*1.2, x(Idx_peak+2), inf,  max_peak, max(x), inf]);
myfittype = fittype('a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c2)^2)', ...
    'dependent', {'YValues'}, 'independent', {'x'}, ...
    'coefficients', {'a1', 'b1', 'c1', 'a2', 'b2', 'c2'}, 'options', options);
myfit = fit(x', YValues', myfittype, 'StartPoint', [max_peak, 2.5791, 5.825, 0.04862, 10.5029, 5.5], 'maxfunevals', 6000);
peak_UVR = myfit.b1;
width_UVR = myfit.c1;
peak_MVR = myfit.b2;
width_MVR = myfit.c2;
hold(ax_85, 'on');
plot(ax_85, myfit, 'k');
set(ax_85, 'units', 'inches', 'xlim', [-6.5,18.5], 'ylim', [0, 0.06], 'box', 'off', 'color', 'none');
xticks(ax_85, [-6.5:5:18.5])
xticklabels(ax_85, {'-5', '0', '5', '10', '15', '20'});
xlabel(ax_85, 'Fluoresence Intensity (ΔF)', 'fontsize', 12);
ylabel(ax_85, 'number of synapse', 'fontsize', 12);
ax_85.Legend.Visible = 'off';
fg_85_title = 'Intensity Distribution';

% histogram of probability
fg_86 = figure(86);
set(fg_86, set_fg{:}, 'position',  [10, 18, 3, 2]);
ax_86 = axes(fg_86);
his = histogram(ax_86, Pro, [0.0475:0.005:1], 'Normalization', 'probability');
bar(ax_86, his.BinEdges(1:end-1) + his.BinWidth/2, his.Values, 'facecolor', 'k', 'edgecolor', 'none', 'FaceAlpha', 0.2, 'barwidth', 0.8);
set(ax_86, set_ax{:}, 'xlim', [0.04,0.25], 'ylim', [0, 0.15]);
xlabel(ax_86, 'probability', 'fontsize', 12);
ylabel(ax_86, 'norm.count', 'fontsize', 12);

Pr_ratio_low = nnz(Pro ==0.05)/numel(Pro);
fg_87 = figure(87);
set(fg_87, set_fg{:}, 'position', [13.5, 18, 1, 2]);
ax_87 = axes(fg_87);
bar(ax_87, 1, Pr_ratio_low);
set(ax_87, set_ax{:});
title(ax_87, 'Pr ratio low');
ylabel(ax_87, 'ratio', 'fontsize', 12);

fg_88 = figure(88);
set(fg_88, set_fg{:}, 'position', [15, 18, 1, 2]);
ax_88 = axes(fg_88);
bar(ax_88, 1, mean([temp_Exp.Dish.Pro_ratio_low], 'omitnan'), 'facecolor', 'k', 'edgecolor', 'none', 'facealpha', 0.2, 'barwidth', 0.8);
hold(ax_88, 'on');
errorbar(ax_88, 1, mean([temp_Exp.Dish.Pro_ratio_low], 'omitnan'), SEM([temp_Exp.Dish.Pro_ratio_low], 'omitnan'), 'k', 'capsize', 5);
set(ax_88, set_ax{:});
title(ax_88, 'Pr ratio low');
ylabel(ax_88, 'ratio', 'fontsize', 12);
fg_88_title = 'Pr ratio low';

% low and high Pr intensity histogram.
if threshold_5th == 0.9
    Int_low_Pr = [];
    Int_high_Pr = [];
    for i = 1:length(temp_Exp.Period)
        if temp_Exp.Period(i).NEvent == threhold_Pro_low                                  % low probability
            Int_low_Pr = [Int_low_Pr, temp_Exp.Period(i).trace.Int5th_roi_traceF];
        else
            Int_high_Pr = [Int_high_Pr, temp_Exp.Period(i).trace.Int5th_roi_traceF];
        end
    end
    % low Pr synapse intensity
    fg_89_1 = figure(891);
    set(fg_89_1, set_fg{:}, 'position', [16.5, 18, 3, 2]);
    ax_89_1 = axes(fg_89_1);
    his = histogram(ax_89_1, Int_low_Pr, [-5:0.5:20], 'Normalization', 'probability');
    YValues_low_Pr = his.Values;
    max_peak = max(YValues_low_Pr);
    x = his.BinEdges(1:end-1)+his.BinWidth/2;
    bar(ax_89_1, x, YValues_low_Pr, 'facecolor', 'k', 'edgecolor', 'none', 'FaceAlpha', 0.2, 'barwidth', 0.8);
    options = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [max_peak*0.8, peak_UVR, width_UVR, 0, peak_MVR, width_MVR], ...
        'Upper', [max_peak*1.2, peak_UVR, width_UVR,  max_peak, peak_MVR, width_MVR]);
    myfittype = fittype('a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c2)^2)', ...
        'dependent', {'YValues'}, 'independent', {'x'}, ...
        'coefficients', {'a1', 'b1', 'c1', 'a2', 'b2', 'c2'}, 'options', options);
    myfit_low = fit(x', YValues_low_Pr', myfittype, 'StartPoint', [max_peak, peak_UVR, width_UVR, 0.04862, peak_MVR, width_MVR], 'maxfunevals', 6000);
    hold(ax_89_1, 'on');
    plot(ax_89_1, myfit_low, 'k');
    set(ax_89_1, set_ax{:});
    xlabel(ax_89_1, 'Intensity', 'fontsize', 12);
    ylabel(ax_89_1, 'norm.count', 'fontsize', 12);
    ax_89_1.Legend.Visible = 'off';

    fg_89_2 = figure(892);
    set(fg_89_2, set_fg{:}, 'position', [20, 18, 3, 2]);
    ax_89_2 = axes(fg_89_2);
    his = histogram(ax_89_2, Int_high_Pr, [-5:0.5:20], 'Normalization', 'probability');
    YValues_high_Pr = his.Values;
    [max_peak, Idx_peak] = max(YValues_high_Pr);
    bar(ax_89_2, x, YValues_high_Pr, 'facecolor', 'r', 'edgecolor', 'none', 'FaceAlpha', 0.2, 'barwidth', 0.8);
    options = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [max_peak*0.8, peak_UVR, width_UVR, 0, peak_MVR, width_MVR], ...
        'Upper', [max_peak*1.2, peak_UVR, width_UVR,  max_peak, peak_MVR, width_MVR]);
    myfittype = fittype('a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c2)^2)', ...
        'dependent', {'YValues'}, 'independent', {'x'}, ...
        'coefficients', {'a1', 'b1', 'c1', 'a2', 'b2', 'c2'}, 'options', options);
    myfit_high = fit(x', YValues_high_Pr', myfittype, 'StartPoint', [max_peak, peak_UVR, width_UVR, 0.04862, peak_MVR, width_MVR], 'maxfunevals', 6000);
    hold(ax_89_2, 'on');
    plot(ax_89_2, myfit_high, 'r');
    set(ax_89_2, set_ax{:});
    xlabel(ax_89_2, 'Intensity', 'fontsize', 12);
    ylabel(ax_89_2, 'norm.count', 'fontsize', 12);
    ax_89_2.Legend.Visible = 'off';

    fg_89 = figure(89);
    set(fg_89, set_fg{:}, 'position', [23.5, 18, 3, 2]);
    ax_89 = axes(fg_89);
    b = bar(ax_89, x, [YValues_low_Pr; YValues_high_Pr], 'facecolor', 'flat', 'facealpha', 0.2, 'edgecolor', 'none');
    b(1).CData = [0, 0, 0];
    b(2).CData = [1, 0, 0];
    hold(ax_89, 'on');
    plot(ax_89, myfit_low, 'k');
    plot(ax_89, myfit_high, 'r');
    set(ax_89, set_ax{:});
    xlabel(ax_89, 'Intensity', 'fontsize', 12);
    ylabel(ax_89, 'norm.count', 'fontsize', 12);
    ax_89.Legend.Visible = 'off';
    fg_89_title = 'Intensity distribution in low vs high Pr';
end

%% simulation
fg_101 = figure(101);
set(fg_101, 'units', 'inches', 'position', [7, 1, 1.75, 1.5], 'color', 'w');
ax_101 = axes(fg_101);
b_101 = bar(ax_101, [Mean_Simulation_Pr_Between, Mean_Simulation_Pr_Fail], 'FaceColor', 'flat', 'CData', [0 0 0; 255, 0, 0]/255, 'BarWidth', 0.6);
set(ax_101, 'unit', 'inch', 'XTicklabel', {'Event', 'Failure'}, 'xlim', [0.5, 2.5], 'Ylim', [0, 0.15], 'YTick', [0: 0.05: 0.15], ...
    'fontsize', 10, 'box', 'off', 'outerposition', [0, 0, 1, 2], 'position', ps_ax_100, 'color', 'none');
title(ax_101, 'Virtual Synapse', 'FontWeight', 'Normal', 'fontsize', 12);
ylabel(ax_101, 'Probability', 'fontsize', 12, 'visible', 'on');
% yticklabels(ax_101, []);
xtickangle(ax_101, 0);
yticklabels(ax_101, ["0.00", "0.05", "0.10", "0.15"]);
hold(ax_101, 'on');
errorbar(ax_101, [Mean_Simulation_Pr_Between, Mean_Simulation_Pr_Fail], [SEM_Simulation_Pr_Between, SEM_Simulation_Pr_Fail], 'k', 'LineStyle', 'none', 'CapSize', 10);
line_ax_101 = refline(ax_101, 0, Mean_Simulation_Pr_Both);
set(line_ax_101, 'linestyle', ':', 'color', 'k');
% ax_101.Position(3) = 0.8;
ax_101.Position(3) = 1.;
ax_101.Position(4) = 0.9;
hold(ax_101, 'off');

%% experiment 
fg_102 = figure(102);
set(fg_102, 'units', 'inches', 'position', [10, 1, 1.75, 1.5], 'color', 'w');
ax_102 = axes(fg_102);
b_102 = bar(ax_102, [Mean_Pr_Between, Mean_Pr_Fail], 'FaceColor', 'flat', 'CData', [0 0 0; 255, 0, 0]/255, 'BarWidth', 0.6);
set(ax_102, 'unit', 'inch', 'XTicklabel', {'Event', 'Failure'}, 'xlim', [0.5, 2.5], 'Ylim', [0, 0.15], 'YTick', [0: 0.05: 0.15], ...
    'fontsize', 10, 'box', 'off', 'outerposition', [0, 0, 1, 2], 'position', ps_ax_100, 'color', 'none');
% title(ax_102, 'MTEP', 'FontWeight', 'Normal', 'fontsize', 12, 'visible', 'on');
ylabel(ax_102, 'Probability', 'fontsize', 12);
hold(ax_102, 'on');
errorbar(ax_102, [Mean_Pr_Between, Mean_Pr_Fail], [SEM_Pr_Between, SEM_Pr_Fail], 'k', 'LineStyle', 'none', 'CapSize', 10);
line_ax_102 = refline(ax_102, 0, Mean_Pr_Both);
xtickangle(ax_102, 0);
yticklabels(ax_102, ["0.00", "0.05", "0.10", "0.15"]);
ax_102.Position(3) = 1.;
ax_102.Position(4) = 0.9;
set(line_ax_102, 'linestyle', ':', 'color', 'k');
hold(ax_102, 'off');
[EF_a, EF_b] = ttest2([Period.Pr_Between], [[Period.Pr_Fail]])

% temp_Simulation_Both_Between = [temp_Simulation_Period_Between, temp_Simulation_Fail_Between];
% fg_11 = figure(11);
% ax_11 = axes(fg_11);
% his_Simulation_Both = histogram(ax_11, All_Simulation_Both, [0.5:30.5], 'Normalization', 'probability');
% his_Simulation_Both_Values = his_Simulation_Both.Values';
% 
% set(ax_11, 'FontSize', 20, 'xlim', [0 31], 'ylim', [0 0.25]);
% xlabel(ax_11, 'Reuse Time (s)', 'FontSize', 20);
% ylabel(ax_11,  'Norm. Count', 'FontSize', 20);
% grid on;
% 
% fg_12 = figure(12);
% ax_12= axes(fg_12);
% his_Simulation_Event = histogram(ax_12, All_Simulation_Between, [0.5:30.5], 'Normalization', 'probability');
% set(ax_12, 'FontSize', 20, 'xlim', [0 31], 'ylim', [0 0.25]);
% xlabel(ax_12, 'Reuse Time (s)', 'FontSize', 20);
% ylabel(ax_12,  'Norm. Count', 'FontSize', 20);
% grid on;
% 
% % temp_Simulation_Fail_Between = temp_Simulation_Fail_Between(1:length(continue_Between));
% fg_13 = figure(13);
% ax_13 = axes(fg_13);
% his_Simulation_fail = histogram(ax_13, All_Simulation_Fail, [0.5:30.5], 'Normalization', 'probability');
% set(ax_13, 'FontSize', 20, 'xlim', [0 31], 'ylim', [0 0.25]);
% xlabel(ax_13, 'Reuse Time (s)', 'FontSize', 20);
% ylabel(ax_13,  'Norm. Count', 'FontSize', 20);
% grid on;

% Probability vs # Burst
Num_Burst = nan(length(temp_Exp.Period), 1);
Num_Burst_Sim = nan(length(temp_Exp.Simulation_Period), 1);
Num_Burst_event = nan(length(temp_Exp.Period), 1);
Num_Burst_event_Sim = nan(length(temp_Exp.Simulation_Period), 1);
for i = 1:length(Num_Burst)
    Num_Burst(i) = sum(temp_Exp.Period(i).eleN_T >= 2);
    Num_Burst_Sim(i) = sum(temp_Exp.Simulation_Period(i).eleN_T >= 2);
    Num_Burst_event(i) = sum(temp_Exp.Period(i).eleN_T(temp_Exp.Period(i).eleN_T ~= 1));
    Num_Burst_event_Sim(i) = sum(temp_Exp.Simulation_Period(i).eleN_T(temp_Exp.Simulation_Period(i).eleN_T ~= 1));
end

Num_Burst_XY = nan(max([temp_Exp.Period.NEvent]), 3);
Num_Burst_XY_Sim = nan(max([temp_Exp.Period.NEvent]), 3);
Num_Burst_event_XY = nan(max([temp_Exp.Period.NEvent]), 3);
Num_Burst_event_XY_Sim = nan(max([temp_Exp.Period.NEvent]), 3);

for i = min([temp_Exp.Period.NEvent]):max([temp_Exp.Period.NEvent])
    Num_Burst_XY(i, 1) = i;
    Num_Burst_XY(i, 2) = mean(Num_Burst(i == [temp_Exp.Period.NEvent]));
    Num_Burst_XY(i, 3) = SEM(Num_Burst(i == [temp_Exp.Period.NEvent]));
    
    Num_Burst_XY_Sim(i, 1) = i;
    Num_Burst_XY_Sim(i, 2) = mean(Num_Burst_Sim(i == [temp_Exp.Period.NEvent]));
    Num_Burst_XY_Sim(i, 3) = SEM(Num_Burst_Sim(i == [temp_Exp.Period.NEvent]));
    
    Num_Burst_event_XY(i, 1) = i;
    Num_Burst_event_XY(i, 2) = mean(Num_Burst_event(i == [temp_Exp.Period.NEvent]));
    Num_Burst_event_XY(i, 3) = SEM(Num_Burst_event(i == [temp_Exp.Period.NEvent]));
    
    Num_Burst_event_XY_Sim(i, 1) = i;
    Num_Burst_event_XY_Sim(i, 2) = mean(Num_Burst_event_Sim(i == [temp_Exp.Period.NEvent]));
    Num_Burst_event_XY_Sim(i, 3) = SEM(Num_Burst_event_Sim(i == [temp_Exp.Period.NEvent]));
end

% fg_11 = figure(11);
% set(fg_11, 'unit', 'inch', 'position', [7, 1, 1.75, 1.5], 'color', 'w');
% ax_11 = axes(fg_11);
% errorbar(ax_11, Num_Burst_XY(1:30, 1), Num_Burst_XY(1:30, 2), Num_Burst_XY(1:30, 3));
% hold(ax_11, 'on');
% errorbar(ax_11, Num_Burst_XY_Sim(1:30, 1), Num_Burst_XY_Sim(1:30, 2), Num_Burst_XY_Sim(1:30, 3));
% set(ax_11, 'unit', 'inch', 'fontsize', 10, 'box', 'off');

% fg_111= figure(111);
% set(fg_111, 'unit', 'inch', 'position', [10, 1, 1.75, 1.5], 'color', 'w');
% ax_111 = axes(fg_111);
% errorbar(ax_111, Num_Burst_event_XY(1:20, 1), Num_Burst_event_XY(1:20, 2), Num_Burst_event_XY(1:20, 3), 'o-k', 'markersize', 5, 'markerfacecolor', 'k');
% hold(ax_111, 'on');
% errorbar(ax_111, Num_Burst_event_XY_Sim(1:20, 1), Num_Burst_event_XY_Sim(1:20, 2), Num_Burst_event_XY_Sim(1:20, 3), 'o-k', 'markersize', 5, 'markerfacecolor', 'w');
% set(ax_111, 'unit', 'inch', 'fontsize', 10, 'box', 'off');
% xlabel(ax_111, '#Event', 'fontsize', 12);
% ylabel(ax_111, '#Event in Burst','fontsize', 12);

% probability all seconds
fg_14 = figure(14);
ax_14 = axes(fg_14);
his_fail = histogram(ax_14, All_Simulation_Event, [0.5:1:200.5], 'Normalization', 'probability');
% his_fail_10Succ = histogram([fail_10Succ], [1:30],'Normalization', 'probability');
set(ax_14, 'FontSize', 20, 'xlim', [0 201]);
xlabel(ax_14, 'Reuse Time (s)', 'FontSize', 20);
ylabel(ax_14,  'Norm. Count', 'FontSize', 20);
grid on;
text(ax_14, 1, his_fail.Values(1)+0.01, num2str(his_fail.Values(1),'%.3f'), 'FontSize', 20);

% histogram of nubmer of event in the burst
fg_20 = figure(20);
set(fg_20, 'unit', 'inch', 'position', [20, 3, 1.75, 1.5], 'color', 'none');
ax_20 = axes(fg_20);
his_eleN_T = histogram(ax_20, temp_Exp.All_eleN_T(temp_Exp.All_eleN_T ~= 1), [2:7], 'Normalization', 'probability');
hold(ax_20, 'on');
his_Sim_eleN_T = histogram(ax_20, temp_Exp.Sim_All_eleN_T(temp_Exp.Sim_All_eleN_T ~= 1), [2:7], 'Normalization', 'probability'); %#ok<*NBRAK2> 
hold(ax_20, 'off');
bar_eleN_T = bar(ax_20, his_eleN_T.BinEdges(1:end-1), [his_eleN_T.Values; his_Sim_eleN_T.Values], 'FaceColor', 'flat');
bar_eleN_T(1).CData = repmat([0, 1, 0], 5, 1);
bar_eleN_T(2).CData = repmat([1, 1, 1], 5, 1);
set(ax_20, 'unit', 'inch','ylim', [0, 1], 'fontsize', 10, 'box', 'off', 'color', 'none');
xlabel(ax_20, '#Events in Burst', 'FontSize', 12);
yticks(ax_20, [0:0.5:1]);
yticklabels(ax_20, ["0.0", "0.5", "1.0"]);
ylabel(ax_20,  'Norm. Count', 'FontSize', 12);
yticks(ax_20, [0:0.5:1]);
% yticklabels(ax_19, {'0', '', '1'});
ax_20.Position([3,4]) = [1.1, 0.95];

fg_20_1 = figure(201);
set(fg_20_1, 'unit', 'inch', 'position', [20, 6, 1.75, 1.5], 'color', 'w');
ax_20_1 = axes(fg_20_1);
his_eleN_T = histogram(ax_20_1, temp_Exp.All_eleN_T, [1:5], 'Normalization', 'probability');
hold(ax_20_1, 'on');
his_Sim_eleN_T = histogram(ax_20_1, temp_Exp.Sim_All_eleN_T(temp_Exp.Sim_All_eleN_T), [1:5], 'Normalization', 'probability');
hold(ax_20_1, 'off');
bar_eleN_T = bar(ax_20_1, his_eleN_T.BinEdges(1:end-1), [his_eleN_T.Values; his_Sim_eleN_T.Values], 'FaceColor', 'flat');
bar_eleN_T(1).CData = repmat([0, 1, 0], 4, 1);
bar_eleN_T(2).CData = repmat([1, 1, 1], 4, 1)
set(ax_20_1, 'unit', 'inch','ylim', [0, 1], 'fontsize', 10, 'box', 'off', 'color', 'w');
xlabel(ax_20_1, '#Events in Burst', 'FontSize', 12);
ylabel(ax_20_1,  'Norm. Count', 'FontSize', 12);
yticks(ax_20_1, [0:0.5:1]);
% yticklabels(ax_19, {'0', '', '1'});
ax_20_1.Position([3,4]) = [1.1, 0.95];
[a, b] = ttest2(temp_Exp.All_eleN_T(temp_Exp.All_eleN_T ~= 1), temp_Exp.Sim_All_eleN_T(temp_Exp.Sim_All_eleN_T ~= 1))

% event number in burst
fg_20 = figure(20);
set(fg_20, 'unit', 'inch', 'position', [21, 6, 1, 1], 'color', 'w');
ax_20 = axes(fg_20);
bar_NuminBurst = bar(ax_20, [temp_Exp.Mean.Pr_Single, temp_Exp.Mean.Pr_NEventBurst; temp_Exp.Simulation_Mean.Pr_Single, temp_Exp.Simulation_Mean.Pr_NEventBurst], 'stacked',...
    'FaceColor', 'Flat', 'FaceAlpha', 0.6)
bar_NuminBurst(1).CData = repmat([1, 0, 0], 2, 1)
bar_NuminBurst(2).CData = repmat([0, 1, 0], 2, 1)
set(ax_20, 'unit', 'inch', 'fontsize', 10, 'xlim', [0.5, 3.5], 'ylim', [0, 1.2], 'box', 'off');
xticklabels(ax_20, {'Real', 'Virtual'});
yticks(ax_20, [0:0.5:1.5]);

% Number of regular Event
% if strcmp(answer_reload, 'Yes')
    
fg_15 = figure(15);
set(fg_15, 'unit', 'inch', 'position', [20, 3, 2, 1.6] , 'color', 'none')
ax_15 = axes(fg_15);
errorbar(ax_15, [5:30], control.mean_TvsEventN(5:30), control.SEM_TvsEventN(5:30), 'color', [0, 0, 0], 'capsize', 2);         % control
hold(ax_15, 'on');
% errorbar(ax_15, [5:30], TvsEventN.mean_TvsEventN(5:30), TvsEventN.SEM_TvsEventN(5:30), 'color', [1, 0, 0], 'capsize', 2);         % control
set(ax_15, 'unit', 'inch', 'FontSize', 10, 'xlim', [0 40], 'ylim', [0 80], 'box', 'off');

line_15 = refline(ax_15, 1, 0);
set(line_15, 'linestyle', ':', 'color', [0.5, 0.5, 0.5])

yticks(ax_15, [0, 40, 80]);
xlabel(ax_15, 'Number of Cluster', 'FontSize', 12);
ylabel(ax_15,  'Number of Event', 'FontSize', 12);
ax_15.Position([3,4]) = [1.1, 0.9];

% title(ax_15, Name_file, 'FontWeight', 'normal');
legend(ax_15, 'unti', 'inch', {'Exp', 'Control'}, 'fontsize', 10, 'box', 'off');
hold(ax_15, 'off');

% histogram of element number of burst
if exist('TvsEventN', 'var')
    fg_16 = figure(16);
    set(fg_16, 'unit', 'inch', 'position', [10, 5, 1.75, 1.5], 'color', 'w');
    ax_16 = axes(fg_16);
    his_eleNinT = histogram(ax_16, TvsEventN.Total_eleN_T, [1:10], 'Normalization', 'probability');
    bar(ax_16, [1:9], his_eleNinT.Values(1:9), 'FaceColor', 'flat', 'CData', [0.5, 0.5, 0.5], 'barwidth', 0.8);
    hold(ax_16, 'on');
    set(ax_16, 'unit', 'inch', 'FontSize', 20, 'xlim', [0 6], 'ylim', [0 1], 'box', 'off', 'color', 'none', 'fontsize', 10);
    xticks(ax_16, [1:1:5]);
    yticks(ax_16, [0, 0.5, 1]);
    xlabel(ax_16, 'Events in Burst', 'FontSize', 12);
    ylabel(ax_16,  'Norm. Count', 'FontSize', 12);
    hold(ax_16, 'off');
    close(fg_16);

    fg_17 = figure(17);
    set(fg_17, 'unit', 'inch', 'position', [15, 1, 2.3, 1.66], 'color', 'none');
    ax_17 = axes(fg_17);
    his_eleNinT_control = histogram(ax_17, control.Total_eleN_T, [1:10], 'Normalization', 'probability');
    control_Total_eleN_T_value = his_eleNinT_control.Values;
    his_eleNinT = histogram(ax_17, TvsEventN.Total_eleN_T, [1:10], 'Normalization', 'probability');
    bar_eleNinT = bar(ax_17, [1:1:9], [control_Total_eleN_T_value; his_eleNinT.Values], 'barwidth', 0.8);
    bar_eleNinT(1).FaceColor = [0, 255, 0]/255;
    bar_eleNinT(2).FaceColor = [0, 0, 255]/255;
    set(ax_17, 'unit', 'inch', 'xlim', [0.5, 5.5], 'ylim', [0, 1], 'fontsize', 10, 'box', 'off');
    xticks(ax_17, [1:1:10]);
    yticks(ax_17, [0, 0.5, 1]);
    xlabel(ax_17, 'Events in Burst', 'FontSize', 12);
    ylabel(ax_17,  'Norm. Count', 'FontSize', 12);
    % title(ax_17, 'Events in Burst', 'FontSize', 12, 'fontweight', 'normal'); % EGTA-AM 25 μM
    % ax_17.Position(3) = 1.66; %figure 2
    ax_17.Position(3) = 1.4; % figure 3
    ax_17.Position(4) = 0.9;
    close(fg_17);
    [T_a, T_b] = ttest2(control.Total_eleN_T, TvsEventN.Total_eleN_T);
    display(T_a);
    [ks_a, ks_b] = kstest2(control.Total_eleN_T, TvsEventN.Total_eleN_T);
    display(ks_a);

    fg_170 = figure(170);
    set(fg_170, 'unit', 'inch', 'position', [15, 9, 1, 1], 'color', 'none');
    ax_170 = axes(fg_170);
    bar_170 = bar(ax_170, [mean(control.Total_eleN_T), mean(TvsEventN.Total_eleN_T)], 'FaceColor', 'flat', 'CData', [0 255 0; 0, 0, 255]/255, 'BarWidth', 0.6);
    set(ax_170, 'unit', 'inch', 'xlim', [0.5, 2.5], 'Ylim', [0, 1.5], 'XTicklabel', [], 'YTick', [0, 1.5], ...
        'fontsize', 10, 'box', 'off', 'outerposition', [0, 0, 1, 2], 'position', ps_ax_100, 'color', 'w'); % 'XTicklabel', {'Cont', '25 μM'},
    % ylabel(ax_170, 'Events in Burst (≥1)', 'fontsize', 12);
    hold(ax_170, 'on');
    errorbar(ax_170, [mean(control.Total_eleN_T), mean(TvsEventN.Total_eleN_T)], [SEM(control.Total_eleN_T), SEM(TvsEventN.Total_eleN_T)],...
        'k', 'LineStyle', 'none', 'CapSize', 4);
    xtickangle(ax_170, 90);
    ax_170.Position([1:4]) = [0.3, 0.1, 0.27, 0.56];
    hold(ax_170, 'off');
    close(fg_170);

    % Event Numbe in Burst exacpt one event burst
    % remove one event burst

    R_Cont_Total_eleN_T = control.Total_eleN_T(control.Total_eleN_T ~= 1);
    R_Exp_Total_eleN_T = TvsEventN.Total_eleN_T(TvsEventN.Total_eleN_T ~= 1);
    [a_T b_T] = ttest2(R_Cont_Total_eleN_T, R_Exp_Total_eleN_T)
    fg_18 = figure(18);
    set(fg_18, 'unit', 'inch', 'position', [18, 1, 1.75, 1.5], 'color', 'w');
    ax_18 = axes(fg_18);
    his_R_eleNinT_control = histogram(ax_18, R_Cont_Total_eleN_T, [2:7], 'Normalization', 'probability');
    R_control_Total_eleN_T_value = his_R_eleNinT_control.Values;
    his_R_eleNinT = histogram(ax_18, R_Exp_Total_eleN_T, [2:7], 'Normalization', 'probability');
    bar_R_eleNinT = bar(ax_18, [2:6], [R_control_Total_eleN_T_value; his_R_eleNinT.Values], 'barwidth', 0.8);
    bar_R_eleNinT(1).FaceColor = [0, 255, 0]/255;
    bar_R_eleNinT(2).FaceColor = [0, 0, 255]/255;
    set(ax_18, 'unit', 'inch', 'xlim', [1.5, 6.5], 'ylim', [0, 1], 'fontsize', 10, 'box', 'off', 'color', 'none');
    yticks(ax_18, [0, 0.5, 1]);
    yticklabels(ax_18, ["0.0", "0.5", "1.0"]);
    xlabel(ax_18, '#Events in Burst', 'FontSize', 12);
    ylabel(ax_18,  'Norm. Count', 'FontSize', 12);
    % title(ax_17, 'Events in Burst', 'FontSize', 12, 'fontweight', 'normal'); % EGTA-AM 25 μM
    ax_18.Position(3) = 1.1;
    ax_18.Position(4) = 0.95;
    [a, b] = kstest2(R_Cont_Total_eleN_T, R_Exp_Total_eleN_T)
    mean(R_Exp_Total_eleN_T)
    SEM(R_Exp_Total_eleN_T)
    length(R_Exp_Total_eleN_T)

    fg_180 = figure(180);
    set(fg_180, 'unit', 'inch', 'position', [20, 3, 1, 1], 'color', 'w');
    ax_180 = axes(fg_180);
    bar_180 = bar(ax_180, [mean(R_Cont_Total_eleN_T), mean(R_Exp_Total_eleN_T)], 'FaceColor', 'flat', 'CData', [0 255 0; 0, 0, 255]/255, 'BarWidth', 0.6);
    set(ax_180, 'unit', 'inch', 'xlim', [0.5, 2.5], 'Ylim', [0, 3], 'XTicklabel', [], 'YTick', [0, 3], ...
        'fontsize', 10, 'box', 'off', 'outerposition', [0, 0, 1, 2], 'position', ps_ax_100, 'color', 'w'); % 'XTicklabel', {'Cont', '25 μM'},
    % ylabel(ax_170, 'Events in Burst (≥1)', 'fontsize', 12);
    hold(ax_180, 'on');
    errorbar(ax_180, [mean(R_Cont_Total_eleN_T), mean(R_Exp_Total_eleN_T)], [SEM(R_Cont_Total_eleN_T), SEM(R_Exp_Total_eleN_T)],...
        'k', 'LineStyle', 'none', 'CapSize', 4);
    xtickangle(ax_180, 90);
    ax_180.Position([1:4]) = [0.3, 0.1, 0.27, 0.56];
    hold(ax_180, 'off');
end

%Failure
% fg_19 = figure(19);
% set(fg_19, 'unit', 'inch', 'position', [20, 5, 2.75, 2.1], 'color', 'w');
% ax_19 = axes(fg_19);
% bar_19 = bar(ax_19, [0.0977, 0.0708, 0.0737, 0.0764], 'FaceColor', 'flat', 'BarWidth', 0.6);
% set(ax_19, 'unit', 'inch', 'xlim', [0.5, 4.5], 'Ylim', [0, 0.15], 'fontsize', 10, 'box', 'off', 'color', 'w'); 
% hold(ax_19, 'on');
% errorbar(ax_19, [0.0977, 0.0708, 0.0737, 0.0764], ...
%     [0.0017, 0.0023, 0.0024, 0.0026], 'k', 'LineStyle', 'none', 'CapSize', 12);
% title(ax_19, 'Failure', 'fontsize', 12', 'fontweight', 'normal');
% xticklabels(ax_19, {'Neuron+Astros', 'Absolute Failure', 'Neuron only', 'ACM'})
% xtickangle(ax_19, 30);
% yticks(ax_19, [0:0.05:0.15]);
% yticklabels(ax_19, {'0.00', '0.05', '0.10', '0.15'});
% ylabel(ax_19, 'Probability', 'fontsize', 12)
% ax_19.Position([3, 4]) = [1.9, 0.95];
% bar_19.CData(1, :) = [1, 0, 0];
% bar_19.CData(2, :) = [192, 0, 0]/255;
% bar_19.CData(3, :) = [255, 204, 204]/255;
% bar_19.CData(4, :) = [255, 204, 204]/255;
% hold(ax_19, 'off');

% single trace
fg_192 = figure(192);
set(fg_192, 'units', 'inches', 'position', [26, 18, 2.3, 2.2], 'color', 'w');
ax_192 = axes(fg_192);
plot(ax_192, [-timeStartEstimulation:timeframe:1/Frequency - timeframe], mean(roi_traceF, 2), 'k');
hold(ax_192, 'on');
errorbar(ax_192, [-timeStartEstimulation:timeframe:1/Frequency - timeframe], mean(roi_traceF, 2), SEM(roi_traceF, 2), 'k.', 'markersize', 1, 'capsize', 2);
set(ax_192, 'units', 'inches', 'box', 'off', 'color', 'none');
xlabel(ax_192, 'time (s)', 'fontsize', 12);
ylabel(ax_192, 'fluorescence', 'fontsize', 12);
fg_192_title = 'single trace';
title(ax_192, fg_192_title);
ax_192.Title.Visible = 'off';

fg_193 = figure(193);
set(fg_193, set_fg{:}, 'position', [28.5, 18, 2.3, 2.2]);
ax_193 = axes(fg_193);
plot(ax_193, (timeframe * (1-interEstimulation):timeframe:1/Frequency- timeframe), mean(N_roi_traceF, 2), '.-k');
hold(ax_193, 'on');
errorbar(ax_193, (timeframe * (1-interEstimulation):timeframe:1/Frequency- timeframe), mean(N_roi_traceF, 2), SEM(N_roi_traceF, 2), 'k.', 'markersize', 1, 'capsize', 2);
set(ax_193, set_ax{:});
xlabel(ax_193, 'time (s)', 'fontsize', 12);
ylabel(ax_193, 'norm.F', 'fontsize', 12);
fg_193_title = 'single trace(norm)';
title(ax_193, fg_193_title);

% MARKER
if isfield(temp_Exp.Period, 'Dist2Marker')||isfield(temp_Exp.Period, 'Dist2Marker_before')||isfield(temp_Exp.Period, 'Dist2Marker_after')
    Marker_threshold = inputdlg('Enter Marker distance threshold', 'Marker Threshold', [1, 35], {'0.5'});
    Marker_threshold = str2num(Marker_threshold{1});
    Pr = [temp_Exp.Period.Prob];
    if isfield(temp_Exp.Period, 'Dist2Marker')
        Dist2Marker = [Period.Dist2Marker];
        Dist2Marker_near_idx = Dist2Marker <= Marker_threshold;
        Dist2Marker_far_idx = ~Dist2Marker_near_idx;
    end
    if isfield(temp_Exp.Period, 'Dist2Marker_before')
        Dist2Marker_before = [Period.Dist2Marker_before];
        Dist2Marker_before_near_idx = Dist2Marker_before <= Marker_threshold;
        Pr_before_near = Pr(Dist2Marker_before_near_idx);
        Pr_before_far = Pr(~Dist2Marker_before_near_idx);
        fg_211 = figure(211);
        set(fg_211, 'unit', 'inch', 'position', [26, 15, 2, 2], 'color', 'w');
        ax_211 = axes(fg_211);
        b_211 = bar(ax_211, [1, 2], [mean(Pr_before_near), mean(Pr_before_far)], 'facecolor', 'flat', 'edgecolor', 'none', 'barwidth', 0.8);
        b_211.CData = [255, 192, 0;112, 173, 71] / 255;
        hold(ax_211, 'on');
        errorbar(ax_211, [1, 2], [mean(Pr_before_near), mean(Pr_before_far)], [SEM(Pr_before_near), SEM(Pr_before_far)], 'k.', 'markersize', 1, 'capsize', 20);
        set(ax_211, 'units', 'inches', 'fontsize', 10, 'xlim', [0.5, 2.5], 'ylim', [0, 0.1001], 'box', 'off', 'color', 'none');
        yticks(ax_211, [0:0.05:0.1]);
        yticklabels(ax_211, {'0.00', '0.05', '0.10'});
        ylabel(ax_211, 'probability', 'fontsize', 12);
        ax_211.XTick = [];
        [~, b_211] = ttest2(Pr_before_near, Pr_before_far);
        display(['t-test for Pr of Marekr: ', num2str(b_211)]);

    end
    if isfield(temp_Exp.Period, 'Dist2Marker_after')
        Dist2Marker_after = [Period.Dist2Marker_after];
        Dist2Marker_after_near_idx = Dist2Marker_after <= Marker_threshold;
        roi_traceF = [Period.roi_traceF] - [Period.Background];
        traceF_after_near = roi_traceF(Dist2Marker_after_near_idx);
        traceF_after_far = roi_traceF(~Dist2Marker_after_near_idx);
    end
end

%% Saving data
answer = questdlg('Would you like to overlap to original data?');                         % ask whether save this data
if strcmp(answer, 'Yes')
    eval([Name_Exp_file, ' = temp_Exp;']);
    cd(path);
    eval(['save(''', Name_Exp_file_mat, ''', ''', Name_Exp_file, ''');']);
elseif strcmp(answer, 'Cancel')
    return
end
% save figures
cd('H:\My Drive\making paper\figure');
split_Name = strsplit(Name_Exp_file, '_');
save_fig_name = strjoin(split_Name(1:end-1), '_');  
for i = [7, 41, 83, 85, 88, 89, 100, 192, 193]
    if i == 89 && ~exist('fg_89_title', 'var')
        continue
    end
    eval(['save_fg_', num2str(i), ' = erase(strjoin({save_fig_name, fg_', num2str(i), '_title}, "_"), "\");']);
    eval(['savefig(fg_', num2str(i), ', save_fg_', num2str(i), ');']);
end

%% MVR
if strcmp(answer_MVR, 'Yes')
    fprintf('%s File is running MVR\n', Name_Exp_file_mat);
    JY_MVR(temp_Exp, Name_Exp_file_mat);
end

