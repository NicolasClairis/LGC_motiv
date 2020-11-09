% function [asc2,eye_signal,eye_signal_SEM] = read_eyelink_2(run_filename)
% reads eye_data and outputs graph with pupil area size at the time around
% of each event (use for each run)
%
% See also: read_eyelink_1_asc_kk.m

if run_filename == 0
    filename = input('Nom du fichier oeil sXyyMMrZ (X:subject,yyMM:date,Z:run)(.asc)?','s'); % eyedata
elseif run_filename ~= 0 % for read_eyelink_across_runs() function
    filename = run_filename;
end

% parameters that can vary
SD_or_SEM = 0; % to show SEM = 0, SD = 1 for graphs
raw_or_mean = 1; % to show raw = 1, mean = 0
TR = 2010; % TR from particular experiment
diff_TR = 20; % tolerance for TR variation

TRmin = TR - diff_TR;
TRmax = TR + diff_TR;

asc = read_eyelink_1_asc_kk([filename,'.asc']);


asc2.header = asc.header;
asc2.msg = asc.msg;
asc2.efix = cell2mat(struct2cell(asc.efix)'); %contains eye R/L, start, end.
asc2.pupil = asc.pupil; % check if pupil AREA or else
asc2.esacc = cell2mat(struct2cell(asc.esacc)'); %contains eye R/L, start, end.

% load TTL, check if they are ok (close to TR)
TTL = cell2mat(struct2cell(asc.button)');
asc2.TTL = TTL(TTL(:,3) == 1,1);
asc2.edata = asc.dat; %contains timestamp, x and y movements, pupil size
clear('asc');

% check if TTL ok compared to TR
TTLok = diff(asc2.TTL);
if (sum(TTLok < TRmin) > 0 || sum(TTLok > TRmax) > 0)
    disp('probleme avec TTL/TR : verifier les donnees');
    return;
    % in case need to exclude some of the first TTL:
% elseif strcmp(filename,'s10102r1a') == 1
%     asc2.TTL = asc2.TTL(5:end); % exclude 4 first TTL (bug)
%     disp('4 first bad TTL excluded, TTL ok');
elseif strcmp(filename,'MS_s01r4') == 1 %|| strcmp(filename,'s20802r4') == 1 || strcmp(filename,'s42402r1') == 1
    asc2.TTL = asc2.TTL(3:end); % exclude 2 first TTL (bug)
    disp('2 first bad TTL excluded, TTL ok');
else
    disp('TTL ok');
end


% calibrate all times to the first TTL
asc2.efix(:,2) = asc2.efix(:,2) - asc2.TTL(1);
asc2.efix(:,3) = asc2.efix(:,3) - asc2.TTL(1);
asc2.esacc(:,2) = asc2.esacc(:,2) - asc2.TTL(1);
asc2.esacc(:,3) = asc2.esacc(:,3) - asc2.TTL(1);
asc2.edata(1,:) = asc2.edata(1,:) - asc2.TTL(1);
% % remove eye-blinks from edata:
% asc2.edata = asc2.edata(:,asc2.edata(4,:)~=0);

asc2.TTL = asc2.TTL(:) - asc2.TTL(1);
asc2.TTL = asc2.TTL/1000;


%% load onsets
date_exp = datestr([filename(5:6),'/',filename(3:4),'/2016'],'dd-mmm-yyyy');
if size(filename,2) == 8
    load([date_exp,'_onsets_filtered_neutralFB_run',filename(8),'.mat']);
elseif size(filename,2) > 8
    load([date_exp,'_onsets_filtered_neutralFB_run',filename(8:9),'.mat']);
end
nberEvents = numel(names); % ER1/2/3/SC/INC/FC/Repeat
events_onsets = onsets; % ER1/2/3/SC/INC/FC/Repeat

eye_signal = zeros(nberEvents,1202); % mean signal for each symbol/run
if SD_or_SEM == 0
    eye_signal_SEM = zeros(nberEvents,1202);
elseif SD_or_SEM == 1
    eye_signal_SD = zeros(nberEvents,1202);
end
eye_signal_bins_ER123 = zeros(4,6,3); % mean signal (4 occurrences/run-lines) per bin (6bins-columns) ER1/2/3 (3 matrices)
eye_signal_bins_SC = zeros(15,6); % mean signal (4 occurrences/run) per bin (6bins)

figure
col = {'r','k','b','m','c','g','y'};
nberEvents = 4; % check only ER1/2/3/SC and not INC/FC/Repeat
events = zeros(1,nberEvents);
for i = 1:nberEvents %nberEvents % ER1/2/3/SC/INC/FC/Repeat
    % remove first line (0) from onsets (made for fMRI)
    events_onsets{1,i} = events_onsets{1,i}(2:end,1)*1000;
    
    %% pupil size reaction
    % look at each occurrence
    eye_signal_j = zeros(numel(events_onsets{:,i}),1202);
    for j = 1:numel(events_onsets{:,i})
        i_lower  = find(asc2.edata(1,:) <= events_onsets{1,i}(j,1),1,'last');
        i_higher = find(asc2.edata(1,:) >= events_onsets{1,i}(j,1),1,'first');
        eye_signal_j(j,:) = asc2.edata(4,i_lower-200:i_higher+1000);
        if raw_or_mean == 1
            x = -200:1:1001;
            events(i) = plot(x,asc2.edata(4,i_lower-200:i_higher+1000),col{i});
            hold on;
        end
        % bins
        for bin = 1:6 % 6 bins (pas de 200ms)
            bin_signal = eye_signal_j(j,(bin-1)*200+1:bin*200);
            bin_signal = bin_signal(1,bin_signal~=0); % remove blinks/eye-closure
            if i < 4 % ER but not the SC
                eye_signal_bins_ER123(j,bin,i) = mean(bin_signal);
            elseif i == 4 % SC
                eye_signal_bins_SC(j,bin) = mean(bin_signal);
            elseif j > 4 % INC/FC/R
                % a penser si besoin (en particulier pour comparer FC et
                % Repeat)
            end
        end
    end
    
    % mean +/- SEM/SD
    eye_signal(i,:) = mean(eye_signal_j,1);
    lineProps.width = 5;
    lineProps.edgestyle = ':';
    lineProps.col = col(i);
    x = -200:1:1001;
    if SD_or_SEM == 0 && raw_or_mean == 0
        eye_signal_SEM(i,:) = grpstats(eye_signal_j,[],'sem');
        %         errorbar(eye_signal(i,:),eye_signal_SEM(i,:),col{i},'MarkerSize',12); % Mean +- SEM
        mseb(x,eye_signal(i,:),eye_signal_SEM(i,:),lineProps,1); % Mean +- SEM
    elseif SD_or_SEM == 1 && raw_or_mean == 0
        eye_signal_SD(i,:) = std(eye_signal_j,0,1);
        %         errorbar(eye_signal(i,:),eye_signal_SD(i,:),col{i},'MarkerSize',12); % Mean +- SD
        mseb(x,eye_signal(i,:),eye_signal_SD(i,:),lineProps,1); % Mean +- SD
    end
    hold on
    
end
ylim([0 4500]);
xlim([-225 1025]);
ylabel('pupil area (mm3)');
xlabel('time around event (ms)');
% mark onset:
x=zeros(1,451);
y=0:10:4500;
onset0 = plot(x,y,'k--');
if raw_or_mean == 0
    title(['Average pupil area around each event occurrence +/- SEM - ',filename(1:2),' Run ',filename(8:end)]);
    legend('INC-COR','COR-INC','Neutral','SC','0 = event onset','Location','northeast'); % 'ER1','ER2','ER3','SC'
    legend('boxoff');
elseif raw_or_mean == 1
    title(['Raws pupil area around each event occurrence - ',filename(1:2),' Run ',filename(8:end)]);
    legend([events(1:nberEvents) onset0],'INC-COR','COR-INC','Neutral','SC','0 = event onset','Location','northeast'); % 'ER1','ER2','ER3','SC','onset 0'
    legend('boxoff');
end

