function [] = read_eyelink_bins()

% reads eye data and outputs a graph showing the evolution for each symbol
% in function of the run

filename = input('Nom du fichier/quel sujet? sXyyMM (X:subject,yyMM:date)(.asc)','s');

%% load subjects DATA

x = dir([filename,'r*','*.asc']); % takes all the runs for the selected subject
logfile = struct2cell(x);
NR = size(logfile,2);  % Nber of Runs
% NRstr = num2str(NR);

run_filename = cell(1,NR);
nber_ER_samples = zeros(1,4); % 4 events in line ER1/2/3/SC

for i = 1:NR
    %% extract run eye bins data
    if size(logfile{1,i},2) == 12 % sX-yy-MM-r(i).asc
        run_filename{i} = logfile{1,i}(1:8);
    elseif size(logfile{1,i},2) == 13 % sX-yy-MM-r(i+/-?)a/b/c/d.asc
        run_filename{i} = logfile{1,i}(1:9);
    else
        numscan = num2str(i);
        run_filename{i} = input(['nom du fichier au nom trop long? (',numscan,'e scan)'],'s');
    end
    [~,~,~,eye_signal_bins_ER123,eye_signal_bins_SC] = read_eyelink_2(run_filename{i});
    
    %% extract only lines with eye data
    mean_ER = zeros(4,6); % 4 events (ER1/2/3/SC) and 6 bins
    ER_sem = zeros(4,6);
    for j = 1:4 % loop through events ER1/2/3/SC
        raw_ER = [];
        if j < 4
            raw_eye = eye_signal_bins_ER123(:,:,j);
        elseif j == 4
            raw_eye = eye_signal_bins_SC(:,:);
        end
        max_raw_eye = size(raw_eye,1); % nber of lines = nber of occurrences
        for line = 1:max_raw_eye
            if sum(raw_eye(line,:)==0) == 0 % exclude null lines (=baseline with blink, cf read_eyelink_2)
                % if blink after the onset of the ER, full trial not deleted,
                % but the part with the blink is excluded (see
                % read_eye_link_2)
                raw_ER = [raw_ER; raw_eye(line,:)];
                nber_ER_samples(j) = nber_ER_samples(j) + 1;
            end
        end
        % save raws
        ER_id = num2str(j);
        save([filename(1:2),'_',filename(3:6),'16_raw_eye_data_around_ER_run',run_filename{i},'_ER',ER_id],'raw_ER');
        
        % mean each bin
        for bin = 1:6
            if isempty(raw_ER) == 0
                mean_ER(j,bin) = mean(raw_ER(:,bin));
                ER_sem(j,bin) = grpstats(raw_ER(:,bin),[],'sem');
            end
        end
    end
    % save all files because using read-eyelink functions takes time
    save([filename(1:2),'_',filename(3:6),'16_raw_eye_data_around_ER_run',run_filename{i}],'nber_ER_samples','mean_ER','ER_sem');
end
