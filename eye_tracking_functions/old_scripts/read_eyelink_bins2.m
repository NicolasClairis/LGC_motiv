% function [] = read_eyelink_bins2()

filename = input('Nom du fichier/quel sujet? sXyyMM (X:subject,yyMM:date)(.asc)','s');

%% load subjects DATA

x = dir([filename,'r*','*.asc']); % takes all the runs for the selected subject
logfile = struct2cell(x);
NR = size(logfile,2);  % Nber of Runs
% NRstr = num2str(NR);
run_filename = cell(1,NR);
Scan = cell(1,NR);

eye_bins = zeros(3,6,NR); % 3 lines: ER1/2/SC; 6 bins; NR (nber of runs) matrices

% figure properties
col = {'r','k','b','m','c','g','y'};
event = {'INC-COR','COR-INC','SC'}; %,'INC','First Correct','Correct Repeat'};
lineProps.width = 5;
lineProps.edgestyle = ':';
% x = 1:6; % for mseb if SEM used mostly
for i = 1:NR
    % number of the run (especially for the case of run divided in many
    % parts 1a/1b/etc.)
    if size(logfile{1,i},2) == 12 % sX-yy-MM-r(i).asc
        run_filename{i} = logfile{1,i}(1:8);
    elseif size(logfile{1,i},2) == 13 % sX-yy-MM-r(i)a/b/c/d.asc
        run_filename{i} = logfile{1,i}(1:9);
    else
        numscan = num2str(i);
        run_filename{i} = input(['nom du fichier au nom trop long? (',numscan,'e scan)'],'s');
    end
    
    load([filename(1:2),'_',filename(3:6),'16_raw_eye_data_around_ER_run',run_filename{i}],'nber_ER_samples','mean_ER','ER_sem');
    for j = 1:4
        mean_ER(j,:) = mean_ER(j,:) - mean_ER(1,:); % soustrait baseline
        for bin = 1:6
            mean_ER(j,bin) = mean_ER(j,bin) - mean_ER(3,bin); % soustrait neutre pour chaque bin correspondant
        end
    end
    eye_bins(1,:,i) = mean_ER(1,:); % ER1
    eye_bins(2,:,i) = mean_ER(2,:); % ER2
    eye_bins(3,:,i) = mean_ER(4,:); % SC
end


for j = 1:3 % 1 graph/event: 7 events ER1/2/3/SC/INC/FC/R
    figure
    for i=1:NR % 1 curve/run
%         lineProps.col = col(i); % 1color/run
        %         mseb(x,eye_bins(i,:),ER_sem(i,:),lineProps,1);
        plot(eye_bins(j,:,i),col{i});
        hold on
        Scan{i} = ['Run',run_filename{i}(8:end)];
    end
    title(['Pupil area size around ',event{j},' - ',filename(1:2),' ',filename(3:4),'/',filename(5:6),'/2016']);
    legend(Scan{1:end});
    legend('boxoff');
    ylabel('Mean pupil area size (mm3) +/- SEM per run');
    xlabel('6 bins - 200ms each - bin1=before onset of ER');
    xlim([0 7]);
    ylim([-6500 6500]);
end

% 
% 
% 
% for i= 1:15
% eye_signal_SC_test = eye_signal_bins_SC(i,:) - eye_signal_bins_SC(i,1);
% plot(eye_signal_SC_test);
% end
% 
% 
% for j=1:3
%     figure
%     for i= 1:4
%         eye_signal_ER123_test = eye_signal_bins_ER123(i,:,j) - eye_signal_bins_ER123(i,1,j);
%         plot(eye_signal_ER123_test);
%         hold on
%     end
% end