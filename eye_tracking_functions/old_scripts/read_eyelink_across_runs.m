function [eye_signal,eye_signal_SEM] = read_eyelink_across_runs()

% reads eye data and outputs a graph showing the evolution for each symbol
% in function of the run

filename = input('Nom du fichier/quel sujet? sXyyMM (X:subject,yyMM:date)(.asc)','s');

%% load subjects DATA

x = dir([filename,'r*','*.asc']); % takes all the runs for the selected subject
logfile = struct2cell(x);
NR = size(logfile,2);  % Nber of Runs
NRstr = num2str(NR);

eye_signal = zeros(7,1202,NR);
eye_signal_SEM = zeros(7,1202,NR);
% 7 lines: ER1/2/3/SC/INC/FC/R
% 1202 columns: raw signal
% NR 3D matrices: 1 matrix/run
run_filename = cell(1,NR);
Scan = cell(1,NR);

for i = 1:NR
    if size(logfile{1,i},2) == 12 % sX-yy-MM-r(i).asc
        run_filename{i} = logfile{1,i}(1:8);
    elseif size(logfile{1,i},2) == 13 % sX-yy-MM-r(i)a/b/c/d.asc
        run_filename{i} = logfile{1,i}(1:9);
    else
        numscan = num2str(i);
        run_filename{i} = input(['nom du fichier au nom trop long? (',numscan,'e scan)'],'s');
    end
    [asc2,eye_signal(:,:,i),eye_signal_SEM(:,:,i)] = read_eyelink_2(run_filename{i});
end

save([filename(1:2),'_',filename(3:6),'16_',NRstr,'runs_raw_eye_data_around_ER'],'eye_signal','eye_signal_SEM');
save([filename(1:2),'_',filename(3:6),'16_',NRstr,'runs_raw_eye_data'],'asc2');

col = {'r','k','b','m','c','g','y'};
event = {'INC-COR','COR-INC','Neutral','SC','INC','First Correct','Correct Repeat'};
lineProps.width = 5;
lineProps.edgestyle = ':';
x=-200:1:1001;

for i = 1:7 % 1 graph/event: 7 events ER1/2/3/SC/INC/FC/R
    figure
    for j=1:NR % 1 curve/run
        lineProps.col = col(j); % 1color/run
        mseb(x,eye_signal(i,:,j),eye_signal_SEM(i,:,j),lineProps,1);
        hold on
        Scan{j} = ['Run',run_filename{j}(8:end)];
    end
    title(['Pupil area size around ',event{i},' - ',filename(1:2),' ',filename(3:4),'/',filename(5:6),'/2016']);
    legend(Scan{1:end});
    legend('boxoff');
    ylabel('Mean pupil area size (mm3) +/- SEM per run');
    xlabel('time (200ms before onset of the event + 1s duration of the event');
    xlim([-225 1025]);
    ylim([0 4500]);
end

