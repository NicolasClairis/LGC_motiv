
root = 'B:\resultats\';
% which study to check
which_study = input(sprintf('Which study? multiseq (0) or MBBjune2016 (1)? \n'));
if which_study == 0
    cd([root,'multiseq_april_may2016']);
elseif which_study == 1
    cd([root,'MBB_june2016']);
end
% which subject
sub_name = input(sprintf('subject filename? \n'),'s');
cd([sub_name,'\eye_data\'])
% which run
runname = input(sprintf('which run? \n'),'s');
% extract subject number for naming the file
if which_study == 1
    subid = sub_name(2:3);
elseif which_study == 0
    if strcmp(sub_name(3),'_')
        subid = sub_name(2);
    elseif strcmp(sub_name(4),'_')
        subid = sub_name(2:3);
    end
end
% different name depending on study and number
if which_study == 0
    if length(subid) < 2
        filename = ['MSs_s',subid,'r',runname,'.asc'];
    elseif length(subid) == 2
        filename = ['MSs_',subid,'r',runname,'.asc'];
    end
elseif which_study == 1
    filename = ['MS_s',subid,'r',runname,'.asc'];
end
asc = read_eyelink_1_asc_kk(filename);

figure
set(get(handle(gcf),'JavaFrame'),'Maximized',1); % maximize window size
hold on
xlim([0 1024]);
ylim([0 768]);

for k = 1:size(asc.dat,2)
    plot(asc.dat(2,k),asc.dat(3,k));
    drawnow;
end