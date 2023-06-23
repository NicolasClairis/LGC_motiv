function[MadeOrNot] = checkIfRunMade(behaviordir, subid, nrun, nbTotalSessions)
% [MadeOrNot] = checkIfRunMade(behaviordir, subid, nrun, nbTotalSessions)
% function to check whether this run number has already been used for this
% subject or not
%
% INPUT:
% behaviordir: folder with behavioral data for this subject
% subid: subject id number as a string
% nrun: run number
% nbTotalSessions: total number of runs to be made for this subject (all
% tasks considered)
%
% OUTPUT
% MadeOrNot: 1 if run already made, 0 if not
%
% Used for taskGripRP, taskMentalRP and taskLearning75 but easy to adapt
% for any other task.
% made by Nicolas Clairis - february 2017


% check all runs already made (= results file created)
iExistCheck = 1;
run_made_or_not = zeros(1,nbTotalSessions);
taskNames = {'learning','gripRP','mentalRP'};
while iExistCheck <= nbTotalSessions
    sessnber = num2str(iExistCheck);
    for task = 1:length(taskNames) % check in all taskNames whether this number of run was used or not
        if exist([behaviordir filesep 'MBB_battery_',taskNames{task},'_onsets_sub',subid,'_sess',sessnber,'.mat'],'file') == 2
            run_made_or_not(iExistCheck) = 1;
        end
    end
    iExistCheck = iExistCheck + 1;
end

% check if the run entered in input here was made or not
MadeOrNot = 0;
if nrun <= max(find(run_made_or_not == 1))
    runname = num2str(nrun);
    warning(['Erreur dans le numero du run ? Vous avez indique run ',runname,' mais un fichier existe deja avec ce numero de run']);
    MadeOrNot = 1;
end

end
