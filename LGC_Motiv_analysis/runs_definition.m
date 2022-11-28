function[runs, n_runs] = runs_definition(study_nm, sub_nm, condition)
% [runs, n_runs] = runs_definition(study_nm, sub_nm, condition)
%runs_definition will define the number and order of the physical and
%mental effort task for the study and subject entered in input.
%
% INPUTS
% study_nm: study name
%
% sub_nm: subject name
%
% condition: in some rare cases we have behavior, but not fMRI data. This
% variable allows to take into account this so that the corresponding runs
% are still included in the behavioral analysis but not in the fMRI
% analysis.
% 'behavior': behavioral analysis
% 'behavior_no_sat': behavioral analysis without the runs where subject
% saturated
% 'fMRI': fMRI analysis
% 'fMRI_no_move': fMRI analysis without runs where too much movement
%
% OUTPUTS
% runs: structure with number of runs for each task and also the order of
% the runs
%
% n_runs: number of runs to include

%% extract main structure of the task for behavior and for fMRI
switch study_nm
    case 'fMRI_pilots'
        switch sub_nm
            case 'pilot_s1'
                runs.tasks = {'Ep','Em'};
            case 'pilot_s2'
                runs.tasks = {'Ep'};
            case 'pilot_s3'
                runs.tasks = {'Em','Ep'};
        end
    case 'study1'
        switch sub_nm
            case {'040'} % only performed 1 run of each task
                runs.tasks = {'Ep','Em'};
            case {'001','002','003','004','005','008','009',...
                    '012','013','015','018','019',...
                    '027',...
                    '030','036','038','039',...
                    '042','046','047','048',...
                    '050','053',...
                    '060','064','065','069',...
                    '072','076',...
                    '080','083','085','087',...
                    '090','093','094','097'}
                runs.tasks = {'Ep','Em','Ep','Em'};
            case {'011','017',...
                    '020','021','022','024','029',...
                    '032','034','035',...
                    '043','044','045','049',...
                    '052','054','055','056','058','059',...
                    '061','062','068',...
                    '071','073','074','075','078','079',...
                    '081','082','086','088',...
                    '091','095','099',...
                    '100'}
                runs.tasks = {'Em','Ep','Em','Ep'};
        end
    case 'study2'
        error('need to attribute runs');
    otherwise
        error('study not recognized');
end

%% remove runs that were problematic
% (either behavioral saturation or runs with too much movement in the fMRI)

switch study_nm
    case 'study1'
        
        % by default include all sessions
        runs.runsToKeep = 1:4;
        runs.runsToIgnore = [];
        
        %% remove subjects where behavior and fMRI could not be performed => remove runs independently of the condition
        switch sub_nm
            case {'030','049'}
                runs.runsToKeep = [];
                runs.runsToIgnore = 1:4;
            case '040' % fMRI crashed during run 3 and subject was already stressing a lot
                % => avoid including this run
                runs.runsToKeep = [1,2];
                runs.runsToIgnore = 3:4;
        end
        
        %% define subject runs to keep depending on condition
        switch condition
            %% removing runs with saturation
            case 'fMRI_noSatTask'
                switch sub_nm
                    case '002'
                        runs.runsToKeep = 3;
                        runs.runsToIgnore = [1,2,4];
                    case '005'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '012'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '017' % first run: fMRI crashed => we have the behavior but not enough trials for fMRI
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '027'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '032'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '043' % first run: fMRI crashed => we have the behavior but not enough trials for fMRI
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '047'
                        runs.runsToKeep = 3;
                        runs.runsToIgnore = [1,2,4];
                    case '048'
                        runs.runsToKeep = [1,3,4];
                        runs.runsToIgnore = 2;
                    case '052'
                        runs.runsToKeep = [2,4];
                        runs.runsToIgnore = [1,3];
                    case '074' % first run: fMRI crashed => we have the behavior but not enough trials for fMRI
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '076'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '095'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '100'
                        runs.runsToKeep = [1,2];
                        runs.runsToIgnore = [3,4];
                end
            %% too much movement cleaning
            case 'fMRI_no_move'
                switch sub_nm
                    %% subjects with too much movement in ALL runs
                    case {'008','022','024'}
                        runs.runsToKeep = [];
                        runs.runsToIgnore = 1:4;
                    %% subjects with some runs with too much movement
                    case '017' % first run: fMRI crashed => we have the behavior but not enough trials for fMRI
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '021'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = 2:4;
                    case '029'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '043' % first run: fMRI crashed => we have the behavior but not enough trials for fMRI
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '044'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '047'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '053'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '054'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '058'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '062'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = 2:4;
                    case '071'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '074' % first run: fMRI crashed => we have the behavior but not enough trials for fMRI
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '076'
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '078'
                        runs.runsToKeep = [2,3];
                        runs.runsToIgnore = [1,4];
                    case '080'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '083'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '087'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = [2,3,4];
                    case '097'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = [2,3,4];
                    case '099'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                end % subject
            case {'fMRI','fMRI_no_move_bis','fMRI_noSatRun'} % all other fMRI cases
                switch sub_nm
                    case {'017','043','074'} % first run: fMRI crashed => we have the behavior but not enough trials for fMRI
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                end % subject 
        end
    otherwise
        error('case not ready yet');
end % study

% update task types depending on the runs to keep
runs.tasks = runs.tasks(runs.runsToKeep);

%% extract number of runs of each task type
runs.nb_runs.Ep = sum(strcmp(runs.tasks,'Ep'));
runs.nb_runs.Em = sum(strcmp(runs.tasks,'Em'));

%% extract number of runs
n_runs = length(runs.tasks);

end % function