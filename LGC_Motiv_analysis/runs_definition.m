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
% See subject_condition.m for the full list.
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
            case {'001','002','003','004','005','008','009',...
                    '012','013','015','018','019',...
                    '027',...
                    '030','036','038','039',...
                    '040','042','046','047','048',...
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
                error([sub_nm,' should not be included under the condition ''',condition,...
                    ''' (tasks not performed).']);
            case '040' % fMRI crashed during run 3 and subject was already stressing a lot
                % => avoid including this run
                runs.runsToKeep = [1,2];
                runs.runsToIgnore = 3:4;
        end
        
        %% define subject runs to keep depending on condition
        switch condition
            case {'behavior'} % nothing to change for those conditions
            %% for all fMRI conditions, need to remove run 1 from those subjects because of fMRI crash
            case {'fMRI','fMRI_noSatRunSub'}
                switch sub_nm
                    case {'017','043','074'} % first run: fMRI crashed => we have the behavior but not enough trials for fMRI
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                end % subject
                %% removing tasks fully saturated
            case {'behavior_noSatTask','behavior_noSatTask_bayesianMdl'}
                switch sub_nm
                    case '027' % Em full saturation
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '047' % Em full saturation
                        error(['subject ',sub_nm,' should not be included (saturated ALL tasks)']);
                    case '052' % Em full saturation
                        runs.runsToKeep = [2,4];
                        runs.runsToIgnore = [1,3];
                    case '069' % Em full saturation
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '076' % Em full saturation
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '095' % Ep full saturation
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                end % subject loop
            case {'fMRI_noSatTask','fMRI_noSatTask_bayesianMdl'}
                switch sub_nm
                    case {'017','043','074'} % first run: fMRI crashed => we have the behavior but not enough trials for fMRI
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '027' % Em full saturation
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '047' % Em full saturation
                        error(['subject ',sub_nm,' should not be included (saturated ALL tasks)']);
                    case '052' % Em full saturation
                        runs.runsToKeep = [2,4];
                        runs.runsToIgnore = [1,3];
                    case '069' % Em full saturation
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '076' % Em full saturation
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '095' % Ep full saturation
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                end % subject loop
                %% removing any run with saturation
            case {'behavior_noSatRun','behavior_noSatRun_bayesianMdl'}
                switch sub_nm
                    case '002'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '004'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '005'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '012'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '022'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '027'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = [2,3,4];
                    case '032'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '038'
                        runs.runsToKeep = [1,3,4];
                        runs.runsToIgnore = 2;
                    case '044'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '047'
                        error(['subject ',sub_nm,' should not be included (saturated ALL tasks)']);
                    case '048'
                        runs.runsToKeep = [1,3,4];
                        runs.runsToIgnore = 2;
                    case '052'
                        runs.runsToKeep = [2,4];
                        runs.runsToIgnore = [1,3];
                    case '054'
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '055'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '058'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '061'
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '062'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '069'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '076'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '081'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '082'
                        runs.runsToKeep = [1,2];
                        runs.runsToIgnore = [3,4];
                    case '083'
                        runs.runsToKeep = [1,2,3];
                        runs.runsToIgnore = 4;
                    case '088'
                        runs.runsToKeep = [1,2,3];
                        runs.runsToIgnore = 4;
                    case '095'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = [2,3,4];
                    case '097'
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '099'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '100'
                        runs.runsToKeep = [1,2];
                        runs.runsToIgnore = [3,4];
                end
                %% respiration
            case 'respiration_and_noSatRun'
                switch sub_nm
                    case '002'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '004'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '005'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '012'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '022'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '027'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = [2,3,4];
                    case '032'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '038'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = 2:4;
                    case '044'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '047'
                        runs.runsToKeep = [];
                        runs.runsToIgnore = 1:4;
                    case '048'
                        runs.runsToKeep = [1,3,4];
                        runs.runsToIgnore = 2;
                    case '050'
                        runs.runsToKeep = 2;
                        runs.runsToIgnore = [1,3,4];
                    case '052'
                        runs.runsToKeep = [2,4];
                        runs.runsToIgnore = [1,3];
                    case '054'
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '055'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '056'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '058'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '061'
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '062'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '069'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '076'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '081'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '082'
                        runs.runsToKeep = [1,2];
                        runs.runsToIgnore = [3,4];
                    case '083'
                        runs.runsToKeep = [1,2,3];
                        runs.runsToIgnore = 4;
                    case '088'
                        runs.runsToKeep = [1,2,3];
                        runs.runsToIgnore = 4;
                    case '095'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = [2,3,4];
                    case '097'
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '099'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '100'
                        runs.runsToKeep = [1,2];
                        runs.runsToIgnore = [3,4];
                end
            case {'fMRI_noSatRun','fMRI_noSatRun_bayesianMdl',...
                    'fMRI_noSatTaskSub_noSatRun', 'fMRI_noSatTaskSub_noMoveSub_noSatRun'}
                switch sub_nm
                    case {'017','043','074'} % first run: fMRI crashed => we have the behavior but not enough trials for fMRI
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '002'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '004'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '005'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '012'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '022'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '027'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = [2,3,4];
                    case '032'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '038'
                        runs.runsToKeep = [1,3,4];
                        runs.runsToIgnore = 2;
                    case '044'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '047'
                        error(['subject ',sub_nm,' should not be included (saturated ALL tasks)']);
                    case '048'
                        runs.runsToKeep = [1,3,4];
                        runs.runsToIgnore = 2;
                    case '052'
                        runs.runsToKeep = [2,4];
                        runs.runsToIgnore = [1,3];
                    case '054'
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '055'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '058'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '061'
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '062'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '069'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '076'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '081'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '082'
                        runs.runsToKeep = [1,2];
                        runs.runsToIgnore = [3,4];
                    case '083'
                        runs.runsToKeep = [1,2,3];
                        runs.runsToIgnore = 4;
                    case '088'
                        runs.runsToKeep = [1,2,3];
                        runs.runsToIgnore = 4;
                    case '095'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = [2,3,4];
                    case '097'
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '099'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '100'
                        runs.runsToKeep = [1,2];
                        runs.runsToIgnore = [3,4];
                end
            case {'fMRI_noSatRun_choiceSplit_Elvl','fMRI_noSatTaskSub_noSatRun_choiceSplit_Elvl'} % no saturation run + if trials are split according to choice made including a high effort regressor
                % => need to remove runs where high or low effort choice always go with the same effort level
                switch sub_nm
                    case {'017','043'} % first run: fMRI crashed => we have the behavior but not enough trials for fMRI
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '074' % first run: fMRI crashed => we have the behavior but not enough trials for fMRI + run 3 saturation
                        runs.runsToKeep = [2,4];
                        runs.runsToIgnore = [1,3];
                        
                    case '001'
                        runs.runsToKeep = [1,3,4];
                        runs.runsToIgnore = 2;
                    case '002'
                        runs.runsToKeep = [1,2];
                        runs.runsToIgnore = [3,4];
                    case '003'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '004'
                        runs.runsToKeep = [1,2];
                        runs.runsToIgnore = [3,4];
                    case '005'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '009'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '012'
                        runs.runsToKeep = 2;
                        runs.runsToIgnore = [1,3,4];
                    case '015'
                        runs.runsToKeep = [3,4];
                        runs.runsToIgnore = [1,2];
                    case '018'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '020'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '022'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '027'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = [2,3,4]; % r3 saturated
                    case '032'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '036'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '038'
                        runs.runsToKeep = [1,3,4];
                        runs.runsToIgnore = 2;
                    case '042'
                        runs.runsToKeep = [1,2,3];
                        runs.runsToIgnore = 4;
                    case '044'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '046'
                        runs.runsToKeep = [1,3,4];
                        runs.runsToIgnore = 2;
                    case '047'
                        runs.runsToKeep = 3;
                        runs.runsToIgnore = [1,2,4];
                    case '048'
                        runs.runsToKeep = [1,3,4];
                        runs.runsToIgnore = 2;
                    case '052'
                        runs.runsToKeep = [2,4];
                        runs.runsToIgnore = [1,3];
                    case '053'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '054'
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '055'
                        runs.runsToKeep = [1,2];
                        runs.runsToIgnore = [3,4];
                    case '058'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '059'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '061'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = 2:4;
                    case '062'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '065'
                        runs.runsToKeep = [1,2,3];
                        runs.runsToIgnore = 4;
                    case '069'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '076'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '078'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '080'
                        runs.runsToKeep = [1,2,3];
                        runs.runsToIgnore = 4;
                    case '081'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '082'
                        runs.runsToKeep = 1:2;
                        runs.runsToIgnore = 3:4;
                    case '083'
                        runs.runsToKeep = [1,2,3];
                        runs.runsToIgnore = 4;
                    case '088'
                        runs.runsToKeep = [1,2,3];
                        runs.runsToIgnore = 4; % saturated run
                    case '095'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = [2,3,4];
                    case '097'
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '099'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '100'
                        runs.runsToKeep = [1,2];
                        runs.runsToIgnore = [3,4];
                end
                
                %% too much movement cleaning
            case 'fMRI_noMove_bis'
                switch sub_nm
                    case {'017','043','074'} % first run: fMRI crashed => we have the behavior but not enough trials for fMRI
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case {'008','022','024'} % subjects with too much movement in ALL runs
                        error([sub_nm,' should not be included under the condition ',condition,...
                            ' (too much movement in all runs).']);
                        % subjects with some runs with too much movement
                    case '021'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = 2:4;
                    case '029'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
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
                        runs.runsToIgnore = 2:4;
                    case '097'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = [2,3,4];
                    case '099'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                end % subject
                %% too much movement cleaning
            case 'fMRI_noMove_ter'
                switch sub_nm
                    case {'017','074'} % first run: fMRI crashed => we have the behavior but not enough trials for fMRI
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                        %% subjects with too much movement in ALL runs
                    case {'008','022','024','087'}
                        error([sub_nm,' should not be included under the condition ',condition,...
                            ' (too much movement in all runs).']);
                        %% subjects with some runs with too much movement
                    case '005'
                        runs.runsToKeep = 2;
                        runs.runsToIgnore = [1,3,4];
                    case '012'
                        runs.runsToKeep = [1,3,4];
                        runs.runsToIgnore = 2;
                    case '018'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '021'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = 2:4;
                    case '029'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '043' % first run: fMRI crashed => we have the behavior but not enough trials for fMRI
                        % + runs 2 and 4 borderline movement
                        runs.runsToKeep = 3;
                        runs.runsToIgnore = [1,2,4];
                    case '044'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '047'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '050'
                        runs.runsToKeep = [1,2,3];
                        runs.runsToIgnore = 4;
                    case '052'
                        runs.runsToKeep = [1,3,4];
                        runs.runsToIgnore = 2;
                    case '053'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '054'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '056'
                        runs.runsToKeep = [1,4];
                        runs.runsToIgnore = [2,3];
                    case '058'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '062'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = 2:4;
                    case '064'
                        runs.runsToKeep = [1,3,4];
                        runs.runsToIgnore = 2;
                    case '065'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '069'
                        runs.runsToKeep = [1,2];
                        runs.runsToIgnore = [3,4];
                    case '071'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '076'
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '078'
                        runs.runsToKeep = [2,3];
                        runs.runsToIgnore = [1,4];
                    case '079'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '080'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '083'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '086'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '090'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '093'
                        runs.runsToKeep = [2,4];
                        runs.runsToIgnore = [1,3];
                    case '094'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '095'
                        runs.runsToKeep = [1,2,3];
                        runs.runsToIgnore = 1;
                    case '097'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = [2,3,4];
                    case '099'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                end % subject
                %% remove saturation tasks AND runs with too movement
            case {'fMRI_noSatTask_noMove_bis','fMRI_noSatTask_noMove_bis_bayesianMdl'}
                switch sub_nm
                    case {'017','043','074'} % first run: fMRI crashed => we have the behavior but not enough trials for fMRI
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case {'008','022','024'} % subjects with too much movement in ALL runs
                        error([sub_nm,' should not be included under the condition ',condition,...
                            ' (too much movement in all runs).']);
                        % subjects with some runs with too much movement
                    case '047' % saturation all tasks
                        error(['subject ',sub_nm,' should not be included (saturated ALL tasks)']);
                    case '021'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = 2:4;
                    case '027' % Em full saturation
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '029'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '044'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '052' % Em full saturation
                        runs.runsToKeep = [2,4];
                        runs.runsToIgnore = [1,3];
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
                    case '069'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '071'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '076'
                        runs.runsToKeep = 3;
                        runs.runsToIgnore = [1,2,4];
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
                    case '095' % Ep full saturation
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '097'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = [2,3,4];
                    case '099'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                end % subject
                %% remove saturation runs AND runs with too movement
            case 'fMRI_noSatRun_noMove_bis'
                switch sub_nm
                    case {'017','043','074'} % first run: fMRI crashed => we have the behavior but not enough trials for fMRI
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case {'008','022','024'} % subjects with too much movement in ALL runs
                        error([sub_nm,' should not be included under the condition ',condition,...
                            ' (too much movement in all runs).']);
                    case {'047','097'} % too much movement OR saturation in all runs
                        error([sub_nm,' should not be included under the condition ',condition,...
                            ' (too much movement or saturation in all runs).']);
                        % subjects with some runs with too much movement or
                        % saturation runs
                    case '002'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '004'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '005'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '012'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '021'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = 2:4;
                    case '027'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = [2,3,4];
                    case '029'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '032'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '038'
                        runs.runsToKeep = [1,3,4];
                        runs.runsToIgnore = 2;
                    case '044'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = [2,3,4];
                    case '048'
                        runs.runsToKeep = [1,3,4];
                        runs.runsToIgnore = 2;
                    case '052'
                        runs.runsToKeep = [2,4];
                        runs.runsToIgnore = [1,3];
                    case '053'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '054'
                        runs.runsToKeep = 3;
                        runs.runsToIgnore = [1,2,4];
                    case '055'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '058'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '061'
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '062'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = 2:4;
                    case '069'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '071'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '076'
                        runs.runsToKeep = 3;
                        runs.runsToIgnore = [1,2,4];
                    case '078'
                        runs.runsToKeep = [2,3];
                        runs.runsToIgnore = [1,4];
                    case '080'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '081'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '082'
                        runs.runsToKeep = [1,2];
                        runs.runsToIgnore = [3,4];
                    case '083'
                        runs.runsToKeep = [1,2];
                        runs.runsToIgnore = [3,4];
                    case '087'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = [2,3,4];
                    case '088'
                        runs.runsToKeep = [1,2,3];
                        runs.runsToIgnore = 4;
                    case '095'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = [2,3,4];
                    case '099'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = [2,3,4];
                    case '100'
                        runs.runsToKeep = [1,2];
                        runs.runsToIgnore = [3,4];
                end % subject
                %% control that some subjects are not included + remove bad runs
            case 'behavior_noSatTaskSub'
                switch sub_nm
                    case {'027','047','052','069','076','095'}
                        error([sub_nm,' should not be included under the condition ',condition]);
                end
            case 'fMRI_noSatTaskSub'
                switch sub_nm
                    case {'017','043','074'} % first run: fMRI crashed => we have the behavior but not enough trials for fMRI
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case {'027','047','052','069','076','095'}
                        error([sub_nm,' should not be included under the condition ',condition]);
                end
            case 'fMRI_noMoveSub'
                switch sub_nm
                    case {'017','043','074'} % first run: fMRI crashed => we have the behavior but not enough trials for fMRI
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case {'008','022','024'}
                        error([sub_nm,' should not be included under the condition ',condition]);
                end
            case 'fMRI_noMoveSub_bis'
                switch sub_nm
                    case {'017','043','074'} % first run: fMRI crashed => we have the behavior but not enough trials for fMRI
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case {'008',...
                            '021','022','024','029',...
                            '044','047',...
                            '053','054','058',...
                            '062',...
                            '071','076','078',...
                            '080','083','087',...
                            '097','099'}
                        error([sub_nm,' should not be included under the condition ',condition]);
                end
            case 'fMRI_noMoveSub_ter'
                switch sub_nm
                    case {'017','074'} % first run: fMRI crashed => we have the behavior but not enough trials for fMRI
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case {'005','008',...
                            '012','018',...
                            '021','022','024','029',...
                            '040','043','044','047',...
                            '050','052','053','054','056','058',...
                            '062','064','065','069',...
                            '071','076','078','079',...
                            '080','083','086','087',...
                            '090','093','094','095','097','099'}
                        error([sub_nm,' should not be included under the condition ',condition]);
                end
            case {'fMRI_noSatTaskSub_noMoveSub','fMRI_noSatTaskSub_noMove_bis_Sub'}
                switch sub_nm
                    case {'017','043','074'} % first run: fMRI crashed => we have the behavior but not enough trials for fMRI
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case {'027','047','052','069','076','095',...
                            '008','022','024'}
                        error([sub_nm,' should not be included under the condition ',condition]);
                end
            case {'pulse'}
                switch sub_nm
                    case {'020','021','024','040','082','083','085','093'} % ignore completely those subjects
                        runs.runsToKeep = [];
                        runs.runsToIgnore = 1:4;
                    case '002'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '005'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '012'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '015'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '017'
                        runs.runsToKeep = [1,2];
                        runs.runsToIgnore = 3:4;
                    case '018'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = 2:4;
                    case '022'
                        runs.runsToKeep = 1:2;
                        runs.runsToIgnore = 3:4;
                    case '034'
                        runs.runsToKeep = [1,3,4];
                        runs.runsToIgnore = 2;
                    case '043'
                        runs.runsToKeep = [1,3,4];
                        runs.runsToIgnore = 2;
                    case '044'
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '047'
                        runs.runsToKeep = [1,4];
                        runs.runsToIgnore = 2:3;
                    case '050'
                        runs.runsToKeep = 2;
                        runs.runsToIgnore = [1,3,4];
                    case '053'
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '075'
                        runs.runsToKeep = [1,2];
                        runs.runsToIgnore = 3:4;
                    case '078'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = 2:4;
                    case '091'
                        runs.runsToKeep = [1,3,4];
                        runs.runsToIgnore = 2;
                    case '097'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '099'
                        runs.runsToKeep = [1,3,4];
                        runs.runsToIgnore = 2;
                    case '100'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                end % subject loop
            case 'pulse_noSatRun'
                switch sub_nm
                    case {'020','021','024','040','047','082','083','085','093'} % ignore completely those subjects
                        runs.runsToKeep = [];
                        runs.runsToIgnore = 1:4;
                    case '002'
                        runs.runsToKeep = 1:2;
                        runs.runsToIgnore = 3:4;
                    case '004'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '005'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '012'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '015'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '017'
                        runs.runsToKeep = [1,2];
                        runs.runsToIgnore = 3:4;
                    case '018'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = 2:4;
                    case '022'
                        runs.runsToKeep = 1:2;
                        runs.runsToIgnore = 3:4;
                    case '027'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = [2,3,4];
                    case '032'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '034'
                        runs.runsToKeep = [1,3,4];
                        runs.runsToIgnore = 2;
                    case '038'
                        runs.runsToKeep = [1,3,4];
                        runs.runsToIgnore = 2;
                    case '043'
                        runs.runsToKeep = [1,3,4];
                        runs.runsToIgnore = 2;
                    case '044'
                        runs.runsToKeep = [2,4];
                        runs.runsToIgnore = [1,3];
                    case '048'
                        runs.runsToKeep = [1,3,4];
                        runs.runsToIgnore = 2;
                    case '050'
                        runs.runsToKeep = 2;
                        runs.runsToIgnore = [1,3,4];
                    case '052'
                        runs.runsToKeep = [2,4];
                        runs.runsToIgnore = [1,3];
                    case '053'
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '054'
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '055'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '058'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '061'
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '062'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '069'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '075'
                        runs.runsToKeep = [1,2];
                        runs.runsToIgnore = 3:4;
                    case '076'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '078'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = 2:4;
                    case '081'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '088'
                        runs.runsToKeep = [1,2,3];
                        runs.runsToIgnore = 4;
                    case '091'
                        runs.runsToKeep = [1,3,4];
                        runs.runsToIgnore = 2;
                    case '095'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = [2,3,4];
                    case '097'
                        runs.runsToKeep = [2,4];
                        runs.runsToIgnore = [1,3];
                    case '099'
                        runs.runsToKeep = [1,4];
                        runs.runsToIgnore = 2:3;
                    case '100'
                        runs.runsToKeep = 1:2;
                        runs.runsToIgnore = 3:4;
                end % subject loop
            case 'fMRI_noSatTaskSub_noMoveSub_noSatRun_noMoveRun' % remove saturated runs + movement runs
                switch sub_nm
                    case {'017','043','074'} % first run: fMRI crashed => we have the behavior but not enough trials for fMRI
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '047'
                        error(['subject ',sub_nm,' should not be included (saturated ALL tasks)']);
                    case {'008','022','024'}
                        error(['subject ',sub_nm,' should not be included (moved in ALL runs)']);
                    case '097'
                        error(['subject ',sub_nm,' should be removed (either moved or saturated every run']);
                    case '002'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '004'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '005'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '012'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '027'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = [2,3,4];
                    case '032'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '038'
                        runs.runsToKeep = [1,3,4];
                        runs.runsToIgnore = 2;
                    case '044'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = [2,3,4];
                    case '048'
                        runs.runsToKeep = [1,3,4];
                        runs.runsToIgnore = 2;
                    case '052'
                        runs.runsToKeep = [2,4];
                        runs.runsToIgnore = [1,3];
                    case '053'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '054'
                        runs.runsToKeep = 3;
                        runs.runsToIgnore = [1,2,4];
                    case '055'
                        runs.runsToKeep = 1:3;
                        runs.runsToIgnore = 4;
                    case '058'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '061'
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '062'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '065'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '069'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '071'
                        runs.runsToKeep = [1,3];
                        runs.runsToIgnore = [2,4];
                    case '076'
                        runs.runsToKeep = 3;
                        runs.runsToIgnore = [1,2,4];
                    case '078'
                        runs.runsToKeep = 3;
                        runs.runsToIgnore = [1,2,4];
                    case '080'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '081'
                        runs.runsToKeep = [1,2,4];
                        runs.runsToIgnore = 3;
                    case '082'
                        runs.runsToKeep = [1,2];
                        runs.runsToIgnore = [3,4];
                    case '083'
                        runs.runsToKeep = [1,2];
                        runs.runsToIgnore = 3:4;
                    case '087'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = 2:4;
                    case '088'
                        runs.runsToKeep = [1,2,3];
                        runs.runsToIgnore = 4;
                    case '093'
                        runs.runsToKeep = 2:4;
                        runs.runsToIgnore = 1;
                    case '095'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = [2,3,4];
                    case '099'
                        runs.runsToKeep = 1;
                        runs.runsToIgnore = 2:4;
                    case '100'
                        runs.runsToKeep = [1,2];
                        runs.runsToIgnore = [3,4];
                end
            otherwise % condition not ready
                error(['condition ',condition,' is not ready yet for runs_definition.m']);
        end
    otherwise
        error('case not ready yet');
end % study

%% control that things were well defined
for iR = 1:4
    if ismember(iR, runs.runsToKeep) && ismember(iR, runs.runsToIgnore)
        error(['run ',num2str(iR),' is in runs.runsToKeep AND in runs.runsToIgnore']);
    elseif ~ismember(iR, runs.runsToKeep) && ~ismember(iR, runs.runsToIgnore)
        error(['run ',num2str(iR),' is NOT in runs.runsToKeep NOR in runs.runsToIgnore => has to be included']);
    end % filter
end % run loop

%% extract index for each task of run kept
Ep_runsToKeep = strcmp(runs.tasks,'Ep').*ismember(1:4,runs.runsToKeep) == 1;
Em_runsToKeep = strcmp(runs.tasks,'Em').*ismember(1:4,runs.runsToKeep) == 1;
runs.Ep.runsToKeep = find(Ep_runsToKeep);
runs.Em.runsToKeep = find(Em_runsToKeep);

%% update task types depending on the runs to keep
runs.tasks = runs.tasks(runs.runsToKeep);

%% extract number of runs of each task type
runs.nb_runs.Ep = sum(strcmp(runs.tasks,'Ep'));
runs.nb_runs.Em = sum(strcmp(runs.tasks,'Em'));

%% extract number of runs
n_runs = length(runs.tasks);

end % function