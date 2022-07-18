function[prm] = prm_extraction(subject_id)
% [prm] = prm_extraction(subject_id)
% prm_extraction will extract the behavioural parameters for the list of
% subjects entered in input. You will need to select whether you want to
% use a bayesian or a more simple model approach and which model number you
% want to use.
%
% INPUTS
% subject_id: list of subject identification numbers
%
% OUTPUTS
% prm: structure with parameters extracted for each participant

%% main parameters
NS = length(subject_id);

%% bayesian across tasks or simple model each task separately?
listPossibleConditions = {'bayesian','simple'};
mdlType_idx = listdlg('promptstring','Which model type?',...
    'ListString',listPossibleConditions);
mdlType = listPossibleConditions{mdlType_idx};

%% which model number to use?
switch mdlType
    case 'bayesian'
        listPossibleConditions = {'3'};
    case 'simple'
        listPossibleConditions = {'1','2','3','4'};
end
mdlN_idx = listdlg('promptstring','Which model number?',...
    'ListString',listPossibleConditions);
mdlN = listPossibleConditions{mdlN_idx};

%% extract parameters of the selected model
switch mdlType
    case 'bayesian'
        error('go in correct path');
        bayesian_mdl = getfield(load('behavioral_prm_tmp.mat',...
            ['bayesian_mdl',mdlN]),...
            ['bayesian_mdl',mdlN]);
        [prm.kR, prm.kP,...
            prm.kEp, prm.kFp,...
            prm.kEm, prm.kFm] = deal(NaN(1,NS));
        for iS = 1:NS
            sub_nm = subject_id{iS};
            % extract parameters
            sub_idx = strcmp(bayesian_mdl.subject_id, sub_nm);
            if sum(sub_idx == 1)
                prm.kR(iS) = bayesian_mdl.kR(sub_idx);
                prm.kP(iS) = bayesian_mdl.kP(sub_idx);
                prm.kEp(iS) = bayesian_mdl.kEp(sub_idx);
                prm.kFp(iS) = bayesian_mdl.kFp(sub_idx);
                prm.kEm(iS) = bayesian_mdl.kEm(sub_idx);
                prm.kFm(iS) = bayesian_mdl.kFm(sub_idx);
            end % filter if subject extracted by Arthur
        end % subject list
    case 'simple'
        %% perform behavioral model
        figDispGroup = 0;
        dispMoneyOrLevels = 'levels';
        [betas_fullList, pvalues_fullList] = logitfit_choices_group(figDispGroup, dispMoneyOrLevels);
        
        mdl_nm = ['mdl_',mdlN];
        % tasks
        tasks = {'Ep','Em'};
        nTasks = length(tasks);
        tasks_ids = {'p','m'};
        % initialize vars of interest
        for iTask = 1:nTasks
            task_id = tasks_ids{iTask};
            switch mdlN
                case '1'
                    [prm.(['kM',task_id]),...
                        prm.(['kE',task_id])] = deal(NaN(1,NS));
                case '2'
                    [prm.(['kR',task_id]), prm.(['kP',task_id]),...
                        prm.(['kE',task_id])] = deal(NaN(1,NS));
                case '3'
                    [prm.(['kM',task_id]),...
                        prm.(['kE',task_id]),...
                        prm.(['kF',task_id])] = deal(NaN(1,NS));
                case '4'
                    [prm.(['kR',task_id]), prm.(['kP',task_id]),...
                        prm.(['kE',task_id]),...
                        prm.(['kF',task_id])] = deal(NaN(1,NS));
            end % model
        end % task loop
        for iS = 1:NS
            sub_nm = subject_id{iS};
            for iTask = 1:nTasks
                task_nm = tasks{iTask};
                task_id = tasks_ids{iTask};
                sub_idx = strcmp(betas_fullList.subList,sub_nm);
                % models pooling R+P or splitting them
                switch mdlN
                    case {'1','3'}
                        prm.(['kM',task_id])(iS) = betas_fullList.(task_nm).(mdl_nm).kMoney(sub_idx);
                    case {'2','4'}
                        prm.(['kR',task_id])(iS) = betas_fullList.(task_nm).(mdl_nm).kR(sub_idx);
                        prm.(['kP',task_id])(iS) = betas_fullList.(task_nm).(mdl_nm).kP(sub_idx);
                end
                prm.(['kE',task_id])(iS) = betas_fullList.(task_nm).(mdl_nm).kEffort(sub_idx);
                % models with vs models without fatigue
                switch mdlN
                    case {'3','4'}
                        prm.(['kF',task_id])(iS) = betas_fullList.(task_nm).(mdl_nm).kFatigue(sub_idx);
                end
            end % task loop
        end % subject loop
end % bayesian/behavioral

% store corresponding subject id
prm.CID = subject_id;

end % function