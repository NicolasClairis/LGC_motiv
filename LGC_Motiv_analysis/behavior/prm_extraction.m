function[prm, mdlType, mdlN] = prm_extraction(subject_id, mdlType, mdlN)
% [prm, mdlType, mdlN] = prm_extraction(subject_id, mdlType, mdlN)
% prm_extraction will extract the behavioural parameters for the list of
% subjects entered in input. You will need to select whether you want to
% use a bayesian or a more simple model approach and which model number you
% want to use.
%
% INPUTS
% subject_id: list of subject identification numbers
%
% mdlType: 'bayesian' or 'simple' model (will be asked if not defined in
% the inputs)
%
% mdlN: string with model number to be used (will be asked if not defined in
% the inputs)
%
% OUTPUTS
% prm: structure with parameters extracted for each participant
%
% mdlType: 'bayesian' or 'simple' model (in case was asked by the script
% allows to remember which one was selected)
%
% mdlN: string with model number to be used (in case was asked by the script
% allows to remember which one was selected)

%% working directories
whichPc = 'lab'; % 'lab'/'home'
switch whichPc
    case 'lab'
        bayesian_root = fullfile('C:','Users','clairis','Desktop','GitHub',...
            'LGC_motiv','LGC_Motiv_analysis','behavior');
        computerRoot = [fullfile('E:'),filesep];
    case 'home'
        error('paths to be updated');
    otherwise
        error(['paths not ready yet for ',whichPc, 'pc']);
end

%% main parameters
NS = length(subject_id);

%% bayesian across tasks or simple model each task separately?
if ~exist('mdlType','var') || isempty(mdlType)
    listPossibleModels = {'bayesian','simple'};
    mdlType_idx = listdlg('promptstring','Which model type?',...
        'ListString',listPossibleModels);
    mdlType = listPossibleModels{mdlType_idx};
end

%% which model number to use?
if ~exist('mdlN','var') || isempty(mdlN)
    switch mdlType
        case 'bayesian'
            listPossibleModelNumbers = {'3'};
        case 'simple'
            listPossibleModelNumbers = {'1','2','3','4'};
    end
    mdlN_idx = listdlg('promptstring','Which model number?',...
        'ListString',listPossibleModelNumbers);
    mdlN = listPossibleModelNumbers{mdlN_idx};
end

%% extract parameters of the selected model
switch mdlType
    case 'bayesian'
        bayesian_mdl = getfield(load([bayesian_root,filesep,...
            'behavioral_prm_tmp.mat'],...
            ['bayesian_mdl',mdlN]),...
            ['bayesian_mdl',mdlN]);
        [prm.kR, prm.kP,...
            prm.kEp, prm.kFp,...
            prm.kEm, prm.kFm] = deal(NaN(1,NS));
        parameter_names = fieldnames(prm);
        for iS = 1:NS
            sub_nm = subject_id{iS};
            % extract parameters
            sub_idx = strcmp(bayesian_mdl.subject_id, sub_nm);
            if sum(sub_idx == 1)
                for iBPrm = 1:length(parameter_names)
                    bayesian_prm_nm = parameter_names{iBPrm};
                    prm.(bayesian_prm_nm)(iS) = bayesian_mdl.(bayesian_prm_nm)(sub_idx);
                end
            end % filter if subject extracted by Arthur
        end % subject list
    case 'simple'
        %% perform behavioral model
        figDispGroup = 0;
        dispMoneyOrLevels = 'levels';
        [betas_fullList, pvalues_fullList] = logitfit_choices_group(subject_id, computerRoot, figDispGroup, dispMoneyOrLevels);
        
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
                otherwise
                    error(['model ',mdlN,' not ready']);
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