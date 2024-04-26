% script to rebuild NV and p(choice = high E) based on parameters for
% model 4

%% general parameters
study_nm = 'study1';
nTrialsPerRun = 54;
nTotalRuns = 4;
condition = 'behavior';

%% working directories
whichPc = 'lab';
switch whichPc
    case 'lab'
        gitRoot = fullfile('C:','Users','clairis','Desktop');
    case 'home'
        gitRoot = fullfile('C:','Users','Loco','Documents');
end
gitPath = [gitRoot, filesep, 'GitHub',filesep,...
    'LGC_motiv',filesep,'LGC_Motiv_results',filesep,...
    study_nm, filesep,...
    'bayesian_modeling', filesep];
root = ['E:',filesep,study_nm,filesep];

%%
bayesian_deltaNV = getfield(load([gitPath, 'bayesian_deltaNV_data.mat'],'bayesian_deltaNV'),'bayesian_deltaNV');
bayesian_pChoice = getfield(load([gitPath, 'bayesian_pChoice_data.mat'],'bayesian_pChoice'),'bayesian_pChoice');

%% load parameters
behavioral_prm_tmp = load([gitPath,'behavioral_prm.mat']);
bayesian_mdl = behavioral_prm_tmp.bayesian_mdl;
bhv_prm = bayesian_mdl.mdl_4;
subject_id = bayesian_mdl.mdl_4.subject_id;
NS = length(subject_id);
bayesian_deltaNV.mdl_4.subject_id = subject_id;
bayesian_pChoice.mdl_4.subject_id = subject_id;

%% loop over subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_nm_bis = ['CID',sub_nm];
    sub_bhv_folder = [root, sub_nm_bis,filesep, 'behavior',filesep];
    % initialize data
    for iR = 1:nTotalRuns
        [bayesian_deltaNV.mdl_4.(sub_nm_bis).(['run',num2str(iR)]),...
            bayesian_pChoice.mdl_4.(sub_nm_bis).(['run',num2str(iR)])] = deal(NaN(nTrialsPerRun,1));
    end % run loop
    
    [runs, n_runs] = runs_definition(study_nm, sub_nm, condition);
    
    % extract parameter
    kR = bhv_prm.kR(iS);
    kP = bhv_prm.kP(iS);
    kEp = bhv_prm.kEp(iS);
    kEm = bhv_prm.kEm(iS);
    kFp = bhv_prm.kFp(iS);
    kLm = bhv_prm.kLm(iS);
    kBiasM = bhv_prm.kBiasM(iS);
    
    for iR = 1:n_runs
        jR = runs.runsToKeep(iR);
        run_nm = num2str(jR);
        run_task_nm = runs.tasks{iR};
        switch run_task_nm
            case 'Ep'
                task_fullName = 'physical';
            case 'Em'
                task_fullName = 'mental';
        end % task name
        % load data
        sub_loadStruct = load([sub_bhv_folder,...
            sub_nm_bis,'_session',num2str(iR),...
            '_',task_fullName,'_task.mat']);
        choice_options = sub_loadStruct.choice_opt;
        dMoney_tmp = abs(choice_options.monetary_amount.left - choice_options.monetary_amount.right);
        dE_tmp = abs(choice_options.E.left - choice_options.E.right);
        RP_tmp = strcmp(choice_options.R_or_P,'R');
        [deltaR_money] = extract_deltaR_money(sub_bhv_folder, sub_nm, run_nm, task_fullName);
        [deltaP_money] = extract_deltaP_money(sub_bhv_folder, sub_nm, run_nm, task_fullName);
        Ep_task = strcmp(run_task_nm,'Ep');
        Em_task = strcmp(run_task_nm,'Em');
        switch run_task_nm
            case 'Em'
                [~,...
                    ~,...
                    prevEfficacy,...
                    ~] = extract_mental_previous_efficacy(sub_bhv_folder, sub_nm, run_nm, task_fullName); % Arthur model: used n.correct/Time (with time starting after 2 first)
                Fp = 0;
            case 'Ep'
                [Fp] = extract_physical_fatigue(sub_bhv_folder, sub_nm, run_nm, task_fullName);
                prevEfficacy = 0;
        end
        
        dNV_tmp = (kR.*deltaR_money + kP.*deltaP_money) +...
            -(Ep_task.*dE_tmp.*(kEp + kFp.*Fp) +...
            Em_task.*dE_tmp.*(kEm - kLm.*prevEfficacy)) +...
            kBiasM;
        error('Need to check which fatigue measure Arthur used, now values are weird after trial 1 for Ep');
        pChoice_tmp = 1./(1 + exp(-dNV_tmp));
        bayesian_deltaNV.mdl_4.(sub_nm_bis).(['run',num2str(iR)])(:) = dNV_tmp';
        bayesian_pChoice.mdl_4.(sub_nm_bis).(['run',num2str(iR)])(:) = pChoice_tmp';
    end % run loop
    
end % subject loop


% %% make one global save
% save([gitPath, 'bayesian_deltaNV_data.mat'],'bayesian_deltaNV');
% save([gitPath, 'bayesian_pChoice_data.mat'],'bayesian_pChoice');