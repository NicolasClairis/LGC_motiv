% script to compute behavioral sensitivities

%% clean workspace before starting
sca;
clearvars;
close all;

%% working directories
% launch within the folder where scripts are stored or will not work

main_folder                 = [pwd]; % you have to be sure that you are in the correct path when you launch the script
results_folder              = ['D:\LGC_motiv\LGC_Motiv_results\fMRI_pilots'];
% BioPac_folder               = [main_folder, 'BioPac_functions' filesep];
% pics_folder                 = [main_task_folder, 'Coin_PNG', filesep];
Matlab_DIY_functions_folder = [main_folder, 'Matlab_DIY_functions', filesep];

% add personal functions (needed for PTB opening at least)
addpath(genpath(main_folder));
% addpath(BioPac_folder);
addpath(results_folder);

%% set run parameters

n_trialsPerSession = 48;
f_fname = [];
g_fname = @g_observation;
n_G_prm = 6;
n_F_prm = 0;
n_hiddenStates = 0;

%% compute for subjects/pilots

    all_files = dir('*task.mat');
    nb_files = size(all_files,1);

    % extract ID of subjects (since we don't start at 1)
 for i_files = 1:length(all_files)
     %save the names for later
    all_files_name(i_files) = {all_files(i_files).name};
    name_tmp = all_files(i_files).name;
    % find each unique subject ID
    if name_tmp(6) ~= '1' || name_tmp(6) ~= '2' || name_tmp(6) ~= '3' || name_tmp(6) ~= '4' || name_tmp(6) ~= '5' || name_tmp(6) ~= '6' || name_tmp(6) ~= '7' || name_tmp(6) ~= '8' || name_tmp(6) ~= '9'
    ID_tmp(i_files) = str2double(name_tmp(5));
    else
        ID_tmp(i_files) = str2double(name_tmp(5:6));
    end
 end
%find all the ID of participants and how many sessions they did
unique_ID = unique(ID_tmp)
nb_sessions_per_subject = histc(ID_tmp(:), unique_ID)

    for i_files = 1:length(unique_ID)
        % empty and prepare another var, in case it's size was not the same between subject (miss a session)
        var = NaN(5,48);
        % prepare also array of all the responses from participants, and it's increment
        all_choices = [];
        i_choice = 1;
        
        % select files for 1 specific subject
        for m_files = 1:length(all_files_name)
            name_tmp =all_files_name{1,m_files};
            selected_name_tmp = strcat('s',num2str(unique_ID(i_files)),'_');
            ID_position = strfind(name_tmp,selected_name_tmp);
            if strcmp(name_tmp(ID_position:ID_position+2),selected_name_tmp)
                selected_name(m_files) = convertCharsToStrings(name_tmp);
            end
        end
        
        % for each session
        for i_session = 1:nb_sessions_per_subject(i_files)
            session_name = char(selected_name(1,i_session));
            physical_effort_ID = strfind(session_name,'physical');
            mental_effort_ID = strfind(session_name,'mental');
            % load the appropriate folder
            load(session_name)
            
            if strcmp(session_name(physical_effort_ID:physical_effort_ID+7),'physical')
                
                % prepare the choices with the chosen values
                choicePhysical = physicalPerf.choice;
                for n_choice = 1:length(choicePhysical)
                    if choicePhysical(i_choice) == -2
                        all_choices(n_choice) = 0.125;
                    elseif choicePhysical(n_choice) == -1
                        all_choices(i_choice) = 0.375;
                    elseif choicePhysical(n_choice) == 1
                        all_choices(i_choice) = 0.625;
                    elseif choicePhysical(n_choice) == 2
                        all_choices(i_choice) = 0.875;
                    end
                    i_choice = i_choice +1;
                end
                
                % prepare all the var informations
                R_or_P = physicalPerf.choiceOptions.R_or_P;
                deltaR = physicalPerf.choiceOptions.monetary_amount.left - physicalPerf.choiceOptions.monetary_amount.right;
                deltaE = physicalPerf.choiceOptions.E.left - physicalPerf.choiceOptions.E.right;
                RP_trials = NaN(1,48);
                Ep_or_Em_trials = ones(n_trialsPerSession,1)';
                for iTrial = 1:n_trialsPerSession
                    switch R_or_P{iTrial}
                        case 'P'
                            RP_trials(iTrial) = 0;
                        case 'R'
                            RP_trials(iTrial) = 1;
                    end
                end
                
                % prepare the matrix of our variables, for each sessions
                var(:,:,i_session) = [deltaR; deltaE; RP_trials; Ep_or_Em_trials;(1:n_trialsPerSession)-1];
                
            elseif strcmp(session_name(mental_effort_ID:mental_effort_ID+5),'mental')

                % prepare the choices with the chosen values
                choiceMental_1 = mentalE_perf.choice;
                for n_choice = 1:length(choiceMental_1)
                    if choiceMental_1(n_choice) == -2
                        all_choices(i_choice) = 0.125;
                    elseif choiceMental_1(n_choice) == -1
                        all_choices(i_choice) = 0.375;
                    elseif choiceMental_1(n_choice) == 1
                        all_choices(i_choice) = 0.625;
                    elseif choiceMental_1(n_choice) == 2
                        all_choices(i_choice) = 0.875;
                    end
                    i_choice = i_choice+1;
                end
                
                % prepare all the var informations
                R_or_P = mentalE_perf.choiceOptions.R_or_P;
                deltaR = mentalE_perf.choiceOptions.monetary_amount.left - mentalE_perf.choiceOptions.monetary_amount.right;
                deltaE = mentalE_perf.choiceOptions.E.left - mentalE_perf.choiceOptions.E.right;
                RP_trials = NaN(1,48);
                Ep_or_Em_trials = ones(n_trialsPerSession,1)'-1;
                for iTrial = 1:n_trialsPerSession
                    switch R_or_P{iTrial}
                        case 'P'
                            RP_trials(iTrial) = 0;
                        case 'R'
                            RP_trials(iTrial) = 1;
                    end
                end
                
                % prepare the matrix of our variables
                var(:,:,i_session) = [deltaR; deltaE; RP_trials; Ep_or_Em_trials;(1:n_trialsPerSession)-1];        
            end
            
        end
        % reshape it as one big session
        var = reshape(permute(var,[2 3 1]),size(var,2)*nb_sessions_per_subject(i_files),[])';
        
        % select options for the VBA toolbox
            options = struct;
            options.sources.type = 0; % binary data
            options.DisplayWin  = 0; % display figure during inversion
            options.verbose     = 0; % display text during inversion (1) or not (0)
            % define number of parameters to estimate
            options.dim = struct('n', n_hiddenStates,... % number of hidden states
                'n_t', size(var,2),... % number of trials across runs
                'n_theta',n_F_prm,... % number of evolution parameters
                'n_phi',n_G_prm); % number of observation parameters
            options.priors.muPhi = zeros(n_G_prm, 1); %(moyenne des priors sur les paramètres)
            options.priors.SigmaPhi = eye(n_G_prm)/2; % (matrice de covariance des paramètres)
         
            % compute sensitivities and save them
        [posterior,out] = VBA_NLStateSpaceModel(all_choices, var, f_fname, g_fname, options.dim, options);
        sensitivities(:,i_files) = posterior.muPhi;
        
    end

%% plot results
sensitivities

%         %% plot results, works for 4 sessions only
%
%         if length(session) == 1
%             figure()
%             x = 1:n_G_prm;
%
%             % extract posterior for each participant
%             for j = 1:nb_run
%                 data_tmp = posterior_session_run(1,j).muPhi';
%                 data_tmp2(j,:) = data_tmp ;
%             end
%             % compute SEM and SD
%             err_sem = squeeze(std(data_tmp2(:,:))./sqrt(size(data_tmp2(:,:),1)));
%             err = squeeze(std(data_tmp2(:,:)));
%             %     if wanted, compute correlation between ''real'' and posterior data
%             for param_i = 1:n_G_prm
%                 corr_tmp = corrcoef(data_tmp2(:,param_i),parameterPhi(param_i,:)');
%                 corr_mat(1,param_i) =corr_tmp(1,2);
%             end
%
%             % plot the average over N simulation.
%             bar(x,mean(squeeze(data_tmp2(:,:))))
%             mean(squeeze(data_tmp2(:,:)))
%             hold on
%
%             er = errorbar(x,mean(squeeze(data_tmp2(:,:))),err,err);
%             er.Color = [0 0 0];
%             er.LineStyle = 'none';
%             er = errorbar(x,mean(squeeze(data_tmp2(:,:))),err_sem,err_sem);
%             er.Color = [1 0 0];
%             er.LineStyle = 'none';
%             set(gca,'XTick',[1 2 3 4])
%             set(gca,'XTickLabel',{'kR','kP','kEp','kEm'})
%             title(['Mean ',num2str(i_run),' simulations with ',num2str(session),' session of efforts'])
%             ylim
%             xlim=get(gca,'xlim');
%             hold on
%             plot(xlim,[0.3 0.3],'m')
%             plot(xlim,[0.6 0.6],'c')
%             %     ylim([0 0.9])
%
%             figure()
%             % plot correlations between ''real'' and posterior for each nb of sessions
%             for j = 1:nb_run
%
%                 data_tmp = posterior_session_run(1,j).muPhi';
%                 data_tmp2(j,:) = data_tmp ;
%             end
%
%             for n = 1:number_of_parameters
%                 for m = 1:number_of_parameters
%                     corrmat_real_tmp = corrcoef(data_tmp2(:,n)',parameterPhi(m,:));
%                     autocorrmat_tmp = corrcoef(data_tmp2(:,n),data_tmp2(:,m));
%                     corrmat(n,m) = corrmat_real_tmp(1,2);
%                     autocorrmat(n,m) = autocorrmat_tmp(1,2);
%                 end
%             end
%
%             % correlation real-estimated plot
%             imagesc(corrmat)
%             if number_of_parameters == 4
%                 set(gca,'XTick',[1 2 3 4])
%                 set(gca,'XTickLabel',{'kR','kP','kEp','kEm'})
%                 set(gca,'YTick',[1 2 3 4])
%                 set(gca,'YTickLabel',{'kR','kP','kEp','kEm'})
%             elseif number_of_parameters == 6
%                 set(gca,'XTick',[1 2 3 4 5 6])
%                 set(gca,'XTickLabel',{'kR','kP','kEp','kEm','kFp','kFm'})
%                 set(gca,'YTick',[1 2 3 4 5 6])
%                 set(gca,'YTickLabel',{'kR','kP','kEp','kEm','kFp','kFm'})
%             end
%             ylabel('"real" parameter')
%             xlabel('estimated parameter')
%             title([num2str(4),' blocks, correlation: estimated and real sens'])
%             colorbar
%
%
%             % autocorrelation plot
%             figure()
%
%             imagesc(autocorrmat)
%             if number_of_parameters == 4
%                 set(gca,'XTick',[1 2 3 4])
%                 set(gca,'XTickLabel',{'kR','kP','kEp','kEm'})
%                 set(gca,'YTick',[1 2 3 4])
%                 set(gca,'YTickLabel',{'kR','kP','kEp','kEm'})
%             elseif number_of_parameters == 6
%                 set(gca,'XTick',[1 2 3 4 5 6])
%                 set(gca,'XTickLabel',{'kR','kP','kEp','kEm','kFp','kFm'})
%                 set(gca,'YTick',[1 2 3 4 5 6])
%                 set(gca,'YTickLabel',{'kR','kP','kEp','kEm','kFp','kFm'})
%             end
%             ylabel('"real" parameter')
%             xlabel('estimated parameter')
%             title([num2str(4),' blocks, autocorrelation: posterior K sens'])
%             colorbar
%             mean([corrmat(1,1),corrmat(2,2),corrmat(3,3),corrmat(4,4),corrmat(5,5),corrmat(6,6)])
%         end
%     end
% end