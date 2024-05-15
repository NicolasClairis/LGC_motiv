function[muPhi, sigmaPhi] = avg_g_observation_mdl_prm(mdl_n, condition, NS)
% [muPhi, sigmaPhi] = avg_g_observation_mdl_prm(mdl_n, condition, NS)
% avg_g_observation_mdl_prm will output the average and sigma for the Phi
% parameters estimated by the g_observation_mdl function with the
% VBA_toolbox for the selected subjects and model.
%
% INPUTS
% mdl_n: model number
%
% condition: condition
%
% NS: number of subjects
% 
% OUTPUTS
% muPhi: average value for each parameter in the following order:
% kR/kP/kEp/kEm/kBias/kFp/kLm (average posteriors)
%
% sigmaPhi: SD for each parameter in the following order:
% kR/kP/kEp/kEm/kBias/kFp/kLm

%% check inputs
[mdl_n, mdl_n_nm] = which_bayesian_mdl_n(mdl_n);
if ~exist('condition','var') || isempty(condition)
    condition = subject_condition;
end
if ~exist('NS','var') || isempty(NS)
    error('need to provide a value for NS');
end

%% prepare size
[mdl_prm] = computational_mdl_prm(mdl_n);
n_G_prm = mdl_prm.n_G_prm;
[muPhi, sigmaPhi] = deal(NaN(1,n_G_prm));

%% extract data
switch mdl_n
    case 5 % kR, kP, kEp, kEm, kBias, kFp, kLm
        switch condition
            case 'behavior_noSatTaskSub_noSatRun'
                switch NS
                    case 67
                        muPhi(1) = ;
                        muPhi(2) = ;
                        muPhi(3) = ;
                        muPhi(4) = ;
                        muPhi(5) = ;
                        muPhi(6) = ;
                        muPhi(7) = ;
                        sigmaPhi(1) = ;
                        sigmaPhi(2) = ;
                        sigmaPhi(3) = ;
                        sigmaPhi(4) = ;
                        sigmaPhi(5) = ;
                        sigmaPhi(6) = ;
                        sigmaPhi(7) = ;
                    otherwise
                        error(['NS = ',num2str(NS),' not ready yet']);
                end
            case 'fMRI_noSatTaskSub_noSatRun'
                muPhi = [];
                sigmaPhi = [];
            otherwise
                error(['condition ',condition,' not ready yet']);
        end
    otherwise
        error(['model ',mdl_n_nm,' not ready.']);
end

end % function