function[x_regs, reg_names] = RT_GLM_regs(GLMprm, potentialRegs)
% [x_regs, reg_names] = RT_GLM_regs(GLMprm, potentialRegs)
%
% INPUTS
% GLMprm: structure with GLM parameters
%
% potentialRegs: structure with all potential regressors
%
% OUTPUTS
% x_regs: matrix with all the regressors to include in the GLM
%
% reg_names: cell with the name of each regressor included in the order
% they were included in x_regs
%
% See also RT_GLM and which_RT_GLM.

%% extract regressors according to GLM parameters
potentialRegressors = fieldnames(GLMprm.regs);
nPotRegs = length(potentialRegressors);
jRegs = 0;
reg_names = {};
x_regs = [];

%% run constants
for iReg = 1:nPotRegs
    reg_nm = potentialRegressors{iReg};
    switch GLMprm.regs.(reg_nm)
        case 'on'
            switch reg_nm
                case 'run_cstt'
                    % add one constant per run
                    % run 1
                    jRegs = jRegs + 1;
                    reg_names{jRegs} = 'run1_cstt';
                    x_regs(:, jRegs) = potentialRegs.(reg_names{jRegs});
                    % run 2
                    jRegs = jRegs + 1;
                    reg_names{jRegs} = 'run2_cstt';
                    x_regs(:, jRegs) = potentialRegs.(reg_names{jRegs});
                    % run 3
                    jRegs = jRegs + 1;
                    reg_names{jRegs} = 'run3_cstt';
                    x_regs(:, jRegs) = potentialRegs.(reg_names{jRegs});
                    % run 4
                    jRegs = jRegs + 1;
                    reg_names{jRegs} = 'run4_cstt';
                    x_regs(:, jRegs) = potentialRegs.(reg_names{jRegs});
                case 'task_cstt'
                    % mental
                    jRegs = jRegs + 1;
                    reg_names{jRegs} = 'Em_cstt';
                    x_regs(:, jRegs) = potentialRegs.(reg_names{jRegs});
                    % physical
                    jRegs = jRegs + 1;
                    reg_names{jRegs} = 'Ep_cstt';
                    x_regs(:, jRegs) = potentialRegs.(reg_names{jRegs});
                case {'conf','confRtg',...
                        'RP','deltaMoney','deltaR','deltaP',...
                        'deltaEffort','deltaEp','deltaEm',...
                        'trialN','RT_raw_prevTrial',...
                        'physical_Fatigue','mental_Facilitation'} % only 1 regressor
                    jRegs = jRegs + 1;
                    reg_names{jRegs} = reg_nm;
                    x_regs(:, jRegs) = potentialRegs.(reg_nm);
                otherwise
                    error([reg_nm,' regressor not ready']);
            end % which regressor
        case 'off'% nothing to do
        case ''
            error(['GLMprm.regs.',reg_nm,' is empty.']);
        otherwise
            error(['case GLMprm.regs.',reg_nm,' = ',GLMprm.regs.(reg_nm),' not ready yet']);
    end
end % regressor loop

end % function