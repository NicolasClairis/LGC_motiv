function[run_bis] = run_conversion(run_idx)
% [run_bis] = run_conversion(run_idx)
% run_conversion will convert physical (or mental) run index (run_idx) from
% 1-4 range to 1-2 range.
%
% INPUTS
% run_idx: run index
%
% OUTPUTS
% run_bis: run index

nRuns = length(run_idx);
run_bis = NaN(1,nRuns);
for iRun = 1:nRuns
    switch run_idx(iRun)
        case {1,2}
            run_bis(iRun) = 1;
        case {3,4}
            run_bis(iRun) = 2;
    end
end % run loop

end % function