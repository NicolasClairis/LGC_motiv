function[T0] = TTL_wait(dummy_scans, trigger_id)
%[T0] = TTL_wait(dummy_scans, trigger_id)
% function waiting for first TTL before launching the rest of the script
%
% INPUTS
% dummy_scans: number of TTL to wait before starting experiment
%
% trigger_id: trigger for TTL code (should correspond to KbName('t')
%
% OUTPUTS
% T0: time of T0 (ie first TTL sent by fMRI)

next = 0;
% wait dummy_scan number of volumes before starting the task
while next < dummy_scans
    [keyisdown, T0IRM, keycode] = KbCheck;
    
    if keyisdown == 1 && keycode(trigger_id) == 1
        if next == 0
            T0 = T0IRM;
        end
        next = next + 1;
        disp([num2str(next),' TTL received']);
        while keycode(trigger_id) == 1
            [keyisdown, T, keycode] = KbCheck;
        end
    end
end

end % function