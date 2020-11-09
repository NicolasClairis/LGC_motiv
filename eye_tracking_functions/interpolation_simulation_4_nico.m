% This script explores the impact of missing data in a signal
% Author Jules BROCHARD (brochard.jules@gmail.com)

n_data = 1e3;
auto_corr_kernel = ones(n_data/10,1)/10;

hole_prct_range = [0.5,1:99,99.5];
[R2_per_hole_prct_per_simu,core_per_hole_prct_per_simu] = deal(NaN(numel(hole_prct_range),1e2));
for hole_prct_idx = 1:numel(hole_prct_range)
    hole_prct = hole_prct_range(hole_prct_idx);
    simu_R2 = NaN(1e2,1);
    for ii=1:1e2
        signal = conv(normrnd(0,1,n_data,1),auto_corr_kernel,'same');
        interpol_holed_signal = interpolate_signal(poke_hole(hole_prct,signal));
        R2_per_hole_prct_per_simu(hole_prct_idx,ii) = 1- var(interpol_holed_signal-signal)/var(signal);
        core_per_hole_prct_per_simu(hole_prct_idx,ii) = corr(interpol_holed_signal,signal);
    end
end

%% plotting
set(figure,'color','w')
signal = conv(normrnd(0,1,n_data,1),auto_corr_kernel,'same');

subplot(2,3,1)
plot(signal,'LineWidth',2)
hold on
plot(interpolate_signal(poke_hole(10,signal)),'LineWidth',1)
title('10% of holes')
ylabel('Signal')

subplot(2,3,2)
plot(signal,'LineWidth',2)
hold on
plot(interpolate_signal(poke_hole(50,signal)),'LineWidth',1)
title('50% of holes')
ylabel('Signal')

subplot(2,3,3)
plot(signal,'LineWidth',2)
hold on
plot(interpolate_signal(poke_hole(95,signal)),'LineWidth',1)
title('95% of holes')
ylabel('Signal')

subplot(2,3,4:6)
plot(hole_prct_range,mean(R2_per_hole_prct_per_simu,2),'LineWidth',2)
hold on
plot(hole_prct_range,mean(core_per_hole_prct_per_simu,2),'LineWidth',1)
title('Accuracy of linear extrapolation')
xlabel('Percentage of holes')
legend({'R2','Correlation'},'Location','Best')
legend('boxoff')
ylim([0,1])

%% miscellaneous function
function holed_signal = poke_hole(hole_prct,signal)
n_holes = ceil(numel(signal)*hole_prct/100);
holes_idx = randsample(numel(signal),n_holes);

holed_signal = signal;
holed_signal(holes_idx) = NaN;
end

function interpolated_signal=interpolate_signal(signal_with_nan)
% do linear interpolation over nan
nan_idx = find(isnan(signal_with_nan));

interpolated_signal = signal_with_nan;
if isempty(nan_idx)
    disp('Great ! No Nan value !')
else
    idx = nan_idx(1);
    while idx <= numel(signal_with_nan)
        % nothing to do here
        if ~isnan(signal_with_nan(idx))
            interpolated_signal(idx) = signal_with_nan(idx);
            idx = idx + 1;
        else
            % get initial value to interpolate from it
            if idx==1
                previous_value = 0; idx =2;
            else
                previous_value = interpolated_signal(idx-1);
            end
            % get next non nan value to interpolate the signal to it
            next_non_NaN_idx = idx + find(~isnan(signal_with_nan(idx:end)),1) -1;
            if isempty(next_non_NaN_idx)
                next_non_NaN_idx = numel(signal_with_nan);
                signal_with_nan(next_non_NaN_idx) = 0;
            end
            
            % interpolation
            interpolated_signal((idx-1):next_non_NaN_idx) = ...
                linspace(previous_value,signal_with_nan(next_non_NaN_idx),next_non_NaN_idx-idx+2);
            idx = next_non_NaN_idx+1;
        end
    end
end
end