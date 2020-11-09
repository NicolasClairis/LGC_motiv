function[PRF_kernel, PRF_kernel_part1, PRF_kernel_part2] = esti_PRF(Y, PLOT_ME,FLIP_GAMMA_1)
% This script evaluates a pupil response function (PRF) from the mean pupil
% response over several subject and trials. The approximation is composed
% of two gamma distribution shaped functions.
%
% INPUT
% Y: pupil diameter across time (locked on the onset of the stimulus)
%
% PLOT_ME: plot fit (1) or not (0)
%
% FLIP_GAMMA_1: expect a constriction (1) or a dilation (0) of the pupil
% after the onset
%
% OUTPUT
% PRF_kernel: signal approximated with the PRF
%
% PRF_kernel_part1: first gamma estimation with the PRF (part mainly
% driven by the changes on screen (luminance, movement))
%
% PRF_kernel_part2: second gamma estimation with the PRF (part still driven
% by luminance but which may also contain more cognitive factors)
%
% Author: Jules BROCHARD
% Last edit: 24 Novemebr 2017

%% Global setting
% PLOT_ME = 1;
% FLIP_GAMMA_1 = 1;

%% clean Y from eventual NaN values
Y = Y(~isnan(Y));
index_Y = find(~isnan(Y));

%% Load the raw data
% Y = getfield(load('tmp.mat'),'Y6'); % Nicolas CLAIRIS grand mean pupil data
% Y = -getfield(load('grand_mean_response_AW_2017-11-24.mat'),'pupil_response'); % Antonius WIELHLER grand mean data
% Y = Y1;

if PLOT_ME
    set(figure,'color','w');
    clf
    subplot(3,4,[6,7,10,11])
    h_raw=plot(index_Y, Y);
end

%% manual approximation & cleaning & shortcut
flat_Y = conv(Y,1/ceil(numel(Y)/20)*ones(ceil(numel(Y)/20),1));
flat_Y = flat_Y(ceil(numel(Y)/20):(end-ceil(numel(Y)/20)));
[~,ii_min_peak] = min(flat_Y);
[~,ii_max_peak] = max(flat_Y(1:ii_min_peak));
rough_connect_point = (ii_max_peak+ii_min_peak)/2;
pre_connect = Y(1:rough_connect_point);
post_connect = Y((rough_connect_point+1):end);

% /!\ This script used the term gamma, not for the usual gamma function but
% for the familly of function that have the shape of gamma distributions
gamma_fun = @(a,b,x) 1/(b^a*gamma(a)) .* (x).^(a-1) .* exp(-x/b);
simple_gamma = @(a,b,sf,cst,dil_x,x) cst + sf.*gamma_fun(a,b,dil_x*x);
if FLIP_GAMMA_1
    do_i_flip = @fliplr; 
    do_i_flip2 = @(x,connect_point) (connect_point-x);
else
    do_i_flip = @(x)x; 
    do_i_flip2 = @(x,connect_point) x;
end

%% first gamma approximation
f01 = fit((1:numel(pre_connect))',do_i_flip(pre_connect)',simple_gamma,...
    'StartPoint',[2,2,2*std(pre_connect)*10,mean(pre_connect),1/100],...
    'Lower',[0,0,0,-max(abs(pre_connect)),0],...
    'Upper',[Inf,Inf,Inf,max(abs(pre_connect)),Inf]);

if PLOT_ME
    subplot(3,4,[1,2])
    h_raw =plot(pre_connect);
    hold on;
    h_appro = plot(do_i_flip(simple_gamma(f01.a,f01.b,f01.sf,f01.cst,f01.dil_x,1:numel(pre_connect))),...
        'Linewidth',2) ;
    legend([h_raw,h_appro],'Raw data','Approximation')
    legend('boxoff')
    title('First approximation of the first gamma')
    xlabel('time')
    ylabel('Pupil response')
    drawnow
end

%% second gamma approximation
f02 = fit((1:numel(post_connect))',(post_connect)',simple_gamma,...
    'StartPoint',[2,2,-2*std(post_connect)*10,mean(post_connect),1/100],...
    'Lower',[0,0,-Inf,-max(abs(post_connect)),0],...
    'Upper',[Inf,Inf,0,Inf]);

if PLOT_ME
    subplot(3,4,[3,4])
    h_raw = plot(post_connect);
    hold on;
    h_appro = plot((simple_gamma(f02.a,f02.b,f02.sf,f02.cst,f02.dil_x,1:numel(post_connect))),...
        'Linewidth',2);
    legend([h_raw,h_appro],'Raw data','Approximation')
    legend('boxoff')
    title('First approximation of the second gamma')
    xlabel('time')
    ylabel('Pupil response')
    drawnow
end

%% PRF approximation
% The PRF has three terms: a constant, a reversed gamma function and a
% negative gamma function. The reversed gamma function is simply a gamma
% function read evaluated from rigth to left, instead of left to right.
% The PRF has 6 parameters: 2 for the first gamma function, 2 for the
% second gamma function, 1 for the constant, 1 for the abscisse (x)
% dilatation and 1 for the connection point between the two gamma function.

PRF_fun = @(a,b,sf,a2,b2,sf2,cst,dil_x,connect_point,x) max(min(...% to avoid infinity estimation
    (cst ...                                                                constant term
    + sf.*gamma_fun(a,b,dil_x*do_i_flip2(x,connect_point)) .* (x<connect_point)...    reverse gamma function up to a given threshold
    + sf2.*gamma_fun(a2,b2,dil_x*(x-connect_point)) .* (x>=connect_point)... negative gamma function from a given threshold
    ),1e10),-1e10);

fo = fit((1:numel(Y))',(Y)',PRF_fun...
    ,'StartPoint',[ % Starting point of the parameter for the optimisaition
    f01.a,...                 a
    f01.b,...                b
    f01.sf,...     sf
    f02.a,...                 a2
    f02.b,...                 b2
    f02.sf,...    sf2
    mean(Y),...         cst
    max(f01.dil_x,f02.dil_x), ...            dil_x
    rough_connect_point ...              connect_point
    ] ...
    ,'Lower',[ % Lower bound of the parameter search space
    0,...                          a
    0,...                          b
    0,...                          sf
    0,...                          a2
    0,...                          b2
    -Inf,...                       sf2
    -Inf,...                       cst
    0,...                          dil_x
    1,...                          connect_point
    ]...
    ,'Upper',[  % Upper bound of the parameter search space
    Inf,...               a
    Inf,...               b
    Inf,...               sf
    Inf,...               a2
    Inf,...               b2
    0,...                 sf2
    Inf, ...              cst
    Inf, ...               dil_x
    numel(Y) ...          connect_point
    ]);


PRF_kernel_part1 = ...
    fliplr(simple_gamma(f01.a,f01.b,f01.sf,f01.cst,f01.dil_x,1:numel(pre_connect)));
PRF_kernel_part2 = ...
    simple_gamma(f02.a,f02.b,f02.sf,f02.cst,f02.dil_x,1:numel(post_connect));
PRF_kernel = PRF_fun(...
    fo.a,fo.b,fo.sf,fo.a2,fo.b2,fo.sf2,fo.cst,fo.dil_x,fo.connect_point,1:numel(Y));

% set(figure,'color','w'); plot(PRF_kernel_part1); hold on,plot(PRF_kernel_part2); plot(PRF_kernel);

if PLOT_ME
    subplot(3,4,[6,7,10,11])
    hold off
    h_raw=plot(index_Y,Y);
    hold on
    h_appro=plot(index_Y,PRF_kernel...
        ,'Linewidth',2);
    legend([h_raw,h_appro],'Raw data','Approximation')
    legend('boxoff','position','best')
    title('Pupil Response Function estimation')
    xlabel('time')
    ylabel('Pupil response')
end

end