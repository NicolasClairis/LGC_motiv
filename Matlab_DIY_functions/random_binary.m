function[binary_rand] = random_binary()
%[binary_rand] = random_binary()
% random_binary will transform rand function output into a binary (0/1)
% variable for when you need one.
%
% OUTPUTS
% binary_rand: number randomly equal to 0 or 1

%% define first trial task type
% (0) odd/even
% (1) lower/higher than 5

% avoid ambiguity with 0.5 number
n_random = 0.5;
while n_random == 0.5
    n_random = rand;
end

%% binarize rand output
if n_random > 0.5
    binary_rand = 0;
elseif n_random < 0.5
    binary_rand = 1;
end

end % function