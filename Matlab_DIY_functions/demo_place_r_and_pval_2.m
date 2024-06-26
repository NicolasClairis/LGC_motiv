% various demonstrations of how place_r_and_pval_2 works

% case where two correlations go in opposite ways
figure;
x = 1:10;
y1 = 1:10;
y2 = 10:-1:1;
lWidth = 2;
plot(x, y1,'b','LineWidth', lWidth);
hold on;
plot(x, y2,'r','LineWidth', lWidth);

% correlation 1
r_corr1 = 0.25;
pval1 = 0.025;
col1 = 'b';
% correlation 2
r_corr2 = -0.45;
pval2 = 0.000001;
col2 = 'r';

place_r_and_pval_2(r_corr1, pval1, col1,...
    r_corr2, pval2, col2);


% case where two correlations go in same way
figure;
x = 1:10;
y1 = 1:10;
y2 = 1:2:20;
lWidth = 2;
plot(x, y1,'b','LineWidth', lWidth);
hold on;
plot(x, y2,'r','LineWidth', lWidth);

% correlation 1
r_corr1 = 0.25;
pval1 = 0.025;
col1 = 'b';
% correlation 2
r_corr2 = 0.45;
pval2 = 0.000001;
col2 = 'r';

place_r_and_pval_2(r_corr1, pval1, col1,...
    r_corr2, pval2, col2);