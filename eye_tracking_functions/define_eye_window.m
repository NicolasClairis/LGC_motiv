
% define areas of interest on the screen
x_screen = 1024;
y_screen = 768;
% word size parameters
x = x_screen/2;
y = y_screen/2;


%% c1D parameters
left_area_screen = [0, 0, x_screen/2-1, y_screen]; % [left_x0 top_y0 right_x1 bottom_y1]
right_area_screen = [x_screen/2+1, 0, x_screen, y_screen];
left_choice1D = [100 300 400 500]; % [left_x0 top_y0 right_x1 bottom_y1]
right_choice1D = [600 300 900 500];
figure;xlim([0 x_screen]);ylim([0 y_screen])
rectangle('Position',[left_choice1D(1),left_choice1D(2),left_choice1D(3)-left_choice1D(1),left_choice1D(4)-left_choice1D(2)],'EdgeColor','r')
rectangle('Position',[right_choice1D(1),right_choice1D(2),right_choice1D(3)-right_choice1D(1),right_choice1D(4)-right_choice1D(2)],'EdgeColor','b')

%% cRE parameters
ha = 45;
houinon = 67;
woui = 139;
wnon = 166;
E_cRE = [(x - 200 - woui/2), (y + ha/2 - 90), (x + 200 + wnon/3), (y - ha/2 + 30)]; % [left_x0 top_y0 right_x1 bottom_y1]
R_cRE = [(x - 200 - woui/2), (y + ha/2 + 30), (x + 200 + wnon/3), (y - ha/2 + 150)];
YES_cRE = [(x - 200 - woui/2), (y + 200), (x - 200 + woui/3), (y + 200 + houinon)];
NO_cRE = [(x + 200 - wnon/2), (y + 200), (x + 200 + wnon/3), (y + 200 + houinon)];

figure;xlim([0 x_screen]);ylim([0 y_screen])
rectangle('Position',[R_cRE(1),R_cRE(2),R_cRE(3)-R_cRE(1),R_cRE(4)-R_cRE(2)],'EdgeColor','r')
rectangle('Position',[E_cRE(1),E_cRE(2),E_cRE(3)-E_cRE(1),E_cRE(4)-E_cRE(2)],'EdgeColor','b')
rectangle('Position',[YES_cRE(1),YES_cRE(2),YES_cRE(3)-YES_cRE(1),YES_cRE(4)-YES_cRE(2)],'EdgeColor','b')
rectangle('Position',[NO_cRE(1),NO_cRE(2),NO_cRE(3)-NO_cRE(1),NO_cRE(4)-NO_cRE(2)],'EdgeColor','b')