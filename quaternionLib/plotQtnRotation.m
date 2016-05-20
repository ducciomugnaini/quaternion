
%% --------------------------------------------------------------- preamble

clear; clc; format compact; clf; close all;

RedColor = [255 46  18 ]/255;
GrnColor = [58  149 72 ]/255;
BluColor = [107 165 231]/255;


%% ------------------------------------------------------------------------ parameter settings

fileName = 'xyzQtnRotation.txt';

%% ------------------------------------------------------------------------ reading XYZCoos
disp('>> reading XYZCoos');

xyzCoos = coosReader.Coo3DReader(fileName);

figure(1);
hold on;
axis equal;
set(gcf,'color','w');
view(14, -7);

grid on; grid minor;
box on;

% ---- plot xyz reference frame
coosReader.plotRef(0.3, 1);

plot3([0 5],[0 5],[0 5],'-o', 'Color', RedColor);   % axis
plot3(1,0,1,'o', 'Color', GrnColor);                % point

for i = 1:size(xyzCoos{1},1)
    plot3(xyzCoos{1}(i,1), xyzCoos{1}(i,2), xyzCoos{1}(i,3),'.', 'Color', BluColor);
    pause(0.1);
end

% plot3(xyzCoos{1}(:,1), xyzCoos{1}(:,2), xyzCoos{1}(:,3),'-', 'Color', BluColor);

hold off;






