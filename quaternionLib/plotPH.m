
%% ------------------------------------------------------------------------ preamble

clear; clc; format compact; clf; close all;

%% ------------------------------------------------------------------------ visual settings

% 0 - demo setting (video, no frame etc.)
% 1 - developer setting
% 2 - multiple ph setting
setting = 0;
switch setting
    case 0
        showBox     = 0;
        showMovie   = 1;
        showGrid    = 0;
        hideFigAxis = 1;
        
        showCP      = 1;
        showCurves  = 1;
        showFSFrame = 1;
        showTI      = 1;
        show_DMs    = 0;
    case 1
        showBox     = 1;
        showMovie   = 0;
        showGrid    = 1;
        hideFigAxis = 0;
        
        showCP      = 1;
        showCurves  = 1;
        showFSFrame = 1;
        showTI      = 1;
        show_DMs    = 1;
    case 2
        showBox     = 1;
        showMovie   = 0;
        showGrid    = 1;
        hideFigAxis = 0;
        
        showCP      = 0;
        showCurves  = 1;
        showFSFrame = 0;
        showTI      = 0;
        show_DMs    = 0;
end

CColor = [0 160 177]/255;
Tcolor = [255 46  18 ]/255;
Pcolor = [58  149 72 ]/255;
Bcolor = [107 165 231]/255;

%% ------------------------------------------------------------------------ parameter settings

% phFileName = 'ph1';               % !! => setting = {0 | 1}
% phFileName = 'pGHelical';         % !! => setting = {0 | 1}
phFileName = 'phHermite';         % !! => setting = {0 | 1}
% phFileName = 'phHermite_thtVar';  % !! => setting = {  2  }

%% ------------------------------------------------------------------------ reading XYZCoos
disp('>> reading XYZCoos');

cpCoos  = coosReader.Coo3DReader([phFileName '_sph5_cp.txt']);
xyzCoos = coosReader.Coo3DReader([phFileName '_sph5.txt']);

if showFSFrame
    APPCoos = coosReader.Coo3DReader([phFileName '_FSF_PHAPPCoo.txt']);
    TCoos   = coosReader.Coo3DReader([phFileName '_FSF_TCoo.txt']);
    PCoos   = coosReader.Coo3DReader([phFileName '_FSF_PCoo.txt']);
    BCoos   = coosReader.Coo3DReader([phFileName '_FSF_BCoo.txt']);
end

%% ------------------------------------------------------------------------ plot spatial PH quitic
disp('>> plot spatial PH quitic');

mainFig = figure(1);
title('Spatial Quintic PH Curve');
set(mainFig,'units','normalized');
if(setting == 0) set(mainFig,'outerposition',[0 0 1 1]); end
view(-20,33);

hold on;

if(showBox)     box on;                   end
if(showGrid)    grid on; grid minor;      end
if(hideFigAxis) set(gca,'visible','off'); end   % gca - get current axis
set(gcf,'color','w');                           % gcf - get current figure

axis equal;
movegui(figure(1),'northwest');
view(-20,33);

coosReader.plotRef(1, 1);

for i = 1:size(cpCoos,2)-1
    if showCP
        plot3(cpCoos{i}(:,1), cpCoos{i}(:,2), cpCoos{i}(:,3), '-o', 'Color', CColor);
    end
    if showCurves
        plot3(xyzCoos{i}(:,1),xyzCoos{i}(:,2),xyzCoos{i}(:,3),'-k','LineWidth',2);
    end
end

if showFSFrame
    %% -------------------------------------------------------------------- plot Frenet-Serret Frame
    disp('>> plot Frenet-Serret Frame');
    
    for i = 1:size(APPCoos{1},1)
        % manually
        % plot3([APPCoos{1}(i,1) TCoos{1}(i,1)],[APPCoos{1}(i,2) TCoos{1}(i,2)],[APPCoos{1}(i,3) TCoos{1}(i,3)],'-', 'Color', Tcolor);
        % plot3([APPCoos{1}(i,1) PCoos{1}(i,1)],[APPCoos{1}(i,2) PCoos{1}(i,2)],[APPCoos{1}(i,3) PCoos{1}(i,3)],'-', 'Color', Pcolor);
        % plot3([APPCoos{1}(i,1) BCoos{1}(i,1)],[APPCoos{1}(i,2) BCoos{1}(i,2)],[APPCoos{1}(i,3) BCoos{1}(i,3)],'-', 'Color', Bcolor);
        
        p1 = APPCoos{1}(i,:);
        
        p2 = TCoos{1}(i,:);
        dp = p2-p1;
        quiver3(p1(1),p1(2),p1(3),dp(1),dp(2),dp(3),0, 'MaxHeadSize',0.5,'Color', Tcolor);
        
        p2 = PCoos{1}(i,:);
        dp = p2-p1;
        quiver3(p1(1),p1(2),p1(3),dp(1),dp(2),dp(3),0, 'MaxHeadSize',0.5,'Color', Pcolor);
        
        p2 = BCoos{1}(i,:);
        dp = p2-p1;
        quiver3(p1(1),p1(2),p1(3),dp(1),dp(2),dp(3),0, 'MaxHeadSize',0.5,'Color', Bcolor);
        
        if(showMovie)
            view(-20+(i/2),33);
            pause(0.1);
        end
    end
    
end

hold off;

if(showTI)
    %% ------------------------------------------------------------------- plot tangent indicatrix
    disp('>> plot tangent indicatrix');
    
    TIcpCoos    = coosReader.Coo3DReader([phFileName '_sph5_tiCp.txt']);
    TIxyzCoos   = coosReader.Coo3DReader([phFileName '_sph5_ti.txt']);
    TIAxisCoos  = coosReader.Coo3DReader([phFileName '_sph5_ax.txt']);
    
    figure(2);
    title('Tangent Indicatrix');
    movegui(figure(2),'northeast');
    hold on
    set(gcf,'color','w');
    grid on; grid minor;
    axis equal;
    view(-145,30);
    
    % ---- plot xyz reference frame
    coosReader.plotRef(0.3, 2);
    
    % ---- plot unit sphere
    r = 1;
    [x,y,z] = sphere(50);
    x0 = 0; y0 = 0; z0 = 0;
    x = x*r + x0;
    y = y*r + y0;
    z = z*r + z0;
    lightGrey = 0.8*[1 1 1];
    surface(x,y,z, 'FaceColor', 'none', 'EdgeColor', lightGrey)
    
    % ---- plot axis of helical ph quintic
    plot3([0 TIAxisCoos{1}(1)], [0 TIAxisCoos{1}(2)], [0 TIAxisCoos{1}(3)]);
    
    % ---- plot tantent indicatrix on sphere
    plot3(TIxyzCoos{1}(:,1),TIxyzCoos{1}(:,2),TIxyzCoos{1}(:,3), 'Color', Bcolor);
    
    % ---- plot origin-tangent indicatrix vector
    plot3([0 TIxyzCoos{1}(1,1)], [0 TIxyzCoos{1}(1,2)], [0 TIxyzCoos{1}(1,3)],'-', 'Color', Pcolor);
    for i=20:20:size(TIxyzCoos{1},1)
        
        plot3([0 TIxyzCoos{1}(i,1)], [0 TIxyzCoos{1}(i,2)], [0 TIxyzCoos{1}(i,3)],'-k');
        
        if(showMovie)
            pause(0.1);
            view(-145-(i/(20*2)),30);
        end
    end
    plot3([0 TIxyzCoos{1}(end,1)], [0 TIxyzCoos{1}(end,2)], [0 TIxyzCoos{1}(end,3)],'-', 'Color', Tcolor);
    
    hold off;
    
end

if(showTI)
    %% ------------------------------------------------------------------- plot differential measures
    disp('>> plot differential measures');
    
    Kpp_Coos = coosReader.Coo3DReader([phFileName '_DM_Kpp.txt']);
    Tau_Coos = coosReader.Coo3DReader([phFileName '_DM_Tau.txt']);
    Omg_Coos = coosReader.Coo3DReader([phFileName '_DM_Omg.txt']);
    
    figure(100);
    title('Curvature Kappa(S(t)/L)');
    xlabel('$\frac{S(t)}{L}$','Interpreter','LaTex')
    ylabel('$\kappa(\frac{S(t)}{L})$','Interpreter','LaTex')
    movegui(figure(100),'southwest');
    hold on
    set(gcf,'color','w');
    grid on; grid minor;
    plot(Kpp_Coos{1}(:,1), Kpp_Coos{1}(:,2),'-','Color',Bcolor);
    hold off;
    
    figure(101);
    title('Torsion Tau(S(t)/L)');
    xlabel('$\frac{S(t)}{L}$','Interpreter','LaTex')
    ylabel('$\tau(\frac{S(t)}{L})$','Interpreter','LaTex');
    movegui(figure(101),'south');
    hold on
    set(gcf,'color','w');
    grid on; grid minor;
    plot(Tau_Coos{1}(:,1), Tau_Coos{1}(:,2),'-','Color',Bcolor);
    hold off;
    
    figure(102);
    title('Rate of Rotation Omega(S(t)/L)');
    xlabel('$\frac{S(t)}{L}$','Interpreter','LaTex')
    ylabel('$\omega(\frac{S(t)}{L})$','Interpreter','LaTex')
    movegui(figure(102),'southeast');
    hold on
    set(gcf,'color','w');
    grid on; grid minor;
    plot(Omg_Coos{1}(:,1), Omg_Coos{1}(:,2),'-','Color',Bcolor);
    hold off;

end

disp('>> done!');





