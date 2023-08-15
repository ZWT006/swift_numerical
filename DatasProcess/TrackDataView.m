%***************************************
%Author: Wentao Zhang
%Date: 2023-8
%E-mail: zwt190315@163.com
%Description: track datas plot picture
%Problem: 
%***************************************
%%%%%%%%%%%%% xy popsition
% load("F:\MATLABWorkSpace\MotionPlan\kinodynamicpath\Alldata.mat")
close all
clear

addpath("F:\MATLABWorkSpace\MotionPlan\kinodynamicpath\PolyOpt\")   % for sdfMap
addpath("F:\MATLABWorkSpace\MotionPlan\kinodynamicpath\LazyPRM\")   % for AngleDelta()

% RealDataAddress = "F:\HUST_MS\1-Aggressive Navigation\test\realCSV\01_18_OVAL_04_ok.csv";
% OfflDataAddress = "E:\datas\Swift\Offline\TRAJ_DATA_map18_OVAL.csv";
% SimulDataAddress = "F:\HUST_MS\1-Aggressive Navigation\test\simulation_traj\18_OVAL_10_ok.csv";
% RealDataAddress = "F:\HUST_MS\1-Aggressive Navigation\test\realCSV\05_18_PRM_04_ok.csv";
% OfflDataAddress = "E:\datas\Swift\Offline\TRAJ_DATA_map18_PRM.csv";
% SimulDataAddress = "F:\HUST_MS\1-Aggressive Navigation\test\simulation_traj\18_PRM_10_ok.csv";
% 
% WaypointAddress = "F:\MATLABWorkSpace\MotionPlan\kinodynamicpath\map\pathnode18.mat";
% MapAddress = "F:\MATLABWorkSpace\MotionPlan\kinodynamicpath\map\map18.png";

RealDataAddress = "F:\HUST_MS\1-Aggressive Navigation\test\realCSV\08_17_OVAL_04_ok.csv";
OfflDataAddress = "E:\datas\Swift\Offline\TRAJ_DATA_map17_OVAL.csv";
SimulDataAddress = "F:\HUST_MS\1-Aggressive Navigation\test\simulation_traj\17_OVAL_10_ok.csv";
% RealDataAddress = "F:\HUST_MS\1-Aggressive Navigation\test\realCSV\07_17_PRM_04_ok.csv";
% OfflDataAddress = "E:\datas\Swift\Offline\TRAJ_DATA_map17_PRM.csv";
% SimulDataAddress = "F:\HUST_MS\1-Aggressive Navigation\test\simulation_traj\17_OVAL_10_ok.csv";


WaypointAddress = "F:\MATLABWorkSpace\MotionPlan\kinodynamicpath\map\pathnode17.mat";
MapAddress = "F:\MATLABWorkSpace\MotionPlan\kinodynamicpath\map\map17.png";

RealTraj = readmatrix(RealDataAddress);
OfflTraj = readmatrix(OfflDataAddress);
SimuTraj = readmatrix(SimulDataAddress);
MapImage = imread(MapAddress);
load(WaypointAddress);

real_dt = 0.1;
offl_dt = 0.01;
simu_dt = 0.1;

OX_n = OfflTraj(:,1);
OY_n = OfflTraj(:,2);
OQ_n = OfflTraj(:,3);
OXv_n   = OfflTraj(:,4);
OYv_n   = OfflTraj(:,5);
OQv_n   = OfflTraj(:,6);
OXa_n   = OfflTraj(:,7);
OYa_n   = OfflTraj(:,8);
OQa_n   = OfflTraj(:,9);

X_bias = OX_n(1);
Y_bias = OY_n(1);

RX_n = RealTraj(:,1) + X_bias;
RY_n = RealTraj(:,2) + Y_bias;
RZ_n = RealTraj(:,3);
RRoll_n     = RealTraj(:,4);
RPitch_n    = RealTraj(:,5);
RYaw_n      = RealTraj(:,6);

path_q = RYaw_n;
path_m =path_q;
for idx=2:length(path_q)
%     path_m(idx)=AngleDelta(path_q(idx-1),path_q(idx));
    path_m(idx)=path_m(idx-1)+AngleDelta(path_q(idx-1),path_q(idx));
end
RYaw_n = path_m;


SX_n = SimuTraj(:,1) + X_bias;
SY_n = SimuTraj(:,2) + Y_bias;
SZ_n = SimuTraj(:,3);
SRoll_n     = SimuTraj(:,4);
SPitch_n    = SimuTraj(:,5);
SYaw_n      = SimuTraj(:,6);

path_q = SYaw_n;
path_m =path_q;
for idx=2:length(path_q)
%     path_m(idx)=AngleDelta(path_q(idx-1),path_q(idx));
    path_m(idx)=path_m(idx-1)+AngleDelta(path_q(idx-1),path_q(idx));
end
SYaw_n = path_m;

RATION = 100;

FoitSizeNum = 28;
GRID_SWITCH = 'on';
BOX_SWITCH = 'on';
%% xy position
fpnum = 1;

fpmap =  figure (fpnum);
clf % 清空图表内容
fpnum = fpnum + 1;

picsdfmap = sdfMap(MapImage);

hold on
picsdfmap.showSDFMap(fpmap);
pt_real = plot(RX_n*RATION , RY_n*RATION , 'b-' , 'LineWidth',1);
pt_simu = plot(SX_n*RATION , SY_n*RATION , 'g-' , 'LineWidth',1);
pt_offl = plot(OX_n*RATION , OY_n*RATION , 'r-' , 'LineWidth',1);

pt_way = scatter(path(:, 1), path(:, 2),'ok');

plot(OX_n(1)*RATION , OY_n(1)*RATION,'d','Color','y','MarkerSize',10,'MarkerEdgeColor','y','MarkerFaceColor','y');
plot(OX_n(end)*RATION , OY_n(end)*RATION,'h','Color','g','MarkerSize',12,'MarkerEdgeColor','g','MarkerFaceColor','g');

legend([pt_real,pt_offl,pt_simu,pt_way],'real world','target','simulation','waypoints','FontSize',FoitSizeNum);
% legend([pt_real,pt_offl],'real world','target','FontSize',FoitSizeNum);

txt = xlabel('$p_x$[m]');
set(txt, 'Interpreter', 'latex');
txt = ylabel('$p_y$[m]');
set(txt, 'Interpreter', 'latex');
grid off
box(BOX_SWITCH)
hold off
axis equal
xlim([0 picsdfmap.rows]);
ylim([0 picsdfmap.cols]);
set(gcf,'Position', [100, 100, 100+picsdfmap.rows, 100+picsdfmap.cols]);
if (picsdfmap.rows == 800)
xticks([0 100 200 300 400 500 600 700 800]);
xticklabels({'0','1','2','3','4','5','6','7','8'});
yticks([0 100 200 300 400 500 600 700 800]);
yticklabels({'0','1','2','3','4','5','6','7','8'});
elseif (picsdfmap.rows == 1200)
xticks([0 100 200 300 400 500 600 700 800 900 1000 1100 1200]);
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11','12'});
yticks([0 100 200 300 400 500 600 700 800]);
yticklabels({'0','1','2','3','4','5','6','7','8'});
end

set(gca,'FontName','Times New Roman','FontSize',FoitSizeNum); %,'LineWidth',1.5

% %%%%%%%%%%%%% angular
%% angular
figure (fpnum)
clf % 清空图表内容
fpnum = fpnum + 1;

% yaw angular position
hold on
Toffl = length(OQ_n) * offl_dt;
Treal = length(RYaw_n) * real_dt;
Tsimu = length(SYaw_n) * simu_dt;

T = max([Toffl,Treal,Tsimu]);

toffl = offl_dt:offl_dt:Toffl;
treal = real_dt:real_dt:Treal;
tsimu = simu_dt:simu_dt:Tsimu;

pt_offl = plot(toffl,rad2deg(OQ_n), 'r-');
pt_real = plot(treal,rad2deg(RYaw_n), 'b-');
pt_simu = plot(tsimu,rad2deg(SYaw_n), 'g-');
plot(toffl(end),rad2deg(OQ_n(end)),'p','Color','k','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k');
plot(treal(end),rad2deg(RYaw_n(end)),'^','Color','k','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');
plot(tsimu(end),rad2deg(SYaw_n(end)),'s','Color','k','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');
hold off
txt = ylabel('$p_{\theta}$[deg]');
set(txt, 'Interpreter', 'latex');
ylim([-200 200]);
xlim([0 ceil(T)+5]);
legend([pt_offl,pt_real,pt_simu],{'target','real world','simulation'}, 'AutoUpdate', 'off');
set(gca,'FontName','Times New Roman','FontSize',FoitSizeNum);
box(BOX_SWITCH)
grid(GRID_SWITCH)

hold off
set(gcf,'Position', [100, 100, 1300, 700]);
%%%<------------------------------------------------------------------->%%%
%% roll pitch acceleration
figure (fpnum)
clf % 清空图表内容
fpnum = fpnum + 1;

subplot(2,1,1) % pitch angular
hold on
pt_offl = line([0,ceil(T)+5],[0,0],'linestyle','-','color','r');
pt_real = plot(treal,rad2deg(wrapToPi(RPitch_n)), 'b-');
pt_simu = plot(tsimu,rad2deg(wrapToPi(SPitch_n)), 'g-');
plot(treal(end),rad2deg(RPitch_n(end)),'^','Color','k','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');
plot(tsimu(end),rad2deg(SPitch_n(end)),'s','Color','k','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');
hold off
txt = ylabel('$pitch$[deg]');
set(txt, 'Interpreter', 'latex');
ylim([-10 10]);
xlim([0 ceil(T)+5]);
legend([pt_offl,pt_real,pt_simu],{'target','real world','simulation'}, 'AutoUpdate', 'off');
set(gca,'FontName','Times New Roman','FontSize',FoitSizeNum);
% set(gca,'xticklabel',[])
box(BOX_SWITCH)
grid(GRID_SWITCH)

subplot(2,1,2) % roll angular
hold on
pt_offl = line([0,ceil(T)+5],[0,0],'linestyle','-','color','r');
pt_real = plot(treal,rad2deg(wrapToPi(RRoll_n)), 'b-');
pt_simu = plot(tsimu,rad2deg(wrapToPi(SRoll_n)), 'g-');
plot(treal(end),rad2deg(RRoll_n(end)),'^','Color','k','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');
plot(tsimu(end),rad2deg(SRoll_n(end)),'s','Color','k','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');
hold off
txt = ylabel('$roll$[deg]');
set(txt, 'Interpreter', 'latex');
ylim([-10 10]);
xlim([0 ceil(T)+5]);
legend([pt_offl,pt_real,pt_simu],{'target','real world','simulation'}, 'AutoUpdate', 'off');
set(gca,'FontName','Times New Roman','FontSize',FoitSizeNum);
box(BOX_SWITCH)
grid(GRID_SWITCH)

hold off
set(gcf,'Position', [100, 100, 1300, 700]);


%%%<------------------------------------------------------------------->%%%
%% xy velocity
% figure (fpnum)
% clf % 清空图表内容
% fpnum = fpnum + 1;
% 
% subx = subplot(2,1,1);
% 
% xvelythold = max([abs(min(OXv_n)),abs(max(OXv_n)),abs(min(SXv)),abs(max(PRMX_ddn)), ...
%          abs(min(Y_ddn)),abs(max(Y_ddn)),abs(min(X_ddn)),abs(max(X_ddn)),pa_max]);
% 
% hold on
% pt_xopt = plot(toffl,X_ddn, 'r-');
% pt_xser = plot(treal,PRMX_ddn, 'r--');
% pt_boun = line([0,ceil(T)+2],[pa_max,pa_max],'linestyle','--','color','k','LineWidth',1);
% line([0,ceil(T)+2],[-pa_max,-pa_max],'linestyle','--','color','k','LineWidth',1);
% plot(toffl(end),X_ddn(end),'p','Color','k','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k');
% plot(treal(end),PRMX_ddn(end),'^','Color','k','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');
% hold off
% 
% txt = ylabel('$a_x$[m/s$^2$]');
% set(txt, 'Interpreter', 'latex');
% ylim([-xvelythold-4 xvelythold+4]);
% legend([pt_xopt,pt_xser,pt_boun],{'optimal','search','bound'},'AutoUpdate', 'off'); %,'FontSize',FoitSizeNum ,'AutoUpdate', 'off'
% % legend('Location','northwestoutside');  %east
% % legend('boxoff');
% set(gca,'FontName','Times New Roman','FontSize',FoitSizeNum); %,'LineWidth',1.5
% % set(gca,'xtick',[],'xticklabel',[])
% set(gca,'xticklabel',[])
% box(BOX_SWITCH)
% grid(GRID_SWITCH)
% 
% suby = subplot(2,1,2);
% hold on
% pt_yopt = plot(toffl,Y_ddn, 'b-');
% pt_yser = plot(treal,PRMY_ddn, 'b--');
% pt_boun = line([0,ceil(T)+2],[pa_max,pa_max],'linestyle','--','color','k','LineWidth',1);
% line([0,ceil(T)+2],[-pa_max,-pa_max],'linestyle','--','color','k','LineWidth',1);
% plot(toffl(end),Y_ddn(end),'p','Color','k','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k');
% plot(treal(end),PRMY_ddn(end),'^','Color','k','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');
% hold off
% txt = ylabel('$a_y$[m/s$^2$]');
% set(txt, 'Interpreter', 'latex');
% xlabel("time [s]",'FontSize',12);
% xlim([0 ceil(T)+2]);
% ylim([-xvelythold-4 xvelythold+4]);
% legend([pt_yopt,pt_yser,pt_boun],{'optimal','search','bound'}, 'AutoUpdate', 'off');
% % legend('Location','northwestoutside');
% % legend('boxoff');
% box(BOX_SWITCH)
% grid(GRID_SWITCH)
% set(gca,'FontName','Times New Roman','FontSize',FoitSizeNum); %,'LineWidth',1.5
% set(gcf,'Position', [100, 100, 1300, 700]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT_DEBUG
