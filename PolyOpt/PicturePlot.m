%***************************************
%Author: Wentao Zhang
%Date: 2023-6
%E-mail: zwt190315@163.com
%Description: plot picture
%Problem: 
%***************************************
%%%%%%%%%%%%% xy popsition
load("F:\MATLABWorkSpace\MotionPlan\kinodynamicpath\Alldata.mat")
close all
fpnum = 9;
FoitSizeNum = 16;
GRID_SWITCH = 'on';
BOX_SWITCH = 'on';
%% xy position
fpmap =  figure (fpnum);
clf % 清空图表内容
fpnum = fpnum + 1;

picsdfmap = sdfMap(map);

hold on
picsdfmap.showSDFMap(fpmap);
pt_opt = plot(X_n, Y_n , 'm-' , 'LineWidth',1);
pt_way = scatter(path(1:size(path, 1), 1)*RATION, path(1:size(path, 1), 2)*RATION,'ok');
pt_ser = plot(QPX_n,QPY_n,'m--' , 'LineWidth',1);
plot(QPX_n(1),QPY_n(1),'d','Color','y','MarkerSize',10,'MarkerEdgeColor','y','MarkerFaceColor','y');
plot(QPX_n(end),QPY_n(end),'h','Color','g','MarkerSize',12,'MarkerEdgeColor','g','MarkerFaceColor','g');
legend([pt_opt,pt_ser,pt_way],'optimal','search','waypoints','FontSize',FoitSizeNum);

txt = xlabel('$p_x$[m]');
set(txt, 'Interpreter', 'latex');
txt = ylabel('$p_y$[m]');
set(txt, 'Interpreter', 'latex');
grid off
box(BOX_SWITCH)
hold off
axis equal
xlim([0 sdfmap.rows]);
ylim([0 sdfmap.cols]);
set(gcf,'Position', [100, 100, 100+sdfmap.rows, 100+sdfmap.cols]);
if (sdfmap.rows == 800)
xticks([0 100 200 300 400 500 600 700 800]);
xticklabels({'0','1','2','3','4','5','6','7','8'});
yticks([0 100 200 300 400 500 600 700 800]);
yticklabels({'0','1','2','3','4','5','6','7','8'});
elseif (sdfmap.rows == 1200)
xticks([0 100 200 300 400 500 600 700 800 900 1000 1100 1200]);
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11','12'});
yticks([0 100 200 300 400 500 600 700 800]);
yticklabels({'0','1','2','3','4','5','6','7','8'});
end

set(gca,'FontName','Times New Roman','FontSize',FoitSizeNum); %,'LineWidth',1.5


%%%%%%%%%%%%% xy velocity
%% xy velocity
figure (fpnum)
clf % 清空图表内容
fpnum = fpnum + 1;
topt=0:tstep:(k-1)*tstep;
topt = topt + Tv0;
tser=0:tstep:(QP_k-1)*tstep;
ythold = max([abs(min(QPY_dn)),abs(max(QPY_dn)),abs(min(QPX_dn)),abs(max(QPX_dn)), ...
         abs(min(Y_dn)),abs(max(Y_dn)),abs(min(X_dn)),abs(max(X_dn)),pv_max]);

subx = subplot(2,1,1);
hold on
pt_xopt = plot(topt,X_dn, 'r-');
pt_xser = plot(tser,QPX_dn, 'r--');
pt_boun = line([0,ceil(T)+2],[pv_max,pv_max],'linestyle','--','color','k','LineWidth',1);
line([0,ceil(T)+2],[-pv_max,-pv_max],'linestyle','--','color','k','LineWidth',1);
plot(topt(end),X_dn(end),'p','Color','k','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k');
plot(tser(end),QPX_dn(end),'^','Color','k','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');
hold off

txt = ylabel('$v_x$[m/s]');
set(txt, 'Interpreter', 'latex');
ylim([-ythold-1 ythold+1]);
legend([pt_xopt,pt_xser,pt_boun],{'optimal','search','bound'},'AutoUpdate', 'off'); %,'FontSize',FoitSizeNum ,'AutoUpdate', 'off'
% legend('Location','northwestoutside');  %east
% legend('boxoff');
set(gca,'FontName','Times New Roman','FontSize',FoitSizeNum); %,'LineWidth',1.5
% set(gca,'xtick',[],'xticklabel',[])
set(gca,'xticklabel',[])
box(BOX_SWITCH)
grid(GRID_SWITCH)

suby = subplot(2,1,2);
hold on
pt_yopt = plot(topt,Y_dn, 'b-');
pt_yser = plot(tser,QPY_dn, 'b--');
pt_boun = line([0,ceil(T)+2],[pv_max,pv_max],'linestyle','--','color','k','LineWidth',1);
line([0,ceil(T)+2],[-pv_max,-pv_max],'linestyle','--','color','k','LineWidth',1);
plot(topt(end),Y_dn(end),'p','Color','k','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k');
plot(tser(end),QPY_dn(end),'^','Color','k','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');
hold off
txt = ylabel('$v_y$[m/s]');
set(txt, 'Interpreter', 'latex');
xlabel("time [s]",'FontSize',12);
xlim([0 ceil(T)+2]);
ylim([-ythold-1 ythold+1]);
legend([pt_yopt,pt_yser,pt_boun],{'optimal','search','bound'}, 'AutoUpdate', 'off');
% legend('Location','northwestoutside');
% legend('boxoff');
box(BOX_SWITCH)
grid(GRID_SWITCH)
set(gca,'FontName','Times New Roman','FontSize',FoitSizeNum); %,'LineWidth',1.5
set(gcf,'Position', [100, 100, 1300, 700]);

%%%%%%%%%%%%% angular
%% angular
figure (fpnum)
clf % 清空图表内容
fpnum = fpnum + 1;

subplot(3,1,1) % angular position
hold on
pt_qopt = plot(topt,rad2deg(Q_n), 'm-');
pt_qser = plot(tser,rad2deg(QPQ_n), 'm--');
plot(topt(end),rad2deg(Q_n(end)),'p','Color','k','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k');
plot(tser(end),rad2deg(QPQ_n(end)),'^','Color','k','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');
hold off
txt = ylabel('$p_{\theta}$[deg]');
set(txt, 'Interpreter', 'latex');
ylim([-200 200]);
legend([pt_qopt,pt_qser],{'optimal','search'}, 'AutoUpdate', 'off');
set(gca,'FontName','Times New Roman','FontSize',FoitSizeNum);
set(gca,'xticklabel',[])
box(BOX_SWITCH)
grid(GRID_SWITCH)


subplot(3,1,2) % angular velocity
hold on
pt_qopt = plot(topt,Q_dn, 'm-');
pt_qser = plot(tser,QPQ_dn, 'm--');
pt_boun = line([0,ceil(T)+2],[wv_max,wv_max],'linestyle','--','color','k','LineWidth',1);
line([0,ceil(T)+2],[-wv_max,-wv_max],'linestyle','--','color','k','LineWidth',1);
plot(topt(end),Q_dn(end),'p','Color','k','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k');
plot(tser(end),QPQ_dn(end),'^','Color','k','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');
hold off
txt = ylabel('$\omega_{\theta}$[rad/s]');
set(txt, 'Interpreter', 'latex');
ythold = max([abs(min(QPQ_dn)),abs(max(QPQ_dn)),abs(min(Q_dn)),abs(max(Q_dn)),wv_max]);
ylim([-ythold-1 ythold+1]);
legend([pt_qopt,pt_qser,pt_boun],{'optimal','search','bound'}, 'AutoUpdate', 'off');
set(gca,'FontName','Times New Roman','FontSize',FoitSizeNum);
set(gca,'xticklabel',[])
box(BOX_SWITCH)
grid(GRID_SWITCH)

subplot(3,1,3) % angular accelaration
hold on
pt_qopt = plot(topt,Q_ddn, 'm-');
pt_qser = plot(tser,QPQ_ddn, 'm--');
pt_boun = line([0,ceil(T)+2],[wa_max,wa_max],'linestyle','--','color','k','LineWidth',1);
line([0,ceil(T)+2],[-wa_max,-wa_max],'linestyle','--','color','k','LineWidth',1);
plot(topt(end),Q_ddn(end),'p','Color','k','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k');
plot(tser(end),QPQ_ddn(end),'^','Color','k','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');
hold off
txt = ylabel('$a_{\theta}$[rad/s${}^2$]'); %[rad/s^2]
set(txt, 'Interpreter', 'latex');
xlabel("time [s]");
xlim([0 ceil(T)+2]);
ythold = max([abs(min(QPQ_ddn)),abs(max(QPQ_ddn)),abs(min(Q_ddn)),abs(max(Q_ddn)),wa_max]);
ylim([-ythold-4 ythold+4]);
legend([pt_qopt,pt_qser,pt_boun],{'optimal','search','bound'}, 'AutoUpdate', 'off');
set(gca,'FontName','Times New Roman','FontSize',FoitSizeNum);
set(gca,'xticklabel',[])
box(BOX_SWITCH)
grid(GRID_SWITCH)

hold off
set(gcf,'Position', [100, 100, 1300, 700]);


%%%<------------------------------------------------------------------->%%%
%% xy acceleration
figure (fpnum)
clf % 清空图表内容
fpnum = fpnum + 1;

topt=0:tstep:(k-1)*tstep;
topt = topt + Tv0;
tser=0:tstep:(QP_k-1)*tstep;
ythold = max([abs(min(QPY_ddn)),abs(max(QPY_ddn)),abs(min(QPX_ddn)),abs(max(QPX_ddn)), ...
         abs(min(Y_ddn)),abs(max(Y_ddn)),abs(min(X_ddn)),abs(max(X_ddn)),pa_max]);

subx = subplot(2,1,1);
hold on
pt_xopt = plot(topt,X_ddn, 'r-');
pt_xser = plot(tser,QPX_ddn, 'r--');
pt_boun = line([0,ceil(T)+2],[pa_max,pa_max],'linestyle','--','color','k','LineWidth',1);
line([0,ceil(T)+2],[-pa_max,-pa_max],'linestyle','--','color','k','LineWidth',1);
plot(topt(end),X_ddn(end),'p','Color','k','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k');
plot(tser(end),QPX_ddn(end),'^','Color','k','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');
hold off

txt = ylabel('$a_x$[m/s$^2$]');
set(txt, 'Interpreter', 'latex');
ylim([-ythold-4 ythold+4]);
legend([pt_xopt,pt_xser,pt_boun],{'optimal','search','bound'},'AutoUpdate', 'off'); %,'FontSize',FoitSizeNum ,'AutoUpdate', 'off'
% legend('Location','northwestoutside');  %east
% legend('boxoff');
set(gca,'FontName','Times New Roman','FontSize',FoitSizeNum); %,'LineWidth',1.5
% set(gca,'xtick',[],'xticklabel',[])
set(gca,'xticklabel',[])
box(BOX_SWITCH)
grid(GRID_SWITCH)

suby = subplot(2,1,2);
hold on
pt_yopt = plot(topt,Y_ddn, 'b-');
pt_yser = plot(tser,QPY_ddn, 'b--');
pt_boun = line([0,ceil(T)+2],[pa_max,pa_max],'linestyle','--','color','k','LineWidth',1);
line([0,ceil(T)+2],[-pa_max,-pa_max],'linestyle','--','color','k','LineWidth',1);
plot(topt(end),Y_ddn(end),'p','Color','k','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k');
plot(tser(end),QPY_ddn(end),'^','Color','k','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k');
hold off
txt = ylabel('$a_y$[m/s$^2$]');
set(txt, 'Interpreter', 'latex');
xlabel("time [s]",'FontSize',12);
xlim([0 ceil(T)+2]);
ylim([-ythold-4 ythold+4]);
legend([pt_yopt,pt_yser,pt_boun],{'optimal','search','bound'}, 'AutoUpdate', 'off');
% legend('Location','northwestoutside');
% legend('boxoff');
box(BOX_SWITCH)
grid(GRID_SWITCH)
set(gca,'FontName','Times New Roman','FontSize',FoitSizeNum); %,'LineWidth',1.5
set(gcf,'Position', [100, 100, 1300, 700]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT_DEBUG
