%%%%% 保存 Trajectory 为 .csv 文件


TRAJ_DATA = [PRMX_n,PRMY_n,PRMQ_n,PRMX_dn,PRMY_dn,PRMQ_dn,PRMX_ddn,PRMY_ddn,PRMQ_ddn];
filename = "E:\datas\Swift\Trajectory\TRAJ_DATA_map18_PRM.csv";
TRAJ_DATA = round(TRAJ_DATA,6);
writematrix(TRAJ_DATA,filename);

TRAJ_DATA = [X_n',Y_n',Q_n',X_dn',Y_dn',Q_dn',X_ddn',Y_ddn',Q_ddn'];
filename = "E:\datas\Swift\Trajectory\TRAJ_DATA_map18_OVAL.csv";
TRAJ_DATA(:,1) = TRAJ_DATA(:,1)/100;
TRAJ_DATA(:,2) = TRAJ_DATA(:,2)/100;
TRAJ_DATA = round(TRAJ_DATA,6);
writematrix(TRAJ_DATA,filename);

TRAJ_DATA = [QPX_n',QPY_n',QPQ_n',QPX_dn',QPY_dn',QPQ_dn',QPX_ddn',QPY_ddn',QPQ_ddn'];
filename = "E:\datas\Swift\Trajectory\TRAJ_DATA_map18_OSQP.csv";
TRAJ_DATA(:,1) = TRAJ_DATA(:,1)/100;
TRAJ_DATA(:,2) = TRAJ_DATA(:,2)/100;
TRAJ_DATA = round(TRAJ_DATA,6);
writematrix(TRAJ_DATA,filename);

TRAJ_DATA = [X_n',Y_n',Q_n',X_dn',Y_dn',Q_dn',X_ddn',Y_ddn',Q_ddn'];
TRAJ_DATA(:,1) = TRAJ_DATA(:,1)/100;
TRAJ_DATA(:,2) = TRAJ_DATA(:,2)/100;
filename = "E:\datas\Swift\Trajectory\TRAJ_DATA_map18_OPT.csv";
TRAJ_DATA = round(TRAJ_DATA,6);
writematrix(TRAJ_DATA,filename);

% 设定符               字段分隔符
% 
% ','     'comma'     逗号。这是默认行为。
% 
% ' '     'space'     空格
% 
% '\t'    'tab'       制表符
% 
% ';'     'semi'      分号
% 
% '|'     'bar'       垂直条
% 
% eg: 'Delimiter','space'

% TRAJ_DATA_READ = readmatrix("TRAJ_DATA.csv");