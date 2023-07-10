%%%%% 保存 Trajectory 为 .csv 文件

TRAJ_DATA = [X_n',Y_n',Q_n',X_dn',Y_dn',Q_dn'];
filename = "E:\datas\Swaft\Trajectory\TRAJ_DATA_SEMICIRCLE.csv";
TRAJ_DATA = round(TRAJ_DATA,4);
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