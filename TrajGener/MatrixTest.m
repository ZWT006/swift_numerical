ts = [0.7757,0.2853,0.6976];
MatQ = getQ(3,7,4,ts);
filename = "E:\datas\Swift\Optimization\MATQ.csv";
TRAJ_DATA = MatQ;
TRAJ_DATA = round(TRAJ_DATA,4);
writematrix(TRAJ_DATA,filename);