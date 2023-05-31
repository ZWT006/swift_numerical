close all
clear
clc

px0 = -0.557813 ;
py0 = 1.06781   ;
pz0 = 0.0812813 ;
pxf = -2.30443  ;
pyf = 3.01808   ;
pzf = 0         ;
vx0 = -1.075    ;
vy0 = 1.475     ;
vz0 = 0.1275    ;
vxf = 0         ;
vyf = 0         ;
vzf = 0         ;
deltapx = pxf-vx0*T-px0;
deltapy = pyf-vy0*T-py0;
deltapz = pzf-vz0*T-pz0;
deltavx = vxf-vx0;
deltavy = vyf-vy0;
deltavz = vzf-vz0;
% A = [[-12/T^3,0,0,6/T^2,0,0];...
%     [0,-12/T^3,0,0,6/T^2,0];...
%     [0,0,-12/T^3,0,0,6/T^2];...
%     [6/T^2,0,0,-2/T,0,0];...
%     [0,6/T^2,0,0,-2/T,0];...
%     [0,0,6/T^2,0,0,-2/T]]
alpha1=(12*(px0 - pxf + T*vx0))/T^3 - (6*(vx0 - vxf))/T^2;
alpha2=(12*(py0 - pyf + T*vy0))/T^3 - (6*(vy0 - vyf))/T^2;
alpha3=(12*(pz0 - pzf + T*vz0))/T^3 - (6*(vz0 - vzf))/T^2;
beta1 =(2*(vx0 - vxf))/T - (6*(px0 - pxf + T*vx0))/T^2;
beta2 =(2*(vy0 - vyf))/T - (6*(py0 - pyf + T*vy0))/T^2;
beta3 =(2*(vz0 - vzf))/T - (6*(pz0 - pzf + T*vz0))/T^2;

%%%%%%%%%%%%%
%T^4
% + (- 4*vx0^2 - 4*vx0*vxf - 4*vxf^2 - 4*vy0^2 - 4*vy0*vyf - 4*vyf^2 - 4*vz0^2 - 4*vz0*vzf - 4*vzf^2)*T^2
% + (24*pxf*vx0 - 24*px0*vxf - 24*px0*vx0 + 24*pxf*vxf - 24*py0*vy0 - 24*py0*vyf + 24*pyf*vy0 + 24*pyf*vyf - 24*pz0*vz0 - 24*pz0*vzf + 24*pzf*vz0 + 24*pzf*vzf)*T
% - 36*px0^2 + 72*px0*pxf - 36*pxf^2 - 36*py0^2 + 72*py0*pyf - 36*pyf^2 - 36*pz0^2 + 72*pz0*pzf - 36*pzf^2
%%%%%%%%%%%%%

poly=[1,0,...
    (- 4*vx0^2 - 4*vx0*vxf - 4*vxf^2 - 4*vy0^2 - 4*vy0*vyf - 4*vyf^2 - 4*vz0^2 - 4*vz0*vzf - 4*vzf^2),...
    (24*pxf*vx0 - 24*px0*vxf - 24*px0*vx0 + 24*pxf*vxf - 24*py0*vy0 - 24*py0*vyf + 24*pyf*vy0 + 24*pyf*vyf - 24*pz0*vz0 - 24*pz0*vzf + 24*pzf*vz0 + 24*pzf*vzf),...
    - 36*px0^2 + 72*px0*pxf - 36*pxf^2 - 36*py0^2 + 72*py0*pyf - 36*pyf^2 - 36*pz0^2 + 72*pz0*pzf - 36*pzf^2];
roots(poly)
%%T^4 + %%
t=0:0.02:5;
costJ=zeros(length(t),1);
for idx=1:length(t)
    T = t(idx);
    costJ(idx)=((4*vx0^2 + 4*vx0*vxf + 4*vxf^2 + 4*vy0^2 + 4*vy0*vyf + 4*vyf^2 + 4*vz0^2 + 4*vz0*vzf + 4*vzf^2)*T^2 ...
        + (12*px0*vx0 + 12*px0*vxf - 12*pxf*vx0 - 12*pxf*vxf + 12*py0*vy0 + 12*py0*vyf - 12*pyf*vy0 - 12*pyf*vyf +...
        12*pz0*vz0 + 12*pz0*vzf - 12*pzf*vz0 - 12*pzf*vzf)*T + 12*px0^2 - 24*px0*pxf + 12*pxf^2 + 12*py0^2 ...
        - 24*py0*pyf + 12*pyf^2 + 12*pz0^2 - 24*pz0*pzf + 12*pzf^2)/T^3;
end
figure
hold on
plot(t,costJ);
grid on
hold off
timet=1;
idxx = timet/0.02;
fprintf("T = %2.2f; cost = %6.2f",t(idxx),costJ(idxx))

%% test p,v free final
% close all
save('path.mat','path')
clear
clc
po=3;pf=7;
vo=2;
t=0:0.02:10;
% T=(pf-po)/vo;
T=3*(pf-po)/vo;
% T = 6;

deltap=po-pf+vo*T;
alpha = 3*deltap/T^3;
J=3*deltap^2/T^3;

%(-3*vo^2)*T^2 + (12*pf*vo - 12*po*vo)*T - 9*pf^2 + 18*pf*po - 9*po^2
poly = [(-3*vo^2),(12*pf*vo - 12*po*vo), - 9*pf^2 + 18*pf*po - 9*po^2];
solt=roots(poly);

c3=alpha/6;
c2=-alpha/2*T;
c1=vo;
c0=po;

pst=c3*t.^3+c2*t.^2+c1*t+c0;
vst=3*c3*t.^2+2*c2*t.^1+c1;
ast=6*c3*t+2*c2;

figure
hold on
plot(t,pst,'g-');
plot(t,vst,'b-');
plot(t,ast,'r-');
plot([T,T],[0,pf],'k-');
legend('p','v','a','goal');
grid on
%% waypoint 计算 obvp MiniAccInputAcc
% save('path.mat','path')
close all
clear
clc
load('path.mat');

[path_length,~]=size(path);

% 行序列翻转
path=flip(path);

statex=zeros(path_length,4);
statey=zeros(path_length,4);
stateq=zeros(path_length,4);
Tarray=zeros(path_length,1);

statev = zeros(path_length,3);
% 初始速度和加速度都是零
statev(1,:) = [0,0,0];
statea(1,:) = [0,0,0];

Vel_factor = 2;
W_factor   = 1;
Racc       = 1;
dt         = 0.01;
RATIO       = 100;

angle = @(x1,y1,x2,y2) atan2((y2-y1),(x2-x1));

angleset=zeros(path_length,1);

for idx = 2:path_length
%     path(idx+1,3)=(angle(path(idx,1),path(idx,2),path(idx+1,1),path(idx+1,2))+angle(path(idx+1,1),path(idx+1,2),path(idx+2,1),path(idx+2,2)))/2;
    angleset(idx)=angle(path(idx-1,1),path(idx-1,2),path(idx,1),path(idx,2));
end
for idx = 1:path_length-2
    path(idx+1,3)=(angle(path(idx,1),path(idx,2),path(idx+1,1),path(idx+1,2))+angle(path(idx+1,1),path(idx+1,2),path(idx+2,1),path(idx+2,2)))/2;
end

statex(:,1) = path(:,1);
statey(:,1) = path(:,2);
stateq(:,1) = path(:,3);

fp = figure(1);
grid on
axis equal
hold on
title('xy position')
fv = figure(2);
grid on
hold on
title('xy velocity')
fa = figure(3);
grid on
hold on
title('xy acceleration')
fqp = figure(4);
grid on
hold on
title('theta position')
fqv = figure(5);
grid on
hold on
title('theta velocity')
fqa = figure(6);
grid on
hold on
title('theta acceleration')

figure(fp);
% hold on
% grid on
% axis equal
[path_length,~] = size(path);
for idx=1:path_length-1
    plot([path(idx,1),path(idx+1,1)],[path(idx,2),path(idx+1,2)],'c--');%,'LineWidth',2
end
detl = 20;
for idx=1:path_length
    plot([path(idx,1),path(idx,1)+cos(path(idx,3))*detl],[path(idx,2),path(idx,2)+sin(path(idx,3))*detl],'-','Color','b','LineWidth',1);
%     plot([path(idx,1),path(idx,1)+cos(angleset(idx))*detl],[path(idx,2),path(idx,2)+sin(angleset(idx))*detl],'-','Color','b','LineWidth',1);
end


%%%%%%%%第一次启动速度不够,所以时间加倍
obvp=SampleOBVP(Vel_factor/2,W_factor/2,Racc,dt);
idx=1;
spi = path(idx,:)/RATIO;
spf = path(idx+1,:)/RATIO;
spi(3)=spi(3)*RATIO;
spf(3)=spf(3)*RATIO;

[obvp,xt,yt,qt,svf] = obvp.SolveMiniAccInputAcc(spi,statev(idx,:),spf);
statev(idx+1,:) = svf;
Tarray(idx) = obvp.refer_T;
figure(fp);
plot(xt(:,1)*RATIO,yt(:,1)*RATIO,'-','Color','r');
figure(fv);
tf = sum(Tarray(1:idx));
ti = tf-Tarray(idx);
t = ti:dt:tf+dt;
plot(t,xt(:,2),'-','Color','b');
plot(t,yt(:,2),'-','Color','r');
figure(fa);
plot(t,xt(:,3),'-','Color','b');
plot(t,yt(:,3),'-','Color','r');
figure(fqp);
plot(t,rad2deg(qt(:,1)),'-','Color','r');
figure(fqv)
plot(t,qt(:,2),'-','Color','b');
figure(fqa);
plot(t,qt(:,3),'-','Color','b');

obvp=SampleOBVP(Vel_factor,W_factor,Racc,dt);
for idx=2:path_length-1
    spi = path(idx,:)/RATIO;
    spf = path(idx+1,:)/RATIO;
    spi(3)=spi(3)*RATIO;
    spf(3)=spf(3)*RATIO;
    [obvp,xt,yt,qt,svf] = obvp.SolveMiniAccInputAcc(spi,statev(idx,:),spf);
    statev(idx+1,:) = svf;
    Tarray(idx) = obvp.refer_T;
    figure(fp);
    plot(xt(:,1)*RATIO,yt(:,1)*RATIO,'-','Color','r');
    figure(fv);
    tf = sum(Tarray(1:idx));
    ti = tf-Tarray(idx);
    t = ti:dt:tf+dt;
    plot(t,xt(:,2),'-','Color','b');
    plot(t,yt(:,2),'-','Color','r');
    figure(fa);
    plot(t,xt(:,3),'-','Color','b');
    plot(t,yt(:,3),'-','Color','r');
    figure(fqp);
    plot(t,rad2deg(qt(:,1)),'-','Color','r');
    figure(fqv)
    plot(t,qt(:,2),'-','Color','b');
    figure(fqa);
    plot(t,qt(:,3),'-','Color','b');
end



hold off
fprintf("OVER !!! \n");

% ang=angle(50,50,97,58.8);
% fprintf("rad = %2.4f; deg = %2.4f ! \n",ang,rad2deg(ang));

% ang=angle(0,0,1,sqrt(3));
% fprintf("rad = %2.4f; deg = %2.4f ! \n",ang,rad2deg(ang));
% ang=angle(0,0,-1,sqrt(3));
% fprintf("rad = %2.4f; deg = %2.4f ! \n",ang,rad2deg(ang));
% ang=angle(0,0,-1,-sqrt(3));
% fprintf("rad = %2.4f; deg = %2.4f ! \n",ang,rad2deg(ang));
% ang=angle(0,0,1,-sqrt(3));
% fprintf("rad = %2.4f; deg = %2.4f ! \n",ang,rad2deg(ang));

%% test MiniAccInputJerk
close all
clear
clc
load('path.mat');

[path_length,~]=size(path);

path=flip(path);

statex=zeros(path_length,4);
statey=zeros(path_length,4);
stateq=zeros(path_length,4);
Tarray=zeros(path_length,1);

statev = zeros(path_length,3);
% 初始速度和加速度都是零
statev(1,:) = [0,0,0];
statea(1,:) = [0,0,0];

Vel_factor = 2;
W_factor   = 1;
Racc       = 1;
dt         = 0.01;
RATIO       = 100;

angle = @(x1,y1,x2,y2) atan2((y2-y1),(x2-x1));

angleset=zeros(path_length,1);

for idx = 2:path_length
%     path(idx+1,3)=(angle(path(idx,1),path(idx,2),path(idx+1,1),path(idx+1,2))+angle(path(idx+1,1),path(idx+1,2),path(idx+2,1),path(idx+2,2)))/2;
    angleset(idx)=angle(path(idx-1,1),path(idx-1,2),path(idx,1),path(idx,2));
end
for idx = 1:path_length-2
    path(idx+1,3)=(angle(path(idx,1),path(idx,2),path(idx+1,1),path(idx+1,2))+angle(path(idx+1,1),path(idx+1,2),path(idx+2,1),path(idx+2,2)))/2;
end

statex(:,1) = path(:,1);
statey(:,1) = path(:,2);
stateq(:,1) = path(:,3);

fp = figure(1);
grid on
axis equal
hold on
title('xy position')
fv = figure(2);
grid on
hold on
title('xy velocity')
fa = figure(3);
grid on
hold on
title('xy acceleration')
fqp = figure(4);
grid on
hold on
title('theta position')
fqv = figure(5);
grid on
hold on
title('theta velocity')
fqa = figure(6);
grid on
hold on
title('theta acceleration')
fsj = figure(7);
grid on
hold on
title('xy jerk')

figure(fp);
% hold on
% grid on
% axis equal
[path_length,~] = size(path);
for idx=1:path_length
    plot(path(:,1),path(:,2),'bo');%,'LineWidth',2
end
for idx=1:path_length-1
    plot([path(idx,1),path(idx+1,1)],[path(idx,2),path(idx+1,2)],'c--');%,'LineWidth',2
end
detl = 20;
for idx=1:path_length
    plot([path(idx,1),path(idx,1)+cos(path(idx,3))*detl],[path(idx,2),path(idx,2)+sin(path(idx,3))*detl],'-','Color','b','LineWidth',1);
%     plot([path(idx,1),path(idx,1)+cos(angleset(idx))*detl],[path(idx,2),path(idx,2)+sin(angleset(idx))*detl],'-','Color','b','LineWidth',1);
end


%%%%%%%%第一次启动速度不够,所以时间加倍
obvp=SampleOBVP(Vel_factor/2,W_factor/2,Racc,dt);
idx=1;
spi = path(idx,:)/RATIO;
spf = path(idx+1,:)/RATIO;
spi(3)=spi(3)*RATIO;
spf(3)=spf(3)*RATIO;

[obvp,xt,yt,qt,svf,saf] = obvp.SolveMiniJerkInputJerk(spi,statev(idx,:),statea(idx,:),spf);
statev(idx+1,:) = svf;
statea(idx+1,:) = saf;
Tarray(idx) = obvp.refer_T;
figure(fp);
plot(xt(:,1)*RATIO,yt(:,1)*RATIO,'-','Color','r');
figure(fv);
tf = sum(Tarray(1:idx));
ti = tf-Tarray(idx);
t = ti:dt:tf+dt;
plot(t,xt(:,2),'-','Color','b');
plot(t,yt(:,2),'-','Color','r');
figure(fa);
plot(t,xt(:,3),'-','Color','b');
plot(t,yt(:,3),'-','Color','r');
figure(fqp);
plot(t,rad2deg(qt(:,1)),'-','Color','r');
figure(fqv)
plot(t,qt(:,2),'-','Color','b');
figure(fqa);
plot(t,qt(:,3),'-','Color','b');

figure(fsj)
plot(t,xt(:,4),'-','Color','b');
plot(t,yt(:,4),'-','Color','r');

figure(8)
plot(0:0.02:5,polyval(obvp.xc,0:0.02:5),'-','Color','b');
grid on

obvp=SampleOBVP(Vel_factor,W_factor,Racc,dt);
for idx=2:path_length-1
    spi = path(idx,:)/RATIO;
    spf = path(idx+1,:)/RATIO;
    spi(3)=spi(3)*RATIO;
    spf(3)=spf(3)*RATIO;
    [obvp,xt,yt,qt,svf,saf] = obvp.SolveMiniJerkInputJerk(spi,statev(idx,:),statea(idx,:),spf);
    statev(idx+1,:) = svf;
    statea(idx+1,:) = saf;
    Tarray(idx) = obvp.refer_T;
    figure(fp);
    plot(xt(:,1)*RATIO,yt(:,1)*RATIO,'-','Color','r');
    figure(fv);
    tf = sum(Tarray(1:idx));
    ti = tf-Tarray(idx);
    t = ti:dt:tf+dt;
    plot(t,xt(:,2),'-','Color','b');
    plot(t,yt(:,2),'-','Color','r');
    figure(fa);
    plot(t,xt(:,3),'-','Color','b');
    plot(t,yt(:,3),'-','Color','r');
    figure(fqp);
    plot(t,rad2deg(qt(:,1)),'-','Color','r');
    figure(fqv)
    plot(t,qt(:,2),'-','Color','b');
    figure(fqa);
    plot(t,qt(:,3),'-','Color','b');
    figure(fsj)
    plot(t,xt(:,4),'-','Color','b');
    plot(t,yt(:,4),'-','Color','r');
end



hold off
fprintf("OVER !!! \n");