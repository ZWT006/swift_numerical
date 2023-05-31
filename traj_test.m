clear
close all
clc
%% 1D test
bvp = OBVP;
sti=[0,1,1];
stf=[1,0,0];
T = 1;
bvp = bvp.FixSolve(sti,stf,T);

[pst,vst,ast,jst]= bvp.ft(0.01);

figure
hold on
plot(pst,'k.');
plot(vst,'b.');
plot(ast,'r.');
plot(jst,'y.');
legend('p','v','a','j');
grid on

%% 1D free final
bvp = OBVP;
sti = [0,3,2];
stf = [1,0,0];

bvp = bvp.FreeFinalFreeTimeSolve(sti,stf);

[pst,vst,ast,jst]= bvp.ft(0.002);

figure
hold on
plot(pst,'k.');
plot(vst,'b.');
plot(ast,'r.');
plot(jst,'y.');
legend('p','v','a','j');
grid on

%% 2D+yaw test
UNIT=0.5;

xti=[0,2,0];
yti=[0,5,0];
qti=[deg2rad(0),0,0];

xtf=[7,2,0.5];
ytf=[9,3,1];
qtf=[deg2rad(45),0,0];

x_obvp = OBVP;
y_obvp = OBVP;
q_obvp = OBVP;

T=5;dt=0.01;
[x_obvp,flag] = x_obvp.FreeFinalFreeTimeSolve(xti,xtf);
fprintf("x obvp flag :%b \n",flag);
[y_obvp,flag] = y_obvp.FreeFinalFreeTimeSolve(yti,ytf);
fprintf("y obvp flag :%b \n",flag);
[q_obvp,flag] = q_obvp.FreeFinalFreeTimeSolve(qti,qtf);
fprintf("q obvp flag :%b \n",flag);

[xpst,xvst,xast,xjst]= x_obvp.ft(0.01);
[ypst,yvst,yast,yjst]= y_obvp.ft(0.01);
[qpst,qvst,qast,qjst]= q_obvp.ft(0.01);

t=0:dt:T;
L = length(t);

figure
hold on
plot(xpst,ypst,'k.');
for idx=1:50:floor(L)
    plot(xpst(idx),ypst(idx),'bo');
    plot([xpst(idx),xpst(idx)+UNIT*cos(qpst(idx))],[ypst(idx),ypst(idx)+UNIT*sin(qpst(idx))],'b--');
    plot([xpst(idx),xpst(idx)+UNIT*cos(atan(yvst(idx)/xvst(idx)))],[ypst(idx),ypst(idx)+UNIT*sin(atan(yvst(idx)/xvst(idx)))],'r--');
end
grid on
axis equal

figure
subplot(1,3,1)
hold on
plot(t,xpst,'k.');
plot(t,xvst,'b.');
plot(t,xast,'r.');
plot(t,xjst,'g.');
legend('p','v','a','j');
subtitle('x')
grid on
hold off
subplot(1,3,2)
hold on
plot(t,ypst,'k.');
plot(t,yvst,'b.');
plot(t,yast,'r.');
plot(t,yjst,'g.');
legend('p','v','a','j');
subtitle('y')
grid on
hold off
subplot(1,3,3)
hold on
plot(t,qpst,'k.');
plot(t,qvst,'b.');
plot(t,qast,'r.');
plot(t,qjst,'g.');
legend('p','v','a','j');
subtitle('q')
grid on
hold off