close all
clear
clc
dt=0.01;T=5;
t=0:dt:T;
%% 1维 final state free;final time free
close all
clear
clc
dt=0.01;T=5;
t=0:dt:T;
si=[0,2,3];
sf=[5,0,0];
pi=si(1);vi=si(2);ai=si(3);
pf=sf(1);vf=sf(2);af=sf(3);
pfi = pf - pi;
% ply_J = [-ai^2,-6*ai*vi,-8*(vi^2-ai*delta_p),20*vi*delta_p,-20*delta_p^2];
ply_J1 = [-1/2*ai,-vi,pfi];
ply_J2 = [ai,4*vi,-6*pfi];


% sJ = polyval(ply_J,t);
sJ1 = polyval(ply_J1,t);
sJ2 = polyval(ply_J2,t);

figure
hold on
% plot(t,sJ);
plot(t,sJ1);
plot(t,sJ2);

plot([0,T],[0,0],'b--');
grid on
hold off
% solve_t = [roots(ply_J)];
% solve_t = [roots(ply_J1)];
% solve_t = [solve_t;roots(ply_J2)];

solve_t = [roots(ply_J1);roots(ply_J2)];

realt=[];
L = length(solve_t);
for idx=1:L
    if isreal(solve_t(idx))
        if solve_t(idx)>0
            realt=[realt;solve_t(idx)];
        end
    end
end
if isempty(realt)
    fprintf("optimal T is not found \n");
else
    fprintf("optimal T is found,num = %d \n",length(realt));
end

realt=[0.5,1,1.76,2,5];

L =length(realt);
for idx=1:L
T = realt(idx);

fprintf("optimal T = %6.2f \n",T);
alpha = 20*(pf-pi-vi*T-ai*T^2/2)/T^5;
fprintf("optimal alpha = %6.2f \n",alpha);
costJ = alpha^2*T^4/20*T;
fprintf("optimal J = %6.2f \n",costJ);

obj.c0 = si(1);
obj.c1 = si(2);
obj.c2 = si(3)/2;
obj.c3 = alpha/12*T^2;
obj.c4 = alpha/24*T;
obj.c5 = alpha/120;

pst =obj.c5*t.^5 -obj.c4*t.^4 +obj.c3*t.^3 +obj.c2*t.^2 +obj.c1*t + obj.c0;
vst =5*obj.c5*t.^4 - 4*obj.c4*t.^3 + 3*obj.c3*t.^2 + 2*obj.c2*t + obj.c1 ;
ast =20*obj.c5*t.^3 - 12*obj.c4*t.^2 + 6*obj.c3*t + 2*obj.c2 ;
jst =60*obj.c5*t.^2 - 24*obj.c4*t.^1 + 6*obj.c3;

figure
hold on
plot(t,pst,'k');
plot(t,vst,'b');
plot(t,ast,'r');
plot(t,jst,'g');
plot([T,T],[0,sf(1)],'k-');
legend('p','v','a','j');
grid on
hold off
end

%%
%%%%%%%%%%
% (9*af^2 - 6*af*ao + 9*ao^2)*T^4 + (48*ao*vf - 48*af*vo - 72*af*vf + 72*ao*vo)*T^3 + (192*vf^2 + 336*vf*vo + 192*vo^2 + 120*af*pf - 120*af*po - 120*ao*pf + 120*ao*po)*T^2 + (720*po*vf - 720*pf*vo - 720*pf*vf + 720*po*vo)*T + 720*pf^2 - 1440*pf*po + 720*po^2
% polyn = (- 18*af^2 + 12*af*ao - 18*ao^2)*T^4 + (216*af*vf + 144*af*vo - 144*ao*vf - 216*ao*vo)*T^3 + (- 768*vf^2 - 1344*vf*vo - 768*vo^2 - 480*af*pf + 480*af*po + 480*ao*pf - 480*ao*po)*T^2 + (3600*pf*vf + 3600*pf*vo - 3600*po*vf - 3600*po*vo)*T - 4320*pf^2 + 8640*pf*po - 4320*po^2
% polydn = (- 9*af^2 + 6*af*ao - 9*ao^2)*T^4 + (144*af*vf + 96*af*vo - 96*ao*vf - 144*ao*vo)*T^3 + (- 576*vf^2 - 1008*vf*vo - 576*vo^2 - 360*af*pf + 360*af*po + 360*ao*pf - 360*ao*po)*T^2 + (2880*pf*vf + 2880*pf*vo - 2880*po*vf - 2880*po*vo)*T - 3600*pf^2 + 7200*pf*po - 3600*po^2
% polydn = (- 18*af^2 + 12*af*ao - 18*ao^2)*T^4 + (216*af*vf + 144*af*vo - 144*ao*vf - 216*ao*vo)*T^3 + (- 768*vf^2 - 1344*vf*vo - 768*vo^2 - 480*af*pf + 480*af*po + 480*ao*pf - 480*ao*po)*T^2 + (3600*pf*vf + 3600*pf*vo - 3600*po*vf - 3600*po*vo)*T - 4320*pf^2 + 8640*pf*po - 4320*po^2
%%%%%%%%%%
si=[1,1,1];
sf=[3,2,4];
po=si(1);vo=si(2);ao=si(3);
pi=si(1);vi=si(2);ai=si(3);
pf=sf(1);vf=sf(2);af=sf(3);

%%%%%%%%%%%%%%%%% 1/T
% polyn = [(- 18*af^2 + 12*af*ao - 18*ao^2),...
%      (216*af*vf + 144*af*vo - 144*ao*vf - 216*ao*vo),...
%      (- 768*vf^2 - 1344*vf*vo - 768*vo^2 - 480*af*pf + 480*af*po + 480*ao*pf - 480*ao*po),...
%      (3600*pf*vf + 3600*pf*vo - 3600*po*vf - 3600*po*vo),...
%      - 4320*pf^2 + 8640*pf*po - 4320*po^2];
% solve_t = roots(polyn);
% sJ = polyval(polyn,t);

% polydn = [(- 9*af^2 + 6*af*ao - 9*ao^2),...
%     (144*af*vf + 96*af*vo - 96*ao*vf - 144*ao*vo),...
%     (- 576*vf^2 - 1008*vf*vo - 576*vo^2 - 360*af*pf + 360*af*po + 360*ao*pf - 360*ao*po),...
%     (2880*pf*vf + 2880*pf*vo - 2880*po*vf - 2880*po*vo),...
%     - 3600*pf^2 + 7200*pf*po - 3600*po^2];
% solve_t = roots(polydn);
% sJ = polyval(polydn,t)./t.^6;

% polyn = [(9*af^2 - 6*af*ao + 9*ao^2),...
%     (48*ao*vf - 48*af*vo - 72*af*vf + 72*ao*vo),...
%     (192*vf^2 + 336*vf*vo + 192*vo^2 + 120*af*pf - 120*af*po - 120*ao*pf + 120*ao*po),...
%     (720*po*vf - 720*pf*vo - 720*pf*vf + 720*po*vo),...
%     720*pf^2 - 1440*pf*po + 720*po^2];
% solve_t = roots(polyn);
% sJ = polyval(polyn,t)./t.^5;

%%%%%%%%%%%%%%%%%%%% non 1/T
% polydn = [(- 18*af^2 + 12*af*ao - 18*ao^2),...
%     (216*af*vf + 144*af*vo - 144*ao*vf - 216*ao*vo),...
%     (- 768*vf^2 - 1344*vf*vo - 768*vo^2 - 480*af*pf + 480*af*po + 480*ao*pf - 480*ao*po),...
%     (3600*pf*vf + 3600*pf*vo - 3600*po*vf - 3600*po*vo),...
%     - 4320*pf^2 + 8640*pf*po - 4320*po^2];
% solve_t = roots(polydn);
% sJ = polyval(polydn,t);

%%%%%%%%%%%%%%%%%%% 1/T cost= jerk^2+T
% T^7 + (9*af^2 - 6*af*ao + 9*ao^2)*T^4 + (48*ao*vf - 48*af*vo - 72*af*vf + 72*ao*vo)*T^3 + (192*vf^2 + 336*vf*vo + 192*vo^2 + 120*af*pf - 120*af*po - 120*ao*pf + 120*ao*po)*T^2 + (720*po*vf - 720*pf*vo - 720*pf*vf + 720*po*vo)*T + 720*pf^2 - 1440*pf*po + 720*po^2
polydn = [1,0,0,...
    (9*af^2 - 6*af*ao + 9*ao^2),...
    (48*ao*vf - 48*af*vo - 72*af*vf + 72*ao*vo),...
    (192*vf^2 + 336*vf*vo + 192*vo^2 + 120*af*pf - 120*af*po - 120*ao*pf + 120*ao*po),...
    (720*po*vf - 720*pf*vo - 720*pf*vf + 720*po*vo),...
    720*pf^2 - 1440*pf*po + 720*po^2];
solve_t = roots(polydn);
sJ = polyval(polydn,t);
n=1;
% for T=0:dt:4
% costJ(n) =  ((9*af^2 - 6*af*ao + 9*ao^2)*T^4 + ...
%     (48*ao*vf - 48*af*vo - 72*af*vf + 72*ao*vo)*T^3 +...
%     (192*vf^2 + 336*vf*vo + 192*vo^2 + 120*af*pf - 120*af*po - 120*ao*pf + 120*ao*po)*T^2 +...
%     (720*po*vf - 720*pf*vo - 720*pf*vf + 720*po*vo)*T + ...
%     720*pf^2 - 1440*pf*po + 720*po^2)/T^6 + T;
% n=n+1;
% end
% figure
% hold on
% plot(t,costJ);
% plot([0,T],[0,0],'b--');
% grid on
% hold off

figure
hold on
plot(t,sJ);
plot([0,T],[0,0],'b--');
grid on
hold off



realt=[];
L = length(solve_t);
for idx=1:L
    if isreal(solve_t(idx))
        if solve_t(idx)>0
            realt=[realt;solve_t(idx)];
        end
    end
end
if isempty(realt)
    fprintf("optimal T is not found \n");
else
    fprintf("optimal T is found,num = %d \n",length(realt));
end

% realt=[1,1.76,2];

L =length(realt);
for idx=1:L
T = realt(idx);

fprintf("optimal T = %6.2f \n",T);
alpha = 20*(pf-pi-vi*T-ai*T^2/2)/T^5;
fprintf("optimal alpha = %6.2f \n",alpha);
%%%%% idx=1
% costJ = alpha^2*T^4/20*T;
%%%%% idx=2
% costJ = ((9*af^2 - 6*af*ao + 9*ao^2)*T^4 + ...
%     (48*ao*vf - 48*af*vo - 72*af*vf + 72*ao*vo)*T^3 +...
%     (192*vf^2 + 336*vf*vo + 192*vo^2 + 120*af*pf - 120*af*po - 120*ao*pf + 120*ao*po)*T^2 +...
%     (720*po*vf - 720*pf*vo - 720*pf*vf + 720*po*vo)*T + ...
%     720*pf^2 - 1440*pf*po + 720*po^2)/T^6 + T;
fprintf("optimal J = %6.2f \n",costJ);

obj.c0 = si(1);
obj.c1 = si(2);
obj.c2 = si(3)/2;
obj.c3 = alpha/12*T^2;
obj.c4 = alpha/24*T;
obj.c5 = alpha/120;

pst =obj.c5*t.^5 -obj.c4*t.^4 +obj.c3*t.^3 +obj.c2*t.^2 +obj.c1*t + obj.c0;
vst =5*obj.c5*t.^4 - 4*obj.c4*t.^3 + 3*obj.c3*t.^2 + 2*obj.c2*t + obj.c1 ;
ast =20*obj.c5*t.^3 - 12*obj.c4*t.^2 + 6*obj.c3*t + 2*obj.c2 ;
jst =60*obj.c5*t.^2 - 24*obj.c4*t.^1 + 6*obj.c3;

figure
hold on
plot(t,pst,'k');
plot(t,vst,'b');
plot(t,ast,'r');
plot(t,jst,'g');
plot([T,T],[0,sf(1)],'k-');
legend('p','v','a','j');
grid on
hold off
end
