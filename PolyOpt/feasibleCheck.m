function checkresult = feasibleCheck(coeffs,ts,segpoly)
%FESIABLECHECK Optimization Result Soft(hard turn to soft) Constrain Check
%   在处理Nonlinear Trajectory Optimization 时,为了降低问题复杂度提高解算速度
%   将

sdfmap = segpoly.sdf;
d_th = segpoly.d_th;
n_seg = segpoly.seg;
Dim = segpoly.Dim;
nodr = segpoly.norder;
CHECK_PLOT = segpoly.CHECK_PLOT;
ORIEN_VEL = segpoly.ORIEN_VEL;
VERDIT_VEL = segpoly.VERDIT_VEL;

%##########################################################################

polytraj = segpoly.traj.setCoeffs(coeffs);
polytraj = polytraj.setTarray(ts);
bars = polytraj.bars;

pv_max =  segpoly.pv_max;
pa_max =  segpoly.pa_max;
wv_max =  segpoly.wv_max;
wa_max =  segpoly.wa_max;

dynamic_rate = segpoly.dyna_rate;
bd_pv_max = dynamic_rate * segpoly.pv_max;
bd_pa_max = dynamic_rate * segpoly.pa_max;
bd_wv_max = dynamic_rate * segpoly.wv_max;
bd_wa_max = dynamic_rate * segpoly.wa_max;
bd_oval   = segpoly.oval_rate;
bd_obstacle_dist = segpoly.d_th;

bd_obstacle_count = 0;
bd_dynamic_vel_count = [0,0,0];
bd_dynamic_acc_count = [0,0,0];
bd_oval_count = 0;

obstacle_count = 0;
obstacle_all = 0;
dynamic_vel_count = [0,0,0];
dynamic_acc_count = [0,0,0];
dynamic_all = 0;
oval_count = 0;
oval_all = 0;
oval_check = zeros(sum(segpoly.traj.bars)+n_seg,1);

for idi = 1:n_seg
    Tts = linspace(0,ts(idi),bars(idi)*5+1);
    xc =  coeffs((idi-1)* nodr* Dim + 1:(idi-1)* nodr* Dim +  nodr);
    yc =  coeffs((idi-1)* nodr* Dim +  nodr + 1:(idi-1)* nodr* Dim +  nodr*2);
    qc =  coeffs((idi-1)* nodr* Dim +  nodr*2 + 1:(idi-1)* nodr* Dim +  nodr*3);
    states = getstatexyq(xc,yc,qc,Tts);
    ids = 1;
    idg = bars(idi)*5+1;
    pos(ids:idg,:) = squeeze(states(1,:,:))';
    obstacle_all = obstacle_all + idg;
    for idj = 0:bars(idi)*5 % 比bars的长度+1
        dotpos = pos(idj+1,:);
        [dotdist,~,~] = sdfmap.getDistAndGrad(dotpos(1),dotpos(2));
        if (dotdist <= 0)
            obstacle_count = obstacle_count + 1;
        end
        if (dotdist <= bd_obstacle_dist)
            bd_obstacle_count = bd_obstacle_count + 1;
        end
    end
end
k = 1;
for idi = 1:n_seg
    [pos,vel,acc] = polytraj.getStates(idi);
    dynamic_all = dynamic_all + polytraj.bars(idi) + 1;
    oval_all = oval_all + polytraj.bars(idi) + 1;
    for idj = 0:polytraj.bars(idi) % 比bars的长度+1
        dotpos = pos(idj+1,:);
        dotvel = vel(idj+1,:);
        dotacc = acc(idj+1,:);
        % 实际超限 ########################################################
        if (dotvel(1) > pv_max)
            dynamic_vel_count(1) = dynamic_vel_count(1) + 1;
        end
        if (dotvel(2) > pv_max)
            dynamic_vel_count(2) = dynamic_vel_count(2) + 1;
        end
        if (dotvel(3) > wv_max)
            dynamic_vel_count(3) = dynamic_vel_count(3) + 1;
        end
        if (dotacc(1) > pa_max)
            dynamic_acc_count(1) = dynamic_acc_count(1) + 1;
        end
        if (dotacc(2) > pa_max)
            dynamic_acc_count(2) = dynamic_acc_count(2) + 1;
        end
        if (dotacc(3) > wa_max)
            dynamic_acc_count(3) = dynamic_acc_count(3) + 1;
        end
        % cost极限边界 ####################################################
        if (dotvel(1) > bd_pv_max)
            bd_dynamic_vel_count(1) = bd_dynamic_vel_count(1) + 1;
        end
        if (dotvel(2) > bd_pv_max)
            bd_dynamic_vel_count(2) = bd_dynamic_vel_count(2) + 1;
        end
        if (dotvel(3) > bd_wv_max)
            bd_dynamic_vel_count(3) = bd_dynamic_vel_count(3) + 1;
        end
        if (dotacc(1) > bd_pa_max)
            bd_dynamic_acc_count(1) = bd_dynamic_acc_count(1) + 1;
        end
        if (dotacc(2) > bd_pa_max)
            bd_dynamic_acc_count(2) = bd_dynamic_acc_count(2) + 1;
        end
        if (dotacc(3) > bd_wa_max)
            bd_dynamic_acc_count(3) = bd_dynamic_acc_count(3) + 1;
        end
        % oval constrain ##################################################
        vel_I = vel(idj+1,1:2); % 惯性系下的速度
        yaw = pos(idj+1,3);
        Rmatrix = [cos(yaw),sin(yaw);-sin(yaw),cos(yaw)];
        vel_B = Rmatrix * vel_I'; % 载体系下的速度
        vel_B2 = vel_B(1)^2/ORIEN_VEL^2+vel_B(2)^2/VERDIT_VEL^2;
        oval_check(k) = vel_B2;
        k = k + 1;
        if (vel_B2 > 1)
            oval_count = oval_count + 1;
        end
        if (vel_B2 > bd_oval)
            bd_oval_count = bd_oval_count + 1;
        end
    end
end

checkresult.obstacle_count  = obstacle_count;
checkresult.obstacle_all    = obstacle_all;
checkresult.dynamic_vel_count = dynamic_vel_count;
checkresult.dynamic_acc_count = dynamic_acc_count;
checkresult.dynamic_all     = dynamic_all;
checkresult.oval_count      = oval_count;
checkresult.oval_all        = oval_all;
checkresult.oval_check      = oval_check;

checkresult.bd_obstacle_count  = bd_obstacle_count;
checkresult.bd_dynamic_vel_count = bd_dynamic_vel_count;
checkresult.bd_dynamic_acc_count = bd_dynamic_acc_count;
checkresult.bd_oval_count      = bd_oval_count;

if (CHECK_PLOT)
    fprintf("obsCheck: whole count =%5d; near bound count =%5d; unfeasible count =%5d \n",obstacle_all,bd_obstacle_count,obstacle_count);
    fprintf("ovaCheck: whole count =%5d; near bound count =%5d; unfeasible count =%5d \n",oval_all,bd_oval_count,oval_count);
    fprintf("dynCheck: whole count =%5d; \n",dynamic_all);
    fprintf("dynCheck vel unfeasible count : x =%5d; y =%5d; q =%5d; \n",dynamic_vel_count(1),dynamic_vel_count(2),dynamic_vel_count(3));
    fprintf("dynCheck vel near bound count : x =%5d; y =%5d; q =%5d; \n",bd_dynamic_vel_count(1),bd_dynamic_vel_count(2),bd_dynamic_vel_count(3));
    fprintf("dynCheck acc unfeasible count : x =%5d; y =%5d; q =%5d; \n",dynamic_acc_count(1),dynamic_acc_count(2),dynamic_acc_count(3));
    fprintf("dynCheck acc near bound count : x =%5d; y =%5d; q =%5d; \n",bd_dynamic_acc_count(1),bd_dynamic_acc_count(2),bd_dynamic_acc_count(3));
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 根据多项式系数和时间计算[pos,vel,acc]
function state = getstatexyq(xc,yc,qc,t)
% TODO: 这里是固定了阶次,之后可以优化代码
tl = length(t);
state(3,3,tl) = 0;
for i=1:tl
    M = getCoeffCons(t(i),7,3);
    xstate = M*xc;
    ystate = M*yc;
    qstate = M*qc;
    state(:,:,i)=[xstate,ystate,qstate];
end
end
