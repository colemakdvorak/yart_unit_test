addpath_yart
%% Interpolate rotate matrices using Rodrigues' formula
ccc
% Get two poses (SE(3))
p1 = cv([0,0,0]);
p2 = cv([0,3,0]);
R1 = rpy2r(360*rand(1,3)*D2R); % on SO(3), Lie Group
R2 = rpy2r(360*rand(1,3)*D2R);
max_tick = 100;
cnt = 0;
for tick = 1:max_tick % for each tick
    rate = tick/max_tick;
    p_t = (1-rate)*p1 + rate*p2;
    % Interpolate using Rodrigues
    R_t = interp_R(R1,R2,rate); % interpolate
    % Animate
    if mod(tick,5) == 0
        fig_idx = 1;
        fig = set_fig(figure(fig_idx),'pos',[0.6,0.4,0.5,0.55],'view_info',[80,26],...
            'axis_info',[-1.1,+1.1,-1.1,+4.1,-1.1,+1.1],'AXIS_EQUAL',1,'GRID_ON',1,...
            'REMOVE_MENUBAR',1,'USE_DRAGZOOM',1,'SET_CAMLIGHT',1,'SET_MATERIAL','METAL',...
            'SET_AXISLABEL',1,'afs',18);
        plot_T(pr2t(p1,R1),'fig_idx',fig_idx,'subfig_idx',1,'alw',4,'PLOT_AXIS_TIP',1);
        plot_T(pr2t(p2,R2),'fig_idx',fig_idx,'subfig_idx',2,'alw',4,'PLOT_AXIS_TIP',1);
        plot_T(pr2t(p_t,R_t),'fig_idx',fig_idx,'subfig_idx',4,'alw',4,'PLOT_AXIS_TIP',1);
        cnt = cnt + 1;
        plot_T(pr2t(p_t,R_t),'fig_idx',fig_idx,'subfig_idx',4+cnt,'alw',1,...
            'PLOT_AXIS_TIP',1,'atr',0.01);
        plot_title(sprintf('[%d/%d] Interploate SO(3)',tick,max_tick),...
            'fig_idx',fig_idx,'tfs',20,'Interpreter','latex');
        drawnow; if ~ishandle(fig), break; end
    end
end % for tick = 1:max_tick % for each tick

%% Get revolute joint names from root to certain joints of interest
ccc
chain_rig = get_common_rig();
chain_rig = add_joi_to_common_rig(chain_rig,'alpha',0.2,'ADD_ELBOW_GUIDE',1);
chain_rig.joi = get_joi_chain(chain_rig);
joi_names = {'rh','lh'};
rev_joint_names = get_rev_joint_names_route_to_jois(chain_rig,joi_names);

%% Optimization-based smoothing with velocity and acceleration threshold
ccc

% Input data
t_min = 0; t_max = 6; HZ = 50;
n_in = round(HZ*(t_max-t_min));
t_in = linspace(t_min,t_max,n_in)';
x_in = t_in.*sin(2*pi/2*t_in) + 0.0*randn(n_in,1);

% Pre-compute matrices to compute the vel of acc of the trajectory
A_vel = zeros(n_in,n_in);
for i_idx = 2:n_in, A_vel(i_idx,i_idx-1:i_idx) = HZ*[-1,+1]; end
A_acc = zeros(n_in,n_in);
A_acc(2,1:2) = HZ*HZ*[-1,1];
for i_idx = 3:n_in, A_acc(i_idx,i_idx-2:i_idx) = HZ*HZ*[1,-2,1]; end

% Vel and acc of the input data
vel_in = A_vel*x_in;
acc_in = A_acc*x_in;

% Optimize
vel_th = 15;
acc_th = 30;
fun = @(x)( norm(x_in-x,2) );
x0 = x_in;
n_optm_iter = 10000;
opt = optimoptions('fmincon', ...
    'OptimalityTolerance', 0, ...
    'StepTolerance', 0, ...
    'MaxFunctionEvaluations', n_optm_iter,...
    'MaxIterations', n_optm_iter, ...
    'Algorithm','interior-point', ...
    'Display', 'off');
[x_smt,f_val] = fmincon(fun,x0,...
    [A_vel;-A_vel;A_acc;-A_acc],...
    [vel_th*ones(n_in,1);vel_th*ones(n_in,1);acc_th*ones(n_in,1);acc_th*ones(n_in,1)],...
    [],[],[],[],[],opt);

% Vel and acc of the smoothed data
vel_smt = A_vel*x_smt;
acc_smt = A_acc*x_smt;

% Figure 1. Original Trajectory
m = 0.1*(max(x_in)-min(x_in));
axis_info = [t_min,t_max,min(x_in)-m,max(x_in)+m]; axes_info = [0.08,0.12,0.88,0.8];
fig_idx = 1; set_fig(figure(fig_idx),'pos',[0.0,0.7,0.3,0.3],'AXIS_EQUAL',0,...
    'USE_DRAGZOOM',0,'axis_info',axis_info,'ax_str','t','ay_str','x(t)','afs',17,...
    'axes_info',axes_info);
plot(t_in,x_in,'-','color','k','linewidth',2);
title_str = sprintf('Original Trajectory');
plot_title(title_str,'fig_idx',fig_idx,'tfs',16);

% Figure 2. Original Velocity
m = 0.1*(max(vel_in)-min(vel_in));
axis_info = [t_min,t_max,min(vel_in)-m,max(vel_in)+m]; axes_info = [0.08,0.12,0.88,0.8];
fig_idx = 2; set_fig(figure(fig_idx),'pos',[0.3,0.7,0.3,0.3],'AXIS_EQUAL',0,...
    'USE_DRAGZOOM',0,'axis_info',axis_info,'ax_str','t','ay_str','$\dot{x}(t)$','afs',17,...
    'axes_info',axes_info);
plot(t_in,vel_in,'-','color','k','linewidth',2);
plot([t_in(1),t_in(end)],+vel_th*[1,1],'--','color',0.5*[1,1,1],'linewidth',2);
plot([t_in(1),t_in(end)],-vel_th*[1,1],'--','color',0.5*[1,1,1],'linewidth',2);
title_str = sprintf('Original Velocity');
plot_title(title_str,'fig_idx',fig_idx,'tfs',16);

% Figure 3. Original Acceleration
m = 0.1*(max(acc_in)-min(acc_in));
axis_info = [t_min,t_max,min(acc_in)-m,max(acc_in)+m]; axes_info = [0.08,0.12,0.88,0.8];
fig_idx = 3; set_fig(figure(fig_idx),'pos',[0.6,0.7,0.3,0.3],'AXIS_EQUAL',0,...
    'USE_DRAGZOOM',0,'axis_info',axis_info,'ax_str','t','ay_str','$\ddot{x}(t)$','afs',17,...
    'axes_info',axes_info);
plot(t_in,acc_in,'-','color','k','linewidth',2);
plot([t_in(1),t_in(end)],+acc_th*[1,1],'--','color',0.5*[1,1,1],'linewidth',2);
plot([t_in(1),t_in(end)],-acc_th*[1,1],'--','color',0.5*[1,1,1],'linewidth',2);
title_str = sprintf('Original Acceleration');
plot_title(title_str,'fig_idx',fig_idx,'tfs',16);

% Figure 4. Smoothed Trajectory
m = 0.1*(max(x_in)-min(x_in));
axis_info = [t_min,t_max,min(x_in)-m,max(x_in)+m]; axes_info = [0.08,0.12,0.88,0.8];
fig_idx = 4; set_fig(figure(fig_idx),'pos',[0.0,0.35,0.3,0.3],'AXIS_EQUAL',0,...
    'USE_DRAGZOOM',0,'axis_info',axis_info,'ax_str','t','ay_str','x(t)','afs',17,...
    'axes_info',axes_info);
h_x_in = plot(t_in,x_in,'-','color','k','linewidth',2);
h_x_smt = plot(t_in,x_smt,'-','color','b','linewidth',2);
legend([h_x_in,h_x_smt],...
    {'Original Trajectory','Smoothed Trajectory'},...
    'fontsize',12,'fontname','consolas','location','northwest');
title_str = sprintf('Smoothed Trajectory');
plot_title(title_str,'fig_idx',fig_idx,'tfs',16);

% Figure 5. Smoothed Velocity
m = 0.1*(max(vel_in)-min(vel_in));
axis_info = [t_min,t_max,min(vel_in)-m,max(vel_in)+m]; axes_info = [0.08,0.12,0.88,0.8];
fig_idx = 5; set_fig(figure(fig_idx),'pos',[0.3,0.35,0.3,0.3],'AXIS_EQUAL',0,...
    'USE_DRAGZOOM',0,'axis_info',axis_info,'ax_str','t','ay_str','$\dot{x}(t)$','afs',17,...
    'axes_info',axes_info);
h_vel_in = plot(t_in,vel_in,'-','color','k','linewidth',2);
h_vel_smt = plot(t_in,vel_smt,'-','color','b','linewidth',2);
h_vel_th = plot([t_in(1),t_in(end)],+vel_th*[1,1],'--','color',0.5*[1,1,1],'linewidth',2);
plot([t_in(1),t_in(end)],-vel_th*[1,1],'--','color',0.5*[1,1,1],'linewidth',2);
legend([h_vel_in,h_vel_smt,h_vel_th],...
    {'Original Velocity','Smoothed Velocity','Velocity Threshold'},...
    'fontsize',12,'fontname','consolas','location','northwest');
title_str = sprintf('Smoothed Velocity');
plot_title(title_str,'fig_idx',fig_idx,'tfs',16);

% Figure 6. Smoothed Acceleration
m = 0.1*(max(acc_in)-min(acc_in));
axis_info = [t_min,t_max,min(acc_in)-m,max(acc_in)+m]; axes_info = [0.08,0.12,0.88,0.8];
fig_idx = 6; set_fig(figure(fig_idx),'pos',[0.6,0.35,0.3,0.3],'AXIS_EQUAL',0,...
    'USE_DRAGZOOM',0,'axis_info',axis_info,'ax_str','t','ay_str','$\ddot{x}(t)$','afs',17,...
    'axes_info',axes_info);
h_acc_in = plot(t_in,acc_in,'-','color','k','linewidth',2);
h_acc_smt = plot(t_in,acc_smt,'-','color','b','linewidth',2);
h_acc_th = plot([t_in(1),t_in(end)],+acc_th*[1,1],'--','color',0.5*[1,1,1],'linewidth',2);
plot([t_in(1),t_in(end)],-acc_th*[1,1],'--','color',0.5*[1,1,1],'linewidth',2);
legend([h_acc_in,h_acc_smt,h_acc_th],...
    {'Original Acceleration','Smoothed Acceleration','Acceleration Threshold'},...
    'fontsize',12,'fontname','consolas','location','northwest');
title_str = sprintf('Smoothed Acceleration');
plot_title(title_str,'fig_idx',fig_idx,'tfs',16);

%% Define addition and subtraction on 2-Sphere
ccc
% Configuration
USE_UV_ADD      = 1;
USE_SPH         = 1;
USE_QUAT        = 1;
SAVE_VID        = 0;

% Random unit vectors
uv1 = uv(cv([1,1,2]));
uv2 = uv(uv1 + randn(3,1));
uv3 = uv(randn(3,1));

% v4 = v3 + (v2-v1)
% Simple formulation
uv4_simple = uv(uv3+uv2-uv1);
% Spherical formulation
[a1,e1,~] = cart2sph(uv1(1),uv1(2),uv1(3));
[a2,e2,~] = cart2sph(uv2(1),uv2(2),uv2(3));
[a3,e3,~] = cart2sph(uv3(1),uv3(2),uv3(3));
[x,y,z] = sph2cart(a3-a1+a2,e3-e1+e2,1);
uv4_sph = cv([x,y,z]);
% Quaternion formulation
q_1to2 = get_q_uv1_to_uv2(uv1,uv2);
R_1to2 = quat2r(q_1to2); % rotation corresponding to the quaternion
rpy_1to2 = quat2rpy(q_1to2); % Euler angle corresponding to the quaternion
uv4_quat = quat_action(uv3,q_1to2);

% Slerp
L = 100;
g1to2 = zeros(L,3);
for tick = 1:L
    g1to2(tick,:) = rv(slerp(uv1,uv2,tick/L));
end
g3to4_simple = zeros(L,3);
for tick = 1:L
    g3to4_simple(tick,:) = rv(slerp(uv3,uv4_simple,tick/L));
end
g3to4_sph = zeros(100,3);
for tick = 1:L
    g3to4_sph(tick,:) = rv(slerp(uv3,uv4_sph,tick/L));
end
g3to4_quat = zeros(100,3);
for tick = 1:L
    g3to4_quat(tick,:) = rv(slerp(uv3,uv4_quat,tick/L));
end

% Plot
vid_path = sprintf('../vid/unit_test/2-sphere_%s.mp4',datestr(now,'yy_mm_dd_HH_MM_SS'));
vobj = init_vid_record(vid_path,'HZ',50,'SAVE_VID',SAVE_VID);
fig_idx = 1;
fig = set_fig(figure(fig_idx),'pos',[0.0,0.5,0.3,0.5],'view_info',[80,26],...
    'axis_info',1.2*[-1,+1,-1,+1,-1,+1],'AXIS_EQUAL',1,'GRID_ON',1,...
    'REMOVE_MENUBAR',1,'USE_DRAGZOOM',1,'SET_CAMLIGHT',1,'SET_MATERIAL','METAL',...
    'SET_AXISLABEL',1,'afs',18);
plot_T(pr2t('',''),'fig_idx',fig_idx,'PLOT_AXIS',0,...
    'PLOT_SPHERE',1,'sr',1.0,'sfc',0.9*[1,1,1],'sfa',0.1,'sec',0.5*[1,1,1]);
sw = 0.05; tw = 0.1; text_fs = 20; text_p2_offset = 0.08;
plot_arrow_3d(cv([0,0,0]),uv1,'fig_idx',fig_idx,'subfig_idx',1,...
    'alpha',0.5,'color','r','sw',sw,'tw',tw,...
    'text_str','$\mathbf{v}_1$','text_fs',text_fs,'text_p2_offset',text_p2_offset,...
    'interpreter','latex');
plot_arrow_3d(cv([0,0,0]),uv2,'fig_idx',fig_idx,'subfig_idx',2,...
    'alpha',0.5,'color','b','sw',sw,'tw',tw,...
    'text_str','$\mathbf{v}_2$','text_fs',text_fs,'text_p2_offset',text_p2_offset,...
    'interpreter','latex');
plot_traj(g1to2,'fig_idx',fig_idx,'subfig_idx',1,'tlc','b','tlw',2);
plot_arrow_3d(cv([0,0,0]),uv3,'fig_idx',fig_idx,'subfig_idx',3,...
    'alpha',0.5,'color','k','sw',sw,'tw',tw,...
    'text_str','$\mathbf{v}_3$','text_fs',text_fs,'text_p2_offset',text_p2_offset,...
    'interpreter','latex');
if USE_UV_ADD
    plot_traj(g3to4_simple,'fig_idx',fig_idx,'subfig_idx',2,'tlc','m','tlw',3);
end
if USE_SPH
    plot_traj(g3to4_sph,'fig_idx',fig_idx,'subfig_idx',4,'tlc','c','tlw',3);
end
if USE_QUAT
    plot_traj(g3to4_quat,'fig_idx',fig_idx,'subfig_idx',5,'tlc','y','tlw',3);
end
rot_axis = uv(q_1to2(2:4));
plot_plane('fig_idx',fig_idx,'subfig_idx',1,...
    'xmin',-1,'xmax',1,'xres',0.1,'ymin',-1,'ymax',1,'yres',0.1,...
    'plane_normal',rot_axis,'plane_center',cv([0,0,0]),'pfc',0.5*[1,1,1]);
plot_arrow_3d(cv([0,0,0]),rot_axis*0.5,'fig_idx',fig_idx,'subfig_idx',7,...
    'alpha',0.5,'color','k','sw',sw,'tw',tw,...
    'text_str','Rotation Axis','text_fs',text_fs,'text_p2_offset',text_p2_offset,...
    'interpreter','latex');
% Loop
for tick = round(linspace(1,L,20)) % loop
    v4_simple_t = cv(g3to4_simple(tick,:));
    v4_sph_t = cv(g3to4_sph(tick,:));
    v4_quat_t = cv(g3to4_quat(tick,:));
    if USE_UV_ADD
        plot_arrow_3d(cv([0,0,0]),v4_simple_t,'fig_idx',fig_idx,'subfig_idx',4,...
            'alpha',0.5,'color','m','sw',sw,'tw',tw,...
            'text_str','$\mathbf{v}_{4}$-simple','text_fs',text_fs,...
            'text_p2_offset',text_p2_offset,'interpreter','latex');
    end
    if USE_SPH
        plot_arrow_3d(cv([0,0,0]),v4_sph_t,'fig_idx',fig_idx,'subfig_idx',5,...
            'alpha',0.5,'color','c','sw',sw,'tw',tw,...
            'text_str','$\mathbf{v}_{4}$-sph','text_fs',text_fs,...
            'text_p2_offset',text_p2_offset,'interpreter','latex');
    end
    if USE_QUAT
        plot_arrow_3d(cv([0,0,0]),v4_quat_t,'fig_idx',fig_idx,'subfig_idx',6,...
            'alpha',0.5,'color','y','sw',sw,'tw',tw,...
            'text_str','$\mathbf{v}_{4}$-quat','text_fs',text_fs,...
            'text_p2_offset',text_p2_offset,'interpreter','latex');
    end
    title_str = sprintf(['[%d/%d] ',...
        '$ \\mathbf{v}_4 = \\mathbf{v}_3 - \\mathbf{v}_1 + \\mathbf{v}_2$'],...
        tick,L);
    plot_title(title_str,'fig_idx',fig_idx,'tfs',20,'interpreter','latex','ALLOW_UNDERBAR',1);
    drawnow;
    record_vid(vobj,'fig',fig);
end % for tick = 1:L % loop
end_vid_record(vobj);

%% Normal plane on a 2-Sphere
ccc
uv1 = uv(randn(3,1));
uv2 = uv(randn(3,1));
L = 100;
g1to2 = zeros(L,3);
for tick = 1:L, g1to2(tick,:) = rv(slerp(uv1,uv2,tick/L)); end

% Plot
fig_idx = 1;
fig = set_fig(figure(fig_idx),'pos',[0.0,0.5,0.3,0.5],'view_info',[80,26],...
    'axis_info',1.5*[-1,+1,-1,+1,-1,+1],'AXIS_EQUAL',1,'GRID_ON',1,...
    'REMOVE_MENUBAR',1,'USE_DRAGZOOM',1,'SET_CAMLIGHT',1,'SET_MATERIAL','METAL',...
    'SET_AXISLABEL',1,'afs',18);
plot_T(pr2t('',''),'fig_idx',fig_idx,'PLOT_AXIS',0,...
    'PLOT_SPHERE',1,'sr',1.0,'sfc',0.9*[1,1,1],'sfa',0.1,'sec',0.5*[1,1,1]);
sw = 0.04; tw = 0.1; text_fs = 20; text_p2_offset = 0.08;
plot_arrow_3d(cv([0,0,0]),uv1,'fig_idx',fig_idx,'subfig_idx',1,...
    'alpha',0.5,'color','r','sw',sw,'tw',tw,...
    'text_str','$\mathbf{v}_1$','text_fs',text_fs,'text_p2_offset',text_p2_offset,...
    'interpreter','latex');
plot_arrow_3d(cv([0,0,0]),uv2,'fig_idx',fig_idx,'subfig_idx',2,...
    'alpha',0.5,'color','b','sw',sw,'tw',tw,...
    'text_str','$\mathbf{v}_2$','text_fs',text_fs,'text_p2_offset',text_p2_offset,...
    'interpreter','latex');
plot_traj(g1to2,'fig_idx',fig_idx,'subfig_idx',1,'tlc','b','tlw',2);
plot_plane('fig_idx',fig_idx,'subfig_idx',1,...
    'xmin',-1,'xmax',1,'xres',0.1,'ymin',-1,'ymax',1,'yres',0.1,...
    'plane_normal',uv1,'plane_center',uv1,'pfc',0.5*[1,1,1]);
plot_title('Normal Plane',...
    'fig_idx',fig_idx,'tfs',20,'interpreter','latex','ALLOW_UNDERBAR',1);

%% Find a rotation offset to align a specific axis of a rotation matrix.
ccc
% Random rotation matrix
R_init = rpy2r(360*D2R*rand(1,3));
% Specific axis of 'R_init' to be aligned with a target axis
which_axis  = cv([1,0,0]); % first axis: red
target_axis = cv([0,0,1]); % third axis: blue
reg_coef    = 0.01;
fun = @(rpy) (...
    norm(rpy2r(rpy)*R_init*which_axis-target_axis) ...      % original objective
    + reg_coef*(abs(rpy(1))+abs(rpy(2))+abs(rpy(3))) ...    % regularizer
    );
rpy_hat = fminsearch(fun,cv([0,0,0]));
R_hat = rpy2r(rpy_hat)*R_init;
% Interpolate 'R_init' and 'R_hat'
p1 = cv([0,0,0]); p2 = cv([0,2,0]);
R1 = R_init; R2 = R_hat;
R_link = R1'*R2; % rotation matrix which links R1 and R2
w_link = r2w(R_link); % equivalent velocity vector fron the rotation matrix
L = 30;
ts = linspace(0,1,L); % t:[0,1]
% Loop
x_traj = []; y_traj = []; z_traj = [];
for tick = 1:length(ts)
    % Interpolate R1 and R2
    t   = ts(tick);
    p_t = (1-t)*p1 + t*p2;
    R_t = R1*rodrigues(w_link/norm(w_link),norm(w_link)*t); % interpolate with slerp
    % Append the trajectory
    all = 0.75;
    x_traj = cat(1,x_traj, rv(p_t)+all*rv(R_t(:,1)));
    y_traj = cat(1,y_traj, rv(p_t)+all*rv(R_t(:,2)));
    z_traj = cat(1,z_traj, rv(p_t)+all*rv(R_t(:,3)));
    % Plot
    set_fig(figure(1),'pos',[0.5,0.4,0.5,0.55],...
        'view_info',[88,10],'axis_info',[-1,+1,-1,+3,-1,+1],'AXIS_EQUAL',1,'GRID_ON',1,...
        'REMOVE_MENUBAR',1,'USE_DRAGZOOM',1,'SET_CAMLIGHT',1,'SET_MATERIAL','METAL',...
        'SET_AXISLABEL',1,'afs',18);
    plot_T(pr2t('',''),...
        'fig_idx',1,'subfig_idx',1,'PLOT_AXIS',1,'all',1.0,'alw',1,'PLOT_AXIS_TIP',0);
    plot_T(pr2t('',R_init),...
        'fig_idx',1,'subfig_idx',2,'PLOT_AXIS',1,'all',all,'alw',3,'PLOT_AXIS_TIP',1);
    plot_T(pr2t(p2,R_hat),...
        'fig_idx',1,'subfig_idx',3,'PLOT_AXIS',1,'all',all,'alw',3,'PLOT_AXIS_TIP',1);
    plot_T(pr2t(p_t,R_t),...
        'fig_idx',1,'subfig_idx',4,'PLOT_AXIS',1,'all',all,'alw',3,'PLOT_AXIS_TIP',1);
    plot_traj(x_traj,'fig_idx',1,'subfig_idx',2,'tlc','r','tlw',1,'tls','--');
    plot_traj(y_traj,'fig_idx',1,'subfig_idx',3,'tlc','g','tlw',1,'tls','--');
    plot_traj(z_traj,'fig_idx',1,'subfig_idx',4,'tlc','b','tlw',1,'tls','--');
    title_str = sprintf('[%d/%d] Align the local y-axis to the global z-axis',tick,L);
    plot_title(title_str,'tfs',25,'interpreter','latex');
    drawnow;
end

%% Shorstest distance from a line segment to a point (works for both 2D and 3D)
ccc
% Random line segment and point
p_a = rand(3,1);
p_b = rand(3,1);
p_x = rand(3,1);
% Shortest distance from line to point
[min_dist,p_closest] = get_dist_point2line(p_x,p_a,p_b);
% Plot
fig_idx = 1;
marg = 0.1;
set_fig(figure(fig_idx),'pos',[0.0,0.4,0.3,0.5],'view_info',[80,10],...
    'axis_info',[0,1,0,1,0,1]+marg*[-1,1,-1,1,-1,1],'AXIS_EQUAL',1,'GRID_ON',1,...
    'REMOVE_MENUBAR',1,'USE_DRAGZOOM',1,'SET_CAMLIGHT',1,'SET_MATERIAL','METAL',...
    'SET_AXISLABEL',1,'afs',18);
plot_T(pr2t(p_a,''),'fig_idx',fig_idx,'subfig_idx',1,...
    'PLOT_AXIS',0,'PLOT_SPHERE',1,'sr',0.02,'sfc','r','sfa',0.9,...
    'text_str','~$p_a$','text_interp','latex','text_fs',20);
plot_T(pr2t(p_b,''),'fig_idx',fig_idx,'subfig_idx',2,...
    'PLOT_AXIS',0,'PLOT_SPHERE',1,'sr',0.02,'sfc','b','sfa',0.9,...
    'text_str','~$p_b$','text_interp','latex','text_fs',20);
plot_line(p_a,p_b,'fig_idx',fig_idx,'subfig_idx',1,'lc','k','lw',2);
plot_T(pr2t(p_x,''),'fig_idx',fig_idx,'subfig_idx',3,...
    'PLOT_AXIS',0,'PLOT_SPHERE',1,'sr',0.02,'sfc','k','sfa',0.9,...
    'text_str','~$p_x$','text_interp','latex','text_fs',20);
plot_T(pr2t(p_closest,''),'fig_idx',fig_idx,'subfig_idx',4,...
    'PLOT_AXIS',0,'PLOT_SPHERE',1,'sr',0.01,'sfc','k','sfa',0.9);
plot_T(pr2t(p_closest,''),'fig_idx',fig_idx,'subfig_idx',5,...
    'PLOT_AXIS',0,'PLOT_SPHERE',1,'sr',0.03,'sfc','g','sfa',0.3,...
    'text_str','~$p_c$','text_interp','latex','text_fs',20);
plot_line(p_closest,p_x,'fig_idx',fig_idx,'subfig_idx',2,'lc','k','lw',2,'ls','--');
plot_title('Shortest Distance from Line to Point','tfs',20);

%% Interpolation of unit vectors using slerp algorithm.
ccc
max_tick = 20;
colors = linspecer(max_tick);
uv1 = uv(randn(3,1));
uv2 = uv(randn(3,1));
set_fig(figure(1),'pos',[0.5,0.5,0.3,0.5],...
    'view_info',[88,10],'axis_info',1.2*[-1,+1,-1,+1,-1,+1],'AXIS_EQUAL',1,'GRID_ON',0,...
    'REMOVE_MENUBAR',1,'USE_DRAGZOOM',1,'SET_CAMLIGHT',1,'SET_MATERIAL','METAL',...
    'SET_AXISLABEL',1,'afs',18);
plot_T(pr2t('',''),'alw',3);
plot_arrow_3d(cv([0,0,0]),uv1,'subfig_idx',1,'color',colors(1,:));
plot_arrow_3d(cv([0,0,0]),uv2,'subfig_idx',2,'color',colors(end,:));
for tick = 1:max_tick
    t = tick/max_tick;
    uv_t = slerp(uv1,uv2,t);
    color = colors(tick,:);
    plot_arrow_3d(cv([0,0,0]),uv_t,'subfig_idx',2+tick,'color',color);
    title_str = sprintf('[%d/%d] Slerp',tick,max_tick);
    plot_title(title_str,'fig_idx',1,'tfs',20,'Interpreter','Latex');
end

%% SE(3) Pose transformations
% Post-multiplication : local transformation
% Pre-multiplication  : global transformation
ccc
% Configuration
transform_type = 'Global'; % 'Global', 'Local'
rotation_axis  = 'Z'; % 'X', 'Y', 'Z'
% Random homogeneous transformation matrix
T_w = pr2t(cv(-1+2*rand(1,3)),rpy2r(360*rand(1,3)*D2R));
root_traj = []; x_traj = []; y_traj = []; z_traj = [];
max_tick = 360;
for tick = 1:max_tick % for each tick
    % Coordinate transform
    switch rotation_axis
        case 'X'
            T_w2a = pr2t(cv([0,0,0]),rpy2r(tick*[1,0,0]*D2R)); % rotate w.r.t. x-axis
        case 'Y'
            T_w2a = pr2t(cv([0,0,0]),rpy2r(tick*[0,1,0]*D2R)); % rotate w.r.t. y-axis
        case 'Z'
            T_w2a = pr2t(cv([0,0,0]),rpy2r(tick*[0,0,1]*D2R)); % rotate w.r.t. z-axis
    end
    switch lower(transform_type)
        case 'global'
            T_a = T_w2a*T_w; % pre-multiplication (global)
        case 'local'
            T_a = T_w*T_w2a; % post-multiplication (local)
    end
    R_t = t2r(T_a); % get the rotation matrix
    % Append the axes trajectories
    root_traj   = [root_traj; rv(t2p(T_a))];
    x_traj      = [x_traj; rv(t2p(T_a)) + 0.5*rv(R_t(:,1))];
    y_traj      = [y_traj; rv(t2p(T_a)) + 0.5*rv(R_t(:,2))];
    z_traj      = [z_traj; rv(t2p(T_a)) + 0.5*rv(R_t(:,3))];
    % Animate
    if mod(tick,5) == 0 % animate
        fig = set_fig(figure(1),'pos',[0.0,0.5,0.3,0.4],...
            'view_info',[80,26],'axis_info',1.6*[-1,+1,-1,+1,-1,+1],'AXIS_EQUAL',1,'GRID_ON',1,...
            'REMOVE_MENUBAR',1,'USE_DRAGZOOM',1,'SET_CAMLIGHT',1,'SET_MATERIAL','METAL',...
            'SET_AXISLABEL',1,'afs',18,'interpreter','latex','NO_MARGIN',0);
        plot_T(pr2t(cv([0,0,0]),eye(3,3)),'fig_idx',1,'subfig_idx',1,...
            'PLOT_AXIS',1,'all',1.0,'alw',3,'PLOT_SPHERE',0,...
            'text_str','World','text_fs',20,'text_interp','latex'); % world coordinate
        plot_T(T_w,'fig_idx',1,'subfig_idx',2,...
            'PLOT_AXIS',1,'all',0.5,'alw',2,'PLOT_SPHERE',0); % initial coordinate
        plot_T(T_a,'fig_idx',1,'subfig_idx',3,...
            'PLOT_AXIS',1,'all',0.5,'alw',2,'PLOT_AXIS_TIP',1,'atr',0.05,...
            'PLOT_SPHERE',1,'sr',0.03,'sfc','k',...
            'text_str','T','text_fs',20,'text_interp','latex'); % transformed coordinate
        plot_traj(root_traj,'fig_idx',1,'subfig_idx',1,'tlc','k','tlw',3);
        plot_traj(x_traj,'fig_idx',1,'subfig_idx',2,'tlc','r','tlw',1);
        plot_traj(y_traj,'fig_idx',1,'subfig_idx',3,'tlc','g','tlw',1);
        plot_traj(z_traj,'fig_idx',1,'subfig_idx',4,'tlc','b','tlw',1);
        title_str = sprintf('[%d/%d] %s Transform w.r.t. %s axis',...
            tick,max_tick,transform_type,rotation_axis);
        plot_title(title_str,'tfs',20,'interpreter','latex');
        drawnow;
        if ~ishandle(fig), break; end
    end % if mod(tick,5) == 0 % animate
end % for tick = 1:360 % for each tick
fprintf('Done.\n');

%% Voronoi optimistic optimization (VOO)
ccc
% Set grid (domain)
res = 200;
[xgrid,ygrid] = meshgrid(linspace(0,10,res),linspace(0,10,res));
xygrid = [xgrid(:),ygrid(:)];

% Set the cost function using GMM
mu = [8,2;2,8;8,8];
sigma = cat(3,[1.0,0.0;0.0,1.0],[1.0,0.0;0.0,1.0],...
    [1.0,0.0;0.0,1.0]);
probs = [1,1,2]; probs = probs / sum(probs);
gm = gmdistribution(mu,sigma,probs);
f = @(x)((gm.pdf(x))); % this will be our score function

% Run VOO
n_sample = 10000; dim = 2; max_exploit = 100; omega = 0.2;
px = @(n)(10*rand(n,2)); % sampler
voo = init_voo(n_sample,dim,max_exploit,omega,f,px);
tk = init_tk('VOO');
while voo_not_finished(voo)
    voo = one_step_voo(voo);
    tk = print_tk(tk,voo.tick,voo.n_sample);
end
voo = summarize_voo(voo);

% Figure1: Plot the score map
fig1 = figure(1); set_fig(fig1,'pos',[0.0,0.5,0.3,0.4]);
fval = f(xygrid); fmtx = reshape(fval,res,res); % reshape
h = pcolor(xgrid,ygrid,fmtx); colormap(linspecer);
cb = colorbar; cb.FontSize = 12; cb.FontName = 'consolas';
set(h,'EdgeColor','none','FaceAlpha',0.8); set(gcf,'color','w');
plot_title('Score Map','fig_idx',1,'tfs',20);
axis([0,10,0,10]); axis('on', 'image'); axis off; dragzoom;

% Figure2: Plot the results of VOO
fig2 = figure(2); set_fig(fig2,'pos',[0.3,0.5,0.3,0.4]); hold on;
fval = f(xygrid); fmtx = reshape(fval,res,res); % reshape
h = pcolor(xgrid,ygrid,fmtx); colormap(linspecer);
plot(voo.x_list(:,1),voo.x_list(:,2),'kx');
cb = colorbar; cb.FontSize = 12; cb.FontName = 'consolas';
set(h,'EdgeColor','none','FaceAlpha',0.3); set(gcf,'color','w');
axis([0,10,0,10]); axis('on', 'image'); axis off; dragzoom;
plot_title('Samples of VOO','fig_idx',2,'tfs',20);

% Figure3: 2D historgram of the points
fig3 = figure(3); set_fig(fig3,'pos',[0.6,0.5,0.3,0.4]); hold on;
ndhist(voo.x_list,'bins',1.0,'axis',[0,10,0,10]); caxis([0,inf]);
cb = colorbar; cb.FontSize = 12; cb.FontName = 'consolas';
plot_title('Histogram of VOO Samples','fig_idx',3,'tfs',20);
axis off; dragzoom;

%% Sampling with Score Matching using Hierachical-VOO
ccc
rng(1);
% Set grid (domain)
res = 200;
[xgrid,ygrid] = meshgrid(linspace(0,10,res),linspace(0,10,res));
xygrid = [xgrid(:),ygrid(:)];

% Set the cost function using GMM
mu = [8,2;2,8;8,8];
sigma = cat(3,[1.0,0.0;0.0,1.0],[1.0,0.0;0.0,1.0],...
    [1.0,0.0;0.0,1.0]);
probs = [1,1,2]; probs = probs / sum(probs);
gm = gmdistribution(mu,sigma,probs);
f = @(x)((gm.pdf(x))); % this will be our score function
px = @(n)(10*rand(n,2)); % sampler

% Initialize HVOO
n_sample = 10000; dim = 2;
omega = 0.2; % epsilon-greedy
n_anchor = 30; % number of anchors
hvoo = init_hvoo(n_sample,dim,omega,n_anchor,f,px);
tk = init_tk('HVOO');
while hvoo_not_finished(hvoo)
    hvoo = one_step_hvoo(hvoo);
    tk = print_tk(tk,hvoo.tick,hvoo.n_sample);
end
hvoo = summarize_voo(hvoo);

% Figure1: Plot the score map
fig1 = figure(1); set_fig(fig1,'pos',[0.0,0.5,0.3,0.4]);
fval = f(xygrid); fmtx = reshape(fval,res,res); % reshape
h = pcolor(xgrid,ygrid,fmtx); colormap(linspecer);
cb = colorbar; cb.FontSize = 12; cb.FontName = 'consolas';
set(h,'EdgeColor','none','FaceAlpha',0.8); set(gcf,'color','w');
plot_title('Score Map','fig_idx',1,'tfs',20);
axis([0,10,0,10]); axis('on', 'image'); axis off; dragzoom;
drawnow;

% Figure2: Plot the results of H-VOO
fig2 = figure(2); set_fig(fig2,'pos',[0.3,0.5,0.3,0.4]); hold on;
fval = f(xygrid); fmtx = reshape(fval,res,res); % reshape
h = pcolor(xgrid,ygrid,fmtx); colormap(linspecer);
plot(hvoo.x_list(:,1),hvoo.x_list(:,2),'kx');
plot(hvoo.x_anchors(:,1),hvoo.x_anchors(:,2),'ko','MarkerSize',12,'LineWidth',2,...
    'MarkerFaceColor','w');
SHOW_REGION_INFO = 0;
if SHOW_REGION_INFO
    for r_idx = 1:hvoo.k
        str = sprintf('[%d] f:%.2f n:%d',...
            r_idx,hvoo.favg_regions(r_idx),hvoo.n_regions(r_idx));
        text(hvoo.x_anchors(r_idx,1)+0.25,hvoo.x_anchors(r_idx,2),str,...
            'fontsize',12,'BackgroundColor','w');
    end
end
cb = colorbar; cb.FontSize = 12; cb.FontName = 'consolas';
set(h,'EdgeColor','none','FaceAlpha',0.3); set(gcf,'color','w');
plot_title('Samples of H-VOO','fig_idx',2,'tfs',20);
axis([0,10,0,10]); axis('on', 'image'); axis off; dragzoom;

% Figure3: 2D historgram of the points
fig3 = figure(3); set_fig(fig3,'pos',[0.6,0.5,0.3,0.4]); hold on;
[ex,ey,nh] = ndhist(hvoo.x_list,'bins',1.0,'axis',[0,10,0,10]); caxis([0,inf]);
cb = colorbar; cb.FontSize = 12; cb.FontName = 'consolas'; set(gcf,'color','w');
plot_title('Histogram of H-VOO Samples','fig_idx',3,'tfs',20);
axis off; dragzoom;

%% Locally Approximating Determinantal Point Process (LA-DPP)
% ccc
n = 1000;
x = [rand([n/2,2]); 0.2*rand([n/2,2])]; % imbalanced inputs
k = 20;
ladpp_idx = get_sub_idx_ladpp(x,k); % get subset indices using LA-DPP
temp = randperm(size(x,1)); unif_idx = temp(1:k);
x_ladpp = x(ladpp_idx,:);
x_unif = x(unif_idx,:);
fig = figure();
set_fig(fig,'pos',[0.0,0.5,0.3,0.4],'axis_info',[0,1,0,1]);
hr = plot(x(:,1),x(:,2),'ko','MarkerSize',4,'LineWidth',1,'MarkerFaceColor',0.8*[1,1,1]);
hs = plot(x_ladpp(:,1),x_ladpp(:,2),'ro','MarkerSize',15,'LineWidth',3);
hu = plot(x_unif(:,1),x_unif(:,2),'bs','MarkerSize',15,'LineWidth',3);
legend([hr,hs,hu],{'Raw Data','LA-DPP Sampling','Uniform Sampling'},...
    'fontsize',15,'interpreter','Latex');
plot_title('Uniform Sampling and LA-DPP Sampling','tfs',20,'interpreter','latex');

%% Surrounding capsule optimization
ccc

% Load a stl file and optimize the surrounding capsule
stl_path = '../../yet-another-robotics-toolbox/urdf/ur10/visual/shoulder.stl';
% stl_path = '../yet-another-robotics-toolbox/urdf/ur10/visual/forearm.stl';
[~,file_name,~] = fileparts(stl_path);
fv = load_stl(stl_path);
cap = optimize_capsule(fv);

% Plot the loaded mesh and optimized capsule
fig_idx = 1; fig_pos = [0.5,0.5,0.3,0.5];
view_info = [86,10]; axis_info = [-0.3,0.3,-0.3,0.3,-0.15,0.7]; AXIS_EQUAL = 1; GRID_ON = 1;
REMOVE_MENUBAR = 1; USE_DRAGZOOM = 1;
fig = set_fig(figure(fig_idx),...
    'pos',fig_pos,'view_info',view_info,'axis_info',axis_info,'AXIS_EQUAL',AXIS_EQUAL,...
    'GRID_ON',GRID_ON,'REMOVE_MENUBAR',REMOVE_MENUBAR,'USE_DRAGZOOM',USE_DRAGZOOM);
patch('faces',fv.faces,'vertices',fv.vertices,...
    'FaceColor',0.2*[1,1,1],'EdgeColor', 'none','FaceLighting','gouraud',...
    'AmbientStrength',0.5,'FaceAlpha',0.7);
plot_capsule(cap,'fig_idx',fig_idx,'cfc','r','cfa',0.2);
plot_T(pr2t('',''),'fig_idx',fig_idx,'subfig_idx',1);
plot_T(cap.T_offset,'fig_idx',fig_idx,'subfig_idx',2);
plot_title(sprintf('[%s] of UR10',file_name),'tfs',25);

%% Handle self-collision with Inverse Kinematics
ccc
% Configuration
robot_name  = 'thormang_rilab';
rseed       = 0;
SAVE_VID    = 0;
% Get robot with collision
rng(rseed); % fix random seed
urdf_path = sprintf('../../yet-another-robotics-toolbox/urdf/%s/%s_urdf.xml',robot_name,robot_name);
chain_robot = get_chain_from_urdf_with_caching(robot_name,...
    'RE',0,'urdf_path',urdf_path,'cache_folder','../cache');
collision_margin = max(chain_robot.sz.xyz_len)/100;
chain_robot.sc_checks = get_sc_checks(chain_robot,'collision_margin',collision_margin);
while 1 % loop until self-collision occurs
    chain_robot = update_chain_q(chain_robot,chain_robot.rev_joint_names,...
        90*D2R*randn(chain_robot.n_rev_joint,1));
    [SC,~] = check_sc(chain_robot,'collision_margin',0);
    if SC, break; end
end % while 1 % loop until self-collision occurs
% Loop to solve IK
vid_path = sprintf('vid/unit_test/sc_handle_%s_seed%02d.mp4',robot_name,rseed); HZ = 10;
tick = 0; vobj = init_vid_record(vid_path,'HZ',HZ,'SAVE_VID',SAVE_VID);
while SC % loop until collision-free
    tick = tick + 1;
    % Handle self-collision with IK (small update)
    len_offset = chain_robot.sz.xyz_len(3)/20;
    unit_dq_rad = 1*D2R;
    [chain_robot,ik_info,sc_link_pairs] = handle_sc_with_ik(chain_robot,...
        'len_offset',len_offset,'UNIT_DQ_HEURISTIC',1,'unit_dq_rad',unit_dq_rad,...
        'collision_margin',0);
    % Check self-collision again
    [SC,~,min_cap_dist] = check_sc(chain_robot,'collision_margin',0);
    % Animate robot and colliding capsules
    axis_info = '';
    ral = chain_robot.sz.xyz_len(3)/20;
    fig = plot_chain(chain_robot,'fig_idx',1,'subfig_idx',1,'fig_pos',[0.5,0.5,0.3,0.5],...
        'axis_info',axis_info,'PLOT_LINK',0,'PLOT_BOX',1,'bec','k','PLOT_ROTATE_AXIS',0,'ral',ral);
    plot_capsule('','RESET',1); % reset capsule
    plot_ik_targets('RESET',1); % reset IK targets
    for sc_idx = 1:size(sc_link_pairs,1) % plot colliding capsules
        [cap_i,T_i] = get_chain_capsule(chain_robot,sc_link_pairs(sc_idx,1));
        [cap_j,T_j] = get_chain_capsule(chain_robot,sc_link_pairs(sc_idx,2));
        plot_capsule(cap_i,'fig_idx',1,'subfig_idx',2*sc_idx-1,...
            'T',T_i,'cfc','r','cfa',0.1,'cec','k','cea',0.5);
        plot_capsule(cap_j,'fig_idx',1,'subfig_idx',2*sc_idx,...
            'T',T_j,'cfc','r','cfa',0.1,'cec','k','cea',0.5);
    end
    ik_plot_info = get_ik_plot_info_from_ik_info(ik_info); % IK information
    adl = chain_robot.sz.xyz_len(3)/5;
    plot_ik_targets('chain_robot',chain_robot,'ik_plot_info',ik_plot_info,'sr',0.02,...
        'PLOT_ARROW',1,'adl',adl,'adsw',adl/10,'adtw',adl/5);
    title_str = sprintf('[%d] SC:[%d](%.3f)',tick,SC,min_cap_dist);
    plot_title(title_str,'tfs',20);
    if SC == 0
        plot_capsule('','RESET',1); % reset capsule
        plot_ik_targets('RESET',1); % reset IK targets
    end
    drawnow; if ~ishandle(fig), break; end
    record_vid(vobj,'fig',fig);
end % while SC % loop until collision-free
end_vid_record(vobj);
fprintf('Done.\n');

%% 1D interpolation using three methods (nearest, spline, linear methods)
ccc
% Different interpolation methods
x = 0:pi/4:2*pi;
v = sin(x);
x_test = 0:pi/16:2*pi;
v_nn = interp1(x,v,x_test,'nearest');
v_spline = interp1(x,v,x_test,'spline');
v_linear = interp1(x,v,x_test,'linear');
% Plot
ca; figure(1); hold on;
plot(x,v,'ko','markersize',12,'linewidth',2);
hn = plot(x_test,v_nn,'k-','linewidth',2);
hl = plot(x_test,v_linear,'b-','linewidth',2);
hs = plot(x_test,v_spline,'r-','linewidth',2);
xlim([0 2*pi]);
legend([hl,hs],{'Linear','Spline'},'fontsize',15);
plot_title('Different Interpolation Methods');

%% Selective smoothing using nonstationary GRP
ccc
rng(3);

% Sample a path from a GP prior
n_test      = 500;
t_test      = cv(linspace(0,1,n_test));
hyp         = [0.07,1/10];
K           = kernel_levse(t_test,t_test,ones(n_test,1),ones(n_test,1),hyp);
z           = randn(n_test,1);
traj        = chol(K+1e-8*eye(n_test,n_test))'*z;

% Get mu, dmu, and ddmu
[mu_hat,dmu_hat,ddmu_hat] = get_grp_mu_dmu_ddmu(t_test,traj,...
    't_test',t_test,'hyp',[1,1/5],'meas_std',1e-4);

% Get time regions with acceleration greater than 100 and smooth
lev_smts        = [0.0,0.2,0.4,0.6,0.8,1.0];
traj_smts       = cell(1,length(lev_smts));
dmu_hat_smts    = cell(1,length(lev_smts));
ddmu_hat_smts   = cell(1,length(lev_smts));
for i_idx = 1:length(lev_smts)
    idxs        = find(abs(ddmu_hat)>100);
    t_in        = t_test;
    x_in        = traj;
    l_in        = 1.0*ones(size(t_test));
    lev_smt     = lev_smts(i_idx);
    l_in(idxs)  = lev_smt; % important tuning parameter
    t_out       = t_test;
    hyp_mu      = [1,1/1];
    hyp_var     = [1,1/1];
    sig2w       = 1e-12;
    grp1d_smt   = init_grp1d(t_in,x_in,l_in,t_out,hyp_mu,hyp_var,...
        'eps_ru','','meas_noise_std',sig2w,'APPLY_LEVERAGE_TO_MU',1);
    traj_smt    = grp1d_smt.mu_test;
    [mu_hat_smt,dmu_hat_smt,ddmu_hat_smt] = get_grp_mu_dmu_ddmu(t_test,traj_smt,...
        't_test',t_test,'hyp',[1,1/5],'meas_std',1e-4);
    % Append
    traj_smts{i_idx}     = traj_smt;
    dmu_hat_smts{i_idx}  = dmu_hat_smt;
    ddmu_hat_smts{i_idx} = ddmu_hat_smt;
end

% Plot the sampled and fitted trajectory
fig_idx = 1;
fig1 = figure(fig_idx); set_fig(fig1,'pos',[0.0,0.7,0.25,0.3],...
    'USE_DRAGZOOM',0,'ax_str','t','ay_str','x(t)','afs',17);
plot(t_test(idxs),traj(idxs),'o','linewidth',2,'markersize',8,...
    'color','b','markerfacecolor','b');
h_org = plot(t_test,traj,'k-','linewidth',2);
hm = plot(t_test,mu_hat,'r--','linewidth',2);
legend([h_org,hm],{'Original','Fit'},...
    'fontsize',12,'fontname','consolas');
axis([0,1,-inf,inf]); axis equal;
title_str = 'Original Trajectory';
plot_title(title_str,'fig_idx',fig_idx,'tfs',16);

% Plot velocity and acceleration
fig_idx = 2;
fig2 = figure(fig_idx); set_fig(fig2,'pos',[0.25,0.7,0.25,0.3],'AXIS_EQUAL',0,...
    'USE_DRAGZOOM',0,'ax_str','t','ay_str','$\dot{x}(t)$ ~and~ $\ddot{x}(t)$','afs',17);
colors = linspecer(4);
plot(t_test(idxs),ddmu_hat(idxs),'o','linewidth',2,'markersize',8,...
    'color','b','markerfacecolor','b');
hv = plot(t_test,dmu_hat,'-','linewidth',2,'color',colors(1,:));
ha = plot(t_test,ddmu_hat,'-','linewidth',2,'color',colors(2,:));
legend([hv,ha],{'Velocity','Acceleration'},...
    'fontsize',12,'fontname','consolas');
axis([0,1,-inf,inf]);
title_str = 'Velocity and Acceleration';
plot_title(title_str,'fig_idx',fig_idx,'tfs',16);

% Plot the sampled and smoothed trajectory
fig_idx = 3;
fig3 = figure(fig_idx); set_fig(fig3,'pos',[0.0,0.35,0.25,0.3],...
    'USE_DRAGZOOM',0,'ax_str','t','ay_str','x(t)','afs',17);
plot(t_test(idxs),traj(idxs),'o','linewidth',2,'markersize',8,...
    'color','b','markerfacecolor','b');
h_org = plot(t_test,traj,'k-','linewidth',2);
hm = plot(t_test,mu_hat,'r--','linewidth',2);
h_list = [h_org,hm]; strs = {'Original','Fit'};
colors = linspecer(length(lev_smts));
for i_idx = 1:length(lev_smts)
    lev_smt  = lev_smts(i_idx);
    traj_smt = traj_smts{i_idx};
    hs       = plot(t_test,traj_smt,':','linewidth',2,'color',colors(i_idx,:));
    h_list   = [h_list, hs];
    strs{i_idx+2} = sprintf('Smoothed (%.2f)',lev_smt);
end
legend(h_list,strs,'fontsize',12,'fontname','consolas');
axis([0,1,-inf,inf]); axis equal;
title_str = 'Smoothed Trajectory';
plot_title(title_str,'fig_idx',fig_idx,'tfs',16);

% Plot smoothed velocity and acceleration
fig_idx = 4;
fig4 = figure(fig_idx); set_fig(fig4,'pos',[0.25,0.35,0.25,0.3],'AXIS_EQUAL',0,...
    'USE_DRAGZOOM',0,'ax_str','t','ay_str','$\dot{x}(t)$ ~and~ $\ddot{x}(t)$','afs',17);
colors = linspecer(4);
plot(t_test(idxs),ddmu_hat(idxs),'o','linewidth',2,'markersize',8,...
    'color','b','markerfacecolor','b');
hv = plot(t_test,dmu_hat,'-','linewidth',2,'color',colors(1,:));
ha = plot(t_test,ddmu_hat,'-','linewidth',2,'color',colors(2,:));
hv_smt = plot(t_test,dmu_hat_smt,':','linewidth',2,'color',colors(3,:));
ha_smt = plot(t_test,ddmu_hat_smt,':','linewidth',2,'color',colors(4,:));
legend([hv,ha,hv_smt,ha_smt],...
    {'Velocity','Acceleration','Smoothed Velocity','Smoothed Acceleration'},...
    'fontsize',12,'fontname','consolas');
axis([0,1,-inf,inf]);
title_str = 'Velocity and Acceleration';
plot_title(title_str,'fig_idx',fig_idx,'tfs',16);

%% Compute the ZMP of IIWA7 motion
ccc
% Configuration
omega = 120*D2R; % joint velocity [rad/s]
robot_name = 'iiwa7';
urdf_path = sprintf('../../yet-another-robotics-toolbox/urdf/%s/%s_urdf.xml',...
    robot_name,robot_name);
chain = get_chain_from_urdf_with_caching(robot_name,'RE',0,'SKIP_CAPSULE',0,...
    'urdf_path',urdf_path,'cache_folder','../cache');
chain = update_chain_mass_inertia_com(chain,'RE',1,'density',300);
chain.dt = 1e-2; % set time diff
tick = 0;
com_traj = []; zmp_traj = []; iclk = clock; RESET = 1;
while true
    % Update
    tick = tick + 1;
    sec = tick*chain.dt;
    wallsec = etime(clock,iclk); % wall clock time
    q = 90*D2R*sin(omega*sec);
    chain = update_chain_q(chain,{'iiwa7_joint_1','iiwa7_joint_2','iiwa7_joint_3'},[q,q,q],...
        'RESET',RESET);
    RESET = 0;
    % Compute the center of mass (com)
    com = get_chain_com(chain);
    com_ground = get_com_ground(chain);
    % Compute the zero moment point (zmp)
    z_ground = com_ground(3);
    zmp_ground = get_zmp_ground(chain,com,z_ground);
    % Accumulate trajectory
    com_traj = cat(1,com_traj,rv(com_ground));
    zmp_traj = cat(1,zmp_traj,rv(zmp_ground));
    max_len = 100;
    if size(com_traj,1) > max_len
        com_traj = com_traj(end-max_len+1:end,:);
    end
    if size(zmp_traj,1) > max_len
        zmp_traj = zmp_traj(end-max_len+1:end,:);
    end
    % Animate
    % if mod(tick,10) == 0 % plot every 10 tick
    if wallsec < sec % simulation time sync with wall-clock
        fig_idx = 1;
        fig_pos = [0.0,0.5,0.3,0.45];
        fig = plot_chain(chain,'fig_idx',fig_idx,'fig_pos',fig_pos,...
            'axis_info',chain.axis_info,...
            'PLOT_LINK',0,...
            'PLOT_COM',1,'csr',0.02,...
            'PLOT_LINK_V',1,'lvar',0.1,'lvsw',0.02,'lvtw',0.05,...
            'PLOT_LINK_W',1,'lwar',0.05,'lwsw',0.02,'lwtw',0.05...
            );
        [~,~,h_com] = plot_T(p2t(com_ground),'fig_idx',fig_idx,'subfig_idx',1,...
            'PLOT_AXIS',0,'PLOT_SPHERE',1,'sr',0.025,'sfc','r','sfa',0.9);
        [~,~,h_zmp] = plot_T(p2t(zmp_ground),'fig_idx',fig_idx,'subfig_idx',2,...
            'PLOT_AXIS',0,'PLOT_SPHERE',1,'sr',0.025,'sfc','b','sfa',0.9);
        plot_traj(com_traj,'fig_idx',fig_idx,'subfig_idx',1,'tlc','r');
        plot_traj(zmp_traj,'fig_idx',fig_idx,'subfig_idx',2,'tlc','b');
        title_str = sprintf('[%.2f]sec',sec);
        plot_title(title_str,'fig_idx',fig_idx,'tfs',25,'interpreter','latex');
        plot_legend([h_com,h_zmp],{'COM','ZMP'},'fig_idx',fig_idx,'ll','southeast');
        drawnow;
        if ~ishandle(fig), break; end
    end
end

%% Handle external collision handling with capsules
ccc
% Configuration
SAVE_VID = 1;
% Robot
urdf_path = '../../yet-another-robotics-toolbox/urdf/ambidex_labs/ambidex_labs_urdf.xml';
cache_folder = '../cache';
chain_robot = get_ambidex_labs('urdf_path',urdf_path,'cache_folder',cache_folder);
chain_robot = update_chain_q(chain_robot,{'jointR2','jointR4','jointR5','jointL2'},...
    [-50,20,20,45]*D2R);
% Capsule
cap1 = get_capsule_shape('T_offset',pr2t(cv([0.4,-0.9,1.0]),''),...
    'radius',0.4,'height',0.3);
cap2 = get_capsule_shape('T_offset',pr2t(cv([0.4,0.9,1.0]),''),...
    'radius',0.5,'height',0.2);
caps = {cap1,cap2};
vid_path = sprintf('../vid/unit_test/ec_handle_%s.mp4',chain_robot.name); HZ = 10;
vobj = init_vid_record(vid_path,'HZ',HZ,'SAVE_VID',SAVE_VID);
for tick = 1:100
    
    % Handle external collision
    collision_margin = 0; % positive margin -> less collision
    collision_pairs = []; % (robot link index, capsule index)
    collision_dist = inf; % minimum signed distance
    EC = 0; % external collision
    for link_idx = 1:chain_robot.n_link % for all links
        link_i = chain_robot.link(link_idx);
        cap_i = link_i.capsule;
        if isempty(cap_i), continue; end
        joint_idx_i = link_i.joint_idx;
        if isempty(joint_idx_i), continue; end
        T_i = pr2t(chain_robot.joint(joint_idx_i).p,chain_robot.joint(joint_idx_i).R);
        cl_i = get_capsule_line(T_i,cap_i);
        for cap_idx = 1:length(caps) % for all external capsules
            cap_j = caps{cap_idx};
            cl_j = get_capsule_line(pr2t('',''),cap_j);
            line_dist = get_dist_lines(cl_i.p1,cl_i.p2,cl_j.p1,cl_j.p2);
            cap_dist = line_dist - cap_i.radius - cap_j.radius;
            if cap_dist < -collision_margin
                collision_pairs = cat(1,collision_pairs,[link_idx,cap_idx]);
                EC = 1;
            end
            if cap_dist < collision_dist
                collision_dist = cap_dist;
            end
        end % for cap_idx = 1:length(caps) % for all external capsules
    end % for link_idx = 1:chain_robot.n_link % for all links
    
    % Define IK info
    len_offset = 0.05;
    ik_info = init_ik_info(chain_robot,'joint_names_to_ctrl',chain_robot.rev_joint_names,...
        'ik_err_th',0.5,'dq_th',5*D2R);
    for collision_idx = 1:size(collision_pairs,1)
        % Robot link capsule (i)
        link_idx = collision_pairs(collision_idx,1); link_i = chain_robot.link(link_idx);
        cap_i = link_i.capsule; joint_idx_i = link_i.joint_idx;
        joint_name_i = chain_robot.joint(joint_idx_i).name;
        T_i = pr2t(chain_robot.joint(joint_idx_i).p,chain_robot.joint(joint_idx_i).R);
        % External capsule (j)
        cap_idx = collision_pairs(collision_idx,2);
        cap_j = caps{cap_idx}; T_j = pr2t('','');
        cl_i = get_capsule_line(T_i,cap_i); cl_j = get_capsule_line(T_j,cap_j);
        p_mid_i = 0.5*(cl_i.p1+cl_i.p2); p_mid_j = 0.5*(cl_j.p1+cl_j.p2);
        p_mid = 0.5*(p_mid_i+p_mid_j);
        uv_repulse_i = uv(p_mid_i-p_mid); uv_repulse_j = uv(p_mid_j-p_mid);
        
        % Add IK information to handle self-collsion for i-th joint
        if ~isempty(chain_robot.joint(joint_idx_i).childs)
            % If child joints exist, use both childs and parent joint for IK
            for child_idx = chain_robot.joint(joint_idx_i).childs
                joint_name_child = chain_robot.joint_names{child_idx};
                T_child = pr2t(chain_robot.joint(child_idx).p,'');
                ik_info = add_ik_info(ik_info,'joint_name',joint_name_child,...
                    'type','IK_P','weight',1,'coord',T_child + pr2t(uv_repulse_i*len_offset,''));
            end
            ik_info = add_ik_info(ik_info,'joint_name',joint_name_i,...
                'type','IK_P','weight',1,'coord',T_i + pr2t(uv_repulse_i*len_offset,''));
        else
            % Otherwise, add the parent joint
            ik_info = add_ik_info(ik_info,'joint_name',joint_name_i,...
                'type','IK_P','weight',1,'coord',T_i + pr2t(uv_repulse_i*len_offset,''));
        end
    end
    
    % Update robot with IK
    UNIT_DQ_HEURISTIC = 1; unit_dq_rad = 1*D2R;
    if ik_info.n_trgt > 0
        [dq,J_use,ik_err,det_J] = get_dq_from_ik_info(chain_robot,ik_info);
        q = get_q_chain(chain_robot,ik_info.joint_names_to_ctrl);
        % Simple hueristics
        if UNIT_DQ_HEURISTIC
            dq = dq / max(abs(dq));
            dq = trim_scale(dq,unit_dq_rad);
        end
        q = q + dq;
        chain_robot = update_chain_q(chain_robot,ik_info.joint_names_to_ctrl,q,'FK',1,'FV',0);
    end
    
    
    % Plot
    fig_idx = 1; view_info = [51,25]; axis_info = [-1,1,-1.5,1.5,0,2]; axes_info = [0.02,0,0.95,0.9];
    fig = plot_chain(chain_robot,'fig_idx',fig_idx,'subfig_idx',1,'fig_pos',[0.0,0.5,0.3,0.45],...
        'view_info',view_info,...
        'axis_info',axis_info,'AXIS_OFF',0,'axes_info',axes_info,'mfa',0.1,...
        'PLOT_CAPSULE',1,...
        'PLOT_LINK',1,'PLOT_BOX_ADDED',1,'bafa',0.3,'PLOT_ROTATE_AXIS',1,'DISREGARD_JOI_GUIDE',1);
    % Reset
    plot_capsule('','RESET',1); % reset capsule
    plot_ik_targets('RESET',1); % reset IK targets
    % Collision links
    for collision_idx = 1:size(collision_pairs,1)
        link_idx = collision_pairs(collision_idx,1);
        link_i = chain_robot.link(link_idx);
        cap_i = link_i.capsule;
        joint_idx_i = link_i.joint_idx;
        T_i = pr2t(chain_robot.joint(joint_idx_i).p,chain_robot.joint(joint_idx_i).R);
        cap_i = get_capsule_shape('T_offset',cap_i.T_offset,'radius',cap_i.radius+0.01,...
            'height',cap_i.height); % increase capsule radius
        plot_capsule(cap_i,'fig_idx',fig_idx,'subfig_idx',collision_idx,...
            'T',T_i,'cfc','b','cfa',0.4,'cec','none');
    end
    for cap_idx = 1:length(caps)
        plot_capsule(caps{cap_idx},'fig_idx',fig_idx,'subfig_idx',50+cap_idx,'cfc','r','cfa',0.2);
    end
    % Plot IK information
    ik_plot_info = get_ik_plot_info_from_ik_info(ik_info); % IK information
    adl = chain_robot.sz.xyz_len(3)/5;
    plot_ik_targets('chain_robot',chain_robot,'ik_plot_info',ik_plot_info,'sr',0.02,...
        'PLOT_ARROW',1,'adl',adl,'adsw',adl/10,'adtw',adl/5);
    if EC
        title_str = sprintf('[%d] External Collision (%.3f)',tick,collision_dist);
    else
        title_str = sprintf('[%d] Collision Free',tick);
    end
    plot_title(title_str,'tfs',20);
    drawnow; record_vid(vobj,'fig',fig);
    if (EC == 0), break; end
end
end_vid_record(vobj);
fprintf(2,'Done.\n');

%% Compute mu, d_mu, and dd_mu using Gaussian random path
ccc
% Reference data
f_ref = @(x)( cos(x) ); % reference function
t_max = 4*pi;
t_ref = linspace(0,t_max,1000)';
x_ref = f_ref(t_ref);
% Anchor dataset (training data)
n_anchor  = 100;
t_anchor  = linspace(0,t_max,n_anchor)';
noise_var = 1e-2;
x_anchor  = f_ref(t_anchor) + sqrt(noise_var)*randn(size(t_anchor));
l_anchor  = ones(n_anchor,1);
% Gaussian random path
n_test = 1000;
t_test = linspace(0,t_max,n_test)';
hyp = [1,3]; % [gain,len]
[k_test,dk_test,ddk_test] = kernel_levse(t_test,t_anchor,ones(n_test,1),l_anchor,hyp);
K_anchor = kernel_levse(t_anchor,t_anchor,l_anchor,l_anchor,hyp);
% Compute mu, d_mu, and dd_mu
meas_std = 1e-2; % expected noise
mu_test = k_test / (K_anchor+meas_std*eye(n_anchor,n_anchor)) * x_anchor;
dmu_test = dk_test / (K_anchor+meas_std*eye(n_anchor,n_anchor)) * x_anchor;
ddmu_test = ddk_test / (K_anchor+meas_std*eye(n_anchor,n_anchor)) * x_anchor;
% Plot mu, dmu, ddmu
set_fig(figure(1),'pos',[0.5,0.4,0.5,0.3],...
    'view_info','','axis_info',[0,4*pi,-1.5,1.5],'AXIS_EQUAL',0,'GRID_ON',1,...
    'REMOVE_MENUBAR',1,'USE_DRAGZOOM',0,'SET_CAMLIGHT',1,'SET_MATERIAL','METAL',...
    'SET_AXISLABEL',1,'ax_str','t','ay_str','f(t)','afs',18,'interpreter','latex');
h_ref    = plot(t_ref,x_ref,'k-','linewidth',3);
h_mu     = plot(t_test,mu_test,'r:','linewidth',3);
h_dmu    = plot(t_test,dmu_test,'-','linewidth',2,'Color','m');
h_ddmu   = plot(t_test,ddmu_test,'-','linewidth',2,'Color','c');
h_anchor = plot(t_anchor,x_anchor,'bo','linewidth',1,'markersize',11);
legend([h_ref,h_anchor,h_mu,h_dmu,h_ddmu],...
    {'$f(t)$','Observation','$\hat{\mu}(t)$',...
    '$\frac{d}{dt} \hat{\mu}(t)$','$\frac{d^2}{dt^2} \hat{\mu}(t)$'},...
    'fontsize',15,'interpreter','latex','location','NorthEastOutside');
plot_title('Estimated First and Second Derivatives using Gaussian Processes',...
    'interpreter','latex','tfs',20);

%% Find a rotational matrix that aligns two vectors using 'get_r_a_to_b'
ccc
p_a = rand(3,1);
p_b = rand(3,1);
R_a2b = get_r_a_to_b(p_a,p_b);
R_p_a = R_a2b*p_a; % <= this aligns with 'p_b'
% Plot
fig = set_fig(figure(1),'pos',[0.6,0.4,0.3,0.55],...
    'view_info',[80,26],'axis_info',1.1*[-1/10,+1,-1/10,+1,-1/10,+1],'AXIS_EQUAL',1,'GRID_ON',1,...
    'REMOVE_MENUBAR',1,'USE_DRAGZOOM',1,'SET_CAMLIGHT',1,'SET_MATERIAL','METAL',...
    'SET_AXISLABEL',1,'afs',18);
plot_line(cv([0,0,0]),p_a,'fig_idx',1,'subfig_idx',1,'lc','r');
plot_line(cv([0,0,0]),p_b,'fig_idx',1,'subfig_idx',2,'lc','b');
plot_line(cv([0,0,0]),R_p_a,'fig_idx',1,'subfig_idx',3,'lc','c','ls','--');
plot_T(pr2t('',''),'fig_idx',1,'subfig_idx',1,'alc','k'); % {W}
plot_T(pr2t(p_a,''),'fig_idx',1,'subfig_idx',2,'PLOT_AXIS',0,...
    'PLOT_SPHERE',1,'sr',0.03,'sfc','r','sfa',0.5,...
    'text_str','$~\mathbf{p}_a$','text_interp','latex','text_fs',21);
plot_T(pr2t(p_b,''),'fig_idx',1,'subfig_idx',3,'PLOT_AXIS',0,...
    'PLOT_SPHERE',1,'sr',0.03,'sfc','b','sfa',0.5,...
    'text_str','$~\mathbf{p}_b$','text_interp','latex','text_fs',21);
plot_T(pr2t(R_p_a,''),'fig_idx',1,'subfig_idx',4,'PLOT_AXIS',0,...
    'PLOT_SPHERE',1,'sr',0.03,'sfc','c','sfa',0.5,...
    'text_str','$~R\mathbf{p}_a$','text_interp','latex','text_fs',21);

%% Interpolate rotation matrices using GRP
%
% Here, we use three different methods
%  1. Optimal geodesic curve on SO(3)
%  2. GRP smoothing on so(3)
%  3. GPR smoothing on 6D continuous representation presented in [1]
%
% [1]. "On the Continuity of Rotation Representations in Neural Networks"
%
ccc
% Get two rotation matrices
p1 = cv([0,0,0]);
p2 = cv([0,3,0]);
R1 = rpy2r(360*rand(1,3)*D2R); % on SO(3), Lie Group
R2 = rpy2r(360*rand(1,3)*D2R);
w1 = r2w(R1); % on so(3)
w2 = r2w(R2);
g1 = norm(w1);
g2 = norm(w2);
b1 = so3_to_sixd(R1);
b2 = so3_to_sixd(R2);
% Interpolate on so(3), Lie Algebra
max_tick = 100;
t_in_ref = cv(linspace(0,1,max_tick));
w_traj = smooth_with_grp_multi_dim([w1';w2'],'t_test',t_in_ref,'hyp_mu',[1,1]);
% Interpolate using Rodrigues' formula
R1to2 = R1'*R2;
w1to2 = r2w(R1to2);
% Interpolate using 6D representation
b_traj = smooth_with_grp_multi_dim([b1';b2'],'t_test',t_in_ref,'hyp_mu',[1,1]);

% Loop
x_traj_so3 = []; y_traj_so3 = []; z_traj_so3 = [];
x_traj_rod = []; y_traj_rod = []; z_traj_rod = [];
x_traj_sixd = []; y_traj_sixd = []; z_traj_sixd = [];
w_rod_traj = zeros(max_tick,3);
for tick = 1:max_tick
    rate = tick/max_tick;
    p_t = (1-rate)*p1 + rate*p2;
    % Interpolate in so(3) GRP
    w_t = w_traj(tick,:);
    R_t_so3 = rodrigues(uv(w_t),norm(w_t));
    x_traj_so3 = [x_traj_so3; 1.0*rv(p_t+R_t_so3(:,1))];
    y_traj_so3 = [y_traj_so3; 1.0*rv(p_t+R_t_so3(:,2))];
    z_traj_so3 = [z_traj_so3; 1.0*rv(p_t+R_t_so3(:,3))];
    % Interpolate using Rodrigues
    R_t_rod = R1*rodrigues(w1to2/norm(w1to2),norm(w1to2)*(rate)); % interpolate
    x_traj_rod = [x_traj_rod; 1.0*rv(p_t+R_t_rod(:,1))];
    y_traj_rod = [y_traj_rod; 1.0*rv(p_t+R_t_rod(:,2))];
    z_traj_rod = [z_traj_rod; 1.0*rv(p_t+R_t_rod(:,3))];
    w_rod_traj(tick,:) = rv(r2w(R_t_rod));
    % Interpolate using 6D representation
    b_t = b_traj(tick,:);
    R_t_sixd = sixd_to_so3(b_t);
    x_traj_sixd = [x_traj_sixd; 1.0*rv(p_t+R_t_sixd(:,1))];
    y_traj_sixd = [y_traj_sixd; 1.0*rv(p_t+R_t_sixd(:,2))];
    z_traj_sixd = [z_traj_sixd; 1.0*rv(p_t+R_t_sixd(:,3))];
    % Animate
    if mod(tick,5) == 0
        fig_idx = 1;
        fig = set_fig(figure(fig_idx),'pos',[0.6,0.4,0.5,0.55],...
            'view_info',[80,26],'axis_info',[-1,+1,-1,+4,-1,+1],'AXIS_EQUAL',1,'GRID_ON',1,...
            'REMOVE_MENUBAR',1,'USE_DRAGZOOM',1,'SET_CAMLIGHT',1,'SET_MATERIAL','METAL',...
            'SET_AXISLABEL',1,'afs',18);
        plot_T(pr2t(p1,R1),'fig_idx',fig_idx,'subfig_idx',1,'alw',4,'PLOT_AXIS_TIP',1);
        plot_T(pr2t(p2,R2),'fig_idx',fig_idx,'subfig_idx',2,'alw',4,'PLOT_AXIS_TIP',1);
        plot_T(pr2t(p_t,R_t_so3),'fig_idx',fig_idx,'subfig_idx',3,'alw',4,'PLOT_AXIS_TIP',1);
        plot_T(pr2t(p_t,R_t_rod),'fig_idx',fig_idx,'subfig_idx',4,'alw',4,'PLOT_AXIS_TIP',1);
        plot_T(pr2t(p_t,R_t_sixd),'fig_idx',fig_idx,'subfig_idx',5,'alw',4,'PLOT_AXIS_TIP',1);
        h_so3 = plot_traj(x_traj_so3,'fig_idx',fig_idx,'subfig_idx',1,'tlc','r','tlw',1,'tls','--');
        plot_traj(y_traj_so3,'fig_idx',fig_idx,'subfig_idx',2,'tlc','g','tlw',1,'tls','--');
        plot_traj(z_traj_so3,'fig_idx',fig_idx,'subfig_idx',3,'tlc','b','tlw',1,'tls','--');
        h_rod = plot_traj(x_traj_rod,'fig_idx',fig_idx,'subfig_idx',4,'tlc','r','tlw',1,'tls','-');
        plot_traj(y_traj_rod,'fig_idx',fig_idx,'subfig_idx',5,'tlc','g','tlw',1,'tls','-');
        plot_traj(z_traj_rod,'fig_idx',fig_idx,'subfig_idx',6,'tlc','b','tlw',1,'tls','-');
        h_sixd = plot_traj(x_traj_sixd,'fig_idx',fig_idx,'subfig_idx',7,...
            'tlc','r','tlw',2,'tls',':');
        plot_traj(y_traj_sixd,'fig_idx',fig_idx,'subfig_idx',8,'tlc','g','tlw',2,'tls',':');
        plot_traj(z_traj_sixd,'fig_idx',fig_idx,'subfig_idx',9,'tlc','b','tlw',2,'tls',':');
        plot_title(sprintf('[%d/%d]',tick,max_tick),...
            'fig_idx',fig_idx,'tfs',20,'Interpreter','latex');
        plot_legend([h_rod,h_so3,h_sixd],...
            {'Optimal Interpolation','GRP on so(3)','GRP on 6D Representation'},...
            'lfs',15,'interpreter','latex','ll','SouthEast');
        drawnow; if ~ishandle(fig), break; end
    end
end

%% Gaussian Random Path Sampling in 1D
ccc
% Anchor points
t_anchor = [0,1/3,2/3,1.0]';
x_anchor = [0,1.0,-2.0,3.0]';
l_anchor = [1,1-0.1,1-0.5,1]';

% t_anchor = [0,0.5,1.0]';
% x_anchor = [0,1.0,3.0]';
% l_anchor = [1,0.3,1]';

eps_ru   = 0.01;

% Test data
n_test = 1000;
t_test = linspace(0,1,n_test)';

% Hyper parameters
hyp_mu  = [1,1.0]; % [gain,len]
hyp_var = [2,0.2]; % [gain,len]

% Initialize GRP
grp1d = init_grp1d(t_anchor,x_anchor,l_anchor,t_test,hyp_mu,hyp_var,'eps_ru',eps_ru,...
    'APPLY_LEVERAGE_TO_MU',1);

% Sample GRP paths
n_path = 100;
randomness = randn(n_test,n_path);
sampled_trajs = sqrt(grp1d.K_max)*grp1d.chol_K*randomness + grp1d.mu_test;

% Plot
plot_grp1d(grp1d,sampled_trajs);

%% Gaussian Random Path Sampling in 2D
ccc
% Initialize 2D GRP
t_anchor  = cv([0,1]);
xy_anchor = [0,-8; 0,8];
n_test    = 1000;
t_test    = linspace(0,1,n_test)';
hyp_mu    = [0.1,1.0]; % [gain,len]
hyp_var   = [1,0.2]; % [gain,len]
eps_ru    = 0.001;
grp_x = init_grp1d(t_anchor,cv(xy_anchor(:,1)),ones(size(t_anchor)),...
    t_test,hyp_mu,hyp_var,'eps_ru',eps_ru);
hyp_mu    = [0.1,1.0]; % [gain,len]
hyp_var   = [1,0.2]; % [gain,len]
eps_ru    = 0.01;
grp_y = init_grp1d(t_anchor,cv(xy_anchor(:,2)),ones(size(t_anchor)),...
    t_test,hyp_mu,hyp_var,'eps_ru',eps_ru);

% Sample GRP paths in 2D
n_path = 10;
randomness = randn(n_test,n_path);
sampled_x_trajs = sqrt(grp_x.K_max)*grp_x.chol_K*randomness + grp_x.mu_test;
sampled_y_trajs = sqrt(grp_y.K_max)*grp_y.chol_K*randomness + grp_y.mu_test;

% Plot sampled paths
figure(1); hold on; axis equal; grid on;
booth_colors = linspecer(n_path);
for i_idx = 1:n_path
    xy_traj = [sampled_x_trajs(:,i_idx),sampled_y_trajs(:,i_idx)];
    color = booth_colors(i_idx,:);
    plot(xy_traj(:,1),xy_traj(:,2),'-','color',color);
end

%% Determinantal Trajectory Process
ccc

t_min = 0; t_max = 1.0;
n_test = 1000;
t_test = cv(linspace(t_min,t_max,n_test));
l_test = ones(n_test,1);
hyp = [1,1/4]; sig2w = 1e-6;
K = kernel_levse(t_test,t_test,l_test,l_test,hyp);
K_chol = chol(K + sig2w*eye(n_test,n_test))';

% Sample from GP prior
n_traj = 200;
trajs = K_chol * randn(n_test,n_traj);

% Now, select subset using DPP with Hilbert Space norm-ish
idx_temp = round(linspace(1,n_test,50));
x = trajs(idx_temp,:)';
k = 5; % number of subset
C = kernel_levse(t_test(idx_temp),t_test(idx_temp),l_test(idx_temp),l_test(idx_temp),hyp);
C = 0.5*(C+C') + sig2w*eye(size(C));
[sub_idx_dpp,K_dpp] = get_sub_idx_ladpp(x,k,'C',C,'k_gain',20);

% Plot
fig_idx = 1; axis_info = [t_min,t_max,-inf,+inf]; axes_info = [0.08,0.12,0.88,0.8];
fig1 = set_fig(figure(fig_idx),'pos',[0.0,0.7,0.3,0.3],'AXIS_EQUAL',0,...
    'USE_DRAGZOOM',0,'axis_info',axis_info,'ax_str','t','ay_str','x(t)','afs',17,...
    'axes_info',axes_info);
% Plot all sampled trajectories
for i_idx = 1:n_traj
    traj = trajs(:,i_idx);
    h_sample = plot(t_test,traj,'-','color',0.5*[1,1,1],'linewidth',1/2);
end
% Plot random-selected ones
sub_idx_rand = randperm(n_traj,k);
for i_idx = 1:k
    traj = trajs(:,sub_idx_rand(i_idx));
    color = 'b';
    h_rand = plot(t_test,traj,'-','color',color,'linewidth',1);
end
% Plot DPP-selected ones
for i_idx = 1:k
    traj = trajs(:,sub_idx_dpp(i_idx));
    color = 'r';
    h_dpp = plot(t_test,traj,'-','color',color,'linewidth',2);
end
lfc = 'w'; lfs = 13; ll = 'best';
plot_legend([h_sample,h_rand,h_dpp],...
    {'GP-sampled Trajectories.',...
    'Randomly-selected Trajectories',...
    'DPP-selected Trajectories'},...
    'fig_idx',fig_idx,'lfc',lfc,'lfs',lfs,'interpreter','latex','ll',ll);
title_str = sprintf('Determinantal Trajectory Process');
plot_title(title_str,'fig_idx',fig_idx,'tfs',18);

%% Quality-Diversity Optimization 
ccc

rng(1);
% Set grid (domain)
xmin = 0; xmax = 10; nx = 200; ymin = 0; ymax = 10; ny = 200;
[xgrid,ygrid] = meshgrid(linspace(xmin,xmax,nx),linspace(ymin,ymax,ny));
xygrid = [xgrid(:),ygrid(:)];

% Define the cost function using GMM
mu = [7.5,2.5;2.5,7.5;7.5,7.5];
sigma = cat(3,[1.0,0.0;0.0,1.0],[1.0,0.0;0.0,1.0],[1.0,0.0;0.0,1.0]);
probs = [1,1,2]; probs = probs/sum(probs);
gm = gmdistribution(mu,sigma,probs);
f = @(x)(10*(gm.pdf(x))); % this will be our score function

% Random sampling first
n_pool = 3000;
x_pool = rand(n_pool,2)*diag([xmax,ymax]) + [xmin,ymin];
f_pool = f(x_pool);
gain = 1; len = 1;
K_pool = kernel_levse(x_pool,x_pool,ones(n_pool,1),ones(n_pool,1),[gain,len]);

% LA-DPP
k = 20; % number of selections
idxs_ladpp = zeros(k,1); remain_idxs = (1:n_pool)'; idx_ladpp = [];
[~,sel_idx] = max(f_pool); % first index to be the one with the maximum score
remain_idxs(remain_idxs==sel_idx) = [];
idxs_ladpp(1) = sel_idx; % append
for k_idx = 2:k % select (k-1) ones
    n_remain = length(remain_idxs);
    score_values = zeros(n_remain,1);
    similarity_values = zeros(n_remain,1);
    for r_idx = 1:n_remain % for the remaining ones
        remain_idx = remain_idxs(r_idx);
        x_prime = x_pool(remain_idx,:);
        score_r = f(x_prime);
        idxs_ladpp_curr = idxs_ladpp(1:(k_idx-1));
        sim_rs = cv(K_pool(remain_idx,idxs_ladpp_curr).^2) ./ ...
            (sigmoid(whitening(f_pool(idxs_ladpp_curr))));
        sim_r = mean(sim_rs);
        % Append score and similarity values
        score_values(r_idx) = score_r;
        similarity_values(r_idx) = sim_r; 
    end
    acqusition_values = whitening(score_values) - whitening(similarity_values);
    [~,max_idx] = max(acqusition_values); % select the one with the maximum acquisition value
    sel_idx = remain_idxs(max_idx); 
    remain_idxs(remain_idxs==sel_idx) = [];
    idxs_ladpp(k_idx) = sel_idx; % append
end

% Figure 1: score map
fig_idx = 1; axes_info = [0.03,0.08,0.93,0.85]; view_info = [88,16];
set_fig(figure(fig_idx),'pos',[0.0,0.6,0.3,0.45],'axes_info',axes_info);
fval = f(xygrid); fmtx = reshape(fval,ny,nx); % reshape
h = pcolor(xgrid,ygrid,fmtx); colormap(linspecer);
cb = colorbar; cb.FontSize = 12; cb.FontName = 'consolas';
set(h,'EdgeColor','none','FaceAlpha',0.98); set(gcf,'color','w');
axis([xmin,xmax,ymin,ymax]); axis('on', 'image'); axis on; dragzoom;
plot_title('Score Map','fig_idx',fig_idx,'tfs',20);

% Figure 2: DPP sample
fig_idx = 2; axes_info = [0.03,0.08,0.93,0.85];
set_fig(figure(fig_idx),'pos',[0.3,0.6,0.3,0.45],'axes_info',axes_info);
h = pcolor(xgrid,ygrid,fmtx); colormap(linspecer);
cb = colorbar; cb.FontSize = 12; cb.FontName = 'consolas';
set(h,'EdgeColor','none','FaceAlpha',0.4); set(gcf,'color','w');
plot(x_pool(:,1),x_pool(:,2),'.','color','k','markersize',1/2);
plot(x_pool(idxs_ladpp,1),x_pool(idxs_ladpp,2),'o',...
    'color','r','markersize',10,'linewidth',2);
axis([xmin,xmax,ymin,ymax]); axis('on', 'image'); axis on; dragzoom;
plot_title('Samples from k-DPP','fig_idx',fig_idx,'tfs',20);
drawnow;

%%


