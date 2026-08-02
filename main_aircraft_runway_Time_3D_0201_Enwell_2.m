%----------------时域：六自由度飞机-三维跑道（高效优化完整版·已修复bug）-------------------
close all; clc; clear;
timing = struct();
%----------------导入跑道结构计算结果----------------
Run_ans=load('/share/home/u22167/data/HouTianxin/Energy0125/Run_result_50_30_sub_7_2.mat');
result_filename = 'Air_Run_result_120_50_30_sub7_C0w0_Enwell_2.mat'; % 结果保存文件名
length_runway = Run_ans.length_runway;    
width_runway = Run_ans.width_runway;
thickness_runway = Run_ans.thickness_runway; 
xx = Run_ans.xx; yy = Run_ans.yy; zz= Run_ans.zz;

FEM_ff = Run_ans.FEM_ff;
FEM_gcoord = Run_ans.FEM_gcoord;  FEM_nodes = Run_ans.FEM_nodes;
Mr = Run_ans.Mr;  Cr = Run_ans.Cr;  Kr = Run_ans.Kr;  

% 输入计算参数
ve = 120/3.6;               % 滑跑速度(m/s)
dt = 0.0025;                % 时间步长(s)
z0 = 25;                    % 前轮初始位置(m)
x0 = width_runway/2;        % 轮距中心x坐标
max_error = 1e-8;           % 收敛误差
Le = length_runway - z0 - 5;  % 计算长度(m)
tm = Le/ve;                 % 总计算时间(s)
nt = floor(tm/dt);          % 总计算步数

% ========== 关键优化2：向量化不平整度计算 ==========
fprintf('预计算不平整度数据...\n');
shice = load('roughness_C0w0.mat');
Z3D = shice.roughness_data;
if size(Z3D, 2) > 1
    Z3D = Z3D(:, 1);
end
Z2D = Z3D';
x_Z2D = 0 : 0.25 : 0.25*(length(Z2D)-1);
yr_ce = csape(x_Z2D, Z2D);

% 定义飞机参数
rho_a = 1.293;  s_w = 62.29*2;  c_l = 0.61;
Fl = 0.5*rho_a*s_w*c_l*ve^2;  % 升力计算

ms = 78472;   Ix = 3394953;  Iy = 1866711;  mf = 578;  mr = 1150;  ml = 1150;
ksf = 58027;  ksr = 996531;  ksl = 996531;  ktf = 1966300; ktr = 2812900; ktl = 2812900;
csf = 131500; csr = 572500;  csl = 572500;  ctf = 4066;  ctr =4066;  ctl =4066;
lf = 14.6;  lm = 1;  br = 2.8;  bl = 2.8;  g = 9.81;  lc = lf + lm;

% 质量、刚度、阻尼矩阵
Ma = diag([ms, Ix, Iy, mf, mr, ml]); 
Ka = [ksf+ksl+ksr,           ksf*lf-ksl*lm-ksr*lm,        -ksl*bl+ksr*br,      -ksf,     -ksr,    -ksl;
      ksf*lf-ksl*lm-ksr*lm,  ksf*lf^2+ksl*lm^2+ksr*lm^2,  ksl*bl*lm-ksr*br*lm, -ksf*lf,  ksr*lm,  ksl*lm;   
      -ksl*bl+ksr*br,        ksl*lm*bl-ksr*lm*br,         ksl*bl^2+ksr*br^2,     0,      -ksr*br, ksl*bl;
      -ksf,                  -ksf*lf,                     0,                   ksf+ktf,  0,       0;
      -ksr,                  ksr*lm,                      -ksr*br,             0,        ksr+ktr, 0;
      -ksl,                  ksl*lm,                      ksl*bl,              0,        0,       ksl+ktl]; 

Ca = [csf+csl+csr,           csf*lf-csl*lm-csr*lm,        -csl*bl+csr*br,      -csf,     -csr,    -csl;
      csf*lf-csl*lm-csr*lm,  csf*lf^2+csl*lm^2+csr*lm^2,  csl*bl*lm-csr*br*lm, -csf*lf,  csr*lm,  csl*lm;
      -csl*bl+csr*br,        csl*lm*bl-csr*lm*br,         csl*bl^2+csr*br^2,     0,      -csr*br, csl*bl;
      -csf,                  -csf*lf,                     0,                   csf+ctf,  0,       0;
      -csr,                  csr*lm,                      -csr*br,             0,        csr+ctr, 0;
      -csl,                  csl*lm,                      csl*bl,              0,        0,       csl+ctl]; 

% 轮子起始位置与自由度参数
A_a = [x0, thickness_runway, z0; x0-br, thickness_runway, z0-lc; x0+bl, thickness_runway, z0-lc];
N_W = size(A_a,1);  N_air_dof = size(Ma,1);

% 网格参数
FEM_no_dof = 3;  FEM_no_nodeE1 = 8;
Num_element_x = Run_ans.Num_element_x;  
Num_element_y = Run_ans.Num_element_y;  
Num_element_z = Run_ans.Num_element_z;
no_node_x = Num_element_x+1;  
no_node_y = Num_element_y+1;  
no_node_z = Num_element_z+1;
FEM_no_nodeSys = no_node_x*no_node_y*no_node_z;  
FEM_Sys_dof = FEM_no_nodeSys*FEM_no_dof;
FEM_no_element = size(FEM_nodes,1);

% ========== 关键优化3：预计算和向量化节点数据 ==========
fprintf('预计算节点数据...\n');

% 预计算顶面节点与轮迹节点
top_nodes = find(abs(FEM_gcoord(:,2)-thickness_runway) < 1e-3);
w_node = zeros(N_W,1);  w_ydof = zeros(N_W,1);
for k = 1:N_W
    node_mask = abs(FEM_gcoord(:,1)-A_a(k,1)) < 1e-2 & ...
                abs(FEM_gcoord(:,2)-A_a(k,2)) < 1e-2 & ...
                abs(FEM_gcoord(:,3)-A_a(k,3)) < 1e-2;
    w_node(k) = find(node_mask,1);  
    w_ydof(k) = (w_node(k)-1)*FEM_no_dof+2;
end

% 预计算运动路径节点 - 向量化版本
w_move_node = zeros(no_node_z, N_W);
node_coords_x = FEM_gcoord(:,1);
node_coords_y = FEM_gcoord(:,2);
node_coords_z = FEM_gcoord(:,3);  % 预存储z坐标
for k = 1:N_W
    mask = abs(node_coords_x - A_a(k,1)) < 1e-2 & abs(node_coords_y - A_a(k,2)) < 1e-2;
    idx = find(mask);
    if length(idx) > no_node_z
        idx = idx(1:no_node_z);
    end
    w_move_node(1:length(idx), k) = idx;
end

% ========== 关键优化4：内存预分配 ==========
fprintf('预分配内存...\n');

% 使用更紧凑的数据类型
Ur = zeros(FEM_Sys_dof, nt+1 );  
Ua = zeros(N_air_dof, nt+1);
Vr = zeros(FEM_Sys_dof, nt+1);  
Va = zeros(N_air_dof, nt+1);
% ==================== 修复1：删除错误转置' ====================
Ar = zeros(FEM_Sys_dof, nt+1);  
Aa = zeros(N_air_dof, nt+1);
ur_load = zeros(N_W, nt+1);  
vr_load = zeros(N_W, nt+1);
% ==================== 修复2：预分配Ecw_each，避免未定义报错 ====================
Ecw_each = zeros(6, nt+1);

% 预分配迭代变量
Ur_iter = zeros(FEM_Sys_dof, 1);
Vr_iter = zeros(FEM_Sys_dof, 1);
fr = zeros(N_W, 1);


%================ 初始条件计算  ===========================================
fprintf('计算初始静平衡...\n');
tic;
% 强制矩阵对称 (修复数值误差导致的非对称报错)
if norm(Kr - Kr', 'inf') > 1e-5
    fprintf('  - 修复刚度矩阵对称性误差...\n');
    Kr = (Kr + Kr') / 2; 
    Mr = (Mr + Mr') / 2;
    Cr = (Cr + Cr') / 2;
end

Ua(:,1) = Ka \ [Fl-ms*g; 0; 0; -mf*g; -mr*g; -ml*g];
fr0_sparse = FEM_ff; 
fr0_sparse = fr0_sparse + sparse(w_ydof, 1, [ktf*Ua(4,1); ktr*Ua(5,1); ktl*Ua(6,1)], FEM_Sys_dof, 1);

fprintf('  - 计算初始预处理器 (ichol)...\n');
ichol_opts = struct('type', 'ict', 'droptol', 1e-2, 'diagcomp', 0.1);        % droptol=1e-2 是权衡速度和内存的最佳值
try
    L_init = ichol(Kr, ichol_opts); 
catch
    warning('内存紧张，降低预处理器精度...');          % 如果内存还不够，增大 droptol
    ichol_opts.droptol = 5e-2;
    L_init = ichol(Kr, ichol_opts);
end
fprintf('  - 开始初始迭代求解...\n');
[Ur(:,1), flag, relres, iter] = pcg(Kr, fr0_sparse, 1e-5, 2000, L_init, L_init');
if flag == 0
    fprintf('  - 初始平衡完成: 迭代 %d 步, 残差 %.2e\n', iter, relres);
else
    fprintf('  - 警告: 初始 PCG 未完全收敛 (Flag=%d, Res=%.2e)\n', flag, relres);
end
clear L_init;           % 立即释放，腾出内存给主循环
timing.part1 = toc;
fprintf('1.初始条件计算完成，用时: %.2f 秒\n', timing.part1);


% Newmark-β参数
gama = 0.5;  beta = 0.25;
a0 = 1/(beta*dt^2);  a1 = gama/(beta*dt);  a2 = 1/(beta*dt);
a3 = 1/(2*beta)-1;  a4 = gama/beta-1;  a5 = dt/2*(gama/beta-2);
a6 = dt*(1-gama);  a7 = gama*dt;

%================ 飞机、跑道矩阵分解  ======================================
fprintf('准备动力学迭代求解器...\n');
Ke_a = Ka + a0*Ma + a1*Ca;
KKe_a = decomposition(Ke_a, 'lu');  

Ke_r = Kr + a0*Mr + a1*Cr;
if norm(Ke_r - Ke_r', 'inf') > 1e-5
    fprintf('  - 修复等效刚度矩阵对称性误差...\n');
    Ke_r = (Ke_r + Ke_r') / 2; 
end
fprintf('  - 计算主循环预处理器 (ichol)... \n');
L_dyn = ichol(Ke_r, ichol_opts); 
L_dyn_T = L_dyn';  % 预转置

% 向量化预计算所有时间步的不平整度
T = dt*(0:nt);
z1 = z0 + T*ve;  z2 = z1 - lc;  z3 = z2;
yr1 = ppval(yr_ce, z1);  yr2 = ppval(yr_ce, z2);  yr3 = ppval(yr_ce, z3);
dyr1 = gradient(yr1)./gradient(T);  
dyr2 = gradient(yr2)./gradient(T);  
dyr3 = gradient(yr3)./gradient(T);


% ==========  主时间循环 ===================================================
fprintf('开始主循环计算...\n'); 
tic;
backNode_dof_base = (0:3:(FEM_Sys_dof-1))';  % 每个节点的基索引：预计算DOF偏移

for j = 2:nt+1    
    % 提取当前时刻参数
    z_j = [z1(j); z2(j); z3(j)];
    yr_j = [yr1(j); yr2(j); yr3(j)];
    dyr_j = [dyr1(j); dyr2(j); dyr3(j)];
    
    % 快速节点查找
    [backNode, frontNode, N4_lt, N8_rt] = fast_find_nodes(z_j, w_move_node, node_coords_z, N_W);
    
    % Newmark不变量 - 向量化计算
    D_a = Ma*(a0*Ua(:,j-1) + a2*Va(:,j-1) + a3*Aa(:,j-1)) + ...
          Ca*(a1*Ua(:,j-1) + a4*Va(:,j-1) + a5*Aa(:,j-1));
    pred_vec_r = a1*Ur(:,j-1) + a4*Vr(:,j-1) + a5*Ar(:,j-1);
    D_r = Mr*(a0*Ur(:,j-1) + a2*Vr(:,j-1) + a3*Ar(:,j-1)) + Cr*pred_vec_r;
    
    % 迭代求解
    Ur_iter(:) = Ur(:,j-1);  
    Vr_iter(:) = Vr(:,j-1);
    back_dof = backNode_dof_base(backNode) + 2;
    front_dof = backNode_dof_base(frontNode) + 2; 
    
    error = 1;  
    no_iterat = 0;  
    max_iter = 50;    
    while error > max_error && no_iterat < max_iter
        % 轮下位移
        Ur_load = N4_lt .* Ur_iter(back_dof) + N8_rt .* Ur_iter(front_dof);
        Vr_load = N4_lt .* Vr_iter(back_dof) + N8_rt .* Vr_iter(front_dof);
        
        % 飞机荷载与响应 - 使用矩阵分解求解
        fa = [Fl-ms*g; 0; 0;
              ktf*(Ur_load(1)+yr_j(1)) + ctf*(Vr_load(1)+dyr_j(1)) - mf*g;
              ktr*(Ur_load(2)+yr_j(2)) + ctr*(Vr_load(2)+dyr_j(2)) - mr*g;
              ktl*(Ur_load(3)+yr_j(3)) + ctl*(Vr_load(3)+dyr_j(3)) - ml*g];
        Ua1 = KKe_a \ (fa + D_a);  
        Aa1 = a0*(Ua1-Ua(:,j-1)) - a2*Va(:,j-1) - a3*Aa(:,j-1);
        Va1 = Va(:,j-1) + a6*Aa(:,j-1) + a7*Aa1;
        
        % 相互作用力
        fr(1) = ktf*(Ua1(4)-Ur_load(1)-yr_j(1)) + ctf*(Va1(4)-Vr_load(1)-dyr_j(1));
        fr(2) = ktr*(Ua1(5)-Ur_load(2)-yr_j(2)) + ctr*(Va1(5)-Vr_load(2)-dyr_j(2));
        fr(3) = ktl*(Ua1(6)-Ur_load(3)-yr_j(3)) + ctl*(Va1(6)-Vr_load(3)-dyr_j(3));
        
        % 构造右端项 RHS，避免 full转换
        indices = [back_dof; front_dof];
        values = [N4_lt .* fr; N8_rt .* fr];
        RHS_r = FEM_ff + sparse(indices, 1, values, FEM_Sys_dof, 1) + D_r;
        
        % PCG 求解跑道位移：最大迭代500、热启动 x0 = Ur(:,j-1)、容差: 1e-6
        [Ur1, flag] = pcg(Ke_r, RHS_r, 1e-5, 500, L_dyn, L_dyn_T, Ur(:,j-1));        
        Ar1 = a0*(Ur1-Ur(:,j-1)) - a2*Vr(:,j-1) - a3*Ar(:,j-1);
        Vr1 = Vr(:,j-1) + a6*Ar(:,j-1) + a7*Ar1;        
      
        % 收敛判断
        error = norm(Ur1 - Ur_iter) / norm(Ur1);
        Ur_iter = Ur1;
        Vr_iter = Vr1;
        no_iterat = no_iterat + 1;
    end
    % 存储当前时刻结果
    Ua(:,j) = Ua1;  Va(:,j) = Va1;  Aa(:,j) = Aa1;
    Ur(:,j) = Ur1;  Vr(:,j) = Vr1;  Ar(:,j) = Ar1;
    ur_load(:,j) = Ur_load;  vr_load(:,j) = Vr_load;
    fprintf('第 %d 时间点计算结束，逐步时间计算累计耗时: %.4f 秒\n', j, toc); 
end
timing.part2 = toc;
fprintf('2.时间逐步计算完成，用时: %.2f 秒\n', timing.part2);



% ========== 关键优化9：高效能量计算 ==========
fprintf('开始能量计算...\n');

%---------------------（1）优化后跑道能量计算方案2--------------------------
fprintf('计算跑道阻尼耗散能量 (向量化)...\n');
tic;
Ecr_total = 0;
block_size = 500; % 分块防止内存溢出
% ==================== 修复3：循环上限改为nt-1，防止range+1越界 ====================
for batch_start = 1:block_size:nt-1
    idx_end = min(batch_start + block_size - 1, nt-1);        % 同步修改上限
    range = batch_start:idx_end;
    
    % 1. 计算这一块的全局阻尼力 F = C * V
    F_damp_block = Cr * Vr(:, range);
    
    % 2. 计算这一块的位移增量 dx = U(t+1) - U(t)
    dX_block = Ur(:, range+1) - Ur(:, range);
    
    % 3. 全局做功累加 sum(F .* dx)
    Ecr_total = Ecr_total + sum(sum(F_damp_block .* dX_block));
end
timing.part3 = toc;
fprintf('3.跑道能量计算完成，用时: %.2f 秒\n', timing.part3);
fprintf('跑道阻尼耗散的总能量为: %.4e J\n', Ecr_total);

% （2）飞机阻尼耗散总能量 - 向量化版本
fprintf('计算飞机阻尼耗散能量...\n');

% 向量化计算相对位移
u_s1 = Ua(1,:) + lf*Ua(2,:) - Ua(4,:);
u_s2 = Ua(1,:) - lm*Ua(2,:) + br*Ua(3,:) - Ua(5,:);
u_s3 = Ua(1,:) - lm*Ua(2,:) - bl*Ua(3,:) - Ua(6,:);
u_t1 = Ua(4,:) - yr1 - ur_load(1,:);
u_t2 = Ua(5,:) - yr2 - ur_load(2,:);
u_t3 = Ua(6,:) - yr3 - ur_load(3,:);

% 向量化位移增量计算
dx_s1 = [0, diff(u_s1)];
dx_s2 = [0, diff(u_s2)];
dx_s3 = [0, diff(u_s3)];
dx_t1 = [0, diff(u_t1)];
dx_t2 = [0, diff(u_t2)];
dx_t3 = [0, diff(u_t3)];

% 向量化相对速度计算
v_s1 = Va(1,:) + lf*Va(2,:) - Va(4,:);
v_s2 = Va(1,:) - lm*Va(2,:) + br*Va(3,:) - Va(5,:);
v_s3 = Va(1,:) - lm*Va(2,:) - bl*Va(3,:) - Va(6,:);
v_t1 = Va(4,:) - dyr1 - vr_load(1,:);
v_t2 = Va(5,:) - dyr2 - vr_load(2,:);
v_t3 = Va(6,:) - dyr3 - vr_load(3,:);

% 向量化阻尼力计算
Fca_s1 = csf * v_s1;
Fca_s2 = csr * v_s2;
Fca_s3 = csl * v_s3;
Fca_t1 = ctf * v_t1;
Fca_t2 = ctr * v_t2;
Fca_t3 = ctl * v_t3;

% 瞬时能耗功率
Ecw_each(1,:)=abs(Fca_s1.*v_s1);
Ecw_each(2,:)=abs(Fca_s2.*v_s2);
Ecw_each(3,:)=abs(Fca_s3.*v_s3);
Ecw_each(4,:)=abs(Fca_t1.*v_t1);
Ecw_each(5,:)=abs(Fca_t2.*v_t2);
Ecw_each(6,:)=abs(Fca_t3.*v_t3);
Ecw_time=sum(Ecw_each,1);           % 各步长飞机能耗功率

% 向量化能量累加
Eca_total = sum(abs(Fca_s1 .* dx_s1)) + ...
            sum(abs(Fca_s2 .* dx_s2)) + ...
            sum(abs(Fca_s3 .* dx_s3)) + ...
            sum(abs(Fca_t1 .* dx_t1)) + ...
            sum(abs(Fca_t2 .* dx_t2)) + ...
            sum(abs(Fca_t3 .* dx_t3));

fprintf('飞机阻尼耗散的总能量为: %.4e J\n', Eca_total);

% （3）系统总能耗
Ec_total = Ecr_total + Eca_total;
fprintf('系统阻尼耗散的总能量为: %.4e J\n', Ec_total);


% 计算时间分析↓--------------------------------------------------------------------------
timing.total = timing.part1 + timing.part2 + timing.part3;
timing.percentage = [
    (timing.part1 / timing.total) * 100, ...
    (timing.part2 / timing.total) * 100, ...
    (timing.part3 / timing.total) * 100      ];
fprintf('\n========== 性能分析报告 ==========\n');
fprintf('总计算时间: %.2f 秒 (约 %.2f 分钟)\n', timing.total, timing.total/60);
fprintf('各部分耗时详情:\n');
fprintf('  1. 第一部分: %8.2f 秒 (%5.1f%%)\n', timing.part1, timing.percentage(1));
fprintf('  2. 第二部分: %8.2f 秒 (%5.1f%%)\n', timing.part2, timing.percentage(2));
fprintf('  3. 第三部分: %8.2f 秒 (%5.1f%%)\n', timing.part3, timing.percentage(3));
fprintf('==================================\n');
total_time = timing.total; 
% 计算时间分析↑--------------------------------------------------------------------------


% 保存关键结果
save(result_filename, 'Ec_total','Ecr_total','Eca_total','Ecw_time','total_time',...
                       'Ua','Va','Aa','Ur','Vr','Ar',...
                       'T','z1','z2','ve','yr1','yr2',...
                       'xx','yy','zz');

fprintf('优化计算完成！总耗时：%.2f 秒\n', total_time);
fprintf('结果已保存至: %s\n', result_filename);

% ========== 所有函数定义放在文件末尾 ==========
function [backNode, frontNode, N4_lt, N8_rt] = fast_find_nodes(z_j, w_move_node, node_coords_z, N_W)
    backNode = zeros(N_W,1);  
    frontNode = zeros(N_W,1);  
    N4_lt = zeros(N_W,1);  
    N8_rt = zeros(N_W,1);
    
    for jj = 1:N_W
        column_nodes = w_move_node(w_move_node(:,jj) > 0, jj);
        if isempty(column_nodes)
            backNode(jj) = 1; frontNode(jj) = 1;
            N4_lt(jj) = 0.5; N8_rt(jj) = 0.5;
            continue;
        end
        
        % 向量化距离计算
        distances = abs(node_coords_z(column_nodes) - z_j(jj));
        [~, sortedIndices] = sort(distances);
        closest_nodes = column_nodes(sortedIndices(1:min(2, length(sortedIndices))));
        
        if length(closest_nodes) == 1
            backNode(jj) = closest_nodes(1);
            frontNode(jj) = closest_nodes(1);
            N4_lt(jj) = 1.0; 
            N8_rt(jj) = 0.0;
        else
            backNode(jj) = min(closest_nodes);
            frontNode(jj) = max(closest_nodes);
            z_back = node_coords_z(backNode(jj));
            z_front = node_coords_z(frontNode(jj));
            if abs(z_front - z_back) > 1e-10
                z_local = 2*(z_j(jj)-z_back)/(z_front-z_back)-1;
                N4_lt(jj) = 0.5*(1-z_local);
                N8_rt(jj) = 0.5*(1+z_local);
            else
                N4_lt(jj) = 0.5; 
                N8_rt(jj) = 0.5;
            end
        end
    end
end