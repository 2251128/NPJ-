clc;clear;close all
tic

% ===================== 全局参数配置 =====================
ve_list = 120/3.6; 
v_kmh_all = ve_list * 3.6;
n_ve = length(ve_list);

% 核心修改：道面不平整参数调整为2.0到4.0，步长0.2
roughness_levels = 2.20:0.01:2.80;  % 生成2.0,2.2,2.4,...,3.8,4.0
n_rough = length(roughness_levels);

% 结果存储（修改为三维结构，适配多道面参数）
results = struct( ...
    'roughness_level', roughness_levels, ...
    'v_kmh', repmat(v_kmh_all, n_rough, 1), ...
    'Eca_total', zeros(n_rough, n_ve), ...
    'CCI1', zeros(n_rough, n_ve), ...
    'x_static_f', zeros(n_rough, n_ve), ...
    'x_static_r', zeros(n_rough, n_ve), ...
    'lift_ratio', zeros(n_rough, n_ve), ...
    'damping_factor', zeros(n_rough, n_ve),...
    'ksf', zeros(n_rough, n_ve), ...
    'ksr', zeros(n_rough, n_ve) ...
);

% 固定参数
z0 =620;               % 前轮初始位置 (m)
lc = 15.6;              % 轮距相关参数 (m)
L =1120;                % 计算长度 (m)
dt = 8e-5;              % 时间步长 (s)
Le = L - z0;            % 有效计算长度 (m)
g = 9.81;               % 重力加速度 (m/s²)
max_nt = 1e5;           % 最大计算步数

% 飞机参数
rho_a = 1.298;                        % 空气密度：单位：kg/m^3
s_w = 62.29 * 2;                      % 机翼面积：单位：m^2
c_l = 1.769;                           % 升力系数：无单位（1.769   1.36   0.61）
ms =32217.8741; Ix = 3394953; Iy = 1866711;  % 质量和转动惯量
mf =361.2187; mr =1449.9961; ml =1449.9961;        % 轮组质量
% 几何参数
lf =13.4; lm =0.9; br =3.4; bl =3.4;  % 单位：m
% 阻尼参数（完全保留原逻辑，不做任何修改）
csf_base =36531.5017;     % 前缓冲器阻尼基础值
csr_base =159044.1801;     % 主缓冲器阻尼基础值
csl_base =159044.1801;
ctf_base =3950.0025;     % 前轮胎阻尼基础值
ctr_base =4016.0189;     % 主轮胎阻尼基础值
ctl_base =4016.0189;

% 刚度约束
ksf_min = 4.0e4;        % 前缓冲器刚度下限
ksf_max = 2.0e6;        % 前缓冲器刚度上限
ksr_min = 7.0e5;        % 主缓冲器刚度下限
ksr_max = 2.0e7;        % 主缓冲器刚度上限

% 轮胎参数
ktf_k0 =546250.0000;       % 前轮胎基准刚度 (N/m)
ktf_k1 =3.2e3;         % 前轮胎刚度系数
ktr_k0 =2585210.8928;       % 主轮胎基准刚度 (N/m)
ktr_k1 = 3.2e3;         % 主轮胎刚度系数
F_tire_min = 500;       % 轮胎最小载荷 (N)

% ===================== 外层循环：道面不平整参数 =====================
for r_idx = 1:n_rough
    rough_level = roughness_levels(r_idx);
    
    % 格式化逻辑：2.0~4.0步长0.2，统一显示1位小数
    rough_fmt_str = sprintf('%.2f', rough_level);
    rough_fmt_num = str2double(rough_fmt_str); % 用于日志输出的数值
    
    fprintf('============================================\n');
    fprintf('开始计算道面不平整参数：%s\n', rough_fmt_str); % 使用格式化字符串
    fprintf('============================================\n');
    
    % ===================== 道面数据加载 =====================
    try
        % 动态拼接文件名：使用格式化后的字符串（保证文件名匹配）
        filename = sprintf('roughnessresult_%s.mat', rough_fmt_str);
        shice = load(filename);
        if ~isfield(shice, 'roughness_data')
            error('mat文件中缺少关键字段：roughness_data');
        end
        Z3D = shice.roughness_data; 
        if size(Z3D, 2) > 1
            Z3D = Z3D(:, 1);  
        end
        Z2D = Z3D';                             
        Z2D_n = size(Z2D,2); 
        x_Z2D = 0 : 0.25 : 0.25 * (Z2D_n-1);    
        
        % 滤波
        b_low = fir1(8, 0.1, 'low');         
        Z2D_filtered = filtfilt(b_low, 1, Z2D);
        Z2D_filtered = Z2D_filtered * 1.0;    
        yr_ce = csape(x_Z2D, Z2D_filtered);    
        fprintf('✅ 道面数据加载成功（%s），数据长度：%d个点\n', rough_fmt_str, Z2D_n);
        
    catch ME
        warning(['道面数据加载失败（%s）：', ME.message, ' → 使用模拟道面'], rough_fmt_str);
        x_Z2D = 0:0.25:720;
        % 模拟道面振幅随不平整参数变化（使用原始数值计算，保证精度）
        amp = 0.005 * rough_level/4.0; % 按比例调整模拟道面振幅
        Z2D_filtered = amp*sin(2*pi*x_Z2D/80); 
        yr_ce = csape(x_Z2D, Z2D_filtered);
        fprintf('⚠️  模拟道面（%s）：80m波长，%.3fmm振幅\n', rough_fmt_str, amp*1000);
    end

    % ===================== 内层循环：速度遍历计算 =====================
    for i = 1:n_ve
        ve = ve_list(i);
        v_kmh = v_kmh_all(i);
        
        % 时间参数
        tm = Le / ve;
        nt = floor(tm / dt);
        nt = min(nt, max_nt);
        T = dt*(0:1:nt);
        
        % 1. 升力系数（核心修改：确保240km/h升力占比100%）
        % 更合理的升力系数曲线（滑跑阶段）
        if v_kmh < 120
            c_l = 0.61;    % 低速（未加速）
        elseif v_kmh <200
            c_l = 1.36;    % 中速（开始加速）
        else
            c_l =1.36;    % 高速（即将离地，接近真实上限）
        end
        Fl = 0.5 * rho_a * s_w * c_l * ve^2 ;
        lift_ratio = (Fl / (ms*g)) * 100;
        % 限制升力占比不超过100%（物理合理性）
        lift_ratio = min(lift_ratio, 100);
        Fl = min(Fl, ms*g);
        F_sprung_net = max(ms * g - Fl, 0);
        
        % 3. 载荷分配
        F_spring_f = F_sprung_net * (lm / (lf + lm));
        F_spring_r = F_sprung_net * (lf / (lf + lm)) / 2;
        
        % 4. 缓冲器压缩量
        coeffs_f = [0.0005, -0.0292, 16.055, -(7963.7 + F_spring_f)];
        roots_f = roots(coeffs_f);
        valid_roots_f = roots_f(imag(roots_f)==0 & roots_f >= 0 & roots_f < 1000);         
        x_static_f_mm = ifelse(isempty(valid_roots_f), 0, min(valid_roots_f));
        
        coeffs_r = [0.0109, -1.183, 186.94, -(97898 + F_spring_r)];
        roots_r = roots(coeffs_r);
        valid_roots_r = roots_r(imag(roots_r)==0 & roots_r >= 0 & roots_r < 1000);
        x_static_r_mm = ifelse(isempty(valid_roots_r), 0, min(valid_roots_r));
        
        % 5. 缓冲器刚度
        k_tan_f_Nmm = 3*0.0005*x_static_f_mm^2 - 2*0.0292*x_static_f_mm + 16.055;
        ksf = k_tan_f_Nmm * 1000;
        ksf = max(min(ksf, ksf_max), ksf_min);
        
        k_tan_r_Nmm = 3*0.0109*x_static_r_mm^2 - 2*1.183*x_static_r_mm + 186.94;
        ksr = k_tan_r_Nmm * 1000;
        ksr = max(min(ksr, ksr_max), ksr_min);
        ksl = ksr;
        csf = csf_base ;
        csr = csr_base ;
        csl = csl_base ;
        ctf = ctf_base ;
        ctr = ctr_base;
        ctl = ctl_base ;
        
        % 7. 轮胎刚度
        F_tire_f_load = max((ms + mf + mr + ml)*g - Fl, F_tire_min) * (lm/(lf+lm));
        F_tire_r_load = max((ms + mf + mr + ml)*g - Fl, F_tire_min) * (lf/(lf+lm))/2;
        ktf = ktf_k0 + ktf_k1 * sqrt(F_tire_f_load);
        ktr = ktr_k0 + ktr_k1 * sqrt(F_tire_r_load);
        ktl = ktr;
        
        % 输出：使用格式化后的道面参数字符串
        fprintf('------------------------------------------------\n');
        fprintf('道面参数：%s | 速度: %.1f km/h | 升力系数: %.3f | 升力占比: %.1f%%\n', ...
            rough_fmt_str, v_kmh, c_l, lift_ratio);
        fprintf('前起压缩: %.2f mm → 刚度: %.2e N/m\n', x_static_f_mm, ksf);
        fprintf('主起压缩: %.2f mm → 刚度: %.2e N/m\n', x_static_r_mm, ksr);
        
        % 8. 动力学矩阵
        msf = ms*lm/(lf+lm) + mf;   
        msr = 0.5*ms*lf/(lf+lm) + mr; 
        msl = 0.5*ms*lf/(lf+lm) + ml;
        Ma = diag([ms, Ix, Iy, mf, mr, ml]);
        
        Ka = [
            ksf+ksl+ksr,           ksf*lf-ksl*lm-ksr*lm,        -ksl*bl+ksr*br,      -ksf,     -ksr,    -ksl;
            ksf*lf-ksl*lm-ksr*lm,  ksf*lf^2+ksl*lm^2+ksr*lm^2,  ksl*bl*lm-ksr*br*lm, -ksf*lf,  ksr*lm,  ksl*lm;   
            -ksl*bl+ksr*br,        ksl*lm*bl-ksr*lm*br,         ksl*bl^2+ksr*br^2,     0,      -ksr*br, ksl*bl;
            -ksf,                  -ksf*lf,                     0,                   ksf+ktf,  0,       0;
            -ksr,                  ksr*lm,                      -ksr*br,             0,        ksr+ktr, 0;
            -ksl,                  ksl*lm,                      ksl*bl,              0,        0,       ksl+ktl;
        ];
        
        Ca = [
            csf+csl+csr,           csf*lf-csl*lm-csr*lm,        -csl*bl+csr*br,      -csf,     -csr,    -csl;
            csf*lf-csl*lm-csr*lm,  csf*lf^2+csl*lm^2+csr*lm^2,  csl*bl*lm-csr*br*lm, -csf*lf,  csr*lm,  csl*lm;
            -csl*bl+csr*br,        csl*lm*bl-csr*lm*br,         csl*bl^2+csr*br^2,     0,      -csr*br, csl*bl;
            -csf,                  -csf*lf,                     0,                   csf+ctf,  0,       0;
            -csr,                  csr*lm,                      -csr*br,             0,        csr+ctr, 0;
            -csl,                  csl*lm,                      csl*bl,              0,        0,       csl+ctl;
        ];
        N_air_dof = size(Ma,1);
        
        % 9. 道面激励
        z1 = z0 + T*ve;  z2 = z1 - lc;  z3 = z2;
        z1 = min(z1, max(x_Z2D));
        z2 = min(z2, max(x_Z2D));
        z3 = min(z3, max(x_Z2D));
        yr1 = ppval(yr_ce, z1);  
        yr2 = ppval(yr_ce, z2);  
        yr3 = ppval(yr_ce, z3);  
        dyr1 = gradient(yr1)./gradient(T);  
        dyr2 = gradient(yr2)./gradient(T);  
        dyr3 = gradient(yr3)./gradient(T);
        
        % 10. Newmark求解
        Ua = zeros(N_air_dof, nt+1);  
        Va = zeros(N_air_dof, nt+1);  
        Aa = zeros(N_air_dof, nt+1);  
        Fa = zeros(N_air_dof, nt+1);  
        
        Ua(:,1) = (Ka + 1e-6*eye(size(Ka))) \ [Fl - ms*g; 0; 0; -mf*g; -mr*g; -ml*g];
        
        gamma = 0.5; beta = 0.25;
        a0 = 1/(beta*dt^2);  a1 = gamma/(beta*dt);
        a2 = 1/(beta*dt);    a3 = 1/(2*beta) - 1;
        a4 = gamma/beta - 1; a5 = dt/2*(gamma/beta - 2);
        a6 = dt*(1 - gamma); a7 = gamma*dt;
        
        Ke = Ka + a0*Ma + a1*Ca;
        for j = 2:nt+1
            Fa(:,j) = [
                Fl - ms*g;
                0;
                0;
                ktf*yr1(j) + ctf*dyr1(j) - mf*g;
                ktr*yr2(j) + ctr*dyr2(j) - mr*g;
                ktl*yr3(j) + ctl*dyr3(j) - ml*g;
            ];
            Fe = Fa(:,j) + Ma*(a0*Ua(:,j-1)+a2*Va(:,j-1)+a3*Aa(:,j-1)) ...
                 + Ca*(a1*Ua(:,j-1)+a4*Va(:,j-1)+a5*Aa(:,j-1));
            Ua1 = Ke \ Fe;
            Aa1 = a0*(Ua1 - Ua(:,j-1)) - a2*Va(:,j-1) - a3*Aa(:,j-1);
            Va1 = Va(:,j-1) + a6*Aa(:,j-1) + a7*Aa1;
            Ua(:,j) = Ua1;
            Va(:,j) = Va1;
            Aa(:,j) = Aa1;
        end
        
        % 11. 能耗计算
        v_s1 = Va(1,:) + lf*Va(2,:) - Va(4,:);
        v_s2 = Va(1,:) - lm*Va(2,:) + br*Va(3,:) - Va(5,:);
        v_s3 = Va(1,:) - lm*Va(2,:) - bl*Va(3,:) - Va(6,:);
        v_t1 = Va(4,:) - dyr1;
        v_t2 = Va(5,:) - dyr2;
        v_t3 = Va(6,:) - dyr3;
        
        Eca_s1 = csf * sum(v_s1.^2) * dt;
        Eca_s2 = csr * sum(v_s2.^2) * dt;
        Eca_s3 = csl * sum(v_s3.^2) * dt;
        Eca_t1 = ctf * sum(v_t1.^2) * dt;
        Eca_t2 = ctr * sum(v_t2.^2) * dt;
        Eca_t3 = ctl * sum(v_t3.^2) * dt;
        
        Eca_total = Eca_s1 + Eca_s2 + Eca_s3 + Eca_t1 + Eca_t2 + Eca_t3;
        
        % 12. CCI指标
        v_s1_abs_sum = sum(abs(v_s1));
        v_s2_abs_sum = sum(abs(v_s2));
        v_s3_abs_sum = sum(abs(v_s3));
        v_t1_abs_sum = sum(abs(v_t1));
        v_t2_abs_sum = sum(abs(v_t2));
        v_t3_abs_sum = sum(abs(v_t3));
        v_sum = v_s1_abs_sum + v_s2_abs_sum + v_s3_abs_sum + v_t1_abs_sum + v_t2_abs_sum + v_t3_abs_sum;
        CCI1 = 100*v_sum/6 * dt * ve / Le;
        
        
        % 存储结果（修改为二维索引：道面参数idx + 速度idx）
        results.Eca_total(r_idx, i) = Eca_total;
        results.CCI1(r_idx, i) = CCI1;
        results.x_static_f(r_idx, i) = x_static_f_mm;
        results.x_static_r(r_idx, i) = x_static_r_mm;
        results.lift_ratio(r_idx, i) = lift_ratio;
        results.ksf(r_idx, i) = ksf;
        results.ksr(r_idx, i) = ksr;
        
        fprintf('总耗散能量: %.4e N·m (%.2f MJ) | CCI1: %.6f\n', ...
            Eca_total, Eca_total/1e6, CCI1);
    end
end

% ===================== 新增：查找各道面参数下能耗最大值对应的速度 =====================
% 初始化存储数组
max_Eca_value = zeros(1, n_rough);    % 各道面参数下的最大能耗值
max_Eca_velocity = zeros(1, n_rough); % 各道面参数下最大能耗对应的速度

for r_idx = 1:n_rough
    % 找到当前道面参数下能耗最大值的索引
    [max_val, max_idx] = max(results.Eca_total(r_idx, :));
    max_Eca_value(r_idx) = max_val;
    % 根据索引获取对应的速度
    max_Eca_velocity(r_idx) = results.v_kmh(r_idx, max_idx);
end
%% ===================== CCI-能耗关系拟合 =====================

% 提取数据
CCI_data = results.CCI1(:,1);

Energy_data = results.Eca_total(:,1)/1e6;   % 转换为 MJ


% 删除异常值
valid = isfinite(CCI_data) & isfinite(Energy_data);

CCI_data = CCI_data(valid);
Energy_data = Energy_data(valid);



%% 多项式拟合
% 二阶拟合
p = polyfit(CCI_data,Energy_data,2);


% 拟合曲线
CCI_fit = linspace(min(CCI_data),max(CCI_data),200);

Energy_fit = polyval(p,CCI_fit);



% R2计算
Energy_pred = polyval(p,CCI_data);

SS_res = sum((Energy_data-Energy_pred).^2);

SS_tot = sum((Energy_data-mean(Energy_data)).^2);

R2 = 1-SS_res/SS_tot;



%% 绘图

figure('Color','w');

scatter(CCI_data,Energy_data,...
    45,'filled');

hold on

plot(CCI_fit,Energy_fit,...
    'r-',...
    'LineWidth',2);


xlabel('CCI');

ylabel('Vibration dissipated energy (MJ)');

title(sprintf('CCI-energy relationship (R^2 = %.3f)',R2));


grid on;


legend('Simulation results',...
       'Polynomial fitting',...
       'Location','best');



fprintf('\n============================\n');
fprintf('CCI-能耗拟合结果\n');
fprintf('============================\n');

fprintf('拟合公式:\n');

fprintf('E = %.4f CCI^2 + %.4f CCI + %.4f\n',...
    p(1),p(2),p(3));


fprintf('R² = %.4f\n',R2);
% % ===================== 可视化 =====================
% figure('Name','多道面参数 升力占比 100% 结果','Position',[100,100,1200,900]);
% set(gcf, 'Color', 'white');

% 定义颜色数组，适配不同道面参数的曲线颜色
colors = lines(n_rough);

% % 子图1：速度-总耗散能量（多道面参数对比，标注最大值点）
% subplot(3,1,1);
% hold on;
% for r_idx = 1:n_rough
%     rough_level = roughness_levels(r_idx);
%     % 格式化逻辑：统一显示1位小数
%     rough_fmt_str = sprintf('%.1f', rough_level);
% 
%     [v_sort, sort_idx] = sort(results.v_kmh(r_idx, :));
%     Eca_sort = results.Eca_total(r_idx, sort_idx);
%     % 绘制能耗曲线
%     plot(v_sort, Eca_sort/1e6, 'o-', ...
%         'LineWidth',1.2, 'MarkerSize',3, 'Color',colors(r_idx,:), ...
%         'MarkerFaceColor',colors(r_idx,:), 'MarkerEdgeColor','white', ...
%         'DisplayName',sprintf('道面参数 %s', rough_fmt_str)); % 格式化标注
%     % 标注最大值点
%     plot(max_Eca_velocity(r_idx), max_Eca_value(r_idx)/1e6, 's', ...
%         'MarkerSize',6, 'MarkerFaceColor','red', 'MarkerEdgeColor','black', ...
%         'DisplayName',sprintf('%s最大值', rough_fmt_str));
% end
% xlabel('滑跑速度 (km/h)','FontSize',9);
% ylabel('总耗散能量 (MJ)','FontSize',9);
% title('不同道面参数下 滑跑速度 → 总耗散能量（标注最大值点）','FontSize',10, 'FontWeight','bold');
% legend('Location','best','FontSize',8);
% grid on; grid minor;
% set(gca, 'GridAlpha', 0.3);
% hold off;
% 
% % 子图2：压缩量（修改为4.0为例展示，适配新参数范围）
% subplot(3,1,2);
% r_idx_demo = find(roughness_levels == 4.0, 1);
% if isempty(r_idx_demo), r_idx_demo = n_rough; end
% % 格式化演示用道面参数
% rough_level_demo = roughness_levels(r_idx_demo);
% rough_fmt_demo = sprintf('%.1f', rough_level_demo);
% 
% plot(results.v_kmh(r_idx_demo, :), results.x_static_f(r_idx_demo, :), 'ro-', ...
%     'LineWidth',1.2, 'MarkerSize',5, 'DisplayName','前起落架');
% hold on;
% plot(results.v_kmh(r_idx_demo, :), results.x_static_r(r_idx_demo, :), 'bo-', ...
%     'LineWidth',1.2, 'MarkerSize',5, 'DisplayName','主起落架');
% xlabel('滑跑速度 (km/h)','FontSize',11);
% ylabel('静压缩量 (mm)','FontSize',11);
% title(sprintf('道面参数 %s 下的缓冲器静压缩量', rough_fmt_demo), 'FontSize',10, 'FontWeight','bold');
% legend('Location','best','FontSize',10);
% grid on; grid minor; set(gca, 'GridAlpha', 0.3);
% hold off;
% 
% % 子图3：刚度（以4.0为例展示）
% subplot(3,1,3);
% plot(results.v_kmh(r_idx_demo, :), results.ksf(r_idx_demo, :)/1e6, 'go-', ...
%     'LineWidth',1.2, 'MarkerSize',5, 'DisplayName','前缓冲器');
% hold on;
% plot(results.v_kmh(r_idx_demo, :), results.ksr(r_idx_demo, :)/1e6, 'Color',[0.5 0 0.5], 'LineStyle','-', 'Marker','^', ...
%     'LineWidth',1.2, 'MarkerSize',5, 'DisplayName','主缓冲器');
% xlabel('滑跑速度 (km/h)','FontSize',11);
% ylabel('缓冲器刚度 (MN/m)','FontSize',11);
% title(sprintf('道面参数 %s 下的缓冲器刚度', rough_fmt_demo), 'FontSize',10, 'FontWeight','bold');
% legend('Location','best','FontSize',10);
% grid on; grid minor; set(gca, 'GridAlpha', 0.3);
% hold off;
% 
% % ===================== 结果汇总 =====================
% total_time = toc;
% fprintf('============================================\n');
% fprintf('总计算完成！总耗时: %s\n', datestr(total_time/86400, 'HH:MM:SS'));
% fprintf('\n===== 各道面参数下能耗最大值汇总 =====\n');
% 
% % 输出各道面参数的关键结果（包含最大值对应速度）
% for r_idx = 1:n_rough
%     rough_level = roughness_levels(r_idx);
%     rough_fmt_str = sprintf('%.1f', rough_level);
% 
%     fprintf('===== 道面参数 %s 结果 =====\n', rough_fmt_str);
%     fprintf('1. 升力占比：%.1f%% (210km/h)\n', ...
%         results.lift_ratio(r_idx, 1)); % 仅单速度，无需区间
%     fprintf('2. 总耗散能量：峰值%.2f MJ (对应速度: %.1f km/h)\n', ...
%         max_Eca_value(r_idx)/1e6, max_Eca_velocity(r_idx));
% end
% 
% % 额外输出简洁的最大值对应速度汇总表
% fprintf('\n===== 能耗最大值对应速度汇总表 =====\n');
% fprintf('道面参数\t最大能耗(MJ)\t对应速度(km/h)\n');
% fprintf('----------------------------------------\n');
% for r_idx = 1:n_rough
%     rough_level = roughness_levels(r_idx);
%     rough_fmt_str = sprintf('%.1f', rough_level);
% 
%     fprintf('%s\t\t%.2f\t\t%.1f\n', ...
%         rough_fmt_str, max_Eca_value(r_idx)/1e6, max_Eca_velocity(r_idx));
% end

% ===================== 辅助函数 =====================
function val = ifelse(cond, val_true, val_false)
    if cond
        val = val_true;
    else
        val = val_false;
    end
end

function yrx = yrx(x,B0,L0,A0)              
A1 = A0 - L0;  A2 = A0 + L0;
L00 = 2*L0;
yrx = 0.*(x<A1) + (0.5*B0*(1-cos(2*pi*(x-A1)/L00))).*(x>=A1 & x<=A2) + 0.*(x>A2);
end

function dyrx = dyrx(x,B0,L0,A0,ve)          
A1 = A0 - L0;  A2 = A0 + L0;
L00 = 2*L0;
dyrx = 0.*(x<A1) + (B0*pi*ve/L00*sin(2*pi*(x-A1)/L00)).*(x>=A1 & x<=A2) + 0.*(x>A2);
end