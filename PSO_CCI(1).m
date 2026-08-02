% 基于粒子群算法优化标准机型CCI参数（考虑机型占比）
close all; clc; clear;

%%  1. 数据准备 
fprintf('====== 粒子群优化标准机型CCI参数 ======\n\n');

Energy_data = xlsread('Energy_Dissipation_Results.xlsx', 'Eca_total', 'B2:K21');
[n_aircraft, n_runways] = size(Energy_data);
fprintf('机型数量: %d\n', n_aircraft);
fprintf('跑道数量: %d\n', n_runways);


runway_number = 1:1:10;

param_range = xlsread('参数范围.xlsx', 'Sheet1', 'B2:Q3');
param_min = param_range(1, :);
param_max = param_range(2, :);

param_names = {'ms', 'mf','mr','ml', 'ksf', 'ksr','ksl', 'ktf', 'ktr','ktl', 'csf', 'csr', 'csl','ctf', 'ctr', 'ctl'};
n_params = length(param_names);

%%  1.1 读取机型占比 
[~, ~, raw_weight] = xlsread('机型占比.xlsx');
aircraft_weights = zeros(n_aircraft, 1);
for i = 1:n_aircraft
    aircraft_weights(i) = raw_weight{i+1, 2};
end

% 确保权重归一化
aircraft_weights = aircraft_weights / sum(aircraft_weights);

fprintf('\n机型占比（已归一化）:\n');
for i = 1:n_aircraft
    fprintf('  机型%2d: %.4f (%.1f%%)\n', i, aircraft_weights(i), aircraft_weights(i)*100);
end

%%  2. 能量归一化 
Energy_norm = zeros(n_aircraft, n_runways);
for i = 1:n_aircraft
    e_min = min(Energy_data(i, :));
    e_max = max(Energy_data(i, :));
    Energy_norm(i, :) = (Energy_data(i, :) - e_min) / (e_max - e_min + 1e-10);
end
fprintf('\n能量归一化完成\n');

%%  3. PSO参数设置 
n_particles = 100;
max_iter = 100;
w = 0.9;
w_min = 0.4;
c1 = 2.0;
c2 = 2.0;

%% 4. 初始化粒子群
positions = zeros(n_particles, n_params);
for p = 1:n_particles
    for i = 1:n_params
        positions(p, i) = param_min(i) + rand() * (param_max(i) - param_min(i));
    end
end

v_max = (param_max - param_min) * 0.2;
velocities = zeros(n_particles, n_params);
for p = 1:n_particles
    for i = 1:n_params
        velocities(p, i) = -v_max(i) + rand() * 2 * v_max(i);
    end
end

pbest_positions = positions;
pbest_fitness = inf(n_particles, 1);

gbest_position = zeros(1, n_params);
gbest_fitness = inf;
gbest_corr_weighted_mean = 0;
gbest_corr_min = 0;
gbest_corr_list = zeros(n_aircraft, 1);

fitness_history = zeros(max_iter, 1);
corr_weighted_mean_history = zeros(max_iter, 1);
corr_min_history = zeros(max_iter, 1);

%%  5. PSO主循环
fprintf('\n开始PSO优化...\n');

for iter = 1:max_iter
    tic;
    w_current = w - (w - w_min) * iter / max_iter;
    
    for p = 1:n_particles
        params = positions(p, :);
        CCI_vector = zeros(1, n_runways);
        
        % 使用跑道编号加载对应的mat文件
        for runway_idx = 1:n_runways
            rwy_num = runway_number(runway_idx);
            mat_filename = sprintf('roughnessresult_%d.mat', rwy_num);
            CCI_vector(runway_idx) = calculate_CCI(params, mat_filename);
        end
        
        % CCI有效性检查
        cci_min_val = min(CCI_vector);
        cci_max_val = max(CCI_vector);
        cci_range = cci_max_val - cci_min_val;
        cci_std = std(CCI_vector);
        
        if cci_range < 1e-6 || cci_std < 1e-8 || any(isnan(CCI_vector)) || any(isinf(CCI_vector))
            fitness = 10;
            corr_weighted_mean = 0;
            corr_min_val = -1;
            corr_list = -ones(n_aircraft, 1);
        else
            CCI_norm = (CCI_vector - cci_min_val) / cci_range;
            
            % 计算与每个机型的相关系数
            corr_list = zeros(n_aircraft, 1);
            for i = 1:n_aircraft
                if std(CCI_norm) < 1e-10 || std(Energy_norm(i,:)) < 1e-10
                    corr_list(i) = 0;
                else
                    R = corrcoef(CCI_norm, Energy_norm(i, :));
                    if ~isnan(R(1,2)) && ~isinf(R(1,2))
                        corr_list(i) = R(1, 2);
                    else
                        corr_list(i) = 0;
                    end
                end
            end
            
            % 计算加权平均相关系数
            corr_weighted_mean = sum(corr_list .* aircraft_weights);
            corr_min_val = min(corr_list);
            
            % 适应度函数
            fitness = 0.7 * (1 - corr_weighted_mean) + 0.3 * (1 - corr_min_val);
        end
        
        % 更新个体最优
        if fitness < pbest_fitness(p)
            pbest_fitness(p) = fitness;
            pbest_positions(p, :) = positions(p, :);
        end
        
        % 更新全局最优
        if fitness < gbest_fitness
            gbest_fitness = fitness;
            gbest_position = positions(p, :);
            gbest_corr_weighted_mean = corr_weighted_mean;
            gbest_corr_min = corr_min_val;
            gbest_corr_list = corr_list;
        end
    end
    
    % 记录历史
    fitness_history(iter) = gbest_fitness;
    corr_weighted_mean_history(iter) = gbest_corr_weighted_mean;
    corr_min_history(iter) = gbest_corr_min;
    
    % 更新粒子速度和位置
    for p = 1:n_particles
        r1 = rand(1, n_params);
        r2 = rand(1, n_params);
        
        velocities(p, :) = w_current * velocities(p, :) + ...
                          c1 * r1 .* (pbest_positions(p, :) - positions(p, :)) + ...
                          c2 * r2 .* (gbest_position - positions(p, :));
        
        % 速度限制
        velocities(p, :) = max(-v_max, min(v_max, velocities(p, :)));
        
        % 位置更新
        positions(p, :) = positions(p, :) + velocities(p, :);
        
        % 边界处理
        for i = 1:n_params
            if positions(p, i) < param_min(i)
                positions(p, i) = param_min(i);
                velocities(p, i) = -velocities(p, i) * 0.5;
            elseif positions(p, i) > param_max(i)
                positions(p, i) = param_max(i);
                velocities(p, i) = -velocities(p, i) * 0.5;
            end
        end
    end
    
    iter_time = toc;
    if mod(iter, 10) == 0 || iter == 1
        fprintf('迭代%3d: fitness=%.4f, corr_weighted=%.4f, corr_min=%.4f, 耗时%.1fs\n', ...
                iter, gbest_fitness, gbest_corr_weighted_mean, gbest_corr_min, iter_time);
    end
end

fprintf('优化完成!\n');

%%  6. 输出结果 
fprintf('\n====== 最优结果 ======\n');
fprintf('最优适应度: %.6f\n', gbest_fitness);
fprintf('加权平均相关系数: %.4f\n', gbest_corr_weighted_mean);
fprintf('最小相关系数: %.4f\n', gbest_corr_min);

fprintf('\n最优标准机型参数:\n');
for i = 1:n_params
    fprintf('  %s = %.4f\n', param_names{i}, gbest_position(i));
end

%%  7. 绘图 
figure;
subplot(1,3,1);
plot(1:max_iter, fitness_history, 'b-', 'LineWidth', 1.5);
xlabel('迭代次数'); ylabel('适应度'); title('适应度收敛'); grid on;

subplot(1,3,2);
plot(1:max_iter, corr_weighted_mean_history, 'r-', 'LineWidth', 1.5);
xlabel('迭代次数'); ylabel('加权平均相关系数'); title('相关系数变化'); grid on;

subplot(1,3,3);
plot(1:max_iter, corr_min_history, 'g-', 'LineWidth', 1.5);
xlabel('迭代次数'); ylabel('最小相关系数'); title('最小相关系数'); grid on;

saveas(gcf, 'PSO_convergence.png');
