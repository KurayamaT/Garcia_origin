function stable_continuous_walker(flag, z0_input)
% 安定化された連続1.5ストライド歩行シミュレーション
% stable_continuous_walker(1,[0,-.1,0,-.3])

clc
close all
format long

% 引数の処理
if nargin == 0
    flag = 1;
    z0_input = [];
elseif nargin == 1
    z0_input = [];
end

% ウォーカーパラメータの設定
if flag == 1
    % Garcia's simplest walker
    walker.M = 1000; walker.m = 1.0; walker.I = 0.00; walker.l = 1.0; walker.w = 0.0; 
    walker.c = 1.0;  walker.r = 0.3; walker.g = 1.0; walker.gam = 0.009; 
    
    % 固定点
    zstar = [0.200161072169750  -0.199906060087682   0.400322144339512  -0.015805473227965];
else 
    % More General round feet walker
    walker.M = 1.0; walker.m = 0.5; walker.I = 0.02; walker.l = 1.0; walker.w = 0.0; 
    walker.c = 0.5; walker.r = 0.2; walker.g = 1.0; walker.gam = 0.01; 
    
    % 固定点
    zstar = [0.189472782205104  -0.239124222551699   0.378945564410209  -0.053691703909393];
end

% 初期値の設定（固定点近傍から開始）
if isempty(z0_input)
    fprintf('=== 安定化された連続歩行シミュレーション ===\n\n');
    fprintf('初期値を固定点近傍から自動設定します。\n');
    
    % 固定点に小さな摂動を加える
    perturbation = 0.01 * randn(1,4);
    z0_input = zstar + perturbation;
    
    fprintf('初期値: [%.6f, %.6f, %.6f, %.6f]\n', z0_input);
end

% パラメータ設定
num_cycles = 5;
max_angle_deviation = 0.5; % 角度の最大偏差（rad）
max_velocity_deviation = 1.0; % 角速度の最大偏差（rad/s）

% データ保存用変数
all_states = zeros(num_cycles+1, 4);
all_states(1, :) = z0_input;
all_times = [0];
success_flags = zeros(num_cycles, 1);
x_positions = [0];

current_z0 = z0_input;
current_time = 0;
current_x = 0;

fprintf('\n=== シミュレーション開始 ===\n');

for cycle = 1:num_cycles
    fprintf('\n--- サイクル %d/%d ---\n', cycle, num_cycles);
    
    % 状態の妥当性チェック
    if any(abs(current_z0([1,3])) > max_angle_deviation) || ...
       any(abs(current_z0([2,4])) > max_velocity_deviation)
        fprintf('⚠️  警告: 状態が不安定になりました。固定点にリセットします。\n');
        current_z0 = zstar + 0.01 * randn(1,4);
    end
    
    try
        % 1.5ストライドの実行
        [z_final, t_elapsed, x_traveled, success] = execute_one_half_stride_safe(current_z0, walker);
        
        if success
            fprintf('✅ 成功: 最終状態 [%.4f, %.4f, %.4f, %.4f]\n', z_final);
            fprintf('   脚間角度差: %.6f rad\n', abs(z_final(1) - z_final(3)));
            fprintf('   移動距離: %.4f m\n', x_traveled);
            
            all_states(cycle+1, :) = z_final;
            current_time = current_time + t_elapsed;
            all_times = [all_times, current_time];
            current_x = current_x + x_traveled;
            x_positions = [x_positions, current_x];
            success_flags(cycle) = 1;
            
            % 次のサイクルの初期条件を生成（安定化処理付き）
            if cycle < num_cycles
                current_z0 = generate_next_initial_condition(z_final, zstar);
            end
        else
            fprintf('❌ 失敗: 歩行が不安定になりました。\n');
            success_flags(cycle) = 0;
            
            % 固定点から再スタート
            current_z0 = zstar + 0.01 * randn(1,4);
            all_states(cycle+1, :) = current_z0;
        end
        
    catch ME
        fprintf('❌ エラー: %s\n', ME.message);
        success_flags(cycle) = 0;
        
        % 固定点から再スタート
        current_z0 = zstar + 0.01 * randn(1,4);
        all_states(cycle+1, :) = current_z0;
    end
end

% 結果の表示
fprintf('\n=== シミュレーション結果 ===\n');
fprintf('成功したサイクル: %d/%d (%.1f%%)\n', sum(success_flags), num_cycles, 100*sum(success_flags)/num_cycles);
fprintf('総移動距離: %.3f m\n', current_x);

% プロット
plot_results(all_states, all_times, x_positions, success_flags, zstar);

end

%% 補助関数

function [z_final, t_elapsed, x_traveled, success] = execute_one_half_stride_safe(z0, walker)
    % 安全な1.5ストライド実行
    
    success = false;
    max_time = 10; % 最大シミュレーション時間
    
    try
        % 初期エネルギー計算
        TE0 = calculate_energy(z0, walker);
        z0_extended = [z0, TE0, 0, 0, walker.l*cos(z0(1)), 0];
        
        % 1ステップ目
        options = odeset('abstol',1e-10,'reltol',1e-10,'events',@collision, 'MaxStep', 0.01);
        [t1, z1] = ode45(@single_stance, [0, max_time/2], z0_extended, options, walker);
        
        if t1(end) < max_time/2
            % 衝突検出
            z_after_collision = heelstrike(t1(end), z1(end,:), walker);
            
            % 0.5ステップ目
            options2 = odeset('abstol',1e-10,'reltol',1e-10,'events',@swing_catches_stance, 'MaxStep', 0.01);
            [t2, z2] = ode45(@single_stance, [0, max_time/2], z_after_collision, options2, walker);
            
            if t2(end) < max_time/2
                z_final = z2(end, 1:4);
                t_elapsed = t1(end) + t2(end);
                x_traveled = z2(end, 6);
                
                % 成功判定
                if abs(z_final(1) - z_final(3)) < 0.01 && all(abs(z_final) < 2.0)
                    success = true;
                end
            else
                z_final = z2(end, 1:4);
                t_elapsed = t1(end) + t2(end);
                x_traveled = z2(end, 6);
            end
        else
            z_final = z0;
            t_elapsed = max_time;
            x_traveled = 0;
        end
        
    catch
        z_final = z0;
        t_elapsed = 0;
        x_traveled = 0;
    end
end

function z_next = generate_next_initial_condition(z_final, zstar)
    % 次の初期条件を安定的に生成
    
    % 基本的な脚の切り替え
    q1_next = z_final(3);
    q2_next = 2 * q1_next;
    
    % 角度の制限
    max_angle = 0.5;
    q1_next = max(min(q1_next, max_angle), -max_angle);
    q2_next = max(min(q2_next, max_angle), -max_angle);
    
    % 速度は固定点の比率を参考に
    velocity_ratio = zstar(2) / zstar(4);
    u1_next = z_final(4);
    u2_next = u1_next / velocity_ratio;
    
    % 速度の制限
    max_velocity = 1.0;
    u1_next = max(min(u1_next, max_velocity), -max_velocity);
    u2_next = max(min(u2_next, max_velocity), -max_velocity);
    
    z_next = [q1_next, u1_next, q2_next, u2_next];
end

function plot_results(all_states, all_times, x_positions, success_flags, zstar)
    % 結果のプロット
    
    figure('Position', [100, 100, 1200, 800]);
    
    % 角度の時間変化
    subplot(2,3,1)
    plot(all_times, all_states(:,1), 'r-', 'LineWidth', 2)
    hold on
    plot(all_times, all_states(:,3), 'b-', 'LineWidth', 2)
    plot(all_times, zstar(1)*ones(size(all_times)), 'k--')
    plot(all_times, zstar(3)*ones(size(all_times)), 'k--')
    xlabel('時間 (s)')
    ylabel('角度 (rad)')
    title('脚角度の時間変化')
    legend('q1 (スタンス)', 'q2 (スイング)', '固定点', 'Location', 'best')
    grid on
    
    % 角速度の時間変化
    subplot(2,3,2)
    plot(all_times, all_states(:,2), 'r-', 'LineWidth', 2)
    hold on
    plot(all_times, all_states(:,4), 'b-', 'LineWidth', 2)
    xlabel('時間 (s)')
    ylabel('角速度 (rad/s)')
    title('脚角速度の時間変化')
    legend('u1', 'u2', 'Location', 'best')
    grid on
    
    % 位相図 (q1-u1)
    subplot(2,3,3)
    plot(all_states(:,1), all_states(:,2), 'r.-', 'LineWidth', 1.5)
    hold on
    plot(zstar(1), zstar(2), 'k*', 'MarkerSize', 12)
    plot(all_states(1,1), all_states(1,2), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g')
    plot(all_states(end,1), all_states(end,2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r')
    xlabel('q1 (rad)')
    ylabel('u1 (rad/s)')
    title('位相図 (スタンス脚)')
    legend('軌道', '固定点', '開始', '終了', 'Location', 'best')
    grid on
    
    % 位相図 (q2-u2)
    subplot(2,3,4)
    plot(all_states(:,3), all_states(:,4), 'b.-', 'LineWidth', 1.5)
    hold on
    plot(zstar(3), zstar(4), 'k*', 'MarkerSize', 12)
    xlabel('q2 (rad)')
    ylabel('u2 (rad/s)')
    title('位相図 (スイング脚)')
    grid on
    
    % 移動距離
    subplot(2,3,5)
    plot(all_times, x_positions, 'k-', 'LineWidth', 2)
    xlabel('時間 (s)')
    ylabel('移動距離 (m)')
    title('累積移動距離')
    grid on
    
    % 成功率
    subplot(2,3,6)
    bar(1:length(success_flags), success_flags)
    xlabel('サイクル番号')
    ylabel('成功 (1) / 失敗 (0)')
    title('各サイクルの成功/失敗')
    ylim([-0.1, 1.1])
    grid on
    
    sgtitle('安定化された連続歩行の結果')
end

%% 元のコードから必要な関数（省略部分は元のコードと同じ）
function TE = calculate_energy(z, walker)
    q1 = z(1); u1 = z(2); q2 = z(3); u2 = z(4);
    M = walker.M;  m = walker.m; I = walker.I;   
    l = walker.l;  c = walker.c; w = walker.w;   
    r = walker.r;  g = walker.g; gam = walker.gam;

    TE = 1/2*m*(((-l*cos(q1)-r)*u1-u1*(-c*cos(q1)+w*sin(q1)))^2+(-l*sin(q1)*u1+u1*(c*sin(q1)+w*cos(q1)))^2)+1/2*m*(((-l*cos(q1)-r)*u1-(u1-u2)*(-c*cos(q1-q2)+w*sin(q1-q2)))^2+(-l*sin(q1)*u1+(u1-u2)*(c*sin(q1-q2)+w*cos(q1-q2)))^2)+1/2*M*((-l*cos(q1)-r)^2*u1^2+l^2*sin(q1)^2*u1^2)+1/2*I*(u1^2+(u1-u2)^2)+2*m*g*cos(gam)*r+2*m*g*l*cos(gam-q1)-m*g*c*cos(gam-q1)-m*g*w*sin(gam-q1)+2*m*g*sin(gam)*r*q1-m*g*c*cos(gam-q1+q2)-m*g*w*sin(gam-q1+q2)+M*g*cos(gam)*r+M*g*l*cos(gam-q1)+M*g*sin(gam)*r*q1;
end

% 以下、single_stance, collision, swing_catches_stance, heelstrike 関数は元のコードと同じ