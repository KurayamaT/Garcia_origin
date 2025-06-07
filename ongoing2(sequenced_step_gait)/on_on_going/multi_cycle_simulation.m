function multi_cycle_simulation()
    % 複数サイクルパッシブウォーカーシミュレーション
    % passive_walker_physics.mの関数を使用
    
    clear; close all;
    
    % パスに物理計算関数があることを確認
    if ~exist('one_and_half_stride_detailed', 'file')
        error('passive_walker_physics.mが見つかりません。同じフォルダに配置してください。');
    end
    
    % === シミュレーション設定 ===
    fprintf('=== 複数サイクルパッシブウォーカーシミュレーション ===\n');
    fprintf('passive_walker_physics.mの物理モデルを使用します\n\n');
    
    % ウォーカータイプの選択
    fprintf('ウォーカータイプを選択してください:\n');
    fprintf('1. Garcia''s Simplest Walker\n');
    fprintf('2. General Round Feet Walker\n');
    flag = input('選択 (1-2): ');
    if isempty(flag)
        flag = 1;
    end
    
    % ウォーカーパラメータの設定
    if flag == 1
        walker.M = 1000; walker.m = 1.0; walker.I = 0.00; walker.l = 1.0; walker.w = 0.0; 
        walker.c = 1.0;  walker.r = 0.0; walker.g = 1.0; walker.gam = 0.009; 
        zstar = [0.200161072169750; -0.199906060087682; 0.400322144339512; -0.015805473227965];
    else 
        walker.M = 1.0; walker.m = 0.5; walker.I = 0.02; walker.l = 1.0; walker.w = 0.0; 
        walker.c = 0.5; walker.r = 0.2; walker.g = 1.0; walker.gam = 0.01; 
        zstar = [0.189472782205104; -0.239124222551699; 0.378945564410209; -0.053691703909393];
    end
    
    % サイクル数の入力
    n_cycles = input('実行するサイクル数を入力してください (1-10): ');
    if isempty(n_cycles) || n_cycles < 1
        n_cycles = 3;
    end
    
    % 初期値の選択
    z0_original = select_initial_state(flag, zstar);
    
    % アニメーション速度
    anim_speed_factor = select_animation_speed();
    
    % 補正方法の選択
    correction_method = select_correction_method();
    
    % === メインシミュレーションループ ===
    [all_times, all_trajectories, all_extended_states, cycle_data] = ...
        run_simulation_loop(n_cycles, z0_original, walker, correction_method, anim_speed_factor);
    
    % === 全体のアニメーション再生 ===
    fprintf('\n全サイクルの連続アニメーションを再生しますか？ (y/n): ');
    if input('', 's') == 'y'
        animate_all_cycles(all_times, all_extended_states, walker, anim_speed_factor);
    end
    
    % === 結果の可視化 ===
    visualize_results(all_times, all_trajectories, all_extended_states, cycle_data, z0_original);
    
    % === サマリー表示 ===
    display_summary(cycle_data, z0_original);
end

%% 補助関数

function z0 = select_initial_state(flag, zstar)
    fprintf('\n初期値の選択:\n');
    fprintf('1. 固定点から開始（推奨）\n');
    fprintf('2. プリセット初期値\n');
    fprintf('3. カスタム初期値\n');
    init_choice = input('選択 (1-3): ');
    
    switch init_choice
        case 1
            z0 = zstar;
            fprintf('固定点から開始します\n');
        case 2
            if flag == 1
                presets = [
                    0.2, -0.2, 0.4, -0.3;
                    0.1, -0.1, 0.2, -0.2;
                    -0.01, -0.1, -0.01, -0.3;
                ];
            else
                presets = [
                    0.15, -0.2, 0.3, -0.3;
                    0.25, -0.3, 0.4, -0.2;
                    0.18, -0.24, 0.36, -0.06;
                ];
            end
            fprintf('\nプリセット選択:\n');
            for i = 1:size(presets, 1)
                fprintf('%d: [%.2f, %.2f, %.2f, %.2f]\n', i, presets(i,:));
            end
            preset_idx = input('選択: ');
            z0 = presets(preset_idx, :)';
        otherwise
            z0 = input('初期値 [q1; u1; q2; u2] を入力: ');
            if size(z0, 2) > 1
                z0 = z0';
            end
    end
end

function speed_factor = select_animation_speed()
    fprintf('\nアニメーション設定:\n');
    fprintf('1. リアルタイム表示\n');
    fprintf('2. 高速表示（2倍速）\n');
    fprintf('3. スロー表示（0.5倍速）\n');
    anim_speed = input('選択 (1-3): ');
    speed_factors = [1, 2, 0.5];
    speed_factor = speed_factors(anim_speed);
end

function method = select_correction_method()
    fprintf('\n角速度補正方法の選択:\n');
    fprintf('1. エネルギー補償による補正（推奨）\n');
    fprintf('2. 固定倍率による補正\n');
    fprintf('3. 手動補正\n');
    method = input('選択 (1-3): ');
end

function [all_times, all_trajectories, all_extended_states, cycle_data] = ...
    run_simulation_loop(n_cycles, z0_original, walker, correction_method, anim_speed_factor)
    
    all_trajectories = {};
    all_times = {};
    all_extended_states = {};
    cycle_data = [];
    z_current = z0_original;
    cumulative_time = 0;
    
    fprintf('\n=== シミュレーション開始 ===\n');
    
    % アニメーション用のfigure準備
    fig_anim = figure('Name', 'パッシブウォーカー アニメーション', ...
                      'Position', [100 100 1200 600]);
    
    for cycle = 1:n_cycles
        fprintf('\n--- サイクル %d/%d ---\n', cycle, n_cycles);
        
        % 1.5歩のシミュレーション実行
        [z_traj, t_traj, z_final, z_midpoint] = one_and_half_stride_detailed(z_current, walker);
        
        if isempty(z_traj) || size(z_traj, 1) < 2
            fprintf('エラー: サイクル %d で歩行が失敗しました\n', cycle);
            break;
        end
        
        % 時間を累積
        t_traj_adjusted = t_traj + cumulative_time;
        cumulative_time = t_traj_adjusted(end);
        
        % 結果を保存
        all_trajectories{cycle} = z_traj(:, 1:4);
        all_extended_states{cycle} = z_traj;
        all_times{cycle} = t_traj_adjusted;
        
        % サイクル情報を記録
        cycle_info = record_cycle_info(cycle, z_current, z_final, z_midpoint, t_traj, walker);
        
        % リアルタイムアニメーション表示
        animate_cycle(fig_anim, t_traj, z_traj, walker, cycle, anim_speed_factor);
        
        % 次のサイクルのための角速度補正
        if cycle < n_cycles
            [z_current, velocity_correction] = apply_velocity_correction(...
                z_final, z0_original, cycle_info, correction_method);
            cycle_info.velocity_correction = velocity_correction;
        end
        
        cycle_data = [cycle_data, cycle_info];
        pause(0.5);
    end
end

function cycle_info = record_cycle_info(cycle, z_current, z_final, z_midpoint, t_traj, walker)
    cycle_info = struct();
    cycle_info.number = cycle;
    cycle_info.initial_state = z_current;
    cycle_info.final_state = z_final;
    cycle_info.midpoint_state = z_midpoint;
    cycle_info.duration = t_traj(end);
    
    % エネルギー計算
    TE_initial = calculate_energy(z_current, walker);
    TE_final = calculate_energy(z_final, walker);
    cycle_info.energy_initial = TE_initial;
    cycle_info.energy_final = TE_final;
    cycle_info.energy_loss = TE_final - TE_initial;
    
    fprintf('初期状態: [%.4f, %.4f, %.4f, %.4f]\n', z_current);
    fprintf('最終状態: [%.4f, %.4f, %.4f, %.4f]\n', z_final);
    fprintf('エネルギー変化: %.6f\n', cycle_info.energy_loss);
end

function [z_next, correction_norm] = apply_velocity_correction(z_final, z0_original, cycle_info, method)
    fprintf('\n次サイクルの準備...\n');
    
    switch method
        case 1  % エネルギー補償
            if cycle_info.energy_loss < 0
                energy_compensation = sqrt(cycle_info.energy_initial / cycle_info.energy_final);
            else
                energy_compensation = 1;
            end
            z_next = [z_final(1); 
                     z0_original(2) * energy_compensation;
                     z_final(3); 
                     z0_original(4) * energy_compensation];
            
        case 2  % 固定倍率
            correction_factor = 1.1;
            z_next = [z_final(1); 
                     z0_original(2) * correction_factor;
                     z_final(3); 
                     z0_original(4) * correction_factor];
            
        case 3  % 手動補正
            fprintf('現在の角度: q1=%.4f, q2=%.4f\n', z_final(1), z_final(3));
            u1_corrected = input('新しい角速度 u1 を入力: ');
            u2_corrected = input('新しい角速度 u2 を入力: ');
            z_next = [z_final(1); u1_corrected; z_final(3); u2_corrected];
    end
    
    fprintf('補正後の角速度: u1=%.4f, u2=%.4f\n', z_next(2), z_next(4));
    correction_norm = norm([z_next(2) - z_final(2); z_next(4) - z_final(4)]);
end

%% アニメーション関数

function animate_cycle(fig, t, z_all, walker, cycle_num, speed_factor)
    figure(fig);
    clf;
    
    % パラメータ
    l = walker.l;
    r = walker.r;
    
    % アニメーション用データの準備
    fps = 10;
    total_frames = round(t(end)*fps);
    t_anim = linspace(t(1), t(end), total_frames);
    
    % 補間
    z_anim = zeros(total_frames, size(z_all, 2));
    for i = 1:size(z_all, 2)
        z_anim(:,i) = interp1(t, z_all(:,i), t_anim);
    end
    
    % ウィンドウ設定
    if size(z_all, 2) >= 8
        min_xh = min(z_all(:,6)); max_xh = max(z_all(:,6));
    else
        min_xh = -2; max_xh = 2;
    end
    window_xmin = min_xh - 1*l; window_xmax = max_xh + 1*l;
    window_ymin = -0.1; window_ymax = 1.1*(l+r);
    
    % サブプロット1: アニメーション
    subplot(1,2,1);
    axis('equal')
    axis([window_xmin window_xmax window_ymin window_ymax])
    title(sprintf('サイクル %d - 歩行アニメーション', cycle_num));
    xlabel('水平位置');
    ylabel('垂直位置');
    grid on;
    
    % サブプロット2: 状態変数
    subplot(1,2,2);
    title('角度と角速度');
    xlabel('時間 [s]');
    grid on;
    
    % アニメーションループ
    dt_anim = 0.1 / speed_factor;
    
    for i = 1:length(t_anim)
        % 左側：歩行アニメーション
        subplot(1,2,1);
        cla;
        plot([window_xmin window_xmax], [0 0], 'k-', 'LineWidth', 2);
        hold on;
        
        % 現在の状態
        q1 = z_anim(i,1); q2 = z_anim(i,3);
        
        if size(z_anim, 2) >= 8
            xh = z_anim(i,6); yh = z_anim(i,8);
        else
            xh = -l*sin(q1) + t_anim(i) * 0.2;
            yh = l*cos(q1) + r;
        end
        
        % ウォーカーを描画
        draw_walker(q1, q2, xh, yh, walker);
        
        % 右側：状態変数のプロット
        subplot(1,2,2);
        cla;
        
        % これまでの軌跡
        idx = find(t <= t_anim(i));
        plot(t(idx), z_all(idx,1), 'b-', 'LineWidth', 2, 'DisplayName', 'q1');
        plot(t(idx), z_all(idx,3), 'r-', 'LineWidth', 2, 'DisplayName', 'q2');
        plot(t(idx), z_all(idx,2), 'b--', 'LineWidth', 1, 'DisplayName', 'u1');
        plot(t(idx), z_all(idx,4), 'r--', 'LineWidth', 1, 'DisplayName', 'u2');
        
        xlim([t(1) t(end)]);
        ylim([-1 1]);
        legend('Location', 'best');
        
        drawnow;
        pause(dt_anim);
    end
end

function animate_all_cycles(all_times, all_extended_states, walker, speed_factor)
    figure('Name', '全サイクル連続アニメーション', 'Position', [100 100 1200 600]);
    
    % 全データを結合
    t_all = [];
    z_all = [];
    for i = 1:length(all_times)
        t_all = [t_all; all_times{i}];
        z_all = [z_all; all_extended_states{i}];
    end
    
    % アニメーション設定
    l = walker.l;
    r = walker.r;
    
    if size(z_all, 2) >= 8
        min_xh = min(z_all(:,6)); max_xh = max(z_all(:,6));
    else
        min_xh = -2; max_xh = 2*length(all_times);
    end
    window_xmin = min_xh - 1*l; window_xmax = max_xh + 1*l;
    window_ymin = -0.1; window_ymax = 1.1*(l+r);
    
    subplot(1,2,1);
    axis('equal')
    axis([window_xmin window_xmax window_ymin window_ymax])
    title('全サイクル 歩行アニメーション');
    xlabel('水平位置');
    ylabel('垂直位置');
    grid on;
    
    subplot(1,2,2);
    title('全サイクル 状態変数');
    xlabel('時間 [s]');
    grid on;
    
    % アニメーションループ
    fps = 10;
    dt_anim = 0.1 / speed_factor;
    t_anim = t_all(1):dt_anim:t_all(end);
    
    hip_trajectory = [];
    
    for i = 1:length(t_anim)
        current_t = t_anim(i);
        
        % 補間して現在の状態を取得
        current_z = zeros(1, size(z_all, 2));
        for j = 1:size(z_all, 2)
            current_z(j) = interp1(t_all, z_all(:,j), current_t);
        end
        
        % 左側：歩行アニメーション
        subplot(1,2,1);
        cla;
        plot([window_xmin window_xmax], [0 0], 'k-', 'LineWidth', 2);
        hold on;
        
        q1 = current_z(1); q2 = current_z(3);
        
        if size(current_z, 2) >= 8
            xh = current_z(6); yh = current_z(8);
        else
            xh = -l*sin(q1) + current_t * 0.2;
            yh = l*cos(q1) + r;
        end
        
        % 腰の軌跡
        hip_trajectory = [hip_trajectory; xh, yh];
        if size(hip_trajectory, 1) > 1
            plot(hip_trajectory(:,1), hip_trajectory(:,2), 'g--', 'LineWidth', 1);
        end
        
        draw_walker(q1, q2, xh, yh, walker);
        
        % 右側：状態変数
        subplot(1,2,2);
        cla;
        idx = find(t_all <= current_t);
        plot(t_all(idx), z_all(idx,1), 'b-', 'LineWidth', 2);
        plot(t_all(idx), z_all(idx,3), 'r-', 'LineWidth', 2);
        plot(t_all(idx), z_all(idx,2), 'b--', 'LineWidth', 1);
        plot(t_all(idx), z_all(idx,4), 'r--', 'LineWidth', 1);
        
        xlim([t_all(1) t_all(end)]);
        ylim([-1 1]);
        
        drawnow;
        pause(dt_anim);
    end
end

function draw_walker(q1, q2, xh, yh, walker)
    l = walker.l;
    r = walker.r;
    
    % スタンス脚（スタンス脚の足は原点）
    stance_foot_x = xh + l*sin(q1) + r*q1;
    stance_foot_y = 0;
    
    % スイング脚
    swing_foot_x = xh + l*sin(q2) + r*(q2-q1);
    swing_foot_y = yh - l*cos(q2) - r;
    
    % 脚を描画
    plot([stance_foot_x, xh], [stance_foot_y, yh], 'b-', 'LineWidth', 3);
    plot([xh, swing_foot_x], [yh, swing_foot_y], 'r-', 'LineWidth', 3);
    
    % 関節と足を描画
    plot(xh, yh, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
    plot(stance_foot_x, stance_foot_y, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    plot(swing_foot_x, swing_foot_y, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    
    % 進行方向を示す矢印
    quiver(xh, yh+0.2, 0.2, 0, 'k', 'LineWidth', 2);
end

%% 結果の可視化

function visualize_results(all_times, all_trajectories, all_extended_states, cycle_data, z0_original)
    figure('Name', '複数サイクル解析結果', 'Position', [50 50 1400 900]);
    
    colors = lines(length(all_times));
    
    % 角度の時系列
    subplot(3,3,1);
    hold on;
    for i = 1:length(all_times)
        t = all_times{i};
        z = all_trajectories{i};
        plot(t, z(:,1), '-', 'Color', colors(i,:), 'LineWidth', 1.5);
        plot(t, z(:,3), '--', 'Color', colors(i,:), 'LineWidth', 1.5);
    end
    xlabel('時間 [s]');
    ylabel('角度 [rad]');
    title('角度の時間変化');
    grid on;
    legend('q1', 'q2', 'Location', 'best');
    
    % 角速度の時系列
    subplot(3,3,2);
    hold on;
    for i = 1:length(all_times)
        t = all_times{i};
        z = all_trajectories{i};
        plot(t, z(:,2), '-', 'Color', colors(i,:), 'LineWidth', 1.5);
        plot(t, z(:,4), '--', 'Color', colors(i,:), 'LineWidth', 1.5);
    end
    xlabel('時間 [s]');
    ylabel('角速度 [rad/s]');
    title('角速度の時間変化');
    grid on;
    
    % エネルギーの推移
    subplot(3,3,3);
    cycles = [cycle_data.number];
    E_initial = [cycle_data.energy_initial];
    E_final = [cycle_data.energy_final];
    
    plot(cycles, E_initial, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
    hold on;
    plot(cycles, E_final, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
    
    xlabel('サイクル');
    ylabel('エネルギー');
    title('サイクルごとのエネルギー変化');
    legend('初期', '最終', 'Location', 'best');
    grid on;
    
    % 位相空間（q1-u1）
    subplot(3,3,4);
    hold on;
    for i = 1:length(all_trajectories)
        z = all_trajectories{i};
        plot(z(:,1), z(:,2), 'Color', colors(i,:), 'LineWidth', 1.5);
    end
    plot(z0_original(1), z0_original(2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
    xlabel('q1 [rad]');
    ylabel('u1 [rad/s]');
    title('位相空間 (q1-u1)');
    grid on;
    
    % 位相空間（q2-u2）
    subplot(3,3,5);
    hold on;
    for i = 1:length(all_trajectories)
        z = all_trajectories{i};
        plot(z(:,3), z(:,4), 'Color', colors(i,:), 'LineWidth', 1.5);
    end
    plot(z0_original(3), z0_original(4), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
    xlabel('q2 [rad]');
    ylabel('u2 [rad/s]');
    title('位相空間 (q2-u2)');
    grid on;
    
    % エネルギーの時系列
    subplot(3,3,6);
    hold on;
    for i = 1:length(all_extended_states)
        if size(all_extended_states{i}, 2) >= 5
            t = all_times{i};
            z = all_extended_states{i};
            plot(t, z(:,5), 'Color', colors(i,:), 'LineWidth', 1.5);
        end
    end
    xlabel('時間 [s]');
    ylabel('総エネルギー');
    title('エネルギーの時間変化');
    grid on;
    
    % ヒップ軌道
    subplot(3,3,7);
    hold on;
    for i = 1:length(all_extended_states)
        if size(all_extended_states{i}, 2) >= 8
            z = all_extended_states{i};
            plot(z(:,6), z(:,8), 'Color', colors(i,:), 'LineWidth', 1.5);
        end
    end
    xlabel('x [m]');
    ylabel('y [m]');
    title('ヒップ軌道');
    axis equal;
    grid on;
    
    % 角速度補正量
    subplot(3,3,8);
    if length(cycle_data) > 1
        corrections = zeros(length(cycle_data)-1, 1);
        for i = 1:length(cycle_data)-1
            if isfield(cycle_data(i), 'velocity_correction')
                corrections(i) = cycle_data(i).velocity_correction;
            end
        end
        bar(1:length(corrections), corrections);
        xlabel('サイクル');
        ylabel('角速度補正量');
        title('各サイクルでの角速度補正');
        grid on;
    end
    
    % 状態ベクトルのノルム
    subplot(3,3,9);
    hold on;
    for i = 1:length(all_trajectories)
        t = all_times{i};
        z = all_trajectories{i};
        state_norm = sqrt(sum(z.^2, 2));
        plot(t, state_norm, 'Color', colors(i,:), 'LineWidth', 1.5);
    end
    xlabel('時間 [s]');
    ylabel('状態ベクトルのノルム');
    title('状態空間での軌道');
    grid on;
    
    sgtitle('パッシブウォーカー: 複数サイクル解析');
end

function display_summary(cycle_data, z0_original)
    fprintf('\n=== シミュレーションサマリー ===\n');
    fprintf('総サイクル数: %d\n', length(cycle_data));
    
    if length(cycle_data) > 0
        avg_energy_change = mean([cycle_data.energy_loss]);
        fprintf('平均エネルギー変化: %.6f\n', avg_energy_change);
        
        if isfield(cycle_data(1), 'velocity_correction')
            avg_correction = mean([cycle_data(1:end-1).velocity_correction]);
            fprintf('平均角速度補正量: %.6f\n', avg_correction);
        end
        
        % 各サイクルの偏差
        fprintf('\n--- 各サイクルの目標からの偏差 ---\n');
        for i = 1:length(cycle_data)
            deviation = norm(cycle_data(i).final_state - z0_original);
            fprintf('サイクル %d: %.6f\n', i, deviation);
        end
    end
    
    fprintf('\n=== シミュレーション完了 ===\n');
end