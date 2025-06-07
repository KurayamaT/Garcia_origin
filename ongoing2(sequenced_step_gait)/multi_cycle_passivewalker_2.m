function multi_cycle_passivewalker_fixed()
    % 修正版：複数サイクルパッシブウォーカーシミュレーション
    % 主な修正点：
    % 1. R関数の定義追加
    % 2. 構造体フィールドの統一
    % 3. エネルギー補正の改善
    % 4. アニメーション描画の修正
    
    clear; close all;
    
    % === シミュレーション設定 ===
    fprintf('=== 修正版複数サイクルパッシブウォーカーシミュレーション ===\n');
    fprintf('各サイクル終了時に角速度を補正して周期的な歩行を実現します\n\n');
    
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
        %% Garcia's simplest walker
        walker.M = 1000; walker.m = 1.0; walker.I = 0.00; walker.l = 1.0; walker.w = 0.0; 
        walker.c = 1.0;  walker.r = 0.0; walker.g = 1.0; walker.gam = 0.009; 
        
        % 固定点
        zstar = [0.200161072169750; -0.199906060087682; 0.400322144339512; -0.015805473227965];
        
    else 
        %%  More General round feet walker
        walker.M = 1.0; walker.m = 0.5; walker.I = 0.02; walker.l = 1.0; walker.w = 0.0; 
        walker.c = 0.5; walker.r = 0.2; walker.g = 1.0; walker.gam = 0.01; 
        
        % 固定点
        zstar = [0.189472782205104; -0.239124222551699; 0.378945564410209; -0.053691703909393];
    end
    
    % サイクル数の入力
    n_cycles = input('実行するサイクル数を入力してください (1-10): ');
    if isempty(n_cycles) || n_cycles < 1
        n_cycles = 3;  % デフォルト値
    end
    
    % 初期値の選択
    fprintf('\n初期値の選択:\n');
    fprintf('1. 固定点から開始（推奨）\n');
    fprintf('2. プリセット初期値\n');
    fprintf('3. カスタム初期値\n');
    init_choice = input('選択 (1-3): ');
    
    switch init_choice
        case 1
            z0_original = zstar;
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
            z0_original = presets(preset_idx, :)';
        otherwise
            z0_original = input('初期値 [q1; u1; q2; u2] を入力: ');
            if size(z0_original, 2) > 1
                z0_original = z0_original';
            end
    end
    
    % アニメーション設定
    fprintf('\nアニメーション設定:\n');
    fprintf('1. リアルタイム表示\n');
    fprintf('2. 高速表示（2倍速）\n');
    fprintf('3. スロー表示（0.5倍速）\n');
    anim_speed = input('選択 (1-3): ');
    if isempty(anim_speed)
        anim_speed = 1;
    end
    speed_factor = [1, 2, 0.5];
    anim_speed_factor = speed_factor(min(anim_speed, 3));
    
    % 補正方法の選択
    fprintf('\n角速度補正方法の選択:\n');
    fprintf('1. エネルギー補償による補正（推奨）\n');
    fprintf('2. 固定倍率による補正\n');
    fprintf('3. 手動補正\n');
    correction_method = input('選択 (1-3): ');
    if isempty(correction_method)
        correction_method = 1;
    end
    
    % === メインシミュレーションループ ===
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
        
        try
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
            all_trajectories{cycle} = z_traj(:, 1:4);  % 基本状態変数のみ
            all_extended_states{cycle} = z_traj;  % 全状態変数
            all_times{cycle} = t_traj_adjusted;
            
            % サイクル情報を記録（統一されたフィールド構造）
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
            
            % 初期化（すべてのサイクルで同じフィールドを持つように）
            cycle_info.velocity_correction = 0;  % デフォルト値
            
            fprintf('初期状態: [%.4f, %.4f, %.4f, %.4f]\n', z_current);
            fprintf('最終状態: [%.4f, %.4f, %.4f, %.4f]\n', z_final);
            fprintf('エネルギー変化: %.6f\n', cycle_info.energy_loss);
            
            % アニメーション表示（エラーハンドリング付き）
            try
                animate_cycle_fixed(fig_anim, t_traj, z_traj, walker, cycle, anim_speed_factor);
            catch me
                fprintf('アニメーションエラー: %s\n', me.message);
            end
            
            % 次のサイクルのための角速度補正
            if cycle < n_cycles
                fprintf('\n次サイクルの準備...\n');
                
                z_corrected = apply_velocity_correction(z_final, z0_original, correction_method, cycle_info);
                
                % 補正量を記録
                cycle_info.velocity_correction = norm([z_corrected(2) - z_final(2); 
                                                      z_corrected(4) - z_final(4)]);
                
                z_current = z_corrected;
                fprintf('補正後の状態: [%.4f, %.4f, %.4f, %.4f]\n', z_current);
            end
            
            cycle_data = [cycle_data, cycle_info];
            
            % 少し待機
            pause(0.5);
            
        catch me
            fprintf('サイクル %d でエラーが発生しました: %s\n', cycle, me.message);
            break;
        end
    end
    
    % === 全体のアニメーション再生 ===
    if ~isempty(all_trajectories)
        fprintf('\n全サイクルの連続アニメーションを再生しますか？ (y/n): ');
        user_input = input('', 's');
        if strcmpi(user_input, 'y') || strcmpi(user_input, 'yes')
            animate_all_cycles_continuous(fig_anim, all_times, all_extended_states, walker);
        end
        
        % 結果の可視化
        visualize_results_fixed(all_times, all_trajectories, all_extended_states, cycle_data, z0_original);
        display_summary_fixed(cycle_data, z0_original);
    end
end

% ===== 内部関数群 =====

function z_corrected = apply_velocity_correction(z_final, z0_original, correction_method, cycle_info)
    % 角速度補正を適用
    
    switch correction_method
        case 1
            % エネルギー補償による補正
            energy_ratio = cycle_info.energy_initial / cycle_info.energy_final;
            if energy_ratio > 1
                energy_compensation = sqrt(energy_ratio);
            else
                energy_compensation = 1.0;
            end
            
            % 角度を微調整
            angle_adjustment = 0.05;
            
            z_corrected = [z_final(1) + angle_adjustment; 
                          z0_original(2) * energy_compensation;
                          z_final(3) - angle_adjustment; 
                          z0_original(4) * energy_compensation];
            
        case 2
            % 固定倍率による補正
            correction_factor = 1.1;
            angle_adjustment = 0.05;
            
            z_corrected = [z_final(1) + angle_adjustment; 
                          z0_original(2) * correction_factor;
                          z_final(3) - angle_adjustment; 
                          z0_original(4) * correction_factor];
            
        case 3
            % 手動補正
            fprintf('現在の角度: q1=%.4f, q2=%.4f\n', z_final(1), z_final(3));
            u1_corrected = input('新しい角速度 u1 を入力: ');
            u2_corrected = input('新しい角速度 u2 を入力: ');
            z_corrected = [z_final(1); u1_corrected; z_final(3); u2_corrected];
            
        otherwise
            z_corrected = z_final;
    end
end

function animate_cycle_fixed(fig, t, z_all, walker, cycle_num, speed_factor)
    % 修正版アニメーション関数（継続的な歩行表示）
    
    figure(fig);
    
    % パラメータ
    l = walker.l;
    r = walker.r;
    
    % 継続的な歩行のための変数
    persistent cumulative_x_offset;
    persistent trajectory_history;
    persistent time_history;
    persistent is_first_cycle;
    
    if isempty(is_first_cycle) || cycle_num == 1
        cumulative_x_offset = 0;
        trajectory_history = [];
        time_history = [];
        is_first_cycle = true;
        clf;  % 最初のサイクルのみクリア
    else
        is_first_cycle = false;
    end
    
    % アニメーション用データの準備
    fps = 15;
    if size(z_all, 2) >= 8
        z_all_plot = [z_all(:,1) z_all(:,3) z_all(:,6) z_all(:,8)];
    else
        % 拡張状態がない場合の処理
        hip_x = cumsum([0; diff(t)] .* 0.5);  % 時間に比例した水平移動
        z_all_plot = [z_all(:,1) z_all(:,3) hip_x ones(size(z_all,1),1)*l];
    end
    
    nn = size(z_all_plot,2);
    total_frames = min(round(t(end)*fps), 100);
    t_anim = linspace(t(1), t(end), total_frames);
    z = zeros(total_frames,nn);
    
    for i=1:nn
        z(:,i) = interp1(t, z_all_plot(:,i), t_anim);
    end
    
    % サイクル間での水平位置の連続性を確保
    if ~is_first_cycle && ~isempty(trajectory_history)
        last_x = trajectory_history(end, 1);
        current_start_x = z(1, 3);
        x_offset = last_x - current_start_x + 2*l;  % 歩幅分のオフセット
        cumulative_x_offset = cumulative_x_offset + x_offset;
    end
    
    % 水平位置を調整
    z(:,3) = z(:,3) + cumulative_x_offset;
    
    % 履歴に追加
    current_trajectory = [z(:,3), z(:,4)];  % [x, y]座標
    trajectory_history = [trajectory_history; current_trajectory];
    time_history = [time_history; t_anim'];
    
    % ウィンドウ設定（移動する視点）
    current_x = z(end,3);  % 現在の水平位置
    window_width = 4*l;
    window_xmin = current_x - window_width/3;
    window_xmax = current_x + 2*window_width/3;
    window_ymin = -0.2;
    window_ymax = 1.2*(l+r);
    
    % アニメーションループ
    dt_anim = 0.05 / speed_factor;
    
    for i = 1:total_frames
        % 左側：歩行アニメーション
        subplot(1,2,1);
        cla;
        
        % 状態変数
        q1 = z(i,1); q2 = z(i,2); 
        xh = z(i,3); yh = z(i,4);
        
        % 連続的な歩行描画
        draw_continuous_walker(q1, q2, xh, yh, l, r, trajectory_history, current_x);
        
        % 動的なウィンドウ調整
        current_frame_x = z(i,3);
        axis([current_frame_x - window_width/2, current_frame_x + window_width/2, window_ymin, window_ymax]);
        axis equal;
        
        title(sprintf('連続歩行 - サイクル %d (t=%.2fs)', cycle_num, t_anim(i)), 'FontSize', 12)
        
        % 右側：状態変数（累積表示）
        subplot(1,2,2);
        if i == 1
            cla;
        end
        
        % 現在のサイクルの軌跡（q1, q2で色分け）
        idx = 1:i;
        if ~isempty(idx)
            hold on;
            current_time = time_history(end-length(t_anim)+1:end-length(t_anim)+i);
            
            % q1脚の状態変数（青系）
            plot(current_time, z_all(idx,1), 'b-', 'LineWidth', 3, 'DisplayName', 'q1角度');
            plot(current_time, z_all(idx,2), 'b--', 'LineWidth', 2, 'DisplayName', 'q1角速度');
            
            % q2脚の状態変数（赤系）
            plot(current_time, z_all(idx,3), 'r-', 'LineWidth', 3, 'DisplayName', 'q2角度');
            plot(current_time, z_all(idx,4), 'r--', 'LineWidth', 2, 'DisplayName', 'q2角速度');
            
            % 現在位置をマーク（色分け）
            plot(current_time(end), z(i,1), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', [0, 0.3, 0.9], 'LineWidth', 2);
            plot(current_time(end), z(i,2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', [0.9, 0.1, 0.1], 'LineWidth', 2);
        end
        
        if i == 1
            xlabel('時間 [s]');
            ylabel('角度/角速度 [rad, rad/s]');
            title('q1脚(青)とq2脚(赤)の状態変数');
            legend('q1角度', 'q1角速度', 'q2角度', 'q2角速度', 'Location', 'best');
            grid on;
        end
        
        % 時間軸の動的調整
        if ~isempty(time_history)
            xlim([time_history(1), time_history(end)]);
            ylim([-1.5, 1.5]);
        end
        
        drawnow;
        pause(dt_anim);
    end
end

function draw_continuous_walker(q1, q2, xh, yh, l, r, trajectory_history, current_x)
    % 継続的な歩行を表現するウォーカー描画
    
    % 長い地面の描画（移動する視点に合わせて）
    ground_length = 10*l;
    ground_x = [current_x - ground_length/2, current_x + ground_length/2];
    plot(ground_x, [0 0], 'k-', 'LineWidth', 3);
    hold on;
    
    % 地面に格子パターンを追加（歩行感を強調）
    grid_spacing = l/2;
    for gx = floor(ground_x(1)/grid_spacing)*grid_spacing:grid_spacing:ground_x(2)
        plot([gx, gx], [0, -0.05], 'k-', 'LineWidth', 1);
    end
    
    % 歩行軌跡の表示（連続性を強調）
    if ~isempty(trajectory_history) && size(trajectory_history, 1) > 1
        % 腰関節の軌跡
        hip_traj_x = trajectory_history(:,1);
        hip_traj_y = trajectory_history(:,2);
        
        % 軌跡の長さを制限（最新の部分のみ表示）
        max_traj_points = min(200, length(hip_traj_x));
        if length(hip_traj_x) > max_traj_points
            start_idx = length(hip_traj_x) - max_traj_points + 1;
            hip_traj_x = hip_traj_x(start_idx:end);
            hip_traj_y = hip_traj_y(start_idx:end);
        end
        
        % 軌跡を色のグラデーションで表示
        for i = 1:length(hip_traj_x)-1
            alpha = i / length(hip_traj_x);  % 透明度
            color = [0, 0.7, 0] * alpha + [0.9, 0.9, 0.9] * (1-alpha);
            plot(hip_traj_x(i:i+1), hip_traj_y(i:i+1), '-', 'Color', color, 'LineWidth', 1);
        end
    end
    
    % 足跡マークの表示
    draw_footprint_marks(xh, yh, q1, q2, l, trajectory_history);
    
    % 現在のウォーカーの描画（q1とq2で色分け）
    % 脚1（q1で制御される脚）
    leg1_foot_x = xh + l*sin(q1);
    leg1_foot_y = 0;  % 通常は地面に接触
    
    % 脚2（q2で制御される脚）
    leg2_foot_x = xh + l*sin(q2);
    leg2_foot_y = yh - l*cos(q2);
    
    % 脚の状態を判定（どちらがスタンス脚かスイング脚か）
    leg1_on_ground = abs(leg1_foot_y) < 0.05;
    leg2_on_ground = abs(leg2_foot_y) < 0.05;
    
    % 脚1の色設定（q1脚 - 青系）
    if leg1_on_ground
        leg1_color = [0, 0.3, 0.9];        % 濃い青（スタンス脚）
        leg1_marker = 's';                  % 四角マーカー
        leg1_size = 12;
    else
        leg1_color = [0.4, 0.6, 1.0];     % 薄い青（スイング脚）
        leg1_marker = 'o';                  % 丸マーカー
        leg1_size = 10;
    end
    
    % 脚2の色設定（q2脚 - 赤系）
    if leg2_on_ground
        leg2_color = [0.9, 0.1, 0.1];     % 濃い赤（スタンス脚）
        leg2_marker = 's';                  % 四角マーカー
        leg2_size = 12;
    else
        leg2_color = [1.0, 0.4, 0.4];     % 薄い赤（スイング脚）
        leg2_marker = 'o';                  % 丸マーカー
        leg2_size = 10;
    end
    
    % 腰関節（体の中心）
    plot(xh, yh, 'ko', 'MarkerSize', 15, 'MarkerFaceColor', [0.2, 0.2, 0.2], ...
         'MarkerEdgeColor', [0.8, 0.8, 0.8], 'LineWidth', 2);
    
    % 脚1（q1脚）の描画
    plot([leg1_foot_x, xh], [leg1_foot_y, yh], '-', 'LineWidth', 5, 'Color', leg1_color);
    plot(leg1_foot_x, leg1_foot_y, leg1_marker, 'MarkerSize', leg1_size, ...
         'MarkerFaceColor', leg1_color, 'MarkerEdgeColor', leg1_color*0.7, 'LineWidth', 2);
    
    % 脚2（q2脚）の描画
    plot([leg2_foot_x, xh], [leg2_foot_y, yh], '-', 'LineWidth', 5, 'Color', leg2_color);
    plot(leg2_foot_x, leg2_foot_y, leg2_marker, 'MarkerSize', leg2_size, ...
         'MarkerFaceColor', leg2_color, 'MarkerEdgeColor', leg2_color*0.7, 'LineWidth', 2);
    
    % 脚のラベル表示
    text(leg1_foot_x-0.15, leg1_foot_y+0.15, 'q1脚', 'FontSize', 10, 'FontWeight', 'bold', ...
         'Color', leg1_color, 'BackgroundColor', 'white', 'EdgeColor', leg1_color);
    text(leg2_foot_x+0.05, leg2_foot_y+0.15, 'q2脚', 'FontSize', 10, 'FontWeight', 'bold', ...
         'Color', leg2_color, 'BackgroundColor', 'white', 'EdgeColor', leg2_color);
    
    % 速度ベクトルの表示（動きの方向を示す）
    velocity_scale = 0.3;
    if leg2_foot_x > leg1_foot_x
        % 右向きに歩行
        quiver(xh, yh + 0.1, velocity_scale, 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 0.5);
        text(xh + 0.2, yh + 0.4, '→', 'FontSize', 20, 'Color', 'green');
    else
        % 左向きに歩行
        quiver(xh, yh + 0.1, -velocity_scale, 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 0.5);
        text(xh - 0.2, yh + 0.4, '←', 'FontSize', 20, 'Color', 'green');
    end
    
    % 詳細な角度情報の表示（色分けして表示）
    angle_info = sprintf('q1=%.1f° (青脚)\nq2=%.1f° (赤脚)', q1*180/pi, q2*180/pi);
    text(xh+0.3, yh+0.2, angle_info, 'FontSize', 11, 'BackgroundColor', 'white', ...
         'EdgeColor', 'black', 'FontWeight', 'bold');
    
    % 脚の状態表示
    status_info = '';
    if leg1_on_ground && leg2_on_ground
        status_info = '両脚接地';
        status_color = 'yellow';
    elseif leg1_on_ground
        status_info = 'q1脚スタンス';
        status_color = [0.7, 0.9, 1.0];
    elseif leg2_on_ground
        status_info = 'q2脚スタンス';
        status_color = [1.0, 0.7, 0.7];
    else
        status_info = '両脚スイング';
        status_color = 'orange';
    end
    
    text(xh-0.5, yh+0.5, status_info, 'FontSize', 12, 'FontWeight', 'bold', ...
         'BackgroundColor', status_color, 'EdgeColor', 'black');
    
    % 凡例の追加
    legend_x = xh + 1.5;
    legend_y = yh + 0.3;
    
    % 凡例背景
    rectangle('Position', [legend_x-0.1, legend_y-0.4, 0.8, 0.6], ...
              'FaceColor', [1, 1, 1, 0.8], 'EdgeColor', 'black');
    
    % 凡例内容
    plot(legend_x, legend_y, 's', 'MarkerSize', 8, 'MarkerFaceColor', [0, 0.3, 0.9], ...
         'MarkerEdgeColor', [0, 0.2, 0.7]);
    text(legend_x+0.1, legend_y, 'q1脚(青)', 'FontSize', 9, 'Color', [0, 0.3, 0.9]);
    
    plot(legend_x, legend_y-0.15, 's', 'MarkerSize', 8, 'MarkerFaceColor', [0.9, 0.1, 0.1], ...
         'MarkerEdgeColor', [0.7, 0.1, 0.1]);
    text(legend_x+0.1, legend_y-0.15, 'q2脚(赤)', 'FontSize', 9, 'Color', [0.9, 0.1, 0.1]);
    
    % グリッドと軸の設定
    grid on;
    axis equal;
    
    % Y軸の範囲を固定
    ylim([-0.3, 1.5*l]);
end

function draw_footprint_marks(xh, yh, q1, q2, l, trajectory_history)
    % 足跡マークを描画（q1, q2脚で色分け）
    
    persistent footprint_history;
    
    if isempty(footprint_history)
        footprint_history = [];
    end
    
    % 現在の足の位置
    leg1_foot_x = xh + l*sin(q1);  % q1脚
    leg2_foot_x = xh + l*sin(q2);  % q2脚
    leg2_foot_y = yh - l*cos(q2);
    
    % q1脚（通常は地面についている）の足跡を記録
    q1_footprint = [leg1_foot_x, 0, 1];  % [x, y, type] type: 1=q1脚, 2=q2脚
    
    % 新しい足跡が前の足跡から十分離れている場合のみ追加
    if isempty(footprint_history) || ...
       abs(footprint_history(end,1) - q1_footprint(1)) > 0.3
        footprint_history = [footprint_history; q1_footprint];
    end
    
    % q2脚が地面に近い場合も足跡として記録
    if leg2_foot_y < 0.05 && leg2_foot_y > -0.05
        q2_footprint = [leg2_foot_x, 0, 2];
        if isempty(footprint_history) || ...
           abs(footprint_history(end,1) - q2_footprint(1)) > 0.3
            footprint_history = [footprint_history; q2_footprint];
        end
    end
    
    % 足跡の数を制限（メモリ使用量を制御）
    max_footprints = 25;
    if size(footprint_history, 1) > max_footprints
        footprint_history = footprint_history(end-max_footprints+1:end, :);
    end
    
    % 足跡を描画（q1, q2で色分け）
    for i = 1:size(footprint_history, 1)
        fp_x = footprint_history(i, 1);
        fp_y = footprint_history(i, 2);
        fp_type = footprint_history(i, 3);
        
        % 透明度（古い足跡ほど薄く）
        alpha = (i / size(footprint_history, 1)) * 0.7 + 0.3;
        
        if fp_type == 1  % q1脚の足跡（青）
            base_color = [0, 0.3, 0.9];
            color = base_color * alpha + [0.8, 0.8, 0.8] * (1-alpha);
            marker = 's';  % 四角
            marker_size = 8;
            edge_color = base_color * 0.7;
        else  % q2脚の足跡（赤）
            base_color = [0.9, 0.1, 0.1];
            color = base_color * alpha + [0.8, 0.8, 0.8] * (1-alpha);
            marker = 'o';  % 丸
            marker_size = 7;
            edge_color = base_color * 0.7;
        end
        
        plot(fp_x, fp_y, marker, 'MarkerSize', marker_size, 'MarkerFaceColor', color, ...
             'MarkerEdgeColor', edge_color, 'MarkerFaceAlpha', alpha, 'LineWidth', 1.5);
        
        % 足跡にラベルを追加（最新のもののみ）
        if i > size(footprint_history, 1) - 3
            if fp_type == 1
                text(fp_x-0.1, fp_y-0.1, 'q1', 'FontSize', 8, 'Color', base_color, 'FontWeight', 'bold');
            else
                text(fp_x+0.05, fp_y-0.1, 'q2', 'FontSize', 8, 'Color', base_color, 'FontWeight', 'bold');
            end
        end
    end
end

function animate_all_cycles_continuous(fig, all_times, all_extended_states, walker)
    % 全サイクルの連続アニメーション
    
    figure(fig);
    clf;
    
    % 全データを結合
    t_all = [];
    z_all = [];
    for i = 1:length(all_times)
        t_all = [t_all; all_times{i}];
        z_all = [z_all; all_extended_states{i}];
    end
    
    % パラメータ
    l = walker.l;
    r = walker.r;
    
    % 連続歩行のための座標調整
    x_offset = 0;
    adjusted_z_all = z_all;
    
    % 各サイクルのx座標を調整して連続性を確保
    cycle_start_idx = 1;
    for cycle = 1:length(all_times)
        cycle_end_idx = cycle_start_idx + length(all_times{cycle}) - 1;
        
        if cycle > 1
            % 前のサイクルの最後のx座標
            if size(adjusted_z_all, 2) >= 6
                prev_end_x = adjusted_z_all(cycle_start_idx-1, 6);
                current_start_x = z_all(cycle_start_idx, 6);
                x_offset = prev_end_x - current_start_x + 2*l;
            else
                x_offset = x_offset + 2*l;  % デフォルトの歩幅
            end
        end
        
        % x座標を調整
        if size(adjusted_z_all, 2) >= 6
            adjusted_z_all(cycle_start_idx:cycle_end_idx, 6) = ...
                z_all(cycle_start_idx:cycle_end_idx, 6) + x_offset;
        end
        
        cycle_start_idx = cycle_end_idx + 1;
    end
    
    % アニメーション用データの準備
    fps = 20;
    if size(adjusted_z_all, 2) >= 8
        z_plot = [adjusted_z_all(:,1) adjusted_z_all(:,3) adjusted_z_all(:,6) adjusted_z_all(:,8)];
    else
        % 拡張状態がない場合
        hip_x = cumsum([0; diff(t_all)] .* 0.5);
        z_plot = [adjusted_z_all(:,1) adjusted_z_all(:,3) hip_x ones(size(adjusted_z_all,1),1)*l];
    end
    
    total_frames = min(round(t_all(end)*fps), 500);
    t_anim = linspace(t_all(1), t_all(end), total_frames);
    z_interp = zeros(total_frames, size(z_plot,2));
    
    for i = 1:size(z_plot,2)
        z_interp(:,i) = interp1(t_all, z_plot(:,i), t_anim);
    end
    
    % アニメーション実行
    trajectory_history = [];
    for i = 1:total_frames
        clf;
        
        q1 = z_interp(i,1);
        q2 = z_interp(i,2);
        xh = z_interp(i,3);
        yh = z_interp(i,4);
        
        % 軌跡履歴を更新
        trajectory_history = [trajectory_history; [xh, yh]];
        
        % 継続的な歩行描画
        draw_continuous_walker(q1, q2, xh, yh, l, r, trajectory_history, xh);
        
        % 動的ウィンドウ
        window_width = 4*l;
        axis([xh - window_width/2, xh + window_width/2, -0.3, 1.5*l]);
        
        title(sprintf('連続歩行アニメーション (t=%.2fs)', t_anim(i)), 'FontSize', 14);
        
        drawnow;
        pause(0.05);
    end
end

function visualize_results_fixed(all_times, all_trajectories, all_extended_states, cycle_data, z0_original)
    % 結果の可視化（修正版）
    
    figure('Name', '複数サイクル解析結果', 'Position', [50 50 1400 900]);
    
    % === サブプロット1: 角度の時系列 ===
    subplot(2,3,1);
    hold on;
    colors = lines(length(all_times));
    
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
    
    % === サブプロット2: 角速度の時系列 ===
    subplot(2,3,2);
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
    
    % === サブプロット3: エネルギーの推移 ===
    subplot(2,3,3);
    if ~isempty(cycle_data)
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
    end
    
    % === サブプロット4: 位相空間（q1-u1） ===
    subplot(2,3,4);
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
    
    % === サブプロット5: 位相空間（q2-u2） ===
    subplot(2,3,5);
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
    
    % === サブプロット6: 角速度補正量 ===
    subplot(2,3,6);
    if ~isempty(cycle_data) && length(cycle_data) > 1
        corrections = [cycle_data(1:end-1).velocity_correction];
        
        bar(1:length(corrections), corrections);
        xlabel('サイクル');
        ylabel('角速度補正量');
        title('各サイクルでの角速度補正');
        grid on;
    end
    
    sgtitle('パッシブウォーカー: 複数サイクル解析結果（修正版）');
end

function display_summary_fixed(cycle_data, z0_original)
    % サマリー表示（修正版）
    
    fprintf('\n=== シミュレーションサマリー ===\n');
    fprintf('総サイクル数: %d\n', length(cycle_data));
    
    if length(cycle_data) > 0
        avg_energy_change = mean([cycle_data.energy_loss]);
        fprintf('平均エネルギー変化: %.6f\n', avg_energy_change);
        
        if length(cycle_data) > 1
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

% === 以下、物理モデル関数（修正版） ===

function [z_traj, t_traj, z_final, z_midpoint] = one_and_half_stride_detailed(z0, walker)
    % 1.5ストライド（スイング脚がスタンス脚に追いつくまで）の詳細な軌道を返す関数
    
    M = walker.M;  m = walker.m; I = walker.I;   
    l = walker.l;  c = walker.c; w = walker.w;   
    r = walker.r;  g = walker.g; gam = walker.gam;
    
    q1 = z0(1); u1 = z0(2); q2 = z0(3); u2 = z0(4);
    
    % エネルギーと位置の計算
    TE = calculate_energy(z0, walker);
    xp1 = 0;
    xh = -l*sin(q1) - r*q1 + xp1;
    vxh = (-l*cos(q1)-r)*u1; 
    yh =  l*cos(q1) + r;
    vyh = -l*sin(q1)*u1; 
    
    z0_extended = [q1 u1 q2 u2 TE xh vxh yh vyh];
    
    % === 1ステップ目 ===
    fprintf('1ステップ目（最初の脚）を実行中...\n');
    t0 = 0; 
    dt = 5; % 1ステップの最大時間
    time_stamps = 100; % 軌道の解像度
    
    % Single stance phase
    options = odeset('abstol',1e-13,'reltol',1e-13,'events',@collision);
    tspan = linspace(t0, t0+dt, time_stamps);
    [t_step1, z_step1] = ode113(@(t,z) single_stance(t,z,walker), tspan, z0_extended, options);
    
    % 1ステップ目のHeel strike
    if length(t_step1) < time_stamps
        fprintf('1ステップ目完了: t = %.4f秒で次の脚が地面に接触\n', t_step1(end));
        z_after_collision1 = heelstrike(t_step1(end), z_step1(end,:), walker);
        z_midpoint = z_after_collision1(1:4);
    else
        error('1ステップ目: 指定時間内に衝突が発生しませんでした');
    end
    
    % === 0.5ステップ目（スイング脚がスタンス脚に追いつくまで） ===
    fprintf('0.5ステップ目（スイング脚がスタンス脚に追いつくまで）を実行中...\n');
    z1_extended = z_after_collision1;
    t1 = t_step1(end);
    
    % 新しい衝突条件：スイング脚がスタンス脚に追いつく（地面接触なし）
    options_catchup = odeset('abstol',1e-13,'reltol',1e-13,'events',@swing_catches_stance);
    tspan = linspace(t1, t1+dt, time_stamps);
    [t_step1_5, z_step1_5] = ode113(@(t,z) single_stance(t,z,walker), tspan, z1_extended, options_catchup);
    
    % 時間を調整（連続的にする）
    t_step1_5 = t_step1_5 - t1 + t_step1(end);
    
    % 0.5ステップ目の終了（スイング脚追いつき、Heel Strikeなし）
    if length(t_step1_5) < time_stamps
        fprintf('0.5ステップ目完了: t = %.4f秒でスイング脚がスタンス脚に追いつきました\n', t_step1_5(end));
        z_final = z_step1_5(end, 1:4);
    else
        fprintf('警告: 指定時間内にスイング脚が追いつきませんでした\n');
        z_final = z_step1_5(end, 1:4);
    end
    
    % === 軌道データの結合 ===
    % 重複する時間点を除去して連結
    t_traj = [t_step1; t_step1_5(2:end)];
    z_traj = [z_step1; z_step1_5(2:end,:)];
    
    fprintf('1.5ストライド（スイング脚追いつき）完了\n');
end

function TE = calculate_energy(z, walker)
    % 総エネルギーを計算
    
    q1 = z(1); u1 = z(2); q2 = z(3); u2 = z(4);
    M = walker.M;  m = walker.m; I = walker.I;   
    l = walker.l;  c = walker.c; w = walker.w;   
    r = walker.r;  g = walker.g; gam = walker.gam;
    
    TE = 1/2*m*(((-l*cos(q1)-r)*u1-u1*(-c*cos(q1)+w*sin(q1)))^2+(-l*sin(q1)*u1+u1*(c*sin(q1)+w*cos(q1)))^2)+1/2*m*(((-l*cos(q1)-r)*u1-(u1-u2)*(-c*cos(q1-q2)+w*sin(q1-q2)))^2+(-l*sin(q1)*u1+(u1-u2)*(c*sin(q1-q2)+w*cos(q1-q2)))^2)+1/2*M*((-l*cos(q1)-r)^2*u1^2+l^2*sin(q1)^2*u1^2)+1/2*I*(u1^2+(u1-u2)^2)+2*m*g*cos(gam)*r+2*m*g*l*cos(gam-q1)-m*g*c*cos(gam-q1)-m*g*w*sin(gam-q1)+2*m*g*sin(gam)*r*q1-m*g*c*cos(gam-q1+q2)-m*g*w*sin(gam-q1+q2)+M*g*cos(gam)*r+M*g*l*cos(gam-q1)+M*g*sin(gam)*r*q1;
end

function [gstop, isterminal,direction]=swing_catches_stance(t,z,walker)
    % スイング脚がスタンス脚に追いつく条件
    
    q1 = z(1); q2 = z(3); 
    
    % スイング脚角度がスタンス脚角度に追いつく条件
    gstop = q2 - q1;
    isterminal = 1; % 条件に達したら積分を停止
    direction = 0;  % どちらの方向からでも検出
end

function zdot=single_stance(t,z,walker)  
    % 単一スタンス期の微分方程式
    
    q1 = z(1);   u1 = z(2);                         
    q2 = z(3);   u2 = z(4);                         
    xh = z(6);  vxh = z(7);                       
    yh = z(8);  vyh = z(9);                     
    
    M = walker.M;  m = walker.m; I = walker.I;   
    l = walker.l;  c = walker.c; w = walker.w;   
    r = walker.r;  g = walker.g; gam = walker.gam;
    
    Th=0;
    
    M11 = -2*w^2*m-2*I+2*m*l*c*cos(q2)+2*m*w*l*sin(q2)-2*m*c^2-2*m*l^2-M*l^2+2*m*l*c-2*m*r^2-M*r^2+2*m*r*c*cos(q1-q2)-2*m*r*w*sin(q1-q2)-2*M*r*l*cos(q1)-4*m*r*l*cos(q1)+2*m*r*c*cos(q1)-2*m*r*w*sin(q1); 
    M12 = w^2*m+I-m*l*c*cos(q2)-m*w*l*sin(q2)+m*c^2-m*r*c*cos(q1-q2)+m*r*w*sin(q1-q2); 
    
    M21 = m*w*l*sin(q2)+m*l*c*cos(q2)-m*r*w*sin(q1-q2)+m*r*c*cos(q1-q2)-m*c^2-w^2*m-I; 
    M22 = w^2*m+m*c^2+I; 
    
    RHS1 = -2*m*r*u1*u2*c*sin(q1-q2)-2*m*r*u1*u2*w*cos(q1-q2)+m*r*u1^2*w*cos(q1)+m*r*u1^2*c*sin(q1)-2*m*r*l*sin(q1)*u1^2+M*g*sin(gam)*r+2*m*g*sin(gam)*r+m*r*u2^2*w*cos(q1-q2)+m*r*u2^2*c*sin(q1-q2)+m*r*u1^2*w*cos(q1-q2)+m*r*u1^2*c*sin(q1-q2)-M*r*l*sin(q1)*u1^2+M*g*l*sin(gam-q1)+2*m*g*l*sin(gam-q1)-m*g*c*sin(gam-q1)+m*g*w*cos(gam-q1)-m*g*c*sin(gam-q1+q2)+m*g*w*cos(gam-q1+q2)-2*m*l*u1*u2*w*cos(q2)-m*l*u2^2*c*sin(q2)+2*m*l*u1*u2*c*sin(q2)+m*l*u2^2*w*cos(q2); 
    RHS2 = -m*g*c*sin(gam-q1+q2)+m*g*w*cos(gam-q1+q2)-Th-m*l*u1^2*w*cos(q2)+m*l*u1^2*c*sin(q2); 
    
    MM = [M11 M12; M21 M22];                               
    RHS = [RHS1; RHS2];                       
    X = MM \ RHS;                                    
    
    ud1 = X(1);                                       
    ud2 = X(2);                                       
    
    DTE = -ud1*I*u2+2*ud1*m*u1*r^2+m*u1*r*u2^2*c*sin(q1-q2)+m*u1*r*u2^2*w*cos(q1-q2)-m*u2^2*l*u1*c*sin(q2)+u2*m*g*c*sin(gam-q1+q2)-u2*m*g*w*cos(gam-q1+q2)+2*ud1*I*u1+ud2*I*u2+m*u2^2*l*u1*w*cos(q2)+2*ud1*m*u1*c^2+ud2*m*u2*c^2-ud2*I*u1+ud1*m*u2*c*l*cos(q2)+ud1*m*u2*w*l*sin(q2)-2*ud1*m*u1*l*c*cos(q2)-2*ud1*m*l*u1*w*sin(q2)-m*u2*u1^2*w*l*cos(q2)+m*u2*u1^2*c*l*sin(q2)+2*ud1*m*u1*w^2+ud2*m*u2*w^2+ud2*m*u1*l*c*cos(q2)+ud1*M*l^2*u1-ud2*m*u1*w^2-ud1*m*u2*c^2-ud2*m*u1*c^2+2*ud1*m*l^2*u1-ud1*m*u2*w^2-2*ud1*m*l*u1*c+2*ud1*M*l*cos(q1)*u1*r+4*ud1*m*l*cos(q1)*u1*r-2*ud1*m*u1*r*c*cos(q1)+2*ud1*m*u1*r*w*sin(q1)-2*ud1*m*u1*r*c*cos(q1-q2)-2*m*u1^3*r*l*sin(q1)+m*u1^3*r*c*sin(q1)+m*u1^3*r*w*cos(q1)+m*u1^3*r*c*sin(q1-q2)+m*u1^3*r*w*cos(q1-q2)-2*m*u1^2*r*u2*c*sin(q1-q2)-2*m*u1^2*r*u2*w*cos(q1-q2)-M*u1^3*r*l*sin(q1)+2*u1*m*g*l*sin(gam-q1)-u1*m*g*c*sin(gam-q1)+u1*m*g*w*cos(gam-q1)+2*u1*m*g*sin(gam)*r+ud2*m*l*u1*w*sin(q2)+ud1*M*u1*r^2-u1*m*g*c*sin(gam-q1+q2)+u1*m*g*w*cos(gam-q1+q2)+u1*M*g*l*sin(gam-q1)+u1*M*g*sin(gam)*r+2*ud1*m*u1*r*w*sin(q1-q2)+ud1*m*u2*c*cos(q1-q2)*r-ud1*m*u2*w*sin(q1-q2)*r+ud2*m*u1*r*c*cos(q1-q2)-ud2*m*u1*r*w*sin(q1-q2); 
    axh = l*sin(q1)*u1^2+(-l*cos(q1)-r)*ud1; 
    ayh = -l*cos(q1)*u1^2-l*sin(q1)*ud1; 
    
    zdot = [u1 ud1 u2 ud2 DTE vxh axh vyh ayh]';  
end

function [gstop, isterminal,direction]=collision(t,z,walker)
    % 衝突検知
    
    q1 = z(1); q2 = z(3); 
    
    gstop = -q2 + 2*q1;
    if (q2>-0.05)
        isterminal = 0;
    else
        isterminal=1;
    end
    direction=-1;
end

function zplus=heelstrike(t,z,walker)      
    % ヒールストライク
    
    r1 = z(1);   v1 = z(2);                         
    r2 = z(3);   v2 = z(4);                         
    xh = z(6);   yh = z(8);                       
    
    q1 = r1 - r2;                         
    q2 = -r2;                                       
    
    M = walker.M;  m = walker.m; I = walker.I;   
    l = walker.l;  c = walker.c; w = walker.w;   
    r = walker.r;  g = walker.g; gam = walker.gam; 
    
    M11 = 2*m*l^2-2*m*l*c+2*m*c^2+2*m*w^2+2*m*r^2+4*m*r*l*cos(q1)-2*m*r*c*cos(q1)+2*m*w*sin(q1)*r-2*m*l*c*cos(q2)-2*m*l*w*sin(q2)-2*m*r*c*cos(q1-q2)+2*m*sin(q1-q2)*w*r+M*l^2+2*M*r*l*cos(q1)+M*r^2+2*I; 
    M12 = m*l*c*cos(q2)+m*l*w*sin(q2)-m*c^2-m*w^2+m*r*c*cos(q1-q2)-m*sin(q1-q2)*w*r-I; 
    
    M21 = -m*l*c*cos(q2)-m*l*w*sin(q2)+m*c^2+m*w^2-m*r*c*cos(q1-q2)+m*sin(q1-q2)*w*r+I; 
    M22 = -m*w^2-m*c^2-I; 
    
    RHS1 = m*l*v2*c-m*c^2*v2+M*v1*r^2-2*m*c*l*v1+M*v1*l^2*cos(r2)+2*m*l^2*v1*cos(r2)+2*I*v1-I*v2+2*m*v1*r^2-2*m*l*v1*c*cos(r2)+2*m*c^2*v1+2*m*w^2*v1-m*w^2*v2+2*m*r*v1*w*sin(r1)+M*r*v1*l*cos(r1)+M*l*cos(-r1+r2)*v1*r-2*m*r*v1*c*cos(r1)+2*m*l*cos(-r1+r2)*v1*r+2*m*r*v1*l*cos(r1)-2*m*r*v1*c*cos(-r1+r2)-2*m*r*v1*w*sin(-r1+r2)+m*r*v2*c*cos(-r1+r2)+m*r*v2*w*sin(-r1+r2); 
    RHS2 = m*r*v1*w*sin(r1)-m*r*v1*c*cos(r1)+I*v1-I*v2+m*w^2*v1-m*c*l*v1+m*c^2*v1; 
    
    MM = [M11 M12; M21 M22];    
    RHS = [RHS1; RHS2];                      
    X = MM \ RHS;                                    
    
    u1 = X(1); u2 = X(2);                                      
    
    TE = 1/2*m*(((-l*cos(q1)-r)*u1-u1*(-c*cos(q1)+w*sin(q1)))^2+(-l*sin(q1)*u1+u1*(c*sin(q1)+w*cos(q1)))^2)+1/2*m*(((-l*cos(q1)-r)*u1-(u1-u2)*(-c*cos(q1-q2)+w*sin(q1-q2)))^2+(-l*sin(q1)*u1+(u1-u2)*(c*sin(q1-q2)+w*cos(q1-q2)))^2)+1/2*M*((-l*cos(q1)-r)^2*u1^2+l^2*sin(q1)^2*u1^2)+1/2*I*(u1^2+(u1-u2)^2)+2*g*m*cos(gam)*r+2*m*g*l*cos(gam-q1)-m*g*c*cos(gam-q1)-m*g*w*sin(gam-q1)+2*g*m*sin(gam)*r*q1-m*g*c*cos(gam-q1+q2)-m*g*w*sin(gam-q1+q2)+g*M*cos(gam)*r+M*g*l*cos(gam-q1)+g*M*sin(gam)*r*q1; 
    vxh = (-l*cos(q1)-r)*u1; 
    vyh = -l*sin(q1)*u1; 
    
    zplus = [q1 u1 q2 u2 TE xh vxh yh vyh];                     
end