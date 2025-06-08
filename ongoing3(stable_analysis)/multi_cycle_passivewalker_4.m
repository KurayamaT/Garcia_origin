function multi_cycle_passivewalker_corrected()
    % 元の物理モデルを使用した複数サイクルパッシブウォーカーシミュレーション
    % single_stride_passivewalker4.mと完全に同じ計算ロジックを使用
    % アニメーション機能付き
    
    clear; close all;
    
    % === シミュレーション設定 ===
    fprintf('=== 元の物理モデルによる複数サイクルシミュレーション ===\n');
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
        walker.c = 1.0;  walker.r = 0.0; walker.g = 1.0; walker.gam = 0;%0.009; 
        
        % 固定点
        zstar = [0.200161072169750; -0.199906060087682; 0.400322144339512; -0.015805473227965];
        
    else 
        %%  More General round feet walker
        walker.M = 1.0; walker.m = 0.5; walker.I = 0.02; walker.l = 1.0; walker.w = 0.0; 
        walker.c = 0.5; walker.r = 0.2; walker.g = 1.0; walker.gam = 0;%0.01; 
        
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
                    0.1, -0.3, 0.2, -0.3;
                    -0.01, -0.1, -0.01, -0.3;
                ];
            else
                presets = [
                    0.2, -0.2, 0.4, -0.3;
                    0.1, -0.1, 0.1, -0.3;
                    -0.01, -0.1, -0.01, -0.3;
                ];
            end
            fprintf('\nプリセット選択:\n');
            for i = 1:size(presets, 1)
                fprintf('%d: [%.2f, %.2f, %.2f, %.2f]\n', i, presets(i,:));
            end
            preset_idx = input('選択: ');
            if isempty(preset_idx) || preset_idx < 1 || preset_idx > size(presets, 1)
                preset_idx = 1;
            end
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
            % 1.5歩のシミュレーション実行（元の物理モデルを使用）
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
            
            % サイクル情報を記録
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
            
            % リアルタイムアニメーション表示
            try
                animate_cycle_original_fixed(fig_anim, t_traj, z_traj, walker, cycle, anim_speed_factor);
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
            try
                animate_all_cycles_original_fixed(fig_anim, all_times, all_extended_states, walker, anim_speed_factor);
            catch me
                fprintf('連続アニメーションエラー: %s\n', me.message);
            end
        end
        
        % メイン関数の最後に、CSV保存機能の呼び出しを追加
        % 結果の可視化
        visualize_results_original(all_times, all_trajectories, all_extended_states, cycle_data, z0_original);
        display_summary(cycle_data, z0_original);
        
        % CSV保存機能を追加
        save_simulation_data_to_csv(all_times, all_trajectories, all_extended_states, cycle_data, z0_original, walker);

        % === 歩行安定性の判定 ===
        fprintf('\n歩行安定性を解析しますか？ (y/n): ');
        user_input = input('', 's');
        if strcmpi(user_input, 'y') || strcmpi(user_input, 'yes')
            stability_result = analyze_walking_stability(all_trajectories, cycle_data, z0_original, walker);
            
            % 安定性判定結果をCSVに保存
            if exist('save_dir', 'var')
                save_stability_analysis_to_csv(save_dir, stability_result, cycle_data);
            end
        end




    end
end

function save_simulation_data_to_csv(all_times, all_trajectories, all_extended_states, cycle_data, z0_original, walker)
    % シミュレーションデータをCSV形式で保存
    
    fprintf('\nデータをCSV形式で保存しますか？ (y/n): ');
    user_input = input('', 's');
    if ~strcmpi(user_input, 'y') && ~strcmpi(user_input, 'yes')
        return;
    end
    
    % 保存ディレクトリの作成
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    save_dir = sprintf('passive_walker_data_%s', timestamp);
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
    
    fprintf('データを保存中... ディレクトリ: %s\n', save_dir);
    
    % 1. ウォーカーパラメータの保存
    save_walker_parameters(save_dir, walker, z0_original);
    
    % 2. 各サイクルの軌道データの保存
    save_trajectory_data(save_dir, all_times, all_trajectories, all_extended_states);
    
    % 3. サイクル統計情報の保存
    save_cycle_statistics(save_dir, cycle_data);
    
    % 4. 統合された時系列データの保存
    save_combined_time_series(save_dir, all_times, all_trajectories, all_extended_states);
    
    % 5. エネルギー解析データの保存
    save_energy_analysis(save_dir, all_times, all_extended_states, cycle_data);
    
    fprintf('データ保存完了: %s\n', save_dir);
    fprintf('保存されたファイル:\n');
    fprintf('  - walker_parameters.csv (ウォーカーパラメータ)\n');
    fprintf('  - cycle_statistics.csv (サイクル統計)\n');
    fprintf('  - combined_time_series.csv (統合時系列データ)\n');
    fprintf('  - energy_analysis.csv (エネルギー解析)\n');
    fprintf('  - trajectory_cycle_*.csv (各サイクルの軌道データ)\n');
end

function save_walker_parameters(save_dir, walker, z0_original)
    % ウォーカーパラメータと初期条件をCSVで保存
    
    filename = fullfile(save_dir, 'walker_parameters.csv');
    
    % パラメータ名と値のペア
    param_names = {'M', 'm', 'I', 'l', 'c', 'w', 'r', 'g', 'gam'};
    param_values = [walker.M, walker.m, walker.I, walker.l, walker.c, ...
                    walker.w, walker.r, walker.g, walker.gam];
    
    % 初期条件
    initial_names = {'q1_initial', 'u1_initial', 'q2_initial', 'u2_initial'};
    initial_values = z0_original';
    
    % CSVファイルの作成
    fid = fopen(filename, 'w');
    fprintf(fid, 'Parameter,Value,Description\n');
    
    % ウォーカーパラメータ
    descriptions = {'Hip mass [kg]', 'Leg mass [kg]', 'Leg inertia [kg*m^2]', ...
                    'Leg length [m]', 'Leg CoM distance [m]', 'Leg width [m]', ...
                    'Foot radius [m]', 'Gravity [m/s^2]', 'Slope angle [rad]'};
    
    for i = 1:length(param_names)
        fprintf(fid, '%s,%.10f,%s\n', param_names{i}, param_values(i), descriptions{i});
    end
    
    % 初期条件
    init_descriptions = {'Initial stance leg angle [rad]', 'Initial stance leg velocity [rad/s]', ...
                        'Initial swing leg angle [rad]', 'Initial swing leg velocity [rad/s]'};
    
    for i = 1:length(initial_names)
        fprintf(fid, '%s,%.10f,%s\n', initial_names{i}, initial_values(i), init_descriptions{i});
    end
    
    fclose(fid);
end

function save_trajectory_data(save_dir, all_times, all_trajectories, all_extended_states)
    % 各サイクルの軌道データを個別のCSVファイルで保存
    
    for cycle = 1:length(all_times)
        filename = fullfile(save_dir, sprintf('trajectory_cycle_%d.csv', cycle));
        
        t = all_times{cycle};
        z_basic = all_trajectories{cycle};
        
        fid = fopen(filename, 'w');
        
        if size(all_extended_states{cycle}, 2) >= 9
            % 拡張状態がある場合
            z_extended = all_extended_states{cycle};
            fprintf(fid, 'time,q1,u1,q2,u2,total_energy,hip_x,hip_vx,hip_y,hip_vy\n');
            
            for i = 1:length(t)
                fprintf(fid, '%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n', ...
                        t(i), z_extended(i,1), z_extended(i,2), z_extended(i,3), z_extended(i,4), ...
                        z_extended(i,5), z_extended(i,6), z_extended(i,7), z_extended(i,8), z_extended(i,9));
            end
        else
            % 基本状態のみの場合
            fprintf(fid, 'time,q1,u1,q2,u2\n');
            
            for i = 1:length(t)
                fprintf(fid, '%.10f,%.10f,%.10f,%.10f,%.10f\n', ...
                        t(i), z_basic(i,1), z_basic(i,2), z_basic(i,3), z_basic(i,4));
            end
        end
        
        fclose(fid);
    end
end

function save_cycle_statistics(save_dir, cycle_data)
    % サイクル統計情報をCSVで保存
    
    filename = fullfile(save_dir, 'cycle_statistics.csv');
    
    fid = fopen(filename, 'w');
    fprintf(fid, 'cycle,duration,energy_initial,energy_final,energy_loss,');
    fprintf(fid, 'q1_initial,u1_initial,q2_initial,u2_initial,');
    fprintf(fid, 'q1_final,u1_final,q2_final,u2_final,');
    fprintf(fid, 'q1_midpoint,u1_midpoint,q2_midpoint,u2_midpoint,');
    fprintf(fid, 'velocity_correction\n');
    
    for i = 1:length(cycle_data)
        cycle = cycle_data(i);
        
        % 角速度補正量（最後のサイクルは0）
        if isfield(cycle, 'velocity_correction')
            vel_correction = cycle.velocity_correction;
        else
            vel_correction = 0;
        end
        
        fprintf(fid, '%d,%.10f,%.10f,%.10f,%.10f,', ...
                cycle.number, cycle.duration, cycle.energy_initial, ...
                cycle.energy_final, cycle.energy_loss);
        
        fprintf(fid, '%.10f,%.10f,%.10f,%.10f,', ...
                cycle.initial_state(1), cycle.initial_state(2), ...
                cycle.initial_state(3), cycle.initial_state(4));
        
        fprintf(fid, '%.10f,%.10f,%.10f,%.10f,', ...
                cycle.final_state(1), cycle.final_state(2), ...
                cycle.final_state(3), cycle.final_state(4));
        
        fprintf(fid, '%.10f,%.10f,%.10f,%.10f,', ...
                cycle.midpoint_state(1), cycle.midpoint_state(2), ...
                cycle.midpoint_state(3), cycle.midpoint_state(4));
        
        fprintf(fid, '%.10f\n', vel_correction);
    end
    
    fclose(fid);
end

function save_combined_time_series(save_dir, all_times, all_trajectories, all_extended_states)
    % 全サイクルを統合した時系列データをCSVで保存
    
    filename = fullfile(save_dir, 'combined_time_series.csv');
    
    % データの結合
    t_combined = [];
    z_basic_combined = [];
    z_extended_combined = [];
    cycle_labels = [];
    
    for cycle = 1:length(all_times)
        t_combined = [t_combined; all_times{cycle}];
        z_basic_combined = [z_basic_combined; all_trajectories{cycle}];
        
        if ~isempty(all_extended_states{cycle})
            z_extended_combined = [z_extended_combined; all_extended_states{cycle}];
        end
        
        cycle_labels = [cycle_labels; cycle * ones(length(all_times{cycle}), 1)];
    end
    
    fid = fopen(filename, 'w');
    
    if size(z_extended_combined, 2) >= 9
        % 拡張状態がある場合
        fprintf(fid, 'cycle,time,q1,u1,q2,u2,total_energy,hip_x,hip_vx,hip_y,hip_vy\n');
        
        for i = 1:length(t_combined)
            fprintf(fid, '%d,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n', ...
                    cycle_labels(i), t_combined(i), z_extended_combined(i,1), z_extended_combined(i,2), ...
                    z_extended_combined(i,3), z_extended_combined(i,4), z_extended_combined(i,5), ...
                    z_extended_combined(i,6), z_extended_combined(i,7), z_extended_combined(i,8), z_extended_combined(i,9));
        end
    else
        % 基本状態のみの場合
        fprintf(fid, 'cycle,time,q1,u1,q2,u2\n');
        
        for i = 1:length(t_combined)
            fprintf(fid, '%d,%.10f,%.10f,%.10f,%.10f,%.10f\n', ...
                    cycle_labels(i), t_combined(i), z_basic_combined(i,1), z_basic_combined(i,2), ...
                    z_basic_combined(i,3), z_basic_combined(i,4));
        end
    end
    
    fclose(fid);
end

function save_energy_analysis(save_dir, all_times, all_extended_states, cycle_data)
    % エネルギー解析データをCSVで保存
    
    filename = fullfile(save_dir, 'energy_analysis.csv');
    
    fid = fopen(filename, 'w');
    fprintf(fid, 'cycle,time_point,kinetic_energy,potential_energy,total_energy,');
    fprintf(fid, 'energy_rate_of_change,cumulative_energy_loss\n');
    
    cumulative_loss = 0;
    
    for cycle = 1:length(all_times)
        if cycle <= length(cycle_data)
            cumulative_loss = cumulative_loss + cycle_data(cycle).energy_loss;
        end
        
        if size(all_extended_states{cycle}, 2) >= 5
            t = all_times{cycle};
            total_energy = all_extended_states{cycle}(:, 5);
            
            % エネルギー変化率の計算
            energy_rate = [0; diff(total_energy) ./ diff(t)];
            
            for i = 1:length(t)
                % 運動エネルギーと位置エネルギーは総エネルギーから推定
                % （正確な分離には追加計算が必要）
                ke_estimate = total_energy(i) * 0.6;  % 推定値
                pe_estimate = total_energy(i) * 0.4;  % 推定値
                
                fprintf(fid, '%d,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n', ...
                        cycle, t(i), ke_estimate, pe_estimate, total_energy(i), ...
                        energy_rate(i), cumulative_loss);
            end
        end
    end
    
    fclose(fid);
end

function z_corrected = apply_velocity_correction(z_final, z0_original, correction_method, cycle_info)
    % 角速度補正を適用
    
    switch correction_method
        case 1
            % エネルギー補償による補正
            if cycle_info.energy_final > 0 && cycle_info.energy_initial > 0
                energy_ratio = cycle_info.energy_initial / cycle_info.energy_final;
                if energy_ratio > 1
                    energy_compensation = sqrt(energy_ratio);
                else
                    energy_compensation = 1.0;
                end
            else
                energy_compensation = 1.0;
            end
            
            % 角度をわずかに調整
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

function visualize_results_original(all_times, all_trajectories, all_extended_states, cycle_data, z0_original)
    % 結果の可視化
    
    figure('Name', '複数サイクル解析結果', 'Position', [50 50 1400 900]);
    
    % === サブプロット1: 角度の時系列 ===
    subplot(3,3,1);
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
    
    % === サブプロット3: エネルギーの推移 ===
    subplot(3,3,3);
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
    
    % === サブプロット5: 位相空間（q2-u2） ===
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
    
    % === サブプロット6: エネルギーの時系列 ===
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
    
    sgtitle('パッシブウォーカー: 複数サイクル解析（元の物理モデル）');
end

function display_summary(cycle_data, z0_original)
    % サマリー表示
    
    fprintf('\n=== シミュレーションサマリー ===\n');
    fprintf('総サイクル数: %d\n', length(cycle_data));
    
    if length(cycle_data) > 0
        avg_energy_change = mean([cycle_data.energy_loss]);
        fprintf('平均エネルギー変化: %.6f\n', avg_energy_change);
        
        if isfield(cycle_data(1), 'velocity_correction') && length(cycle_data) > 1
            corrections = [cycle_data(1:end-1).velocity_correction];
            avg_correction = mean(corrections);
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

% animate_cycle_original_fixed関数の修正版
function animate_cycle_original_fixed(fig, t, z_all, walker, cycle_num, speed_factor)
    % サイクルごとのアニメーション（物理的な脚の色を一貫して保持）

    if ~isvalid(fig)
        return;
    end

    figure(fig);
    clf;

    % パラメータ
    M = walker.M;  m = walker.m; I = walker.I;   
    l = walker.l;  c = walker.c; w = walker.w;   
    r = walker.r;  g = walker.g; gam = walker.gam;

    % アニメーション用データの準備
    fps = 10;
    if size(z_all, 2) >= 8
        z_all_plot = [z_all(:,1) z_all(:,3) z_all(:,6) z_all(:,8)];
    else
        hip_x = cumsum([0; diff(t)] .* 0.5);
        z_all_plot = [z_all(:,1) z_all(:,3) hip_x ones(size(z_all,1),1)*l];
    end

    nn = size(z_all_plot,2);
    total_frames = round(t(end)*fps);
    if total_frames < 1
        total_frames = 1;
    end
    t_anim = linspace(t(1), t(end), total_frames);
    z = zeros(total_frames,nn);

    for i=1:nn
        try
            z(:,i) = interp1(t, z_all_plot(:,i), t_anim, 'linear', 'extrap');
        catch
            z(:,i) = z_all_plot(end,i) * ones(total_frames, 1);
        end
    end

    % ウィンドウ設定
    min_xh = min(z(:,3)); max_xh = max(z(:,3)); 
    window_xmin = min_xh - 1*l; window_xmax = max_xh + 1*l;
    window_ymin = -0.1; window_ymax = 1.1*(l+r);

    % サブプロット1: アニメーション
    subplot(1,2,1);
    axis('equal')
    axis([window_xmin window_xmax window_ymin window_ymax])
    axis off
    set(gcf,'Color',[1,1,1])

    mm = size(z,1);

    % 軌跡を描画するための準備
    trajectory_x = z(:,3);
    trajectory_y = z(:,4);

    lines_for_feet = 4;
    counter = 2 + 2*lines_for_feet;
    th = 0.25;

    %%% create object for hinge %%%%%
    hingepic=line('xdata',0,'ydata',0, 'marker','.','markersize',[20], ...
              'color','black');

    %%% create trajectory line %%%
    trajpic=line('xdata',[],'ydata',[], 'linewidth', 1, 'color','green', 'linestyle', '--');

    %%%% create object for legs and feet %%%%
    barref = [0 0; 0 -1];
    y = [0;-1];
    O = [0; 0];

    % === 物理的な脚の色を固定 ===
    % 右脚は常に赤、左脚は常に青として固定
    right_leg_color = 'red';
    left_leg_color = 'blue';
    
    % 脚オブジェクトの初期化（色は後で決定）
    barpic(1) = line('xdata',barref(1,:),'ydata',barref(2,:),'linewidth', 3);
    barpic(2) = line('xdata',barref(1,:),'ydata',barref(2,:),'linewidth', 3);

    %%% feet in blue %%%
    for p = 3:counter
        barpic(p)=line('xdata',barref(1,:),'ydata',barref(2,:),'linewidth', 2, 'color','blue');
    end

    %%%% create ramp %%%%
    rampref=[window_xmin window_xmax ; 0 0];
    line('xdata',rampref(1,:),'ydata',rampref(2,:), 'linewidth', 2,'color','black');

    title(sprintf('サイクル %d - 歩行アニメーション（物理脚色固定）', cycle_num), 'FontSize', 14)

    % サブプロット2: 状態変数
    subplot(1,2,2);
    title('角度と角速度（物理脚色対応）');
    xlabel('時間 [s]');
    grid on;
    hold on;

    % アニメーションループ
    dt_anim = 0.1 / speed_factor;

    for i = 1:mm
        try
            % 左側：歩行アニメーション
            subplot(1,2,1);

            q1 = z(i,1); q2 = z(i,2); 
            xh = z(i,3); yh = z(i,4);  

            %%% hinge coordinates
            hinge=[xh; yh];   

            %%% leg coordinates
            A = [q1 -(q2-q1)];

            % === 物理的な脚の判定（角度の大小で判定） ===
            % より右側（大きい角度）にある脚を右脚、左側（小さい角度）にある脚を左脚とする
            if q1 > -(q2-q1)  % q1の方が右側
                q1_color = right_leg_color;  % q1 = 右脚
                q2_color = left_leg_color;   % q2 = 左脚
            else  % q2の方が右側
                q1_color = left_leg_color;   % q1 = 左脚
                q2_color = right_leg_color;  % q2 = 右脚
            end
            
            % 脚の色を更新
            set(barpic(1), 'color', q1_color);
            set(barpic(2), 'color', q2_color);

            for p = 1:2
                bar(:,:,p) = [hinge, hinge] + (l+r)*R(A(p))*barref;
                center(:,:,p) = hinge + l*R(A(p))*y;
            end

            %%%% feet coordinates
            B = [-th -2*th; 0 -th; th 0; 2*th th];

            incr = 3;
            for p=1:2
                for q=1:4
                    C = A(p) + B(q,:);
                    bar(:,:,incr) = [center(:,:,p), center(:,:,p)] + r*R(C(1))*[O,y] + r*R(C(2))*[y, O]; 
                    incr = incr + 1;
                end
            end

            %%% animate now    
            set(hingepic,'xdata',hinge(1),'ydata',hinge(2));
            set(trajpic,'xdata',trajectory_x(1:i),'ydata',trajectory_y(1:i));

            % 脚の位置を更新
            for p=1:2
                set(barpic(p),'xdata',bar(1,:,p),'ydata',bar(2,:,p));
            end

            % 足の部分
            for p=3:counter
                set(barpic(p),'xdata',bar(1,:,p),'ydata',bar(2,:,p));
            end

            % 右側：状態変数のプロット（物理的な脚の色に対応）
            subplot(1,2,2);
            cla;

            idx = find(t <= t_anim(i));
            if ~isempty(idx)
                plot(t(idx), z_all(idx,1), '-', 'Color', q1_color, 'LineWidth', 2, 'DisplayName', 'q1');
                plot(t(idx), z_all(idx,3), '-', 'Color', q2_color, 'LineWidth', 2, 'DisplayName', 'q2');
                plot(t(idx), z_all(idx,2), '--', 'Color', q1_color, 'LineWidth', 1, 'DisplayName', 'u1');
                plot(t(idx), z_all(idx,4), '--', 'Color', q2_color, 'LineWidth', 1, 'DisplayName', 'u2');

                % 現在位置をマーク
                if i <= length(z)
                    plot(t_anim(i), z(i,1), 'o', 'Color', q1_color, 'MarkerSize', 8, 'MarkerFaceColor', q1_color);
                    plot(t_anim(i), z(i,2), 'o', 'Color', q2_color, 'MarkerSize', 8, 'MarkerFaceColor', q2_color);
                end
            end

            xlim([t(1) t(end)]);
            ylim([-1 1]);
            legend('Location', 'best');

            drawnow;
            pause(dt_anim);

        catch me
            fprintf('アニメーションフレーム %d でエラー: %s\n', i, me.message);
            continue;
        end
    end
end

function animate_all_cycles_original_fixed(fig, all_times, all_extended_states, walker, speed_factor)
    % 全サイクルの連続アニメーション（物理的な脚の色を一貫して保持）
    
    if ~isvalid(fig)
        return;
    end
    
    figure(fig);
    clf;
    
    % パラメータ
    M = walker.M;  m = walker.m; I = walker.I;   
    l = walker.l;  c = walker.c; w = walker.w;   
    r = walker.r;  g = walker.g; gam = walker.gam;
    
    % 全データを結合（重複時間点を適切に処理）
    t_all = [];
    z_all = [];
    cycle_labels = [];  % 各時刻がどのサイクルに属するかを記録
    
    for i = 1:length(all_times)
        if i == 1
            t_all = all_times{i};
            z_all = all_extended_states{i};
            cycle_labels = ones(length(all_times{i}), 1);
        else
            t_new = all_times{i};
            z_new = all_extended_states{i};
            
            if t_new(1) <= t_all(end)
                non_duplicate_idx = t_new > t_all(end);
                if any(non_duplicate_idx)
                    t_all = [t_all; t_new(non_duplicate_idx)];
                    z_all = [z_all; z_new(non_duplicate_idx, :)];
                    cycle_labels = [cycle_labels; i * ones(sum(non_duplicate_idx), 1)];
                end
            else
                t_all = [t_all; t_new];
                z_all = [z_all; z_new];
                cycle_labels = [cycle_labels; i * ones(length(t_new), 1)];
            end
        end
    end
    
    if isempty(t_all) || size(z_all, 1) < 2
        fprintf('連続アニメーション用のデータが不十分です\n');
        return;
    end
    
    % 時間データの単調性を確保
    [t_all, unique_idx] = unique(t_all);
    z_all = z_all(unique_idx, :);
    cycle_labels = cycle_labels(unique_idx);
    
    % アニメーション用データの準備
    fps = 10;
    if size(z_all, 2) >= 8
        z_all_plot = [z_all(:,1) z_all(:,3) z_all(:,6) z_all(:,8)];
    else
        hip_x = cumsum([0; diff(t_all)] .* 0.5);
        z_all_plot = [z_all(:,1) z_all(:,3) hip_x ones(size(z_all,1),1)*l];
    end
    
    nn = size(z_all_plot,2);
    total_frames = round(t_all(end)*fps);
    if total_frames < 2
        total_frames = 2;
    end
    t = linspace(t_all(1), t_all(end), total_frames);
    z = zeros(total_frames,nn);
    
    try
        for i=1:nn
            z(:,i) = interp1(t_all, z_all_plot(:,i), t, 'linear', 'extrap');
        end
    catch me
        fprintf('補間エラー: %s\n', me.message);
        return;
    end
    
    % ウィンドウ設定
    min_xh = min(z(:,3)); max_xh = max(z(:,3)); 
    window_xmin = min_xh - 1*l; window_xmax = max_xh + 1*l;
    window_ymin = -0.1; window_ymax = 1.1*(l+r);
    
    subplot(1,2,1);
    axis('equal')
    axis([window_xmin window_xmax window_ymin window_ymax])
    axis off
    set(gcf,'Color',[1,1,1])
    
    mm = size(z,1);
    
    % 軌跡を描画するための準備
    trajectory_x = z(:,3);
    trajectory_y = z(:,4);
    
    lines_for_feet = 4;
    counter = 2 + 2*lines_for_feet;
    th = 0.25;
    
    %%% create object for hinge %%%%%
    hingepic=line('xdata',0,'ydata',0, 'marker','.','markersize',[20], ...
              'color','black');
    
    %%% create trajectory line %%%
    trajpic=line('xdata',[],'ydata',[], 'linewidth', 1, 'color','green', 'linestyle', '--');
       
    %%%% create object for legs and feet %%%%
    barref = [0 0; 0 -1];
    y = [0;-1];
    O = [0; 0];
    
    % === 物理的な脚の色を固定 ===
    right_leg_color = 'red';
    left_leg_color = 'blue';
    
    % 脚オブジェクトの初期化
    for p = 1:2
        barpic(p)=line('xdata',barref(1,:),'ydata',barref(2,:),'linewidth', 3);
    end
        
    %%% feet in blue %%%
    for p = 3:counter
        barpic(p)=line('xdata',barref(1,:),'ydata',barref(2,:),'linewidth', 2, 'color','blue');
    end
    
    %%%% create ramp %%%%
    rampref=[window_xmin window_xmax ; 0 0];
    line('xdata',rampref(1,:),'ydata',rampref(2,:), 'linewidth', 2,'color','black');
    
    title('全サイクル 歩行アニメーション（物理脚色固定）', 'FontSize', 14)
    
    % サブプロット2: 状態変数
    subplot(1,2,2);
    title('全サイクル 状態変数（物理脚色対応）');
    xlabel('時間 [s]');
    grid on;
    hold on;
    
    % アニメーションループ
    dt_anim = 0.1 / speed_factor;
    
    for i = 1:mm
        try
            % 左側：歩行アニメーション
            subplot(1,2,1);
            
            q1 = z(i,1); q2 = z(i,2); 
            xh = z(i,3); yh = z(i,4);  
            
            % === 物理的な脚の判定（角度の大小で判定） ===
            if q1 > -(q2-q1)  % q1の方が右側
                q1_color = right_leg_color;  % q1 = 右脚
                q2_color = left_leg_color;   % q2 = 左脚
            else  % q2の方が右側
                q1_color = left_leg_color;   % q1 = 左脚
                q2_color = right_leg_color;  % q2 = 右脚
            end
            
            % 脚の色を更新
            set(barpic(1), 'color', q1_color);
            set(barpic(2), 'color', q2_color);
            
            %%% hinge coordinates
            hinge=[xh; yh];   
            
            %%% leg coordinates
            A = [q1 -(q2-q1)];
            
            for p = 1:2
                bar(:,:,p) = [hinge, hinge] + (l+r)*R(A(p))*barref;
                center(:,:,p) = hinge + l*R(A(p))*y;
            end

            %%%% feet coordinates
            B = [-th -2*th; 0 -th; th 0; 2*th th];
            
            incr = 3;
            for p=1:2
                for q=1:4
                    C = A(p) + B(q,:);
                    bar(:,:,incr) = [center(:,:,p), center(:,:,p)] + r*R(C(1))*[O,y] + r*R(C(2))*[y, O]; 
                    incr = incr + 1;
                end
            end
                  
            %%% animate now    
            set(hingepic,'xdata',hinge(1),'ydata',hinge(2));
            set(trajpic,'xdata',trajectory_x(1:i),'ydata',trajectory_y(1:i));
            
            % 脚の位置を更新
            for p=1:2
                set(barpic(p),'xdata',bar(1,:,p),'ydata',bar(2,:,p));
            end
            
            % 足の部分
            for p=3:counter
                set(barpic(p),'xdata',bar(1,:,p),'ydata',bar(2,:,p));
            end
                        
            % 右側：状態変数（物理的な脚の色に対応）
            subplot(1,2,2);
            cla;
            
            idx = find(t_all <= t(i));
            if ~isempty(idx)
                plot(t_all(idx), z_all(idx,1), '-', 'Color', q1_color, 'LineWidth', 2, 'DisplayName', 'q1');
                plot(t_all(idx), z_all(idx,3), '-', 'Color', q2_color, 'LineWidth', 2, 'DisplayName', 'q2');
                plot(t_all(idx), z_all(idx,2), '--', 'Color', q1_color, 'LineWidth', 1, 'DisplayName', 'u1');
                plot(t_all(idx), z_all(idx,4), '--', 'Color', q2_color, 'LineWidth', 1, 'DisplayName', 'u2');

                % 現在位置をマーク
                if i <= length(z)
                    plot(t(i), z(i,1), 'o', 'Color', q1_color, 'MarkerSize', 8, 'MarkerFaceColor', q1_color);
                    plot(t(i), z(i,2), 'o', 'Color', q2_color, 'MarkerSize', 8, 'MarkerFaceColor', q2_color);
                end
            end
            
            xlim([t_all(1) t_all(end)]);
            ylim([-1 1]);
            legend('Location', 'best');
            
            drawnow;
            pause(dt_anim);
            
        catch me
            fprintf('連続アニメーションフレーム %d でエラー: %s\n', i, me.message);
            continue;
        end
    end
end



% === 以下、元のコードと完全に同じ物理モデル関数を使用 ===

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
    [t_step1, z_step1] = ode113(@single_stance, tspan, z0_extended, options, walker);
    
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
    [t_step1_5, z_step1_5] = ode113(@single_stance, tspan, z1_extended, options_catchup, walker);
    
    % 時間を調整（連続的にする）
    t_step1_5 = t_step1_5 - t1 + t_step1(end);
    
    % 0.5ステップ目の終了（スイング脚追いつき、Heel Strikeなし）
    if length(t_step1_5) < time_stamps
        fprintf('0.5ステップ目完了: t = %.4f秒でスイング脚がスタンス脚に追いつきました\n', t_step1_5(end));
        fprintf('最終角度: q1=%.4f, q2=%.4f, 差=%.6f\n', z_step1_5(end,1), z_step1_5(end,3), z_step1_5(end,3)-z_step1_5(end,1));
        % Heel Strikeは発生させない（スイング脚が追いついただけ）
        z_final = z_step1_5(end, 1:4);
    else
        fprintf('警告: 指定時間内にスイング脚が追いつきませんでした\n');
        fprintf('最終角度: q1=%.4f, q2=%.4f, 差=%.6f\n', z_step1_5(end,1), z_step1_5(end,3), z_step1_5(end,3)-z_step1_5(end,1));
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
    % q2 - q1 = 0 になったときに停止
    gstop = q2 - q1;
    isterminal = 1; % 条件に達したら積分を停止
    direction = 0;  % どちらの方向からでも検出（0 = 両方向）
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

function rotation = R(A)
    % 回転行列
    rotation = [cos(A) -sin(A); sin(A) cos(A)];
end

%% === 歩行安定性判定機能 ===
% このコードをmulti_cycle_passivewalker_4.mの最後に追加してください

function stability_result = analyze_walking_stability(all_trajectories, cycle_data, z0_original, walker)
    % 歩行の安定性を総合的に判定する関数
    %
    % 入力:
    %   all_trajectories: 各サイクルの軌道データ
    %   cycle_data: 各サイクルの統計情報
    %   z0_original: 初期状態（目標状態）
    %   walker: ウォーカーパラメータ
    %
    % 出力:
    %   stability_result: 安定性判定結果の構造体
    
    stability_result = struct();
    
    % 判定基準の設定
    criteria = struct();
    criteria.max_state_deviation = 0.1;        % 状態偏差の許容値
    criteria.max_energy_variation = 0.05;      % エネルギー変動の許容値（相対値）
    criteria.max_period_variation = 0.1;       % 周期変動の許容値（相対値）
    criteria.min_cycles_for_stability = 3;     % 安定判定に必要な最小サイクル数
    criteria.convergence_threshold = 0.01;     % 収束判定のしきい値
    
    % 1. 基本チェック
    n_cycles = length(cycle_data);
    if n_cycles < criteria.min_cycles_for_stability
        stability_result.is_stable = false;
        stability_result.reason = sprintf('サイクル数不足（%d < %d）', n_cycles, criteria.min_cycles_for_stability);
        display_stability_result(stability_result);
        return;
    end
    
    % 2. 状態変数の周期性チェック
    state_analysis = analyze_state_periodicity(cycle_data, z0_original, criteria);
    stability_result.state_analysis = state_analysis;
    
    % 3. エネルギーの安定性チェック
    energy_analysis = analyze_energy_stability(cycle_data, criteria);
    stability_result.energy_analysis = energy_analysis;
    
    % 4. 歩行周期の安定性チェック
    period_analysis = analyze_period_stability(cycle_data, criteria);
    stability_result.period_analysis = period_analysis;
    
    % 5. 位相空間での軌道の収束性チェック
    convergence_analysis = analyze_phase_space_convergence(all_trajectories, criteria);
    stability_result.convergence_analysis = convergence_analysis;
    
    % 6. 総合判定
    stability_result.is_stable = state_analysis.is_periodic && ...
                                energy_analysis.is_stable && ...
                                period_analysis.is_stable && ...
                                convergence_analysis.is_converging;
    
    if stability_result.is_stable
        stability_result.reason = '全ての安定性基準を満たしています';
        stability_result.stability_score = calculate_stability_score(stability_result);
    else
        reasons = {};
        if ~state_analysis.is_periodic
            reasons{end+1} = '状態変数が周期的でない';
        end
        if ~energy_analysis.is_stable
            reasons{end+1} = 'エネルギーが不安定';
        end
        if ~period_analysis.is_stable
            reasons{end+1} = '歩行周期が不安定';
        end
        if ~convergence_analysis.is_converging
            reasons{end+1} = '位相空間で収束していない';
        end
        stability_result.reason = strjoin(reasons, ', ');
        stability_result.stability_score = 0;
    end
    
    % 結果の表示
    display_stability_result(stability_result);
    
    % 詳細な安定性解析グラフの作成
    create_stability_analysis_plots(stability_result, cycle_data, all_trajectories);
end

function state_analysis = analyze_state_periodicity(cycle_data, z0_original, criteria)
    % 状態変数の周期性を解析
    
    state_analysis = struct();
    n_cycles = length(cycle_data);
    
    % 各サイクルの初期状態と最終状態の偏差を計算
    initial_deviations = zeros(n_cycles, 1);
    final_deviations = zeros(n_cycles, 1);
    
    for i = 1:n_cycles
        initial_deviations(i) = norm(cycle_data(i).initial_state - z0_original);
        final_deviations(i) = norm(cycle_data(i).final_state - z0_original);
    end
    
    % 偏差の傾向を分析
    state_analysis.mean_initial_deviation = mean(initial_deviations);
    state_analysis.mean_final_deviation = mean(final_deviations);
    state_analysis.max_deviation = max([initial_deviations; final_deviations]);
    
    % 偏差が減少傾向にあるかチェック
    if n_cycles >= 3
        p = polyfit(1:n_cycles, final_deviations', 1);
        state_analysis.deviation_trend = p(1);  % 傾き
        state_analysis.is_converging = p(1) < 0;
    else
        state_analysis.deviation_trend = 0;
        state_analysis.is_converging = false;
    end
    
    % 周期性の判定
    state_analysis.is_periodic = state_analysis.max_deviation < criteria.max_state_deviation;
end

function energy_analysis = analyze_energy_stability(cycle_data, criteria)
    % エネルギーの安定性を解析
    
    energy_analysis = struct();
    n_cycles = length(cycle_data);
    
    % エネルギー損失の統計
    energy_losses = [cycle_data.energy_loss];
    energy_analysis.mean_loss = mean(energy_losses);
    energy_analysis.std_loss = std(energy_losses);
    
    % エネルギーレベルの統計
    initial_energies = [cycle_data.energy_initial];
    final_energies = [cycle_data.energy_final];
    
    energy_analysis.mean_energy = mean([initial_energies, final_energies]);
    energy_analysis.energy_variation = std([initial_energies, final_energies]) / energy_analysis.mean_energy;
    
    % エネルギー変動の傾向
    if n_cycles >= 3
        p = polyfit(1:n_cycles, initial_energies, 1);
        energy_analysis.energy_trend = p(1);
    else
        energy_analysis.energy_trend = 0;
    end
    
    % 安定性の判定
    energy_analysis.is_stable = energy_analysis.energy_variation < criteria.max_energy_variation;
end

function period_analysis = analyze_period_stability(cycle_data, criteria)
    % 歩行周期の安定性を解析
    
    period_analysis = struct();
    n_cycles = length(cycle_data);
    
    % 各サイクルの継続時間
    durations = [cycle_data.duration];
    period_analysis.mean_period = mean(durations);
    period_analysis.std_period = std(durations);
    period_analysis.period_variation = period_analysis.std_period / period_analysis.mean_period;
    
    % 周期の変動傾向
    if n_cycles >= 3
        p = polyfit(1:n_cycles, durations, 1);
        period_analysis.period_trend = p(1);
        period_analysis.is_converging = abs(p(1)) < 0.01;  % 周期がほぼ一定
    else
        period_analysis.period_trend = 0;
        period_analysis.is_converging = false;
    end
    
    % 安定性の判定
    period_analysis.is_stable = period_analysis.period_variation < criteria.max_period_variation;
end

function convergence_analysis = analyze_phase_space_convergence(all_trajectories, criteria)
    % 位相空間での軌道の収束性を解析
    
    convergence_analysis = struct();
    n_cycles = length(all_trajectories);
    
    if n_cycles < 2
        convergence_analysis.is_converging = false;
        convergence_analysis.convergence_rate = 0;
        return;
    end
    
    % ポアンカレ断面での点の収束を解析
    poincare_points = zeros(n_cycles, 4);
    for i = 1:n_cycles
        % 各サイクルの中間点を取得（ポアンカレ断面として使用）
        traj = all_trajectories{i};
        mid_idx = round(size(traj, 1) / 2);
        poincare_points(i, :) = traj(mid_idx, :);
    end
    
    % 連続するポアンカレ点間の距離
    distances = zeros(n_cycles-1, 1);
    for i = 1:n_cycles-1
        distances(i) = norm(poincare_points(i+1, :) - poincare_points(i, :));
    end
    
    % 収束率の計算（距離が指数的に減少しているかチェック）
    if length(distances) >= 2
        convergence_analysis.mean_distance = mean(distances);
        convergence_analysis.distance_reduction_rate = distances(end) / distances(1);
        convergence_analysis.is_converging = convergence_analysis.distance_reduction_rate < 1 && ...
                                            distances(end) < criteria.convergence_threshold;
    else
        convergence_analysis.mean_distance = distances(1);
        convergence_analysis.distance_reduction_rate = 1;
        convergence_analysis.is_converging = false;
    end
    
    % フロケ乗数の推定（簡易版）
    if n_cycles >= 3
        % 最後の3サイクルから固有値を推定
        last_points = poincare_points(end-2:end, :);
        deviations = diff(last_points);
        if size(deviations, 1) >= 2
            ratio = norm(deviations(2, :)) / norm(deviations(1, :));
            convergence_analysis.estimated_floquet_multiplier = ratio;
        else
            convergence_analysis.estimated_floquet_multiplier = 1;
        end
    else
        convergence_analysis.estimated_floquet_multiplier = 1;
    end
end

function score = calculate_stability_score(stability_result)
    % 安定性スコアを計算（0-100）
    
    score = 0;
    
    % 状態周期性スコア（25点）
    if stability_result.state_analysis.is_periodic
        deviation_score = max(0, 1 - stability_result.state_analysis.max_deviation / 0.1);
        score = score + 25 * deviation_score;
    end
    
    % エネルギー安定性スコア（25点）
    if stability_result.energy_analysis.is_stable
        energy_score = max(0, 1 - stability_result.energy_analysis.energy_variation / 0.05);
        score = score + 25 * energy_score;
    end
    
    % 周期安定性スコア（25点）
    if stability_result.period_analysis.is_stable
        period_score = max(0, 1 - stability_result.period_analysis.period_variation / 0.1);
        score = score + 25 * period_score;
    end
    
    % 収束性スコア（25点）
    if stability_result.convergence_analysis.is_converging
        convergence_score = max(0, 1 - stability_result.convergence_analysis.distance_reduction_rate);
        score = score + 25 * convergence_score;
    end
end

function display_stability_result(stability_result)
    % 安定性判定結果を表示
    
    fprintf('\n');
    fprintf('==========================================================\n');
    fprintf('              歩行安定性判定結果                          \n');
    fprintf('==========================================================\n');
    
    if stability_result.is_stable
        fprintf('判定: ✓ 安定した歩行\n');
        fprintf('安定性スコア: %.1f/100\n', stability_result.stability_score);
    else
        fprintf('判定: ✗ 不安定な歩行\n');
        fprintf('理由: %s\n', stability_result.reason);
    end
    
    fprintf('\n--- 詳細解析結果 ---\n');
    
    % 状態変数の周期性
    if isfield(stability_result, 'state_analysis')
        sa = stability_result.state_analysis;
        fprintf('\n1. 状態変数の周期性:\n');
        fprintf('   - 最大偏差: %.4f %s\n', sa.max_deviation, ...
            ternary(sa.is_periodic, '✓', '✗'));
        fprintf('   - 平均初期偏差: %.4f\n', sa.mean_initial_deviation);
        fprintf('   - 平均最終偏差: %.4f\n', sa.mean_final_deviation);
        if isfield(sa, 'is_converging')
            fprintf('   - 収束傾向: %s (傾き: %.6f)\n', ...
                ternary(sa.is_converging, 'あり', 'なし'), sa.deviation_trend);
        end
    end
    
    % エネルギーの安定性
    if isfield(stability_result, 'energy_analysis')
        ea = stability_result.energy_analysis;
        fprintf('\n2. エネルギーの安定性:\n');
        fprintf('   - エネルギー変動率: %.2f%% %s\n', ...
            ea.energy_variation * 100, ternary(ea.is_stable, '✓', '✗'));
        fprintf('   - 平均エネルギー損失: %.6f\n', ea.mean_loss);
        fprintf('   - エネルギー変化傾向: %.6f\n', ea.energy_trend);
    end
    
    % 歩行周期の安定性
    if isfield(stability_result, 'period_analysis')
        pa = stability_result.period_analysis;
        fprintf('\n3. 歩行周期の安定性:\n');
        fprintf('   - 周期変動率: %.2f%% %s\n', ...
            pa.period_variation * 100, ternary(pa.is_stable, '✓', '✗'));
        fprintf('   - 平均周期: %.4f秒\n', pa.mean_period);
        fprintf('   - 周期標準偏差: %.4f秒\n', pa.std_period);
    end
    
    % 位相空間での収束性
    if isfield(stability_result, 'convergence_analysis')
        ca = stability_result.convergence_analysis;
        fprintf('\n4. 位相空間での収束性:\n');
        fprintf('   - 収束判定: %s %s\n', ...
            ternary(ca.is_converging, '収束', '発散'), ...
            ternary(ca.is_converging, '✓', '✗'));
        fprintf('   - 距離減少率: %.4f\n', ca.distance_reduction_rate);
        if isfield(ca, 'estimated_floquet_multiplier')
            fprintf('   - 推定フロケ乗数: %.4f\n', ca.estimated_floquet_multiplier);
        end
    end
    
    fprintf('\n==========================================================\n');
end

function result = ternary(condition, true_val, false_val)
    % 三項演算子の代替
    if condition
        result = true_val;
    else
        result = false_val;
    end
end

function create_stability_analysis_plots(stability_result, cycle_data, all_trajectories)
    % 安定性解析の詳細グラフを作成
    
    figure('Name', '歩行安定性解析', 'Position', [100 100 1200 800]);
    
    n_cycles = length(cycle_data);
    cycles = 1:n_cycles;
    
    % サブプロット1: 状態偏差の推移
    subplot(2,3,1);
    initial_deviations = zeros(n_cycles, 1);
    final_deviations = zeros(n_cycles, 1);
    for i = 1:n_cycles
        initial_deviations(i) = norm(cycle_data(i).initial_state - cycle_data(1).initial_state);
        final_deviations(i) = norm(cycle_data(i).final_state - cycle_data(1).initial_state);
    end
    
    plot(cycles, initial_deviations, 'bo-', 'LineWidth', 2, 'DisplayName', '初期状態');
    hold on;
    plot(cycles, final_deviations, 'ro-', 'LineWidth', 2, 'DisplayName', '最終状態');
    
    % トレンドライン
    if n_cycles >= 3
        p = polyfit(cycles, final_deviations', 1);
        trend_line = polyval(p, cycles);
        plot(cycles, trend_line, 'r--', 'LineWidth', 1, 'DisplayName', 'トレンド');
    end
    
    xlabel('サイクル');
    ylabel('状態偏差のノルム');
    title('状態変数の偏差推移');
    legend('Location', 'best');
    grid on;
    
    % サブプロット2: エネルギー推移
    subplot(2,3,2);
    initial_energies = [cycle_data.energy_initial];
    final_energies = [cycle_data.energy_final];
    
    plot(cycles, initial_energies, 'bo-', 'LineWidth', 2, 'DisplayName', '初期エネルギー');
    hold on;
    plot(cycles, final_energies, 'ro-', 'LineWidth', 2, 'DisplayName', '最終エネルギー');
    
    % 平均線
    mean_energy = mean([initial_energies, final_energies]);
    plot([1, n_cycles], [mean_energy, mean_energy], 'g--', 'LineWidth', 1, 'DisplayName', '平均');
    
    xlabel('サイクル');
    ylabel('総エネルギー');
    title('エネルギーレベルの推移');
    legend('Location', 'best');
    grid on;
    
    % サブプロット3: 歩行周期
    subplot(2,3,3);
    durations = [cycle_data.duration];
    
    plot(cycles, durations, 'ko-', 'LineWidth', 2);
    hold on;
    
    % 平均と標準偏差
    mean_duration = mean(durations);
    std_duration = std(durations);
    plot([1, n_cycles], [mean_duration, mean_duration], 'b--', 'LineWidth', 1);
    fill([1, n_cycles, n_cycles, 1], ...
         [mean_duration-std_duration, mean_duration-std_duration, ...
          mean_duration+std_duration, mean_duration+std_duration], ...
         'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    xlabel('サイクル');
    ylabel('周期 [秒]');
    title('歩行周期の変動');
    grid on;
    
    % サブプロット4: ポアンカレ断面
    subplot(2,3,4);
    colors = jet(n_cycles);
    
    for i = 1:n_cycles
        traj = all_trajectories{i};
        mid_idx = round(size(traj, 1) / 2);
        if i == 1
            plot(traj(mid_idx, 1), traj(mid_idx, 2), 'o', ...
                'Color', colors(i,:), 'MarkerSize', 10, 'MarkerFaceColor', colors(i,:));
        else
            plot(traj(mid_idx, 1), traj(mid_idx, 2), 'o', ...
                'Color', colors(i,:), 'MarkerSize', 8, 'MarkerFaceColor', colors(i,:));
        end
        hold on;
    end
    
    xlabel('q1 [rad]');
    ylabel('u1 [rad/s]');
    title('ポアンカレ断面 (q1-u1)');
    grid on;
    
    % サブプロット5: 収束性指標
    subplot(2,3,5);
    if n_cycles >= 2
        distances = zeros(n_cycles-1, 1);
        for i = 1:n_cycles-1
            traj1 = all_trajectories{i};
            traj2 = all_trajectories{i+1};
            mid_idx1 = round(size(traj1, 1) / 2);
            mid_idx2 = round(size(traj2, 1) / 2);
            distances(i) = norm(traj2(mid_idx2, :) - traj1(mid_idx1, :));
        end
        
        semilogy(1:n_cycles-1, distances, 'mo-', 'LineWidth', 2);
        xlabel('サイクル間');
        ylabel('ポアンカレ点間距離（対数）');
        title('収束性の指標');
        grid on;
    end
    
    % サブプロット6: 安定性サマリー
    subplot(2,3,6);
    axis off;
    
    % テキストで結果を表示
    text_y = 0.9;
    text(0.1, text_y, '安定性判定結果', 'FontSize', 14, 'FontWeight', 'bold');
    text_y = text_y - 0.15;
    
    if stability_result.is_stable
        text(0.1, text_y, sprintf('総合判定: 安定 (スコア: %.1f/100)', ...
            stability_result.stability_score), 'FontSize', 12, 'Color', 'green');
    else
        text(0.1, text_y, '総合判定: 不安定', 'FontSize', 12, 'Color', 'red');
    end
    
    text_y = text_y - 0.1;
    text(0.1, text_y, sprintf('理由: %s', stability_result.reason), ...
        'FontSize', 10, 'Interpreter', 'none');
    
    % 各項目の判定
    text_y = text_y - 0.15;
    items = {'状態周期性', 'エネルギー安定性', '周期安定性', '位相空間収束性'};
    checks = [stability_result.state_analysis.is_periodic, ...
              stability_result.energy_analysis.is_stable, ...
              stability_result.period_analysis.is_stable, ...
              stability_result.convergence_analysis.is_converging];
    
    for i = 1:length(items)
        if checks(i)
            marker = '✓';
            color = 'green';
        else
            marker = '✗';
            color = 'red';
        end
        text(0.1, text_y, sprintf('%s %s', marker, items{i}), ...
            'FontSize', 11, 'Color', color);
        text_y = text_y - 0.08;
    end
    
    sgtitle('パッシブウォーカー 歩行安定性解析', 'FontSize', 16);
end

%% メイン関数の最後に追加するコード
% multi_cycle_passivewalker_corrected() の最後、
% save_simulation_data_to_csv の呼び出しの後に以下を追加:


%% 安定性解析結果をCSVに保存する関数
function save_stability_analysis_to_csv(save_dir, stability_result, cycle_data)
    % 安定性解析結果をCSVファイルに保存
    
    filename = fullfile(save_dir, 'stability_analysis.csv');
    
    fid = fopen(filename, 'w');
    fprintf(fid, 'Analysis Item,Value,Status\n');
    
    % 総合判定
    fprintf(fid, 'Overall Stability,%s,%s\n', ...
        ternary(stability_result.is_stable, 'Stable', 'Unstable'), ...
        stability_result.reason);
    
    if stability_result.is_stable
        fprintf(fid, 'Stability Score,%.1f,/100\n', stability_result.stability_score);
    end
    
    % 状態変数の周期性
    sa = stability_result.state_analysis;
    fprintf(fid, 'State Periodicity,%s,Max deviation: %.4f\n', ...
        ternary(sa.is_periodic, 'Periodic', 'Not periodic'), sa.max_deviation);
    fprintf(fid, 'Mean Initial Deviation,%.6f,\n', sa.mean_initial_deviation);
    fprintf(fid, 'Mean Final Deviation,%.6f,\n', sa.mean_final_deviation);
    if isfield(sa, 'deviation_trend')
        fprintf(fid, 'Deviation Trend,%.6f,%s\n', sa.deviation_trend, ...
            ternary(sa.is_converging, 'Converging', 'Diverging'));
    end
    
    % エネルギーの安定性
    ea = stability_result.energy_analysis;
    fprintf(fid, 'Energy Stability,%s,Variation: %.2f%%\n', ...
        ternary(ea.is_stable, 'Stable', 'Unstable'), ea.energy_variation * 100);
    fprintf(fid, 'Mean Energy Loss,%.6f,\n', ea.mean_loss);
    fprintf(fid, 'Energy Trend,%.6f,\n', ea.energy_trend);
    
    % 歩行周期の安定性
    pa = stability_result.period_analysis;
    fprintf(fid, 'Period Stability,%s,Variation: %.2f%%\n', ...
        ternary(pa.is_stable, 'Stable', 'Unstable'), pa.period_variation * 100);
    fprintf(fid, 'Mean Period,%.4f,seconds\n', pa.mean_period);
    fprintf(fid, 'Period Std Dev,%.4f,seconds\n', pa.std_period);
    
    % 位相空間での収束性
    ca = stability_result.convergence_analysis;
    fprintf(fid, 'Phase Space Convergence,%s,Distance reduction: %.4f\n', ...
        ternary(ca.is_converging, 'Converging', 'Not converging'), ...
        ca.distance_reduction_rate);
    if isfield(ca, 'estimated_floquet_multiplier')
        fprintf(fid, 'Estimated Floquet Multiplier,%.4f,\n', ...
            ca.estimated_floquet_multiplier);
    end
    
    fclose(fid);
    
    fprintf('安定性解析結果を保存しました: %s\n', filename);
end