function multi_cycle_passivewalker_original()
    % 元の物理モデルを使用した複数サイクルパッシブウォーカーシミュレーション
    % single_stride_passivewalker2.mと完全に同じ計算ロジックを使用
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
    speed_factor = [1, 2, 0.5];
    anim_speed_factor = speed_factor(anim_speed);
    
    % 補正方法の選択
    fprintf('\n角速度補正方法の選択:\n');
    fprintf('1. エネルギー補償による補正（推奨）\n');
    fprintf('2. 固定倍率による補正\n');
    fprintf('3. 手動補正\n');
    correction_method = input('選択 (1-3): ');
    
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
        
        % 1.5歩のシミュレーション実行（single_strideと同じ関数を使用）
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
        animate_cycle_original(fig_anim, t_traj, z_traj, walker, cycle, anim_speed_factor);
        
        % 次のサイクルのための角速度補正
        if cycle < n_cycles
            fprintf('\n次サイクルの準備...\n');
            
            switch correction_method
                case 1
                    % エネルギー補償による補正
                    % エネルギーが増加した場合、角速度を調整する必要がある
                    if cycle_info.energy_loss > 0  % エネルギーが増加した場合
                        % 角速度を減少させる
                        energy_compensation = sqrt(TE_initial / TE_final);
                    else
                        % エネルギーが減少した場合は、元の角速度を維持
                        energy_compensation = 1.0;
                    end
                    
                    % 角度をわずかに調整して、次のステップで衝突が起きやすくする
                    angle_adjustment = 0.05;  % 小さな角度調整
                    
                    z_current = [z_final(1) + angle_adjustment; 
                                z0_original(2) * energy_compensation;
                                z_final(3) - angle_adjustment; 
                                z0_original(4) * energy_compensation];
                    
                case 2
                    % 固定倍率による補正
                    correction_factor = 1.1;
                    angle_adjustment = 0.05;
                    
                    z_current = [z_final(1) + angle_adjustment; 
                                z0_original(2) * correction_factor;
                                z_final(3) - angle_adjustment; 
                                z0_original(4) * correction_factor];
                    
                case 3
                    % 手動補正
                    fprintf('現在の角度: q1=%.4f, q2=%.4f\n', z_final(1), z_final(3));
                    u1_corrected = input('新しい角速度 u1 を入力: ');
                    u2_corrected = input('新しい角速度 u2 を入力: ');
                    z_current = [z_final(1); u1_corrected; z_final(3); u2_corrected];
            end
            
            fprintf('補正後の状態: [%.4f, %.4f, %.4f, %.4f]\n', z_current);
            
            cycle_info.velocity_correction = norm([z_current(2) - z_final(2); 
                                                  z_current(4) - z_final(4)]);
        end
        
        cycle_data = [cycle_data, cycle_info];
        
        % 少し待機
        pause(0.5);
    end
    
    % === 全体のアニメーション再生 ===
    fprintf('\n全サイクルの連続アニメーションを再生しますか？ (y/n): ');
    if input('', 's') == 'y'
        animate_all_cycles_original(fig_anim, all_times, all_extended_states, walker, anim_speed_factor);
    end
    
    % === 結果の可視化 ===
    visualize_results_original(all_times, all_trajectories, all_extended_states, cycle_data, z0_original);
    
    % === サマリー表示 ===
    display_summary(cycle_data, z0_original);
end

function animate_cycle_original(fig, t, z_all, walker, cycle_num, speed_factor)
    % single_stride_passivewalker2.mと同じアニメーション形式を使用
    
    figure(fig);
    clf;
    
    % パラメータ
    M = walker.M;  m = walker.m; I = walker.I;   
    l = walker.l;  c = walker.c; w = walker.w;   
    r = walker.r;  g = walker.g; gam = walker.gam;
    
    % アニメーション用データの準備
    fps = 10;
    z_all_plot = [z_all(:,1) z_all(:,3) z_all(:,6) z_all(:,8)];
    nn = size(z_all_plot,2);
    total_frames = round(t(end)*fps);
    t_anim = linspace(t(1), t(end), total_frames);
    z = zeros(total_frames,nn);
    for i=1:nn
        z(:,i) = interp1(t, z_all_plot(:,i), t_anim);
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
    
    %%%% legs in red %%%      
    for p = 1:2
        barpic(p)=line('xdata',barref(1,:),'ydata',barref(2,:),'linewidth', 3, 'color','red');
    end
    
    %%% feet in blue %%%
    for p = 3:counter
        barpic(p)=line('xdata',barref(1,:),'ydata',barref(2,:),'linewidth', 2, 'color','blue');
    end
    
    %%%% create ramp %%%%
    rampref=[window_xmin window_xmax ; 0 0];
    line('xdata',rampref(1,:),'ydata',rampref(2,:), 'linewidth', 2,'color','black');
    
    title(sprintf('サイクル %d - 歩行アニメーション', cycle_num), 'FontSize', 14)
    
    % サブプロット2: 状態変数
    subplot(1,2,2);
    title('角度と角速度');
    xlabel('時間 [s]');
    grid on;
    hold on;
    
    % アニメーションループ
    dt_anim = 0.1 / speed_factor;
    
    for i = 1:mm
        % 左側：歩行アニメーション
        subplot(1,2,1);
        
        % z_all_plot = [z_all(:,1) z_all(:,3) z_all(:,6) z_all(:,8)]
        % つまり z(:,1)=q1, z(:,2)=q2, z(:,3)=xh, z(:,4)=yh
        q1 = z(i,1); q2 = z(i,2); 
        xh = z(i,3); yh = z(i,4);  
        
        %%% hinge coordinates
        hinge=[xh; yh];   
        
        %%% leg coordinates
        % A(1) = q1 (スタンス脚の角度)
        % A(2) = -(q2-q1) (スイング脚の相対角度の負値)
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
        
        for p=1:counter
            set(barpic(p),'xdata',bar(1,:,p),'ydata',bar(2,:,p));
        end
        
        % 右側：状態変数のプロット
        subplot(1,2,2);
        cla;
        
        % これまでの軌跡
        idx = find(t <= t_anim(i));
        plot(t(idx), z_all(idx,1), 'b-', 'LineWidth', 2, 'DisplayName', 'q1');
        plot(t(idx), z_all(idx,3), 'r-', 'LineWidth', 2, 'DisplayName', 'q2');
        plot(t(idx), z_all(idx,2), 'b--', 'LineWidth', 1, 'DisplayName', 'u1');
        plot(t(idx), z_all(idx,4), 'r--', 'LineWidth', 1, 'DisplayName', 'u2');
        
        % 現在位置をマーク
        plot(t_anim(i), z(i,1), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
        plot(t_anim(i), z(i,2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        
        xlim([t(1) t(end)]);
        ylim([-1 1]);
        legend('Location', 'best');
        
        drawnow;
        pause(dt_anim);
    end
end

function animate_all_cycles_original(fig, all_times, all_extended_states, walker, speed_factor)
    % 全サイクルの連続アニメーション（single_stride_passivewalker2.m形式）
    
    figure(fig);
    clf;
    
    % パラメータ
    M = walker.M;  m = walker.m; I = walker.I;   
    l = walker.l;  c = walker.c; w = walker.w;   
    r = walker.r;  g = walker.g; gam = walker.gam;
    
    % 全データを結合
    t_all = [];
    z_all = [];
    for i = 1:length(all_times)
        t_all = [t_all; all_times{i}];
        z_all = [z_all; all_extended_states{i}];
    end
    
    % アニメーション用データの準備
    fps = 10;
    z_all_plot = [z_all(:,1) z_all(:,3) z_all(:,6) z_all(:,8)];
    nn = size(z_all_plot,2);
    total_frames = round(t_all(end)*fps);
    t = linspace(t_all(1), t_all(end), total_frames);
    z = zeros(total_frames,nn);
    for i=1:nn
        z(:,i) = interp1(t_all, z_all_plot(:,i), t);
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
    
    %%%% legs in red %%%      
    for p = 1:2
        barpic(p)=line('xdata',barref(1,:),'ydata',barref(2,:),'linewidth', 3, 'color','red');
    end
    
    %%% feet in blue %%%
    for p = 3:counter
        barpic(p)=line('xdata',barref(1,:),'ydata',barref(2,:),'linewidth', 2, 'color','blue');
    end
    
    %%%% create ramp %%%%
    rampref=[window_xmin window_xmax ; 0 0];
    line('xdata',rampref(1,:),'ydata',rampref(2,:), 'linewidth', 2,'color','black');
    
    title('全サイクル 歩行アニメーション', 'FontSize', 14)
    
    % サブプロット2: 状態変数
    subplot(1,2,2);
    title('全サイクル 状態変数');
    xlabel('時間 [s]');
    grid on;
    hold on;
    
    % アニメーションループ
    dt_anim = 0.1 / speed_factor;
    
    for i = 1:mm
        % 左側：歩行アニメーション
        subplot(1,2,1);
        
        % z_all_plot = [z_all(:,1) z_all(:,3) z_all(:,6) z_all(:,8)]
        % つまり z(:,1)=q1, z(:,2)=q2, z(:,3)=xh, z(:,4)=yh
        q1 = z(i,1); q2 = z(i,2); 
        xh = z(i,3); yh = z(i,4);  
        
        %%% hinge coordinates
        hinge=[xh; yh];   
        
        %%% leg coordinates
        % A(1) = q1 (スタンス脚の角度)
        % A(2) = -(q2-q1) (スイング脚の相対角度の負値)
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
        
        for p=1:counter
            set(barpic(p),'xdata',bar(1,:,p),'ydata',bar(2,:,p));
        end
        
        % 右側：状態変数
        subplot(1,2,2);
        cla;
        idx = find(t_all <= t(i));
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

function visualize_results_original(all_times, all_trajectories, all_extended_states, cycle_data, z0_original)
    % 結果の可視化（元のコードのスタイルを参考に）
    
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
    
    % === サブプロット6: エネルギーの時系列（拡張状態から） ===
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
    
    % === サブプロット7: ヒップ軌道（もし利用可能なら） ===
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
    
    % === サブプロット8: 角速度補正量 ===
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
    
    % === サブプロット9: 状態ベクトルのノルム ===
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
    
    sgtitle('パッシブウォーカー: 複数サイクル解析（元の物理モデル）');
end

function display_summary(cycle_data, z0_original)
    % サマリー表示
    
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

% === 以下、single_stride_passivewalker2.mと完全に同じ物理モデル関数を使用 ===

%===================================================================
function [z_traj, t_traj, z_final, z_midpoint] = one_and_half_stride_detailed(z0, walker)
%===================================================================
% 1.5ストライド（スイング脚がスタンス脚に追いつくまで）の詳細な軌道を返す関数
% single_stride_passivewalker2.mの実装と完全に一致

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

%===================================================================
function TE = calculate_energy(z, walker)
%===================================================================
% 総エネルギーを計算

q1 = z(1); u1 = z(2); q2 = z(3); u2 = z(4);
M = walker.M;  m = walker.m; I = walker.I;   
l = walker.l;  c = walker.c; w = walker.w;   
r = walker.r;  g = walker.g; gam = walker.gam;

TE = 1/2*m*(((-l*cos(q1)-r)*u1-u1*(-c*cos(q1)+w*sin(q1)))^2+(-l*sin(q1)*u1+u1*(c*sin(q1)+w*cos(q1)))^2)+1/2*m*(((-l*cos(q1)-r)*u1-(u1-u2)*(-c*cos(q1-q2)+w*sin(q1-q2)))^2+(-l*sin(q1)*u1+(u1-u2)*(c*sin(q1-q2)+w*cos(q1-q2)))^2)+1/2*M*((-l*cos(q1)-r)^2*u1^2+l^2*sin(q1)^2*u1^2)+1/2*I*(u1^2+(u1-u2)^2)+2*m*g*cos(gam)*r+2*m*g*l*cos(gam-q1)-m*g*c*cos(gam-q1)-m*g*w*sin(gam-q1)+2*m*g*sin(gam)*r*q1-m*g*c*cos(gam-q1+q2)-m*g*w*sin(gam-q1+q2)+M*g*cos(gam)*r+M*g*l*cos(gam-q1)+M*g*sin(gam)*r*q1;

end

%===================================================================
function [gstop, isterminal,direction]=swing_catches_stance(t,z,walker)
%===================================================================
% スイング脚がスタンス脚に追いつく条件

q1 = z(1); q2 = z(3); 

% スイング脚角度がスタンス脚角度に追いつく条件
% q2 - q1 = 0 になったときに停止
gstop = q2 - q1;
isterminal = 1; % 条件に達したら積分を停止
direction = 0;  % どちらの方向からでも検出（0 = 両方向）

% デバッグ用
fprintf('Debug: t=%.4f, q1=%.4f, q2=%.4f, diff=%.6f\n', t, q1, q2, gstop);

end

%===================================================================
function zdot=single_stance(t,z,walker)  
%===================================================================

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

%===================================================================
function [gstop, isterminal,direction]=collision(t,z,walker)
%===================================================================

q1 = z(1); q2 = z(3); 

gstop = -q2 + 2*q1;
if (q2>-0.05)
    isterminal = 0;
else
    isterminal=1;
end
direction=-1;

end

%===================================================================
function zplus=heelstrike(t,z,walker)      
%===================================================================

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

%===================================================================
function rotation = R(A)
%===================================================================
rotation = [cos(A) -sin(A); sin(A) cos(A)];
end