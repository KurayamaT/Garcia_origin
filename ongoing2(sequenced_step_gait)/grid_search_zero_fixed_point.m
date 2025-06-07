function grid_search_zero_fixed_point(flag)
% GRID_SEARCH_ZERO_FIXED_POINT - q1=0, q2=0の固定点を持つ連続歩行のための初期値探索
% 使い方:
% grid_search_zero_fixed_point(1)  % Garcia's Simplest Walker
% grid_search_zero_fixed_point(2)  % General Round Feet Walker

clc
close all
format long

% 引数の処理
if nargin == 0
    flag = 1; % デフォルトはGarcia's Simplest Walker
end

% ウォーカーパラメータの設定
if flag == 1
    %% Garcia's simplest walker
    walker.M = 1000; walker.m = 1.0; walker.I = 0.00; walker.l = 1.0; walker.w = 0.0; 
    walker.c = 1.0;  walker.r = 0.3; walker.g = 1.0; walker.gam = 0.009; 
    fprintf('=== Garcia''s Simplest Walker - 固定点探索 (q1=q2=0) ===\n\n');
else 
    %%  More General round feet walker
    walker.M = 1.0; walker.m = 0.5; walker.I = 0.02; walker.l = 1.0; walker.w = 0.0; 
    walker.c = 0.5; walker.r = 0.2; walker.g = 1.0; walker.gam = 0.01; 
    fprintf('=== General Round Feet Walker - 固定点探索 (q1=q2=0) ===\n\n');
end

% グリッドサーチパラメータ
fprintf('グリッドサーチパラメータを設定してください:\n');
fprintf('探索モードを選択:\n');
fprintf('1. 角度固定モード (q1, q2を固定値に設定)\n');
fprintf('2. 角度範囲モード (q1, q2も範囲で探索)\n');
search_mode = input('モード選択 (1 or 2) [1]: ');
if isempty(search_mode), search_mode = 1; end

if search_mode == 1
    % 角度固定モード
    fprintf('\n固定角度値を設定:\n');
    default_q1_fixed = 0;
    default_q2_fixed = 0;
    
    q1_fixed = input(sprintf('q1の固定値 [%g]: ', default_q1_fixed));
    if isempty(q1_fixed), q1_fixed = default_q1_fixed; end
    
    q2_fixed = input(sprintf('q2の固定値 [%g]: ', default_q2_fixed));
    if isempty(q2_fixed), q2_fixed = default_q2_fixed; end
    
    fprintf('\n速度範囲の設定:\n');
    fprintf('推奨範囲: u1 ∈ [-1.0, 0.0], u2 ∈ [-1.0, 0.0]\n');
    
    % デフォルト値の設定
    default_u1_min = -0.5;
    default_u1_max = -0.1;
    default_u2_min = -0.8;
    default_u2_max = -0.1;
    default_n_points = 20;
    
    % ユーザー入力
    u1_min = input(sprintf('u1の最小値 [%g]: ', default_u1_min));
    if isempty(u1_min), u1_min = default_u1_min; end
    
    u1_max = input(sprintf('u1の最大値 [%g]: ', default_u1_max));
    if isempty(u1_max), u1_max = default_u1_max; end
    
    u2_min = input(sprintf('u2の最小値 [%g]: ', default_u2_min));
    if isempty(u2_min), u2_min = default_u2_min; end
    
    u2_max = input(sprintf('u2の最大値 [%g]: ', default_u2_max));
    if isempty(u2_max), u2_max = default_u2_max; end
    
    n_points_u = input(sprintf('速度の分割数 [%d]: ', default_n_points));
    if isempty(n_points_u), n_points_u = default_n_points; end
    
    % グリッドの作成
    u1_range = linspace(u1_min, u1_max, n_points_u);
    u2_range = linspace(u2_min, u2_max, n_points_u);
    q1_range = q1_fixed;
    q2_range = q2_fixed;
    
else
    % 角度範囲モード
    fprintf('\n角度範囲の設定:\n');
    fprintf('推奨範囲: q1, q2 ∈ [-0.3, 0.3]\n');
    
    % 角度のデフォルト値
    default_q1_min = -0.2;
    default_q1_max = 0.2;
    default_q2_min = -0.2;
    default_q2_max = 0.2;
    
    q1_min = input(sprintf('q1の最小値 [%g]: ', default_q1_min));
    if isempty(q1_min), q1_min = default_q1_min; end
    
    q1_max = input(sprintf('q1の最大値 [%g]: ', default_q1_max));
    if isempty(q1_max), q1_max = default_q1_max; end
    
    q2_min = input(sprintf('q2の最小値 [%g]: ', default_q2_min));
    if isempty(q2_min), q2_min = default_q2_min; end
    
    q2_max = input(sprintf('q2の最大値 [%g]: ', default_q2_max));
    if isempty(q2_max), q2_max = default_q2_max; end
    
    n_points_q = input(sprintf('角度の分割数 [5]: '));
    if isempty(n_points_q), n_points_q = 5; end
    
    fprintf('\n速度範囲の設定:\n');
    fprintf('推奨範囲: u1 ∈ [-1.0, 0.0], u2 ∈ [-1.0, 0.0]\n');
    
    % 速度のデフォルト値
    default_u1_min = -0.5;
    default_u1_max = -0.1;
    default_u2_min = -0.8;
    default_u2_max = -0.1;
    
    u1_min = input(sprintf('u1の最小値 [%g]: ', default_u1_min));
    if isempty(u1_min), u1_min = default_u1_min; end
    
    u1_max = input(sprintf('u1の最大値 [%g]: ', default_u1_max));
    if isempty(u1_max), u1_max = default_u1_max; end
    
    u2_min = input(sprintf('u2の最小値 [%g]: ', default_u2_min));
    if isempty(u2_min), u2_min = default_u2_min; end
    
    u2_max = input(sprintf('u2の最大値 [%g]: ', default_u2_max));
    if isempty(u2_max), u2_max = default_u2_max; end
    
    n_points_u = input(sprintf('速度の分割数 [10]: '));
    if isempty(n_points_u), n_points_u = 10; end
    
    % グリッドの作成
    q1_range = linspace(q1_min, q1_max, n_points_q);
    q2_range = linspace(q2_min, q2_max, n_points_q);
    u1_range = linspace(u1_min, u1_max, n_points_u);
    u2_range = linspace(u2_min, u2_max, n_points_u);
end

% 結果保存用の配列
results = struct('q1', [], 'q2', [], 'u1', [], 'u2', [], ...
                'q1_final', [], 'q2_final', [], 'u1_final', [], 'u2_final', [], ...
                'error', [], 'stable', [], 'energy_loss', [], 'converged', []);

% 進捗表示の準備
if search_mode == 1
    total_tests = n_points_u * n_points_u;
    fprintf('\n総テスト数: %d (q1=%.3f, q2=%.3f固定)\n', total_tests, q1_fixed, q2_fixed);
else
    total_tests = length(q1_range) * length(q2_range) * length(u1_range) * length(u2_range);
    fprintf('\n総テスト数: %d (4次元グリッドサーチ)\n', total_tests);
end
fprintf('グリッドサーチ開始...\n\n');

% タイマー開始
tic;

% グリッドサーチ実行
test_count = 0;
stable_count = 0;
best_error = inf;
best_result = [];

% 結果保存用
all_results = zeros(total_tests, 12); % q1, q2, u1, u2, q1_f, q2_f, u1_f, u2_f, error, stable, energy_loss, converged

if search_mode == 1
    % 角度固定モード（2次元グリッドサーチ）
    for i = 1:n_points_u
        for j = 1:n_points_u
            test_count = test_count + 1;
            
            % 初期値
            q1_init = q1_fixed;
            q2_init = q2_fixed;
            u1_init = u1_range(i);
            u2_init = u2_range(j);
            z0 = [q1_init, u1_init, q2_init, u2_init];
            
            % 進捗表示
            if mod(test_count, 10) == 0
                fprintf('進捗: %d/%d (%.1f%%) - 安定歩行: %d\n', ...
                        test_count, total_tests, 100*test_count/total_tests, stable_count);
            end
            
            try
                % 複数ストライドのシミュレーション（収束性チェック）
                [is_stable, final_state, error, energy_loss, n_steps] = test_periodic_walk(z0, walker);
                
                % 結果の保存
                all_results(test_count, :) = [q1_init, q2_init, u1_init, u2_init, ...
                                              final_state(1), final_state(3), ...
                                              final_state(2), final_state(4), ...
                                              error, is_stable, energy_loss, n_steps];
                
                % 安定な歩行が見つかった場合
                if is_stable && error < best_error
                    best_error = error;
                    best_result = all_results(test_count, :);
                    stable_count = stable_count + 1;
                    
                    fprintf('\n*** 新しい最良解を発見! ***\n');
                    fprintf('初期値: q1=%.6f, q2=%.6f, u1=%.6f, u2=%.6f\n', q1_init, q2_init, u1_init, u2_init);
                    fprintf('最終状態: q1=%.6f, q2=%.6f\n', final_state(1), final_state(3));
                    fprintf('誤差: %.6e\n', error);
                    fprintf('エネルギー損失/ステップ: %.6f\n\n', energy_loss);
                elseif is_stable
                    stable_count = stable_count + 1;
                end
                
            catch ME
                % エラーが発生した場合（発散など）
                all_results(test_count, :) = [q1_init, q2_init, u1_init, u2_init, ...
                                              NaN, NaN, NaN, NaN, inf, 0, NaN, 0];
            end
        end
    end
else
    % 角度範囲モード（4次元グリッドサーチ）
    for i1 = 1:length(q1_range)
        for i2 = 1:length(q2_range)
            for j1 = 1:length(u1_range)
                for j2 = 1:length(u2_range)
                    test_count = test_count + 1;
                    
                    % 初期値
                    q1_init = q1_range(i1);
                    q2_init = q2_range(i2);
                    u1_init = u1_range(j1);
                    u2_init = u2_range(j2);
                    z0 = [q1_init, u1_init, q2_init, u2_init];
                    
                    % 進捗表示
                    if mod(test_count, 50) == 0
                        fprintf('進捗: %d/%d (%.1f%%) - 安定歩行: %d\n', ...
                                test_count, total_tests, 100*test_count/total_tests, stable_count);
                    end
                    
                    try
                        % 複数ストライドのシミュレーション（収束性チェック）
                        [is_stable, final_state, error, energy_loss, n_steps] = test_periodic_walk(z0, walker);
                        
                        % 結果の保存
                        all_results(test_count, :) = [q1_init, q2_init, u1_init, u2_init, ...
                                                      final_state(1), final_state(3), ...
                                                      final_state(2), final_state(4), ...
                                                      error, is_stable, energy_loss, n_steps];
                        
                        % 安定な歩行が見つかった場合
                        if is_stable && error < best_error
                            best_error = error;
                            best_result = all_results(test_count, :);
                            stable_count = stable_count + 1;
                            
                            fprintf('\n*** 新しい最良解を発見! ***\n');
                            fprintf('初期値: q1=%.6f, q2=%.6f, u1=%.6f, u2=%.6f\n', q1_init, q2_init, u1_init, u2_init);
                            fprintf('最終状態: q1=%.6f, q2=%.6f\n', final_state(1), final_state(3));
                            fprintf('誤差: %.6e\n', error);
                            fprintf('エネルギー損失/ステップ: %.6f\n\n', energy_loss);
                        elseif is_stable
                            stable_count = stable_count + 1;
                        end
                        
                    catch ME
                        % エラーが発生した場合（発散など）
                        all_results(test_count, :) = [q1_init, q2_init, u1_init, u2_init, ...
                                                      NaN, NaN, NaN, NaN, inf, 0, NaN, 0];
                    end
                end
            end
        end
    end
end

% タイマー終了
elapsed_time = toc;

% 結果の解析と表示
fprintf('\n=== グリッドサーチ完了 ===\n');
fprintf('総実行時間: %.2f秒\n', elapsed_time);
fprintf('安定歩行の数: %d/%d (%.1f%%)\n', stable_count, total_tests, 100*stable_count/total_tests);

% 有効な結果のみを抽出
valid_idx = all_results(:,10) == 1; % stable == 1
valid_results = all_results(valid_idx, :);

if ~isempty(valid_results)
    % エラーでソート
    [~, sort_idx] = sort(valid_results(:,9));
    valid_results = valid_results(sort_idx, :);
    
    % 上位10個の結果を表示
    fprintf('\n=== 上位10個の安定歩行 ===\n');
    fprintf('順位 |   q1    |   q2    |   u1    |   u2    | q1_final | q2_final |   誤差   | エネルギー損失\n');
    fprintf('-----|---------|---------|---------|---------|----------|----------|----------|---------------\n');
    
    n_display = min(10, size(valid_results, 1));
    for k = 1:n_display
        fprintf('%4d | %7.4f | %7.4f | %7.4f | %7.4f | %8.5f | %8.5f | %.2e | %.6f\n', ...
                k, valid_results(k,1), valid_results(k,2), valid_results(k,3), valid_results(k,4), ...
                valid_results(k,5), valid_results(k,6), ...
                valid_results(k,9), valid_results(k,11));
    end
    
    % 最良解の詳細表示
    if ~isempty(best_result)
        fprintf('\n=== 最良解の詳細 ===\n');
        fprintf('初期値: [%.6f, %.6f, %.6f, %.6f]\n', best_result(1), best_result(3), best_result(2), best_result(4));
        fprintf('最終状態: [%.6f, %.6f, %.6f, %.6f]\n', ...
                best_result(5), best_result(7), best_result(6), best_result(8));
        fprintf('目標との誤差: %.6e\n', best_result(9));
        fprintf('エネルギー損失/ステップ: %.6f\n', best_result(11));
        fprintf('収束ステップ数: %d\n', best_result(12));
        
        % 最良解でのアニメーション表示オプション
        fprintf('\n最良解のアニメーションを表示しますか？ (y/n): ');
        show_anim = input('', 's');
        if strcmpi(show_anim, 'y') || strcmpi(show_anim, 'yes')
            z0_best = [best_result(1), best_result(3), best_result(2), best_result(4)];
            simulate_and_animate(z0_best, walker);
        end
    end
    
    % 結果のプロット
    if search_mode == 1
        plot_grid_search_results_2d(all_results, u1_range, u2_range);
    else
        plot_grid_search_results_4d(all_results, q1_range, q2_range, u1_range, u2_range);
    end
    
    % 結果をファイルに保存
    save_results = input('\n結果をファイルに保存しますか？ (y/n): ', 's');
    if strcmpi(save_results, 'y') || strcmpi(save_results, 'yes')
        filename = sprintf('grid_search_results_flag%d_%s.mat', flag, datestr(now, 'yyyymmdd_HHMMSS'));
        save(filename, 'all_results', 'u1_range', 'u2_range', 'walker', 'best_result');
        fprintf('結果を %s に保存しました。\n', filename);
    end
    
else
    fprintf('\n警告: 安定な歩行が見つかりませんでした。\n');
    fprintf('グリッドの範囲を変更して再試行してください。\n');
end

end

%===================================================================
function [is_stable, final_state, error, avg_energy_loss, n_steps] = test_periodic_walk(z0, walker)
%===================================================================
% 周期的な歩行をテストし、安定性を評価する

max_steps = 50; % 最大ステップ数
tolerance = 1e-4; % 収束判定の許容誤差
n_check = 5; % 収束チェックするステップ数

% 歩行履歴の保存
history = zeros(max_steps, 4);
energy_history = zeros(max_steps, 1);

% 初期エネルギー
energy_prev = calculate_energy(z0, walker);
total_energy_loss = 0;

% 現在の状態
z_current = z0;

% 目標状態（初期角度と同じ）
target_q1 = z0(1);
target_q2 = z0(3);

for step = 1:max_steps
    % 1ステップ実行
    try
        [z_next, energy_loss] = one_step_walk(z_current, walker);
        
        % エネルギー計算
        energy_current = calculate_energy(z_next, walker);
        step_energy_loss = energy_prev - energy_current;
        total_energy_loss = total_energy_loss + abs(step_energy_loss);
        energy_prev = energy_current;
        
        % 履歴保存
        history(step, :) = z_next;
        energy_history(step) = energy_current;
        
        % 発散チェック
        if any(abs(z_next) > 10) || any(isnan(z_next)) || any(isinf(z_next))
            is_stable = false;
            final_state = z_next;
            error = inf;
            avg_energy_loss = NaN;
            n_steps = step;
            return;
        end
        
        % 収束チェック（n_checkステップごと）
        if step >= n_check
            % 目標状態との誤差
            error_q1 = abs(z_next(1) - target_q1);
            error_q2 = abs(z_next(3) - target_q2);
            
            % 速度の変動チェック
            recent_u1 = history(step-n_check+1:step, 2);
            recent_u2 = history(step-n_check+1:step, 4);
            u1_variation = std(recent_u1);
            u2_variation = std(recent_u2);
            
            % 収束判定
            if error_q1 < tolerance && error_q2 < tolerance && ...
               u1_variation < tolerance && u2_variation < tolerance
                is_stable = true;
                final_state = z_next;
                error = sqrt(error_q1^2 + error_q2^2);
                avg_energy_loss = total_energy_loss / step;
                n_steps = step;
                return;
            end
        end
        
        % 次のステップへ
        z_current = z_next;
        
    catch ME
        % エラー発生（数値的な問題など）
        is_stable = false;
        final_state = z_current;
        error = inf;
        avg_energy_loss = NaN;
        n_steps = step;
        return;
    end
end

% 最大ステップ数に達した場合
final_state = z_current;
error_q1 = abs(final_state(1) - target_q1);
error_q2 = abs(final_state(3) - target_q2);
error = sqrt(error_q1^2 + error_q2^2);
avg_energy_loss = total_energy_loss / max_steps;
n_steps = max_steps;

% 安定性の判定（誤差が小さければ安定とみなす）
is_stable = (error < tolerance * 10);

end

%===================================================================
function [z_after, energy_loss] = one_step_walk(z0, walker)
%===================================================================
% 1ステップの歩行をシミュレート

% 初期エネルギー
energy_before = calculate_energy(z0, walker);

% 拡張状態ベクトルの準備
q1 = z0(1); u1 = z0(2); q2 = z0(3); u2 = z0(4);
xh = 0; % ヒップのx座標（簡略化のため0とする）
vxh = (-walker.l*cos(q1)-walker.r)*u1;
yh = walker.l*cos(q1) + walker.r;
vyh = -walker.l*sin(q1)*u1;
TE = energy_before;

z0_extended = [q1 u1 q2 u2 TE xh vxh yh vyh];

% Single stance phase
options = odeset('abstol',1e-13,'reltol',1e-13,'events',@collision);
tspan = [0 5];
[t, z] = ode113(@single_stance, tspan, z0_extended, options, walker);

% Heel strike
if length(t) < 100
    z_after_extended = heelstrike(t(end), z(end,:), walker);
    z_after = z_after_extended(1:4);
    
    % エネルギー損失
    energy_after = calculate_energy(z_after, walker);
    energy_loss = energy_before - energy_after;
else
    error('ステップが完了しませんでした');
end

end

%===================================================================
function plot_grid_search_results_2d(all_results, u1_range, u2_range)
%===================================================================
% 2次元グリッドサーチの結果をプロット（q1,q2固定モード）

figure('Name', 'グリッドサーチ結果（2D）', 'Position', [100, 100, 1200, 800]);

% データの整形
n_u1 = length(u1_range);
n_u2 = length(u2_range);
error_grid = reshape(all_results(:,9), n_u2, n_u1)';
stable_grid = reshape(all_results(:,10), n_u2, n_u1)';

% 対数スケールに変換（エラーが大きい場合の可視化のため）
error_grid_log = log10(error_grid + 1e-10);
error_grid_log(isinf(error_grid_log)) = 10; % inf を大きな値に置換

% エラーマップ
subplot(2,2,1)
imagesc(u2_range, u1_range, error_grid_log)
colorbar
xlabel('u2 初期値')
ylabel('u1 初期値')
title('誤差マップ (log10スケール)')
axis xy
colormap(flipud(hot))

% 安定性マップ
subplot(2,2,2)
imagesc(u2_range, u1_range, stable_grid)
colorbar
xlabel('u2 初期値')
ylabel('u1 初期値')
title('安定性マップ (1=安定, 0=不安定)')
axis xy
colormap([1 0.8 0.8; 0.8 1 0.8])

% 最終角度 q1
subplot(2,2,3)
q1_final_grid = reshape(all_results(:,5), n_u2, n_u1)';
imagesc(u2_range, u1_range, q1_final_grid)
colorbar
xlabel('u2 初期値')
ylabel('u1 初期値')
title('最終 q1 値')
axis xy
colormap(jet)

% 最終角度 q2
subplot(2,2,4)
q2_final_grid = reshape(all_results(:,6), n_u2, n_u1)';
imagesc(u2_range, u1_range, q2_final_grid)
colorbar
xlabel('u2 初期値')
ylabel('u1 初期値')
title('最終 q2 値')
axis xy
colormap(jet)

q1_fixed = all_results(1,1);
q2_fixed = all_results(1,2);
sgtitle(sprintf('グリッドサーチ結果の可視化 (q1=%.3f, q2=%.3f固定)', q1_fixed, q2_fixed))

% 3Dプロット
figure('Name', '3Dエラーサーフェス', 'Position', [200, 200, 800, 600]);
[U1, U2] = meshgrid(u1_range, u2_range);
surf(U1, U2, error_grid_log')
xlabel('u1 初期値')
ylabel('u2 初期値')
zlabel('log10(誤差)')
title(sprintf('誤差の3Dサーフェス (q1=%.3f, q2=%.3f固定)', q1_fixed, q2_fixed))
colorbar
shading interp
view(45, 30)

end

%===================================================================
function plot_grid_search_results_4d(all_results, q1_range, q2_range, u1_range, u2_range)
%===================================================================
% 4次元グリッドサーチの結果をプロット（スライス表示）

% 安定な結果のみを抽出
stable_idx = all_results(:,10) == 1;
stable_results = all_results(stable_idx, :);

if isempty(stable_results)
    fprintf('安定な結果が見つかりませんでした。\n');
    return;
end

figure('Name', 'グリッドサーチ結果（4D）', 'Position', [100, 100, 1400, 900]);

% 1. 初期角度空間での安定領域
subplot(2,3,1)
scatter(stable_results(:,1), stable_results(:,2), 50, log10(stable_results(:,9)+1e-10), 'filled')
colorbar
xlabel('q1 初期値')
ylabel('q2 初期値')
title('安定領域 (q1-q2空間、色:log10誤差)')
grid on

% 2. 初期速度空間での安定領域
subplot(2,3,2)
scatter(stable_results(:,3), stable_results(:,4), 50, log10(stable_results(:,9)+1e-10), 'filled')
colorbar
xlabel('u1 初期値')
ylabel('u2 初期値')
title('安定領域 (u1-u2空間、色:log10誤差)')
grid on

% 3. q1-u1 断面
subplot(2,3,3)
scatter(stable_results(:,1), stable_results(:,3), 50, log10(stable_results(:,9)+1e-10), 'filled')
colorbar
xlabel('q1 初期値')
ylabel('u1 初期値')
title('安定領域 (q1-u1断面、色:log10誤差)')
grid on

% 4. 最終状態の分布
subplot(2,3,4)
scatter(stable_results(:,5), stable_results(:,6), 50, stable_results(:,11), 'filled')
colorbar
xlabel('q1 最終値')
ylabel('q2 最終値')
title('最終状態分布（色:エネルギー損失）')
grid on

% 5. 収束速度の分布
subplot(2,3,5)
scatter(stable_results(:,1), stable_results(:,2), 50, stable_results(:,12), 'filled')
colorbar
xlabel('q1 初期値')
ylabel('q2 初期値')
title('収束ステップ数')
grid on

% 6. エラーのヒストグラム
subplot(2,3,6)
histogram(log10(stable_results(:,9)+1e-10), 20)
xlabel('log10(誤差)')
ylabel('頻度')
title('誤差分布')
grid on

sgtitle('4次元グリッドサーチ結果の可視化')

% 追加の3Dプロット
figure('Name', '3D安定領域', 'Position', [300, 200, 800, 600]);
scatter3(stable_results(:,1), stable_results(:,2), stable_results(:,3), ...
         50, log10(stable_results(:,9)+1e-10), 'filled')
xlabel('q1 初期値')
ylabel('q2 初期値')
zlabel('u1 初期値')
title('3D安定領域（色:log10誤差）')
colorbar
view(45, 30)
grid on

end

%===================================================================
function simulate_and_animate(z0, walker)
%===================================================================
% 選択された初期値で歩行をシミュレートしてアニメーション表示

fprintf('\n選択された初期値での歩行シミュレーション...\n');

% 複数ステップのシミュレーション
n_steps = 10;
z_current = z0;

figure('Name', '連続歩行アニメーション', 'Position', [300, 100, 800, 600]);

for step = 1:n_steps
    fprintf('ステップ %d/%d を実行中...\n', step, n_steps);
    
    % 1ステップ実行
    try
        % 拡張状態ベクトルの準備
        q1 = z_current(1); u1 = z_current(2); 
        q2 = z_current(3); u2 = z_current(4);
        xh = 0;
        vxh = (-walker.l*cos(q1)-walker.r)*u1;
        yh = walker.l*cos(q1) + walker.r;
        vyh = -walker.l*sin(q1)*u1;
        TE = calculate_energy(z_current, walker);
        
        z_extended = [q1 u1 q2 u2 TE xh vxh yh vyh];
        
        % Single stance phase
        options = odeset('abstol',1e-13,'reltol',1e-13,'events',@collision);
        tspan = [0 5];
        [t, z] = ode113(@single_stance, tspan, z_extended, options, walker);
        
        % アニメーション表示
        animate_step(t, z, walker, step);
        
        % Heel strike
        if length(t) < 100
            z_after_extended = heelstrike(t(end), z(end,:), walker);
            z_current = z_after_extended(1:4);
            
            fprintf('  状態: [%.4f, %.4f, %.4f, %.4f]\n', z_current);
        else
            fprintf('  警告: ステップが完了しませんでした\n');
            break;
        end
        
    catch ME
        fprintf('  エラー: %s\n', ME.message);
        break;
    end
    
    pause(0.5);
end

fprintf('シミュレーション完了\n');

end

%===================================================================
function animate_step(t, z, walker, step_num)
%===================================================================
% 1ステップのアニメーション（簡略版）

% 最後のフレームのみ表示
q1 = z(end,1); q2 = z(end,3);
xh = z(end,6); yh = z(end,8);

clf
hold on

% 地面
plot([-2, 2], [0, 0], 'k-', 'LineWidth', 2);

% ヒップ
plot(xh, yh, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

% スタンス脚
x1 = xh + walker.l*sin(q1);
y1 = yh - walker.l*cos(q1);
plot([xh, x1], [yh, y1], 'b-', 'LineWidth', 3);
plot(x1, 0, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

% スイング脚
x2 = xh + walker.l*sin(q2);
y2 = yh - walker.l*cos(q2);
plot([xh, x2], [yh, y2], 'g-', 'LineWidth', 3);

xlim([-1, 1]);
ylim([-0.2, 1.2]);
axis equal
title(sprintf('ステップ %d - q1=%.3f, q2=%.3f', step_num, q1, q2));
drawnow

end

%=== 以下、元のコードから必要な関数をコピー ===

%===================================================================
function TE = calculate_energy(z, walker)
%===================================================================
q1 = z(1); u1 = z(2); q2 = z(3); u2 = z(4);
M = walker.M;  m = walker.m; I = walker.I;   
l = walker.l;  c = walker.c; w = walker.w;   
r = walker.r;  g = walker.g; gam = walker.gam;

TE = 1/2*m*(((-l*cos(q1)-r)*u1-u1*(-c*cos(q1)+w*sin(q1)))^2+(-l*sin(q1)*u1+u1*(c*sin(q1)+w*cos(q1)))^2)+1/2*m*(((-l*cos(q1)-r)*u1-(u1-u2)*(-c*cos(q1-q2)+w*sin(q1-q2)))^2+(-l*sin(q1)*u1+(u1-u2)*(c*sin(q1-q2)+w*cos(q1-q2)))^2)+1/2*M*((-l*cos(q1)-r)^2*u1^2+l^2*sin(q1)^2*u1^2)+1/2*I*(u1^2+(u1-u2)^2)+2*m*g*cos(gam)*r+2*m*g*l*cos(gam-q1)-m*g*c*cos(gam-q1)-m*g*w*sin(gam-q1)+2*m*g*sin(gam)*r*q1-m*g*c*cos(gam-q1+q2)-m*g*w*sin(gam-q1+q2)+M*g*cos(gam)*r+M*g*l*cos(gam-q1)+M*g*sin(gam)*r*q1;
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