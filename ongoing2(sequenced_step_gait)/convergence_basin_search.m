function convergence_basin_search(params)
% 実際に収束する初期値を探索する関数（並列処理強化版 + 計算時間推定機能付き）
% params: 探索パラメータ構造体（grid_runnerと同じ形式）

    fprintf('\n=== 収束保証付き初期値探索 ===\n');
    
    % Walker設定
    walker.M = 1000; walker.m = 1.0; walker.I = 0.00; walker.l = 1.0; walker.w = 0.0; 
    walker.c = 1.0;  walker.r = 0.3; walker.g = 1.0; walker.gam = 0.009;
    
    % 探索範囲の詳細表示
    fprintf('\n📋 探索範囲:\n');
    fprintf('q1 (スタンス脚角度):  %.3f から %.3f まで %.3f 刻み (%d点)\n', ...
            params.q1_min, params.q1_max, params.q1_step, length(params.q1_range));
    fprintf('u1 (スタンス脚角速度): %.3f から %.3f まで %.3f 刻み (%d点)\n', ...
            params.u1_min, params.u1_max, params.u1_step, length(params.u1_range));
    fprintf('q2 (スイング脚角度):  %.3f から %.3f まで %.3f 刻み (%d点)\n', ...
            params.q2_min, params.q2_max, params.q2_step, length(params.q2_range));
    fprintf('u2 (スイング脚角速度): %.3f から %.3f まで %.3f 刻み (%d点)\n', ...
            params.u2_min, params.u2_max, params.u2_step, length(params.u2_range));
    fprintf('総探索数: %d\n', params.total_combinations);
    
    % 探索範囲の実際の値も表示
    fprintf('\n🔍 実際の探索値:\n');
    fprintf('q1: ['); fprintf('%.3f ', params.q1_range); fprintf(']\n');
    fprintf('u1: ['); fprintf('%.3f ', params.u1_range); fprintf(']\n');
    fprintf('q2: ['); fprintf('%.3f ', params.q2_range); fprintf(']\n');
    fprintf('u2: ['); fprintf('%.3f ', params.u2_range); fprintf(']\n');
    
    % 全組み合わせ作成
    [Q1, U1, Q2, U2] = ndgrid(params.q1_range, params.u1_range, ...
                             params.q2_range, params.u2_range);
    all_conditions = [Q1(:), U1(:), Q2(:), U2(:)];
    total = size(all_conditions, 1);
    
    % ⏱️ 計算時間推定
    fprintf('\n⏱️ 計算時間推定中...\n');
    
    % サンプリングで実行時間を推定（最大10個のサンプル）
    sample_size = min(10, total);
    sample_indices = randsample(total, sample_size);
    
    % 収束テストパラメータ
    test_steps = 30;  % 収束判定のための歩数
    convergence_threshold = 1e-3;  % 収束判定閾値
    
    % サンプル計算を実行
    fprintf('サンプル計算実行中（%d個）...', sample_size);
    tic;
    for i = 1:sample_size
        z0 = all_conditions(sample_indices(i), :);
        test_convergence_from_initial(z0, walker, test_steps, convergence_threshold);
    end
    sample_time = toc;
    fprintf(' 完了\n');
    
    % 1個あたりの平均時間を計算
    time_per_condition = sample_time / sample_size;
    
    % 並列処理の確認と設定
    use_parallel = true;
    pool = gcp('nocreate');
    num_workers = 1;
    
    if isempty(pool)
        if total > 50  % 50個以上の場合のみ並列処理を提案
            % 推定時間を表示して判断材料を提供
            estimated_sequential_time = time_per_condition * total;
            fprintf('\n📊 推定実行時間:\n');
            fprintf('  - 逐次処理: %.1f秒 (%.1f分)\n', estimated_sequential_time, estimated_sequential_time/60);
            
            % 利用可能なコア数を取得
            max_workers = feature('numcores');
            typical_workers = min(max_workers - 1, 8);  % 通常は最大8ワーカー
            estimated_parallel_time = time_per_condition * total / typical_workers * 1.2;  % 1.2はオーバーヘッド
            fprintf('  - 並列処理（%dワーカー想定）: %.1f秒 (%.1f分)\n', ...
                    typical_workers, estimated_parallel_time, estimated_parallel_time/60);
            
            answer = input('\n並列処理を使用しますか？ (y/n) [推奨: y]: ', 's');
            if isempty(answer) || strcmpi(answer, 'y')
                fprintf('並列プールを起動中...\n');
                parpool;
                pool = gcp('nocreate');
                num_workers = pool.NumWorkers;
                use_parallel = true;
            else
                use_parallel = false;
            end
        else
            use_parallel = false;
            estimated_time = time_per_condition * total;
            fprintf('\n📊 推定実行時間: %.1f秒\n', estimated_time);
        end
    else
        num_workers = pool.NumWorkers;
        fprintf('既存の並列プール（ワーカー数: %d）を使用\n', num_workers);
        use_parallel = true;
    end
    
    % 最終的な推定時間を計算・表示
    if use_parallel
        estimated_time = time_per_condition * total / num_workers * 1.2;  % オーバーヘッド込み
        fprintf('\n⏱️ 最終推定実行時間（並列%dワーカー）: %.1f秒 (%.1f分)\n', ...
                num_workers, estimated_time, estimated_time/60);
    else
        estimated_time = time_per_condition * total;
        fprintf('\n⏱️ 最終推定実行時間（逐次処理）: %.1f秒 (%.1f分)\n', ...
                estimated_time, estimated_time/60);
    end
    
    % 実績データに基づく補正（220個で20秒の実績から）
    reference_rate = 20 / 220;  % 秒/個
    if abs(time_per_condition - reference_rate) > reference_rate * 0.5
        fprintf('💡 参考: 前回の実績では220個で約20秒でした（%.3f秒/個）\n', reference_rate);
        alternative_estimate = reference_rate * total;
        if use_parallel
            alternative_estimate = alternative_estimate / num_workers * 1.2;
        end
        fprintf('   実績ベースの推定: %.1f秒 (%.1f分)\n', alternative_estimate, alternative_estimate/60);
    end
    
    % 続行確認
    if estimated_time > 300  % 5分以上かかる場合
        fprintf('\n⚠️ 計算に%.1f分以上かかる見込みです。\n', estimated_time/60);
        continue_answer = input('続行しますか？ (y/n) [n]: ', 's');
        if ~strcmpi(continue_answer, 'y')
            fprintf('処理を中止しました。\n');
            return;
        end
    end
    
    % 結果保存用
    convergent_conditions = [];
    convergence_info = [];
    
    % 進捗表示
    if use_parallel
        fprintf('\n🚀 並列処理で収束テスト実行中...\n');
    else
        fprintf('\n🚀 逐次処理で収束テスト実行中...\n');
    end
    
    % 実際の計算開始
    actual_start_time = tic;
    
    if use_parallel
        % 並列処理（改良版）
        fprintf('並列処理中（%d個の初期条件を評価）...\n', total);
        results_cell = cell(total, 1);
        
        % プログレスバー的な表示のため、一定間隔で進捗を表示
        fprintf('進捗: ');
        
        parfor idx = 1:total
            z0 = all_conditions(idx, :);
            result = test_convergence_from_initial(z0, walker, test_steps, convergence_threshold);
            if result.converged
                results_cell{idx} = result;
            end
        end
        
        fprintf('完了！\n');
        fprintf('並列処理完了。結果を集計中...\n');
        
        % 結果集計
        success_count = 0;
        for idx = 1:total
            if ~isempty(results_cell{idx})
                result = results_cell{idx};
                convergent_conditions = [convergent_conditions; result.initial_condition];
                convergence_info = [convergence_info; result];
                success_count = success_count + 1;
                
                fprintf('✅ 収束確認 #%d: [%.3f, %.3f, %.3f, %.3f] → 固定点まで%.4f (収束: %d歩)\n', ...
                        success_count, result.initial_condition, result.final_distance, result.steps_to_converge);
            end
        end
    else
        % 逐次処理（進捗表示付き）
        progress_interval = max(1, floor(total/20));
        last_progress_time = tic;
        
        for idx = 1:total
            if mod(idx, progress_interval) == 0
                elapsed = toc(actual_start_time);
                progress_percent = 100 * idx / total;
                remaining_time = elapsed / idx * (total - idx);
                fprintf('進捗: %d/%d (%.1f%%) - 経過: %.1f秒, 残り推定: %.1f秒\n', ...
                        idx, total, progress_percent, elapsed, remaining_time);
            end
            
            z0 = all_conditions(idx, :);
            result = test_convergence_from_initial(z0, walker, test_steps, convergence_threshold);
            
            if result.converged
                convergent_conditions = [convergent_conditions; result.initial_condition];
                convergence_info = [convergence_info; result];
                
                fprintf('✅ 収束確認 #%d: [%.3f, %.3f, %.3f, %.3f] → 固定点まで%.4f (収束: %d歩)\n', ...
                        size(convergent_conditions, 1), z0, result.final_distance, result.steps_to_converge);
            end
        end
    end
    
    actual_elapsed = toc(actual_start_time);
    
    % 結果まとめ
    num_convergent = size(convergent_conditions, 1);
    fprintf('\n%s\n', repmat('=', 1, 60));
    fprintf('🎉 収束領域探索完了！\n');
    fprintf('実行時間: %.1f秒 (推定: %.1f秒, 誤差: %.1f%%)\n', ...
            actual_elapsed, estimated_time, 100*abs(actual_elapsed-estimated_time)/estimated_time);
    if use_parallel && ~isempty(pool)
        fprintf('並列処理（%dワーカー）で高速化されました\n', pool.NumWorkers);
    end
    fprintf('収束する初期値: %d/%d (%.1f%%)\n', num_convergent, total, 100*num_convergent/total);
    fprintf('%s\n', repmat('=', 1, 60));
    
    if num_convergent > 0
        % 収束速度でソート（早く収束する順）
        [~, sort_idx] = sort([convergence_info.steps_to_converge]);
        sorted_conditions = convergent_conditions(sort_idx, :);
        sorted_info = convergence_info(sort_idx);
        
        % 最良の条件を表示
        fprintf('\n🏆 【最速収束する初期値】\n');
        best = sorted_info(1);
        fprintf('初期値: [%.4f, %.4f, %.4f, %.4f]\n', best.initial_condition);
        fprintf('固定点: [%.4f, %.4f, %.4f, %.4f]\n', best.fixed_point);
        fprintf('収束歩数: %d歩\n', best.steps_to_converge);
        fprintf('最終誤差: %.6f\n', best.final_distance);
        
        % 上位10個を表示
        fprintf('\n📊 【収束が速い初期値TOP10】\n');
        fprintf('No. | q1      | u1      | q2      | u2      | 収束歩数 | 最終誤差\n');
        fprintf('----|---------|---------|---------|---------|----------|----------\n');
        
        display_count = min(10, num_convergent);
        for i = 1:display_count
            info = sorted_info(i);
            fprintf('%3d | %7.4f | %7.4f | %7.4f | %7.4f | %8d | %.2e\n', ...
                    i, info.initial_condition, info.steps_to_converge, info.final_distance);
        end
        
        if num_convergent > 10
            fprintf('... 他 %d 個の収束条件\n', num_convergent - 10);
        end
        
        % ワークスペースに保存
        assignin('base', 'convergent_initials', sorted_conditions);
        assignin('base', 'convergence_details', sorted_info);
        assignin('base', 'best_convergent_initial', best);
        assignin('base', 'convergence_params', params);
        
        fprintf('\n📁 結果がワークスペースに保存されました:\n');
        fprintf('- convergent_initials: 収束する全初期値（収束速度順）\n');
        fprintf('- convergence_details: 収束の詳細情報\n');
        fprintf('- best_convergent_initial: 最速収束する初期値\n');
        fprintf('- convergence_params: 探索パラメータ\n');
        
        % 🎨 可視化
        fprintf('\n🎨 結果を可視化中...\n');
        visualize_convergence_results(sorted_info, all_conditions, params);
        
        % 💾 結果をスプレッドシートに保存
        fprintf('\n💾 結果をスプレッドシートに保存中...\n');
        save_convergence_results_to_spreadsheet(sorted_info, all_conditions, params);
        
    else
        fprintf('\n❌ 収束する初期値が見つかりませんでした。\n');
        fprintf('💡 対策案:\n');
        fprintf('  - 探索範囲を論文の値 [0.2, -0.2, 0.4, -0.3] の近くに設定\n');
        fprintf('  - 刻み幅を細かくする\n');
        fprintf('  - 収束判定の歩数を増やす（現在: %d歩）\n', test_steps);
    end
end

function result = test_convergence_from_initial(z0, walker, max_steps, threshold)
% 指定した初期値から実際に歩行させて収束をテスト
    
    result.initial_condition = z0;
    result.converged = false;
    result.fixed_point = [NaN, NaN, NaN, NaN];
    result.steps_to_converge = NaN;
    result.final_distance = NaN;
    result.convergence_history = [];
    
    try
        % まず固定点を見つける
        options = optimset('TolFun',1e-12,'TolX',1e-12,'Display','off');
        [zstar, fval, exitflag] = fsolve(@(z) fixedpt(z, walker), z0, options);
        
        if exitflag ~= 1 || norm(fval) > 1e-6
            return;  % 固定点が見つからない
        end
        
        % 固定点の安定性をチェック
        J = partialder(@(z) onestep(z, walker), zstar, walker);
        eigenvalues = eig(J);
        if max(abs(eigenvalues)) >= 1
            return;  % 不安定な固定点
        end
        
        result.fixed_point = zstar;
        
        % 実際に初期値から歩行させる
        current_state = z0;
        convergence_history = zeros(max_steps, 1);
        
        for step = 1:max_steps
            % 1歩進める
            [next_state, ~] = onestep(current_state, walker, 1);
            
            % 最後の状態（1歩後）を取得
            if size(next_state, 1) > 1
                current_state = next_state(end, 1:4);
            else
                current_state = next_state;
            end
            
            % 固定点からの距離を計算
            distance = norm(current_state - zstar);
            convergence_history(step) = distance;
            
            % 収束判定
            if distance < threshold
                result.converged = true;
                result.steps_to_converge = step;
                result.final_distance = distance;
                result.convergence_history = convergence_history(1:step);
                return;
            end
            
            % 発散判定（距離が増大し続ける場合）
            if distance > 10 || any(~isfinite(current_state))
                return;  % 発散
            end
        end
        
        % max_steps後も収束しなかった場合
        result.final_distance = distance;
        result.convergence_history = convergence_history;
        
    catch
        % エラーが発生（転倒など）
        return;
    end
end

function visualize_convergence_results(convergence_info, all_conditions, params)
% 収束結果の可視化
    
    figure('Name', '収束する初期値の分布', 'Position', [100, 100, 1200, 800]);
    
    % 収束する初期値を抽出
    conv_initials = vertcat(convergence_info.initial_condition);
    conv_steps = [convergence_info.steps_to_converge];
    
    % 全条件のデータ
    all_q1 = all_conditions(:, 1);
    all_u1 = all_conditions(:, 2);
    all_q2 = all_conditions(:, 3);
    all_u2 = all_conditions(:, 4);
    
    % サブプロット1: 収束速度の分布（色で表示）
    subplot(2, 2, 1);
    % 失敗した点を薄く表示
    scatter3(all_q1, all_u1, all_q2, 15, [0.8 0.8 0.8], 'o', 'MarkerFaceAlpha', 0.3);
    hold on;
    % 収束する点を色分けして表示
    scatter3(conv_initials(:,1), conv_initials(:,2), conv_initials(:,3), ...
             80, conv_steps, 'filled', 'MarkerEdgeColor', 'k');
    xlabel('q1 (スタンス脚角度)');
    ylabel('u1 (スタンス脚角速度)');
    zlabel('q2 (スイング脚角度)');
    title('収束する初期値（色：収束歩数）');
    c = colorbar;
    ylabel(c, '収束歩数');
    grid on;
    view(45, 30);
    % 最速収束点をハイライト
    plot3(conv_initials(1,1), conv_initials(1,2), conv_initials(1,3), ...
          'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'red', 'LineWidth', 2);
    hold off;
    
    % サブプロット2: 収束歩数のヒストグラム
    subplot(2, 2, 2);
    histogram(conv_steps, 20, 'FaceColor', [0.3 0.6 0.9]);
    xlabel('収束までの歩数');
    ylabel('頻度');
    title('収束速度の分布');
    grid on;
    
    % サブプロット3: 初期値と固定点の関係
    subplot(2, 2, 3);
    fixed_points = vertcat(convergence_info.fixed_point);
    initial_distances = zeros(length(convergence_info), 1);
    for i = 1:length(convergence_info)
        initial_distances(i) = norm(convergence_info(i).initial_condition - ...
                                   convergence_info(i).fixed_point);
    end
    scatter(initial_distances, conv_steps, 50, 'filled');
    xlabel('初期値と固定点の距離');
    ylabel('収束歩数');
    title('初期距離 vs 収束速度');
    grid on;
    
    % サブプロット4: 収束履歴の例（最速5つ）
    subplot(2, 2, 4);
    hold on;
    colors = lines(5);
    for i = 1:min(5, length(convergence_info))
        history = convergence_info(i).convergence_history;
        plot(1:length(history), history, 'LineWidth', 2, 'Color', colors(i,:));
    end
    xlabel('歩数');
    ylabel('固定点からの距離');
    title('収束履歴（最速5例）');
    set(gca, 'YScale', 'log');
    grid on;
    legend('1位', '2位', '3位', '4位', '5位', 'Location', 'northeast');
    hold off;
    
    % 統計情報を表示
    fprintf('\n📊 収束統計:\n');
    fprintf('収束条件数: %d\n', length(conv_steps));
    fprintf('平均収束歩数: %.1f歩\n', mean(conv_steps));
    fprintf('最速収束: %d歩\n', min(conv_steps));
    fprintf('最遅収束: %d歩\n', max(conv_steps));
    fprintf('標準偏差: %.1f歩\n', std(conv_steps));
    
    % パラメータ範囲
    fprintf('\n🎯 収束する初期値の範囲:\n');
    fprintf('q1: %.3f ～ %.3f (平均: %.3f)\n', min(conv_initials(:,1)), max(conv_initials(:,1)), mean(conv_initials(:,1)));
    fprintf('u1: %.3f ～ %.3f (平均: %.3f)\n', min(conv_initials(:,2)), max(conv_initials(:,2)), mean(conv_initials(:,2)));
    fprintf('q2: %.3f ～ %.3f (平均: %.3f)\n', min(conv_initials(:,3)), max(conv_initials(:,3)), mean(conv_initials(:,3)));
    fprintf('u2: %.3f ～ %.3f (平均: %.3f)\n', min(conv_initials(:,4)), max(conv_initials(:,4)), mean(conv_initials(:,4)));
    
    fprintf('\n💡 図の見方:\n');
    fprintf('- グレーの点: 収束しない条件\n');
    fprintf('- カラーの点: 収束する条件（色は収束歩数）\n');
    fprintf('- 赤いダイヤ: 最速で収束する条件\n');
    fprintf('- 3D図は回転・ズーム可能です\n');
end

function save_convergence_results_to_spreadsheet(convergence_info, all_conditions, params)
% 収束結果をスプレッドシートに保存
    
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    
    fprintf('\n💾 収束結果をスプレッドシートに保存中...\n');
    
    %% 1. 収束条件の詳細結果を保存
    if length(convergence_info) > 0
        % 収束結果のテーブル作成
        convergence_table = table();
        convergence_table.No = (1:length(convergence_info))';
        
        % 初期値
        initials = vertcat(convergence_info.initial_condition);
        convergence_table.q1_initial = initials(:, 1);
        convergence_table.u1_initial = initials(:, 2);
        convergence_table.q2_initial = initials(:, 3);
        convergence_table.u2_initial = initials(:, 4);
        
        % 固定点
        fixed = vertcat(convergence_info.fixed_point);
        convergence_table.q1_fixed = fixed(:, 1);
        convergence_table.u1_fixed = fixed(:, 2);
        convergence_table.q2_fixed = fixed(:, 3);
        convergence_table.u2_fixed = fixed(:, 4);
        
        % 収束情報
        convergence_table.steps_to_converge = [convergence_info.steps_to_converge]';
        convergence_table.final_distance = [convergence_info.final_distance]';
        
        % 初期距離
        initial_distances = zeros(length(convergence_info), 1);
        for i = 1:length(convergence_info)
            initial_distances(i) = norm(initials(i,:) - fixed(i,:));
        end
        convergence_table.initial_distance = initial_distances;
        
        % 収束速度評価を追加
        speed_rating = cell(length(convergence_info), 1);
        for i = 1:length(convergence_info)
            steps = convergence_info(i).steps_to_converge;
            if steps <= 5
                speed_rating{i} = '超高速';
            elseif steps <= 10
                speed_rating{i} = '高速';
            elseif steps <= 20
                speed_rating{i} = '標準';
            else
                speed_rating{i} = '低速';
            end
        end
        convergence_table.speed_rating = speed_rating;
        
        % 収束結果を保存
        convergence_filename = sprintf('walker_convergence_results_%s.xlsx', timestamp);
        try
            writetable(convergence_table, convergence_filename, 'Sheet', 'Convergence_Results');
            fprintf('✅ 収束条件を保存: %s\n', convergence_filename);
        catch
            % Excelが使えない場合はCSVで保存
            convergence_filename_csv = sprintf('walker_convergence_results_%s.csv', timestamp);
            writetable(convergence_table, convergence_filename_csv);
            fprintf('✅ 収束条件を保存: %s\n', convergence_filename_csv);
        end
    end
    
    %% 2. 統計情報を保存
    % 統計情報のテーブル作成
    stats_table = table();
    
    % 探索パラメータ
    param_names = {'q1_min'; 'q1_max'; 'q1_step'; 'q1_points'; ...
                   'u1_min'; 'u1_max'; 'u1_step'; 'u1_points'; ...
                   'q2_min'; 'q2_max'; 'q2_step'; 'q2_points'; ...
                   'u2_min'; 'u2_max'; 'u2_step'; 'u2_points'; ...
                   'total_combinations'; 'convergent_count'; 'convergence_rate_percent'};
    
    param_values = [params.q1_min; params.q1_max; params.q1_step; length(params.q1_range); ...
                    params.u1_min; params.u1_max; params.u1_step; length(params.u1_range); ...
                    params.q2_min; params.q2_max; params.q2_step; length(params.q2_range); ...
                    params.u2_min; params.u2_max; params.u2_step; length(params.u2_range); ...
                    params.total_combinations; length(convergence_info); ...
                    100 * length(convergence_info) / params.total_combinations];
    
    stats_table.Parameter = param_names;
    stats_table.Value = param_values;
    
    % 収束統計を追加
    if length(convergence_info) > 0
        conv_steps = [convergence_info.steps_to_converge]';
        conv_stats_names = {'best_initial_q1'; 'best_initial_u1'; 'best_initial_q2'; 'best_initial_u2'; ...
                           'best_steps'; 'mean_steps'; 'std_steps'; 'min_steps'; 'max_steps'};
        
        best = convergence_info(1);
        conv_stats_values = [best.initial_condition'; ...
                            best.steps_to_converge; ...
                            mean(conv_steps); std(conv_steps); ...
                            min(conv_steps); max(conv_steps)];
        
        conv_stats_table = table(conv_stats_names, conv_stats_values, ...
                                'VariableNames', {'Statistic', 'Value'});
        
        % 統計を保存
        stats_filename = sprintf('walker_convergence_stats_%s.xlsx', timestamp);
        try
            writetable(stats_table, stats_filename, 'Sheet', 'Search_Parameters');
            writetable(conv_stats_table, stats_filename, 'Sheet', 'Convergence_Statistics');
            fprintf('✅ 統計情報を保存: %s\n', stats_filename);
        catch
            % CSVで保存
            stats_filename_csv = sprintf('walker_convergence_stats_%s.csv', timestamp);
            writetable(stats_table, stats_filename_csv);
            conv_stats_filename_csv = sprintf('walker_convergence_conv_stats_%s.csv', timestamp);
            writetable(conv_stats_table, conv_stats_filename_csv);
            fprintf('✅ 統計情報を保存: %s, %s\n', stats_filename_csv, conv_stats_filename_csv);
        end
    end
    
    %% 3. 要約レポートを作成
    report_filename = sprintf('walker_convergence_report_%s.txt', timestamp);
    fid = fopen(report_filename, 'w');
    
    fprintf(fid, '=== パッシブウォーカー 収束領域探索レポート ===\n');
    fprintf(fid, '実行日時: %s\n\n', datestr(now));
    
    fprintf(fid, '【探索パラメータ】\n');
    fprintf(fid, 'q1 (スタンス脚角度): %.3f ～ %.3f (刻み: %.3f, %d点)\n', ...
            params.q1_min, params.q1_max, params.q1_step, length(params.q1_range));
    fprintf(fid, 'u1 (スタンス脚角速度): %.3f ～ %.3f (刻み: %.3f, %d点)\n', ...
            params.u1_min, params.u1_max, params.u1_step, length(params.u1_range));
    fprintf(fid, 'q2 (スイング脚角度): %.3f ～ %.3f (刻み: %.3f, %d点)\n', ...
            params.q2_min, params.q2_max, params.q2_step, length(params.q2_range));
    fprintf(fid, 'u2 (スイング脚角速度): %.3f ～ %.3f (刻み: %.3f, %d点)\n', ...
            params.u2_min, params.u2_max, params.u2_step, length(params.u2_range));
    fprintf(fid, '総探索数: %d\n\n', params.total_combinations);
    
    fprintf(fid, '【結果サマリー】\n');
    fprintf(fid, '収束する初期値数: %d / %d (%.1f%%)\n', ...
            length(convergence_info), params.total_combinations, ...
            100 * length(convergence_info) / params.total_combinations);
    
    if length(convergence_info) > 0
        fprintf(fid, '\n【最速収束する初期値】\n');
        best = convergence_info(1);
        fprintf(fid, '初期値: [%.6f, %.6f, %.6f, %.6f]\n', best.initial_condition);
        fprintf(fid, '収束歩数: %d\n', best.steps_to_converge);
        fprintf(fid, '最終誤差: %.6e\n', best.final_distance);
        fprintf(fid, '固定点: [%.6f, %.6f, %.6f, %.6f]\n', best.fixed_point);
        
        fprintf(fid, '\n【収束統計】\n');
        conv_steps = [convergence_info.steps_to_converge];
        fprintf(fid, '平均収束歩数: %.1f\n', mean(conv_steps));
        fprintf(fid, '標準偏差: %.1f\n', std(conv_steps));
        fprintf(fid, '最速: %d歩, 最遅: %d歩\n', min(conv_steps), max(conv_steps));
    end
    
    fclose(fid);
    fprintf('✅ レポートを保存: %s\n', report_filename);
    
    fprintf('\n📁 保存完了！以下のファイルが作成されました:\n');
    if length(convergence_info) > 0
        fprintf('  📊 収束結果: %s\n', convergence_filename);
    end
    fprintf('  📈 統計: %s\n', stats_filename);
    fprintf('  📄 レポート: %s\n', report_filename);
end

% 必要な補助関数（grid_runnerからコピー）
function zdiff = fixedpt(z0, walker)
    zdiff = onestep(z0, walker) - z0; 
end

function J = partialder(FUN, z, walker)
    pert = 1e-5;
    n = length(z);
    J = zeros(n, n);
    
    for i = 1:n
        ztemp1 = z; ztemp2 = z;
        ztemp1(i) = ztemp1(i) + pert;
        ztemp2(i) = ztemp2(i) - pert;
        J(:,i) = (FUN(ztemp1) - FUN(ztemp2));
    end
    J = J / (2*pert);
end

function [z, t] = onestep(z0, walker, steps)
    M = walker.M;  m = walker.m; I = walker.I;   
    l = walker.l;  c = walker.c; w = walker.w;   
    r = walker.r;  g = walker.g; gam = walker.gam;

    flag = 1;
    if nargin < 2
        error('need more inputs to onestep');
    elseif nargin < 3
        flag = 0;
        steps = 1;
    end

    q1 = z0(1);
    u1 = z0(2);
    q2 = z0(3);
    u2 = z0(4);

    TE = 1/2*m*(((-l*cos(q1)-r)*u1-u1*(-c*cos(q1)+w*sin(q1)))^2+(-l*sin(q1)*u1+u1*(c*sin(q1)+w*cos(q1)))^2)+1/2*m*(((-l*cos(q1)-r)*u1-(u1-u2)*(-c*cos(q1-q2)+w*sin(q1-q2)))^2+(-l*sin(q1)*u1+(u1-u2)*(c*sin(q1-q2)+w*cos(q1-q2)))^2)+1/2*M*((-l*cos(q1)-r)^2*u1^2+l^2*sin(q1)^2*u1^2)+1/2*I*(u1^2+(u1-u2)^2)+2*m*g*cos(gam)*r+2*m*g*l*cos(gam-q1)-m*g*c*cos(gam-q1)-m*g*w*sin(gam-q1)+2*m*g*sin(gam)*r*q1-m*g*c*cos(gam-q1+q2)-m*g*w*sin(gam-q1+q2)+M*g*cos(gam)*r+M*g*l*cos(gam-q1)+M*g*sin(gam)*r*q1; 
    xp1 = 0;
    xh = -l*sin(q1) - r*q1 + xp1;
    vxh = (-l*cos(q1)-r)*u1; 
    yh =  l*cos(q1) + r;
    vyh = -l*sin(q1)*u1; 
    
    z0 = [q1 u1 q2 u2 TE xh vxh yh vyh];

    t0 = 0; 
    dt = 5;
    time_stamps = 100;
    t_ode = t0;
    z_ode = z0;

    for i = 1:steps
        options = odeset('abstol',1e-13,'reltol',1e-13,'events',@(t,z) collision(t,z,walker));
        tspan = linspace(t0,t0+dt,time_stamps);
        [t_temp, z_temp] = ode113(@(t,z) single_stance(t,z,walker), tspan, z0, options);
        
        zplus = heelstrike(t_temp(end),z_temp(end,:),walker); 
        
        z0 = zplus;
        t0 = t_temp(end);
        
        t_ode = [t_ode; t_temp(2:end)];
        z_ode = [z_ode; z_temp(2:end,:)];
    end

    z = zplus(1:4);

    if flag == 1
       z = z_ode;
       t = t_ode;
    end
end

function zdot = single_stance(t,z,walker)  
    q1 = z(1);   u1 = z(2);                         
    q2 = z(3);   u2 = z(4);                         
    xh = z(6);  vxh = z(7);                       
    yh = z(8);  vyh = z(9);                     

    M = walker.M;  m = walker.m; I = walker.I;   
    l = walker.l;  c = walker.c; w = walker.w;   
    r = walker.r;  g = walker.g; gam = walker.gam;

    Th = 0;               

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

function [gstop, isterminal,direction] = collision(t,z,walker)
    q1 = z(1); q2 = z(3); 

    gstop = -q2 + 2*q1;
    if (q2 > -0.05)
        isterminal = 0;
    else
        isterminal = 1;
    end
    direction = -1;
end

function zplus = heelstrike(t,z,walker)      
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

    u1 = X(1);                                       
    u2 = X(2);                                      

    TE = 1/2*m*(((-l*cos(q1)-r)*u1-u1*(-c*cos(q1)+w*sin(q1)))^2+(-l*sin(q1)*u1+u1*(c*sin(q1)+w*cos(q1)))^2)+1/2*m*(((-l*cos(q1)-r)*u1-(u1-u2)*(-c*cos(q1-q2)+w*sin(q1-q2)))^2+(-l*sin(q1)*u1+(u1-u2)*(c*sin(q1-q2)+w*cos(q1-q2)))^2)+1/2*M*((-l*cos(q1)-r)^2*u1^2+l^2*sin(q1)^2*u1^2)+1/2*I*(u1^2+(u1-u2)^2)+2*g*m*cos(gam)*r+2*m*g*l*cos(gam-q1)-m*g*c*cos(gam-q1)-m*g*w*sin(gam-q1)+2*g*m*sin(gam)*r*q1-m*g*c*cos(gam-q1+q2)-m*g*w*sin(gam-q1+q2)+g*M*cos(gam)*r+M*g*l*cos(gam-q1)+g*M*sin(gam)*r*q1; 
    vxh = (-l*cos(q1)-r)*u1; 
    vyh = -l*sin(q1)*u1; 

    zplus = [q1 u1 q2 u2 TE xh vxh yh vyh];                     
end