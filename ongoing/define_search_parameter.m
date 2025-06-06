function simple_grid_gui()
% シンプルなグリッドサーチGUI（計算時間推定機能付き）

    % Figure作成
    fig = figure('Position', [200 200 500 450], ...
                 'Name', 'パッシブウォーカー グリッドサーチ', ...
                 'NumberTitle', 'off', ...
                 'MenuBar', 'none');
    
    % タイトル
    uicontrol('Style', 'text', ...
              'String', 'パッシブウォーカー初期条件グリッドサーチ', ...
              'Position', [10 410 480 30], ...
              'FontSize', 14, ...
              'FontWeight', 'bold');
    
    % q1設定
    uicontrol('Style', 'text', ...
              'String', 'q1（スタンス脚角度）:', ...
              'Position', [20 370 150 20], ...
              'HorizontalAlignment', 'left', ...
              'FontWeight', 'bold');
    
    uicontrol('Style', 'text', 'String', '最小値:', ...
              'Position', [30 340 60 20], 'HorizontalAlignment', 'right');
    q1_min = uicontrol('Style', 'edit', 'String', '0.0', ...
                       'Position', [100 340 80 25]);
    
    uicontrol('Style', 'text', 'String', '最大値:', ...
              'Position', [190 340 60 20], 'HorizontalAlignment', 'right');
    q1_max = uicontrol('Style', 'edit', 'String', '0.1', ...
                       'Position', [260 340 80 25]);
    
    uicontrol('Style', 'text', 'String', '刻み:', ...
              'Position', [350 340 50 20], 'HorizontalAlignment', 'right');
    q1_step = uicontrol('Style', 'edit', 'String', '0.01', ...
                        'Position', [410 340 60 25]);
    
    % u1設定
    uicontrol('Style', 'text', ...
              'String', 'u1（スタンス脚角速度）:', ...
              'Position', [20 290 150 20], ...
              'HorizontalAlignment', 'left', ...
              'FontWeight', 'bold');
    
    uicontrol('Style', 'text', 'String', '最小値:', ...
              'Position', [30 260 60 20], 'HorizontalAlignment', 'right');
    u1_min = uicontrol('Style', 'edit', 'String', '-0.25', ...
                       'Position', [100 260 80 25]);
    
    uicontrol('Style', 'text', 'String', '最大値:', ...
              'Position', [190 260 60 20], 'HorizontalAlignment', 'right');
    u1_max = uicontrol('Style', 'edit', 'String', '-0.15', ...
                       'Position', [260 260 80 25]);
    
    uicontrol('Style', 'text', 'String', '刻み:', ...
              'Position', [350 260 50 20], 'HorizontalAlignment', 'right');
    u1_step = uicontrol('Style', 'edit', 'String', '0.01', ...
                        'Position', [410 260 60 25]);
    
    % q2設定
    uicontrol('Style', 'text', ...
              'String', 'q2（スイング脚角度）:', ...
              'Position', [20 210 150 20], ...
              'HorizontalAlignment', 'left', ...
              'FontWeight', 'bold');
    
    uicontrol('Style', 'text', 'String', '最小値:', ...
              'Position', [30 180 60 20], 'HorizontalAlignment', 'right');
    q2_min = uicontrol('Style', 'edit', 'String', '0.0', ...
                       'Position', [100 180 80 25]);
    
    uicontrol('Style', 'text', 'String', '最大値:', ...
              'Position', [190 180 60 20], 'HorizontalAlignment', 'right');
    q2_max = uicontrol('Style', 'edit', 'String', '0.45', ...
                       'Position', [260 180 80 25]);
    
    uicontrol('Style', 'text', 'String', '刻み:', ...
              'Position', [350 180 50 20], 'HorizontalAlignment', 'right');
    q2_step = uicontrol('Style', 'edit', 'String', '0.01', ...
                        'Position', [410 180 60 25]);
    
    % u2設定
    uicontrol('Style', 'text', ...
              'String', 'u2（スイング脚角速度）:', ...
              'Position', [20 130 150 20], ...
              'HorizontalAlignment', 'left', ...
              'FontWeight', 'bold');
    
    uicontrol('Style', 'text', 'String', '最小値:', ...
              'Position', [30 100 60 20], 'HorizontalAlignment', 'right');
    u2_min = uicontrol('Style', 'edit', 'String', '-0.35', ...
                       'Position', [100 100 80 25]);
    
    uicontrol('Style', 'text', 'String', '最大値:', ...
              'Position', [190 100 60 20], 'HorizontalAlignment', 'right');
    u2_max = uicontrol('Style', 'edit', 'String', '-0.25', ...
                       'Position', [260 100 80 25]);
    
    uicontrol('Style', 'text', 'String', '刻み:', ...
              'Position', [350 100 50 20], 'HorizontalAlignment', 'right');
    u2_step = uicontrol('Style', 'edit', 'String', '0.01', ...
                        'Position', [410 100 60 25]);
    
    % ボタン
    uicontrol('Style', 'pushbutton', ...
              'String', '設定確認', ...
              'Position', [20 40 80 35], ...
              'FontSize', 10, ...
              'Callback', @check_settings);
    
    uicontrol('Style', 'pushbutton', ...
              'String', '📊 グリッド解析', ...
              'Position', [110 40 100 35], ...
              'FontSize', 10, ...
              'FontWeight', 'bold', ...
              'BackgroundColor', [0.8 0.8 1], ...
              'Callback', @start_grid_analysis);
    
    uicontrol('Style', 'pushbutton', ...
              'String', '🎯 収束探索', ...
              'Position', [220 40 100 35], ...
              'FontSize', 10, ...
              'FontWeight', 'bold', ...
              'BackgroundColor', [0.8 1 0.8], ...
              'Callback', @start_convergence_search);
    
    uicontrol('Style', 'pushbutton', ...
              'String', '閉じる', ...
              'Position', [330 40 80 35], ...
              'FontSize', 10, ...
              'Callback', @(~,~) close(fig));
    
    % コールバック関数
    function check_settings(~, ~)
        params = get_parameters();
        
        fprintf('\n=== 設定確認 ===\n');
        fprintf('q1: %.2f から %.2f まで %.3f 刻み (%d 点)\n', ...
                params.q1_min, params.q1_max, params.q1_step, length(params.q1_range));
        fprintf('u1: %.2f から %.2f まで %.3f 刻み (%d 点)\n', ...
                params.u1_min, params.u1_max, params.u1_step, length(params.u1_range));
        fprintf('q2: %.2f から %.2f まで %.3f 刻み (%d 点)\n', ...
                params.q2_min, params.q2_max, params.q2_step, length(params.q2_range));
        fprintf('u2: %.2f から %.2f まで %.3f 刻み (%d 点)\n', ...
                params.u2_min, params.u2_max, params.u2_step, length(params.u2_range));
        fprintf('総探索数: %d\n', params.total_combinations);
        
        % ⏱️ 計算時間推定を追加
        fprintf('\n⏱️ === 計算時間推定 ===\n');
        
        % 実績データに基づく推定（220個で20秒）
        reference_rate = 20 / 220;  % 秒/個
        
        % グリッド解析の時間推定（固定点計算とヤコビアン計算）
        grid_rate = reference_rate * 0.3;  % グリッド解析は収束テストより高速
        grid_time_sequential = grid_rate * params.total_combinations;
        
        % 収束探索の時間推定
        convergence_rate = reference_rate;  % 実績ベース
        convergence_time_sequential = convergence_rate * params.total_combinations;
        
        % 並列処理での推定（利用可能なコア数を考慮）
        pool = gcp('nocreate');
        if isempty(pool)
            max_workers = feature('numcores');
            typical_workers = min(max_workers - 1, 8);  % 通常は最大8ワーカー
        else
            typical_workers = pool.NumWorkers;
        end
        
        % 並列処理の効率を考慮（オーバーヘッド1.2倍）
        parallel_efficiency = 1.2;
        grid_time_parallel = grid_time_sequential / typical_workers * parallel_efficiency;
        convergence_time_parallel = convergence_time_sequential / typical_workers * parallel_efficiency;
        
        % グリッド解析の推定表示
        fprintf('\n📊 グリッド解析の推定時間:\n');
        fprintf('  - 逐次処理: %.1f秒 (%.1f分)\n', grid_time_sequential, grid_time_sequential/60);
        fprintf('  - 並列処理（%dワーカー想定）: %.1f秒 (%.1f分)\n', ...
                typical_workers, grid_time_parallel, grid_time_parallel/60);
        
        % 収束探索の推定表示
        fprintf('\n🎯 収束探索の推定時間:\n');
        fprintf('  - 逐次処理: %.1f秒 (%.1f分)\n', convergence_time_sequential, convergence_time_sequential/60);
        fprintf('  - 並列処理（%dワーカー想定）: %.1f秒 (%.1f分)\n', ...
                typical_workers, convergence_time_parallel, convergence_time_parallel/60);
        
        % 警告表示
        if convergence_time_parallel > 300  % 5分以上
            fprintf('\n⚠️ 注意: 計算に時間がかかる可能性があります。\n');
            fprintf('   探索範囲を狭めるか、刻み幅を大きくすることを検討してください。\n');
        elseif convergence_time_parallel > 60  % 1分以上
            fprintf('\n💡 ヒント: 並列処理を使用することで計算時間を短縮できます。\n');
        end
        
        % メモリ使用量の推定
        memory_per_condition = 8 * 9 * 2;  % double型 * 変数数 * 係数
        total_memory_mb = params.total_combinations * memory_per_condition / 1024 / 1024;
        fprintf('\n💾 推定メモリ使用量: %.1f MB\n', total_memory_mb);
        
        % ワークスペースに設定を保存
        assignin('base', 'grid_params', params);
        fprintf('\n✅ 設定がワークスペースに保存されました（変数名: grid_params）\n');
    end
    
    function start_grid_analysis(~, ~)
        params = get_parameters();
        
        % 時間推定を含む確認ダイアログ
        if params.total_combinations > 100
            % 簡易推定
            grid_time = params.total_combinations * 0.03;  % グリッド解析の推定
            
            message = sprintf(['探索数: %d\n' ...
                             '推定時間: %.1f秒 (%.1f分)\n\n' ...
                             '続行しますか？'], ...
                             params.total_combinations, grid_time, grid_time/60);
            
            answer = questdlg(message, '確認', 'はい', 'いいえ', 'いいえ');
            if ~strcmp(answer, 'はい')
                return;
            end
        end
        
        fprintf('\n📊 グリッド解析を開始します...\n');
        fprintf('総探索数: %d\n', params.total_combinations);
        
        % grid_runner.m を呼び出し
        try
            grid_runner(params);
        catch ME
            fprintf('エラーが発生しました: %s\n', ME.message);
            fprintf('エラー詳細:\n%s\n', getReport(ME));
        end
    end
    
    function start_convergence_search(~, ~)
        params = get_parameters();
        
        % 時間推定を含む確認ダイアログ
        if params.total_combinations > 100
            % 実績ベースの推定
            convergence_time = params.total_combinations * (20/220);  % 220個で20秒の実績
            
            message = sprintf(['探索数: %d\n' ...
                             '推定時間: %.1f秒 (%.1f分)\n\n' ...
                             '続行しますか？'], ...
                             params.total_combinations, convergence_time, convergence_time/60);
            
            answer = questdlg(message, '確認', 'はい', 'いいえ', 'いいえ');
            if ~strcmp(answer, 'はい')
                return;
            end
        end
        
        fprintf('\n🎯 収束探索を開始します...\n');
        fprintf('総探索数: %d\n', params.total_combinations);
        
        % convergence_basin_search.m を呼び出し
        try
            convergence_basin_search(params);
        catch ME
            fprintf('エラーが発生しました: %s\n', ME.message);
            fprintf('エラー詳細:\n%s\n', getReport(ME));
        end
    end
    
    function params = get_parameters()
        % パラメータを取得
        params.q1_min = str2double(get(q1_min, 'String'));
        params.q1_max = str2double(get(q1_max, 'String'));
        params.q1_step = str2double(get(q1_step, 'String'));
        params.u1_min = str2double(get(u1_min, 'String'));
        params.u1_max = str2double(get(u1_max, 'String'));
        params.u1_step = str2double(get(u1_step, 'String'));
        params.q2_min = str2double(get(q2_min, 'String'));
        params.q2_max = str2double(get(q2_max, 'String'));
        params.q2_step = str2double(get(q2_step, 'String'));
        params.u2_min = str2double(get(u2_min, 'String'));
        params.u2_max = str2double(get(u2_max, 'String'));
        params.u2_step = str2double(get(u2_step, 'String'));
        
        % 範囲作成
        params.q1_range = params.q1_min:params.q1_step:params.q1_max;
        if isempty(params.q1_range), params.q1_range = params.q1_min; end
        
        params.u1_range = params.u1_min:params.u1_step:params.u1_max;
        if isempty(params.u1_range), params.u1_range = params.u1_min; end
        
        params.q2_range = params.q2_min:params.q2_step:params.q2_max;
        if isempty(params.q2_range), params.q2_range = params.q2_min; end
        
        params.u2_range = params.u2_min:params.u2_step:params.u2_max;
        if isempty(params.u2_range), params.u2_range = params.u2_min; end
        
        % 総組み合わせ数
        params.total_combinations = length(params.q1_range) * length(params.u1_range) * ...
                                   length(params.q2_range) * length(params.u2_range);
    end
end