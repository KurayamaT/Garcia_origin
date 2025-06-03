function grid_search_gui_simple()
% 改良版GUI - 設定から計算実行まで一括で行う

    % Figure作成
    fig = figure('Position', [200 200 500 550], ...
                 'Name', '直立歩行グリッドサーチ', ...
                 'NumberTitle', 'off', ...
                 'MenuBar', 'none');
    
    % タイトル
    uicontrol('Style', 'text', ...
              'String', '直立歩行初期条件グリッドサーチ', ...
              'Position', [10 510 480 30], ...
              'FontSize', 14, ...
              'FontWeight', 'bold');
    
    uicontrol('Style', 'text', ...
              'String', '参考値（オリジナル）: q1=0.2, u1=-0.2, q2=0.4, u2=-0.3', ...
              'Position', [10 480 480 20], ...
              'FontSize', 10);
    
    % q1設定
    uicontrol('Style', 'text', ...
              'String', 'q1（スタンス脚角度）:', ...
              'Position', [20 440 200 20], ...
              'HorizontalAlignment', 'left', ...
              'FontWeight', 'bold', ...
              'ForegroundColor', [0.8 0 0]);
    
    uicontrol('Style', 'text', 'String', '最小値:', ...
              'Position', [30 410 60 20], ...
              'HorizontalAlignment', 'right');
    q1_min = uicontrol('Style', 'edit', ...
                       'String', '0.15', ...
                       'Position', [100 410 80 25], ...
                       'BackgroundColor', [1 0.95 0.95]);
    
    uicontrol('Style', 'text', 'String', '最大値:', ...
              'Position', [190 410 60 20], ...
              'HorizontalAlignment', 'right');
    q1_max = uicontrol('Style', 'edit', ...
                       'String', '0.25', ...
                       'Position', [260 410 80 25], ...
                       'BackgroundColor', [1 0.95 0.95]);
    
    uicontrol('Style', 'text', 'String', '刻み:', ...
              'Position', [350 410 50 20], ...
              'HorizontalAlignment', 'right');
    q1_step = uicontrol('Style', 'edit', ...
                        'String', '0.05', ...
                        'Position', [410 410 60 25], ...
                        'BackgroundColor', [1 0.95 0.95]);
    
    % u1設定
    uicontrol('Style', 'text', ...
              'String', 'u1（スタンス脚角速度）:', ...
              'Position', [20 360 200 20], ...
              'HorizontalAlignment', 'left', ...
              'FontWeight', 'bold');
    
    uicontrol('Style', 'text', 'String', '最小値:', ...
              'Position', [30 330 60 20], ...
              'HorizontalAlignment', 'right');
    u1_min = uicontrol('Style', 'edit', ...
                       'String', '-0.25', ...
                       'Position', [100 330 80 25]);
    
    uicontrol('Style', 'text', 'String', '最大値:', ...
              'Position', [190 330 60 20], ...
              'HorizontalAlignment', 'right');
    u1_max = uicontrol('Style', 'edit', ...
                       'String', '-0.15', ...
                       'Position', [260 330 80 25]);
    
    uicontrol('Style', 'text', 'String', '刻み:', ...
              'Position', [350 330 50 20], ...
              'HorizontalAlignment', 'right');
    u1_step = uicontrol('Style', 'edit', ...
                        'String', '0.05', ...
                        'Position', [410 330 60 25]);
    
    % q2設定
    uicontrol('Style', 'text', ...
              'String', 'q2（スイング脚角度）:', ...
              'Position', [20 280 200 20], ...
              'HorizontalAlignment', 'left', ...
              'FontWeight', 'bold', ...
              'ForegroundColor', [0 0 0.8]);
    
    uicontrol('Style', 'text', 'String', '最小値:', ...
              'Position', [30 250 60 20], ...
              'HorizontalAlignment', 'right');
    q2_min = uicontrol('Style', 'edit', ...
                       'String', '0.35', ...
                       'Position', [100 250 80 25], ...
                       'BackgroundColor', [0.95 0.95 1]);
    
    uicontrol('Style', 'text', 'String', '最大値:', ...
              'Position', [190 250 60 20], ...
              'HorizontalAlignment', 'right');
    q2_max = uicontrol('Style', 'edit', ...
                       'String', '0.45', ...
                       'Position', [260 250 80 25], ...
                       'BackgroundColor', [0.95 0.95 1]);
    
    uicontrol('Style', 'text', 'String', '刻み:', ...
              'Position', [350 250 50 20], ...
              'HorizontalAlignment', 'right');
    q2_step = uicontrol('Style', 'edit', ...
                        'String', '0.05', ...
                        'Position', [410 250 60 25], ...
                        'BackgroundColor', [0.95 0.95 1]);
    
    % u2設定
    uicontrol('Style', 'text', ...
              'String', 'u2（スイング脚角速度）:', ...
              'Position', [20 200 200 20], ...
              'HorizontalAlignment', 'left', ...
              'FontWeight', 'bold', ...
              'ForegroundColor', [0 0 0.8]);
    
    uicontrol('Style', 'text', 'String', '最小値:', ...
              'Position', [30 170 60 20], ...
              'HorizontalAlignment', 'right');
    u2_min = uicontrol('Style', 'edit', ...
                       'String', '-0.35', ...
                       'Position', [100 170 80 25], ...
                       'BackgroundColor', [0.95 0.95 1]);
    
    uicontrol('Style', 'text', 'String', '最大値:', ...
              'Position', [190 170 60 20], ...
              'HorizontalAlignment', 'right');
    u2_max = uicontrol('Style', 'edit', ...
                       'String', '-0.25', ...
                       'Position', [260 170 80 25], ...
                       'BackgroundColor', [0.95 0.95 1]);
    
    uicontrol('Style', 'text', 'String', '刻み:', ...
              'Position', [350 170 50 20], ...
              'HorizontalAlignment', 'right');
    u2_step = uicontrol('Style', 'edit', ...
                        'String', '0.05', ...
                        'Position', [410 170 60 25], ...
                        'BackgroundColor', [0.95 0.95 1]);
    
    % 進捗表示エリア
    progress_text = uicontrol('Style', 'text', ...
                              'String', '設定を確認して「計算開始」ボタンを押してください', ...
                              'Position', [10 120 480 30], ...
                              'FontSize', 10, ...
                              'HorizontalAlignment', 'center', ...
                              'BackgroundColor', [0.95 0.95 0.95]);
    
    % ボタン
    uicontrol('Style', 'pushbutton', ...
              'String', '設定確認', ...
              'Position', [50 70 100 35], ...
              'FontSize', 11, ...
              'Callback', @check_settings);
    
    uicontrol('Style', 'pushbutton', ...
              'String', '🚀 計算開始', ...
              'Position', [170 70 120 35], ...
              'FontSize', 11, ...
              'FontWeight', 'bold', ...
              'BackgroundColor', [0.8 1 0.8], ...
              'Callback', @start_calculation);
    
    uicontrol('Style', 'pushbutton', ...
              'String', '閉じる', ...
              'Position', [310 70 100 35], ...
              'FontSize', 11, ...
              'Callback', @(~,~) close(fig));
    
    % コールバック関数
    function check_settings(~, ~)
        [vals, ranges] = get_settings();
        total = length(ranges.q1) * length(ranges.u1) * length(ranges.q2) * length(ranges.u2);
        
        % 表示更新
        info_text = sprintf('総探索数: %d | q1:%d点, u1:%d点, q2:%d点, u2:%d点', ...
                          total, length(ranges.q1), length(ranges.u1), ...
                          length(ranges.q2), length(ranges.u2));
        set(progress_text, 'String', info_text);
        
        fprintf('\n=== 設定確認 ===\n');
        fprintf('q1: %.2f から %.2f まで %.3f 刻み (%d 点)\n', ...
                vals.q1_min, vals.q1_max, vals.q1_step, length(ranges.q1));
        fprintf('u1: %.2f から %.2f まで %.3f 刻み (%d 点)\n', ...
                vals.u1_min, vals.u1_max, vals.u1_step, length(ranges.u1));
        fprintf('q2: %.2f から %.2f まで %.3f 刻み (%d 点)\n', ...
                vals.q2_min, vals.q2_max, vals.q2_step, length(ranges.q2));
        fprintf('u2: %.2f から %.2f まで %.3f 刻み (%d 点)\n', ...
                vals.u2_min, vals.u2_max, vals.u2_step, length(ranges.u2));
        fprintf('総探索数: %d\n', total);
    end
    
    function start_calculation(~, ~)
        % 設定取得
        [vals, ranges] = get_settings();
        total = length(ranges.q1) * length(ranges.u1) * length(ranges.q2) * length(ranges.u2);
        
        % 警告チェック
        if total > 1000
            answer = questdlg(sprintf('探索数が %d と多いです。時間がかかる可能性があります。続行しますか？', total), ...
                            '確認', 'はい', 'いいえ', 'いいえ');
            if ~strcmp(answer, 'はい')
                return;
            end
        end
        
        % 進捗表示更新
        set(progress_text, 'String', '計算開始中... コマンドウィンドウを確認してください');
        
        % ワークスペースに保存
        assignin('base', 'q1_range', ranges.q1);
        assignin('base', 'u1_range', ranges.u1);
        assignin('base', 'q2_range', ranges.q2);
        assignin('base', 'u2_range', ranges.u2);
        
        % GUIを無効化
        set(findobj(fig, 'Type', 'uicontrol'), 'Enable', 'off');
        
        fprintf('\n🚀 グリッドサーチを開始します...\n');
        
        % グリッドサーチ実行
        try
            run_grid_search_calculation(ranges);
            set(progress_text, 'String', '✅ 計算完了！結果をコマンドウィンドウで確認してください');
        catch ME
            set(progress_text, 'String', '❌ 計算エラーが発生しました');
            fprintf('エラー: %s\n', ME.message);
        end
        
        % GUIを再有効化
        set(findobj(fig, 'Type', 'uicontrol'), 'Enable', 'on');
    end
    
    function [vals, ranges] = get_settings()
        % 値を取得
        vals = struct();
        vals.q1_min = str2double(get(q1_min, 'String'));
        vals.q1_max = str2double(get(q1_max, 'String'));
        vals.q1_step = str2double(get(q1_step, 'String'));
        vals.u1_min = str2double(get(u1_min, 'String'));
        vals.u1_max = str2double(get(u1_max, 'String'));
        vals.u1_step = str2double(get(u1_step, 'String'));
        vals.q2_min = str2double(get(q2_min, 'String'));
        vals.q2_max = str2double(get(q2_max, 'String'));
        vals.q2_step = str2double(get(q2_step, 'String'));
        vals.u2_min = str2double(get(u2_min, 'String'));
        vals.u2_max = str2double(get(u2_max, 'String'));
        vals.u2_step = str2double(get(u2_step, 'String'));
        
        % 範囲作成
        ranges = struct();
        ranges.q1 = vals.q1_min:vals.q1_step:vals.q1_max;
        if isempty(ranges.q1), ranges.q1 = vals.q1_min; end
        ranges.u1 = vals.u1_min:vals.u1_step:vals.u1_max;
        ranges.q2 = vals.q2_min:vals.q2_step:vals.q2_max;
        ranges.u2 = vals.u2_min:vals.u2_step:vals.u2_max;
    end
end

%% グリッドサーチ実行関数
function run_grid_search_calculation(ranges)
    % Walker設定
    walker.M = 1000; walker.m = 1.0; walker.I = 0.00; walker.l = 1.0; walker.w = 0.0;
    walker.c = 1.0;  walker.r = 0.3; walker.g = 1.0; walker.gam = 0.009;
    
    % 並列処理の確認
    use_parallel = input('並列処理を使用しますか？ (y/n) [推奨: y]: ', 's');
    if isempty(use_parallel) || strcmpi(use_parallel, 'y')
        pool = gcp('nocreate');
        if isempty(pool)
            fprintf('並列プールを起動中...\n');
            parpool;
        else
            fprintf('既存の並列プール（ワーカー数: %d）を使用\n', pool.NumWorkers);
        end
        use_parallel = true;
    else
        use_parallel = false;
    end
    
    % 全組み合わせ作成
    if length(ranges.q1) == 1
        [U1, Q2, U2] = meshgrid(ranges.u1, ranges.q2, ranges.u2);
        Q1 = ranges.q1 * ones(size(U1));
    else
        [Q1, U1, Q2, U2] = ndgrid(ranges.q1, ranges.u1, ranges.q2, ranges.u2);
    end
    all_z0 = [Q1(:), U1(:), Q2(:), U2(:)];
    total = size(all_z0, 1);
    
    fprintf('\n総探索数: %d\n', total);
    fprintf('探索中... (成功例は即座に表示されます)\n\n');
    
    % 探索実行
    results = cell(total, 1);
    success_count = 0;
    tic;
    
    if use_parallel
        % 並列処理
        parfor idx = 1:total
            z0 = all_z0(idx, :);
            result = evaluate_single(z0, walker);
            
            if result.success
                fprintf('【成功】 q1=%6.3f, u1=%6.3f, q2=%6.3f, u2=%6.3f | λ_max=%6.4f | θ_max=%5.1f°\n', ...
                        z0(1), z0(2), z0(3), z0(4), result.max_eig, result.max_angle*180/pi);
            end
            
            results{idx} = result;
        end
    else
        % 逐次処理
        for idx = 1:total
            if mod(idx, 100) == 0
                fprintf('進捗: %d/%d\n', idx, total);
            end
            
            z0 = all_z0(idx, :);
            result = evaluate_single(z0, walker);
            
            if result.success
                success_count = success_count + 1;
                fprintf('【成功 #%d】 q1=%6.3f, u1=%6.3f, q2=%6.3f, u2=%6.3f | λ_max=%6.4f | θ_max=%5.1f°\n', ...
                        success_count, z0(1), z0(2), z0(3), z0(4), result.max_eig, result.max_angle*180/pi);
            end
            
            results{idx} = result;
        end
    end
    
    elapsed = toc;
    
    % 結果集計
    results_array = [results{:}];
    successful = results_array([results_array.success]);
    
    fprintf('\n🎉 完了！ 時間: %.1f秒\n', elapsed);
    fprintf('成功: %d/%d (%.1f%%)\n', length(successful), total, 100*length(successful)/total);
    
    if ~isempty(successful)
        % 最良の結果を表示
        [~, best_idx] = min([successful.max_eig]);
        best = successful(best_idx);
        fprintf('\n【最良の初期条件】\n');
        fprintf('q1=%.3f, u1=%.3f, q2=%.3f, u2=%.3f\n', ...
                best.q1, best.u1, best.q2, best.u2);
        fprintf('最大固有値: %.4f\n', best.max_eig);
        
        % 成功例のサマリー
        fprintf('\n【成功例のサマリー】\n');
        fprintf('No. | q1      | u1      | q2      | u2      | λ_max    | θ_max\n');
        fprintf('----|---------|---------|---------|---------|----------|--------\n');
        for i = 1:min(10, length(successful))
            s = successful(i);
            fprintf('%3d | %7.3f | %7.3f | %7.3f | %7.3f | %8.4f | %6.1f°\n', ...
                    i, s.q1, s.u1, s.q2, s.u2, s.max_eig, s.max_angle*180/pi);
        end
        
        if length(successful) > 10
            fprintf('... 他 %d 個の成功例\n', length(successful)-10);
        end
        
        % ワークスペースに保存
        assignin('base', 'best_result', best);
    else
        fprintf('\n成功例が見つかりませんでした。範囲を調整して再試行してください。\n');
    end
    
    % 保存
    assignin('base', 'search_results', results_array);
    assignin('base', 'successful_results', successful);
end

%% 個別評価関数
function result = evaluate_single(z0, walker)
    result.q1 = z0(1); 
    result.u1 = z0(2);
    result.q2 = z0(3); 
    result.u2 = z0(4);
    result.success = false;
    result.max_eig = inf;
    result.max_angle = inf;
    
    try
        % 固定点探索
        options = optimset('TolFun',1e-12,'TolX',1e-12,'Display','off');
        [zstar, ~, exitflag] = fsolve(@fixedpt, z0, options, walker);
        
        if exitflag == 1
            result.zstar = zstar;
            
            % 安定性チェック
            J = partialder(@onestep, zstar, walker);
            eigenvalues = eig(J);
            result.max_eig = max(abs(eigenvalues));
            
            if result.max_eig < 1
                % 実際の歩行テスト
                [z, ~] = onestep(z0, walker, 3);  % 3ステップで高速化
                
                % 転倒チェック（修正版）
                result.max_angle = max(abs(z(:,1)));
                result.min_height = 0;  % デフォルト値
                result.success = (result.max_angle < pi/2);  % 角度のみで判定
            end
        end
    catch
        % エラーは無視
    end
end