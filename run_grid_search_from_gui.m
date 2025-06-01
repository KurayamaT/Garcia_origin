function run_grid_search_from_gui()
% GUIで設定した値を使ってグリッドサーチを実行

    % 設定値の確認
    if ~evalin('base', 'exist(''u1_range'', ''var'')') || ...
       ~evalin('base', 'exist(''q2_range'', ''var'')') || ...
       ~evalin('base', 'exist(''u2_range'', ''var'')')
        error('先にgrid_search_gui_simple()を実行して設定してください。');
    end
    
    % 範囲を取得
    u1_range = evalin('base', 'u1_range');
    q2_range = evalin('base', 'q2_range');
    u2_range = evalin('base', 'u2_range');
    
    q1_fixed = 0.0;  % 直立状態に固定
    
    % 確認
    fprintf('\n=== グリッドサーチ開始 ===\n');
    fprintf('q1 = %.1f (固定)\n', q1_fixed);
    fprintf('u1: %d 点\n', length(u1_range));
    fprintf('q2: %d 点\n', length(q2_range));
    fprintf('u2: %d 点\n', length(u2_range));
    fprintf('総探索数: %d\n', length(u1_range)*length(q2_range)*length(u2_range));
    
    % 結果保存用
    results = [];
    successful = [];
    counter = 0;
    total = length(u1_range)*length(q2_range)*length(u2_range);
    
    tic;
    fprintf('\n探索中...\n');
    
    % プログレスバー
    fprintf('[');
    progress_marks = 0;
    
    for i = 1:length(u1_range)
        for j = 1:length(q2_range)
            for k = 1:length(u2_range)
                counter = counter + 1;
                
                % プログレス表示
                new_marks = floor(counter/total * 50);
                if new_marks > progress_marks
                    fprintf(repmat('=', 1, new_marks-progress_marks));
                    progress_marks = new_marks;
                end
                
                % 初期条件
                q1 = q1_fixed;
                u1 = u1_range(i);
                q2 = q2_range(j);
                u2 = u2_range(k);
                
                % ここで実際の評価を行う
                % （元のプログラムの fixedpt, onestep などを呼び出す）
                
                % 仮の結果（実際の評価に置き換える）
                result = struct();
                result.q1 = q1;
                result.u1 = u1;
                result.q2 = q2;
                result.u2 = u2;
                result.success = false;  % 実際の判定に置き換える
                
                results = [results; result];
                
                if result.success
                    successful = [successful; result];
                end
            end
        end
    end
    
    fprintf(']\n');
    elapsed = toc;
    
    % 結果表示
    fprintf('\n=== 結果 ===\n');
    fprintf('探索時間: %.1f 秒\n', elapsed);
    fprintf('成功数: %d / %d\n', length(successful), total);
    
    % ワークスペースに保存
    assignin('base', 'grid_results', results);
    assignin('base', 'successful_results', successful);
    
    fprintf('\n結果は以下の変数に保存されました:\n');
    fprintf('- grid_results: 全結果\n');
    fprintf('- successful_results: 成功した結果\n');
end