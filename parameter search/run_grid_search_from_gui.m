function run_grid_search_from_gui()
% 元のpassive walkerコードの関数を呼び出すグリッドサーチ
% parforによる並列処理で高速化、成功例をリアルタイム表示

    %% 範囲の取得
    if ~evalin('base', 'exist(''u1_range'', ''var'')')
        error('先にgrid_search_gui_simple()を実行してください。');
    end
    
    if evalin('base', 'exist(''q1_range'', ''var'')')
        q1_range = evalin('base', 'q1_range');
    else
        q1_range = 0;
    end
    u1_range = evalin('base', 'u1_range');
    q2_range = evalin('base', 'q2_range');
    u2_range = evalin('base', 'u2_range');
    
    %% Walker設定
    walker.M = 1000; walker.m = 1.0; walker.I = 0.00; walker.l = 1.0; walker.w = 0.0;
    walker.c = 1.0;  walker.r = 0.3; walker.g = 1.0; walker.gam = 0.009;
    
    %% 並列処理の確認
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
    
    %% 全組み合わせ作成
    if length(q1_range) == 1
        [U1, Q2, U2] = meshgrid(u1_range, q2_range, u2_range);
        Q1 = q1_range * ones(size(U1));
    else
        [Q1, U1, Q2, U2] = ndgrid(q1_range, u1_range, q2_range, u2_range);
    end
    all_z0 = [Q1(:), U1(:), Q2(:), U2(:)];
    total = size(all_z0, 1);
    
    fprintf('\n総探索数: %d\n', total);
    fprintf('探索中... (成功例は即座に表示されます)\n\n');
    
    %% 探索実行
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
    
    %% 結果集計
    results_array = [results{:}];
    successful = results_array([results_array.success]);
    
    fprintf('\n完了！ 時間: %.1f秒\n', elapsed);
    fprintf('成功: %d/%d (%.1f%%)\n', length(successful), total, 100*length(successful)/total);
    
    if ~isempty(successful)
        % 最良の結果を表示
        [~, best_idx] = min([successful.max_eig]);
        best = successful(best_idx);
        fprintf('\n【最良の初期条件】\n');
        fprintf('q1=%.3f, u1=%.3f, q2=%.3f, u2=%.3f\n', ...
                best.q1, best.u1, best.q2, best.u2);
        fprintf('最大固有値: %.4f\n', best.max_eig);
        fprintf('固定点: [%.4f, %.4f, %.4f, %.4f]\n', best.zstar);
        
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
    end
    
    % 保存
    assignin('base', 'search_results', results_array);
    assignin('base', 'successful_results', successful);
    if ~isempty(successful)
        assignin('base', 'best_result', best);
    end
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
                [z, ~] = onestep(z0, walker, 5);
                
                % 転倒チェック
                result.max_angle = max(abs(z(:,1)));
                result.min_height = 0;  % デフォルト値
                result.success = (result.max_angle < pi/2);  % 角度のみで判定
            end
        end
    catch
        % エラーは無視
    end
end