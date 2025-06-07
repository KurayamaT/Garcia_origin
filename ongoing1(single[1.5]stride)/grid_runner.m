function grid_runner(params)
% グリッドサーチ実行関数
% params: GUIから渡されるパラメータ構造体

    fprintf('\n=== パッシブウォーカー グリッドサーチ実行 ===\n');
    
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
    
    % 並列処理の確認と設定
    use_parallel = true;
    pool = gcp('nocreate');
    if isempty(pool)
        if total > 50  % 50個以上の場合のみ並列処理を提案
            answer = input('並列処理を使用しますか？ (y/n) [推奨: y]: ', 's');
            if isempty(answer) || strcmpi(answer, 'y')
                fprintf('並列プールを起動中...\n');
                parpool;
                use_parallel = true;
            else
                use_parallel = false;
            end
        else
            use_parallel = false;
        end
    else
        fprintf('既存の並列プール（ワーカー数: %d）を使用\n', pool.NumWorkers);
        use_parallel = true;
    end
    
    % 結果保存用
    results_cell = cell(total, 1);
    
    % 進捗表示
    if use_parallel
        fprintf('\n🚀 並列処理で探索中...\n');
    else
        fprintf('\n🚀 逐次処理で探索中...\n');
    end
    tic;
    
    % 各条件をテスト（並列処理 or 逐次処理）
    if use_parallel
        % 並列処理
        parfor idx = 1:total
            z0 = all_conditions(idx, :);
            [is_success, max_eigenvalue, fixed_point] = evaluate_condition(z0);
            
            if is_success
                results_cell{idx} = [z0, max_eigenvalue, fixed_point];
            end
        end
        
        % 並列処理後の結果集計
        fprintf('並列処理完了。結果を集計中...\n');
        results = [];
        success_count = 0;
        for idx = 1:total
            if ~isempty(results_cell{idx})
                success_count = success_count + 1;
                results = [results; results_cell{idx}];
                z0 = results_cell{idx}(1:4);
                max_eigenvalue = results_cell{idx}(5);
                fprintf('【成功 #%d】 q1=%.3f, u1=%.3f, q2=%.3f, u2=%.3f | 最大固有値=%.4f\n', ...
                        success_count, z0, max_eigenvalue);
            end
        end
    else
        % 逐次処理
        results = [];
        success_count = 0;
        for idx = 1:total
            % 進捗表示
            if mod(idx, max(1, floor(total/20))) == 0
                fprintf('進捗: %d/%d (%.1f%%)\n', idx, total, 100*idx/total);
            end
            
            z0 = all_conditions(idx, :);
            [is_success, max_eigenvalue, fixed_point] = evaluate_condition(z0);
            
            if is_success
                success_count = success_count + 1;
                results = [results; z0, max_eigenvalue, fixed_point];
                
                fprintf('【成功 #%d】 q1=%.3f, u1=%.3f, q2=%.3f, u2=%.3f | 最大固有値=%.4f\n', ...
                        success_count, z0, max_eigenvalue);
            end
        end
    end
    
    elapsed = toc;
    
    % 結果まとめ
    fprintf('\n%s\n', repmat('=', 1, 60));
    fprintf('🎉 グリッドサーチ完了！\n');
    fprintf('実行時間: %.1f秒\n', elapsed);
    if use_parallel && ~isempty(pool)
        fprintf('並列処理（%dワーカー）で高速化されました\n', pool.NumWorkers);
    end
    fprintf('成功条件: %d/%d (%.1f%%)\n', success_count, total, 100*success_count/total);
    fprintf('%s\n', repmat('=', 1, 60));
    
    if success_count > 0
        % 安定性でソート（最大固有値の小さい順）
        [~, sort_idx] = sort(results(:, 5));
        sorted_results = results(sort_idx, :);
        
        % 最良の条件を表示
        best = sorted_results(1, :);
        fprintf('\n🏆 【最も安定な条件】\n');
        fprintf('q1=%.4f, u1=%.4f, q2=%.4f, u2=%.4f\n', best(1:4));
        fprintf('最大固有値: %.6f\n', best(5));
        fprintf('固定点: [%.6f, %.6f, %.6f, %.6f]\n', best(6:9));
        
        % 成功例一覧（上位10個）
        fprintf('\n📊 【成功例一覧】（安定性順）\n');
        fprintf('No. | q1      | u1      | q2      | u2      | 最大固有値\n');
        fprintf('----|---------|---------|---------|---------|----------\n');
        
        display_count = min(10, success_count);
        for i = 1:display_count
            r = sorted_results(i, :);
            fprintf('%3d | %7.4f | %7.4f | %7.4f | %7.4f | %9.6f\n', ...
                    i, r(1:5));
        end
        
        if success_count > 10
            fprintf('... 他 %d 個の成功例\n', success_count - 10);
        end
        
        % ワークスペースに結果を保存
        assignin('base', 'grid_results', sorted_results);
        assignin('base', 'best_condition', best);
        assignin('base', 'search_params', params);
        
        fprintf('\n📁 結果がワークスペースに保存されました:\n');
        fprintf('- grid_results: 全成功例（安定性順）\n');
        fprintf('- best_condition: 最良の条件\n');
        fprintf('- search_params: 探索パラメータ\n');
        
        % 🎨 結果の可視化
        fprintf('\n🎨 結果を可視化中...\n');
        visualize_results_simple(sorted_results, all_conditions, params);
        
        % 💾 結果をスプレッドシートに保存
        fprintf('\n💾 結果をスプレッドシートに保存中...\n');
        save_results_to_spreadsheet(sorted_results, all_conditions, params);
        
    else
        fprintf('\n❌ 安定な条件が見つかりませんでした。\n');
        fprintf('💡 対策案:\n');
        fprintf('  - 探索範囲を広げてみてください\n');
        fprintf('  - 刻み幅を細かくしてみてください\n');
        fprintf('  - 元の論文の条件 [0.2, -0.2, 0.4, -0.3] 周辺を探索してください\n');
    end
end

function visualize_results_simple(results, all_conditions, params)
% シンプルな3D可視化関数
    
    if size(results, 1) == 0
        fprintf('可視化する成功データがありません。\n');
        return;
    end
    
    % 成功した条件のデータ
    success_q1 = results(:, 1);
    success_u1 = results(:, 2);
    success_q2 = results(:, 3);
    success_u2 = results(:, 4);
    success_eigenvalues = results(:, 5);
    
    % 全条件のデータ
    all_q1 = all_conditions(:, 1);
    all_u1 = all_conditions(:, 2);
    all_q2 = all_conditions(:, 3);
    all_u2 = all_conditions(:, 4);
    
    %% 3D散布図を作成
    figure('Name', 'パラメータ空間での安定領域 (3D)', 'Position', [100, 100, 1200, 800]);
    
    % サブプロット1: q1-u1-q2（色: u2）
    subplot(2, 2, 1);
    % 失敗した点を薄く表示
    scatter3(all_q1, all_u1, all_q2, 15, [0.8 0.8 0.8], 'o', 'MarkerFaceAlpha', 0.3);
    hold on;
    % 成功した点を色分けして表示
    scatter3(success_q1, success_u1, success_q2, 80, success_u2, 'filled', 'MarkerEdgeColor', 'k');
    xlabel('q1 (スタンス脚角度)');
    ylabel('u1 (スタンス脚角速度)');
    zlabel('q2 (スイング脚角度)');
    title('q1-u1-q2 (色: u2)');
    colorbar;
    grid on;
    view(45, 30);
    % 最良の点をハイライト
    plot3(success_q1(1), success_u1(1), success_q2(1), 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'red', 'LineWidth', 2);
    hold off;
    
    % サブプロット2: q1-u1-u2（色: q2）
    subplot(2, 2, 2);
    scatter3(all_q1, all_u1, all_u2, 15, [0.8 0.8 0.8], 'o', 'MarkerFaceAlpha', 0.3);
    hold on;
    scatter3(success_q1, success_u1, success_u2, 80, success_q2, 'filled', 'MarkerEdgeColor', 'k');
    xlabel('q1 (スタンス脚角度)');
    ylabel('u1 (スタンス脚角速度)');
    zlabel('u2 (スイング脚角速度)');
    title('q1-u1-u2 (色: q2)');
    colorbar;
    grid on;
    view(45, 30);
    plot3(success_q1(1), success_u1(1), success_u2(1), 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'red', 'LineWidth', 2);
    hold off;
    
    % サブプロット3: q1-q2-u2（色: u1）
    subplot(2, 2, 3);
    scatter3(all_q1, all_q2, all_u2, 15, [0.8 0.8 0.8], 'o', 'MarkerFaceAlpha', 0.3);
    hold on;
    scatter3(success_q1, success_q2, success_u2, 80, success_u1, 'filled', 'MarkerEdgeColor', 'k');
    xlabel('q1 (スタンス脚角度)');
    ylabel('q2 (スイング脚角度)');
    zlabel('u2 (スイング脚角速度)');
    title('q1-q2-u2 (色: u1)');
    colorbar;
    grid on;
    view(45, 30);
    plot3(success_q1(1), success_q2(1), success_u2(1), 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'red', 'LineWidth', 2);
    hold off;
    
    % サブプロット4: 安定性表示（色: 最大固有値）
    subplot(2, 2, 4);
    scatter3(success_q1, success_u1, success_q2, 80, success_eigenvalues, 'filled', 'MarkerEdgeColor', 'k');
    xlabel('q1 (スタンス脚角度)');
    ylabel('u1 (スタンス脚角速度)');
    zlabel('q2 (スイング脚角度)');
    title('安定性 (色: 最大固有値)');
    c = colorbar;
    ylabel(c, '最大固有値');
    grid on;
    view(45, 30);
    plot3(success_q1(1), success_u1(1), success_q2(1), 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'red', 'LineWidth', 2);
    
    % 統計情報を表示
    fprintf('\n📊 可視化統計:\n');
    fprintf('成功条件数: %d\n', length(success_q1));
    fprintf('最良条件 (赤ダイヤ): q1=%.3f, u1=%.3f, q2=%.3f, u2=%.3f\n', ...
            success_q1(1), success_u1(1), success_q2(1), success_u2(1));
    fprintf('最大固有値: %.6f\n', success_eigenvalues(1));
    
    % パラメータ範囲
    fprintf('\n🎯 成功パラメータ範囲:\n');
    fprintf('q1: %.3f ～ %.3f (平均: %.3f)\n', min(success_q1), max(success_q1), mean(success_q1));
    fprintf('u1: %.3f ～ %.3f (平均: %.3f)\n', min(success_u1), max(success_u1), mean(success_u1));
    fprintf('q2: %.3f ～ %.3f (平均: %.3f)\n', min(success_q2), max(success_q2), mean(success_q2));
    fprintf('u2: %.3f ～ %.3f (平均: %.3f)\n', min(success_u2), max(success_u2), mean(success_u2));
    
    fprintf('\n💡 図の見方:\n');
    fprintf('- グレーの点: 失敗した条件\n');
    fprintf('- カラーの点: 成功した条件（色は4番目のパラメータ）\n');
    fprintf('- 赤いダイヤ: 最も安定な条件\n');
    fprintf('- 3D図は回転・ズーム可能です\n');
end
function [is_success, max_eigenvalue, fixed_point] = evaluate_condition(z0)
% 個別条件の評価関数
% pranav_passivewalker_originの関数を使用

    is_success = false;
    max_eigenvalue = inf;
    fixed_point = [NaN, NaN, NaN, NaN];
    
    try
        % Walker設定（pranav_passivewalker_originと同じ）
        walker.M = 1000; walker.m = 1.0; walker.I = 0.00; walker.l = 1.0; walker.w = 0.0; 
        walker.c = 1.0;  walker.r = 0.3; walker.g = 1.0; walker.gam = 0.009; 
        
        % 固定点探索（pranav_passivewalker_originのfsolveと同じ設定）
        options = optimset('TolFun',1e-12,'TolX',1e-12,'Display','off');
        [zstar, fval, exitflag] = fsolve(@(z) fixedpt(z, walker), z0, options);
        
        % 収束チェック
        if exitflag ~= 1
            return;
        end
        
        % 残差チェック
        if norm(fval) > 1e-6
            return;
        end
        
        fixed_point = zstar;
        
        % 安定性解析（エラーハンドリング強化）
        try
            J = partialder(@(z) onestep(z, walker), zstar, walker);
            
            % ヤコビアンの妥当性チェック
            if any(any(~isfinite(J)))
                return;
            end
            
            eigenvalues = eig(J);
            
            % 固有値の妥当性チェック
            if any(~isfinite(eigenvalues))
                return;
            end
            
            max_eigenvalue = max(abs(eigenvalues));
            
            % 安定条件（全ての固有値の絶対値が1未満）
            if max_eigenvalue < 1
                is_success = true;
            end
            
        catch ME
            % 安定性解析でエラーが発生した場合
            % デバッグ用（必要に応じてコメントアウト）
            % fprintf('安定性解析エラー: %s\n', ME.message);
            return;
        end
        
    catch ME
        % 全体でエラーが発生した場合
        % デバッグ用（必要に応じてコメントアウト）
        % fprintf('評価関数エラー: %s\n', ME.message);
        is_success = false;
    end
end

% pranav_passivewalker_originから必要な関数をコピー
function zdiff = fixedpt(z0, walker)
    zdiff = onestep(z0, walker) - z0; 
end

function J = partialder(FUN, z, walker)
    pert = 1e-5;
    n = length(z);
    J = zeros(n, n);
    
    % Using central difference, accuracy quadratic
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
        flag = 0; % send only last state, for root finder and jacobian
        steps = 1;
    end

    q1 = z0(1);
    u1 = z0(2);
    q2 = z0(3);
    u2 = z0(4);

    % Derived variables
    TE = 1/2*m*(((-l*cos(q1)-r)*u1-u1*(-c*cos(q1)+w*sin(q1)))^2+(-l*sin(q1)*u1+u1*(c*sin(q1)+w*cos(q1)))^2)+1/2*m*(((-l*cos(q1)-r)*u1-(u1-u2)*(-c*cos(q1-q2)+w*sin(q1-q2)))^2+(-l*sin(q1)*u1+(u1-u2)*(c*sin(q1-q2)+w*cos(q1-q2)))^2)+1/2*M*((-l*cos(q1)-r)^2*u1^2+l^2*sin(q1)^2*u1^2)+1/2*I*(u1^2+(u1-u2)^2)+2*m*g*cos(gam)*r+2*m*g*l*cos(gam-q1)-m*g*c*cos(gam-q1)-m*g*w*sin(gam-q1)+2*m*g*sin(gam)*r*q1-m*g*c*cos(gam-q1+q2)-m*g*w*sin(gam-q1+q2)+M*g*cos(gam)*r+M*g*l*cos(gam-q1)+M*g*sin(gam)*r*q1; 
    xp1 = 0;
    xh = -l*sin(q1) - r*q1 + xp1;
    vxh = (-l*cos(q1)-r)*u1; 
    yh =  l*cos(q1) + r;
    vyh = -l*sin(q1)*u1; 
    
    z0 = [q1 u1 q2 u2 TE xh vxh yh vyh];

    t0 = 0; 
    dt = 5; % might need to be changed based on time taken for one step
    time_stamps = 100;
    t_ode = t0;
    z_ode = z0;

    for i = 1:steps
        options = odeset('abstol',1e-13,'reltol',1e-13,'events',@(t,z) collision(t,z));
        tspan = linspace(t0,t0+dt,time_stamps);
        [t_temp, z_temp] = ode113(@(t,z) single_stance(t,z,walker), tspan, z0, options);
        
        zplus = heelstrike(t_temp(end),z_temp(end,:),walker); 
        
        z0 = zplus;
        t0 = t_temp(end);
        
        % Ignore time stamps for heelstrike and first integration point
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

    Th = 0;   % external hip torque, if needed               

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

function [gstop, isterminal,direction] = collision(t,z)
    q1 = z(1); q2 = z(3); 

    gstop = -q2 + 2*q1;
    if (q2 > -0.05) % allow legs to pass through for small hip angles
        isterminal = 0;
    else
        isterminal = 1; % ode should terminate is conveyed by 1
    end
    direction = -1; % The t_final can be approached by any direction
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

function save_results_to_spreadsheet(results, all_conditions, params)
% グリッドサーチ結果をスプレッドシートに保存する関数

    % タイムスタンプ付きファイル名を作成
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    
    fprintf('\n💾 結果をスプレッドシートに保存中...\n');
    
    %% 1. 成功条件の詳細結果を保存
    if size(results, 1) > 0
        % 成功結果のテーブル作成
        success_table = table();
        success_table.No = (1:size(results, 1))';
        success_table.q1_stance_angle = results(:, 1);
        success_table.u1_stance_velocity = results(:, 2);
        success_table.q2_swing_angle = results(:, 3);
        success_table.u2_swing_velocity = results(:, 4);
        success_table.max_eigenvalue = results(:, 5);
        success_table.fixed_point_q1 = results(:, 6);
        success_table.fixed_point_u1 = results(:, 7);
        success_table.fixed_point_q2 = results(:, 8);
        success_table.fixed_point_u2 = results(:, 9);
        
        % 安定性評価を追加
        stability_rating = cell(size(results, 1), 1);
        for i = 1:size(results, 1)
            if results(i, 5) < 0.5
                stability_rating{i} = '非常に安定';
            elseif results(i, 5) < 0.7
                stability_rating{i} = '安定';
            elseif results(i, 5) < 0.9
                stability_rating{i} = 'やや安定';
            else
                stability_rating{i} = '弱い安定';
            end
        end
        success_table.stability_rating = stability_rating;
        
        % 成功結果を保存
        success_filename = sprintf('walker_gridsearch_success_%s.xlsx', timestamp);
        try
            writetable(success_table, success_filename, 'Sheet', 'Success_Results');
            fprintf('✅ 成功条件を保存: %s\n', success_filename);
        catch
            % Excelが使えない場合はCSVで保存
            success_filename_csv = sprintf('walker_gridsearch_success_%s.csv', timestamp);
            writetable(success_table, success_filename_csv);
            fprintf('✅ 成功条件を保存: %s\n', success_filename_csv);
        end
    end
    
    %% 2. 全探索結果を保存（成功/失敗含む）
    % 全結果のテーブル作成
    all_table = table();
    all_table.No = (1:size(all_conditions, 1))';
    all_table.q1_stance_angle = all_conditions(:, 1);
    all_table.u1_stance_velocity = all_conditions(:, 2);
    all_table.q2_swing_angle = all_conditions(:, 3);
    all_table.u2_swing_velocity = all_conditions(:, 4);
    
    % 成功/失敗フラグを追加
    success_flag = false(size(all_conditions, 1), 1);
    max_eigenvalue_all = nan(size(all_conditions, 1), 1);
    
    for i = 1:size(all_conditions, 1)
        % この条件が成功結果に含まれているかチェック
        if size(results, 1) > 0
            match_idx = find(abs(results(:,1) - all_conditions(i,1)) < 1e-6 & ...
                           abs(results(:,2) - all_conditions(i,2)) < 1e-6 & ...
                           abs(results(:,3) - all_conditions(i,3)) < 1e-6 & ...
                           abs(results(:,4) - all_conditions(i,4)) < 1e-6, 1);
            if ~isempty(match_idx)
                success_flag(i) = true;
                max_eigenvalue_all(i) = results(match_idx, 5);
            end
        end
    end
    
    all_table.success = success_flag;
    all_table.max_eigenvalue = max_eigenvalue_all;
    
    % 結果ステータスを追加
    result_status = cell(size(all_conditions, 1), 1);
    for i = 1:size(all_conditions, 1)
        if success_flag(i)
            result_status{i} = '成功';
        else
            result_status{i} = '失敗';
        end
    end
    all_table.result_status = result_status;
    
    % 全結果を保存
    all_filename = sprintf('walker_gridsearch_all_%s.xlsx', timestamp);
    try
        writetable(all_table, all_filename, 'Sheet', 'All_Results');
        fprintf('✅ 全探索結果を保存: %s\n', all_filename);
    catch
        % Excelが使えない場合はCSVで保存
        all_filename_csv = sprintf('walker_gridsearch_all_%s.csv', timestamp);
        writetable(all_table, all_filename_csv);
        fprintf('✅ 全探索結果を保存: %s\n', all_filename_csv);
    end
    
    %% 3. 探索パラメータと統計を保存
    % 統計情報のテーブル作成
    stats_table = table();
    
    % 探索パラメータ
    param_names = {'q1_min'; 'q1_max'; 'q1_step'; 'q1_points'; ...
                   'u1_min'; 'u1_max'; 'u1_step'; 'u1_points'; ...
                   'q2_min'; 'q2_max'; 'q2_step'; 'q2_points'; ...
                   'u2_min'; 'u2_max'; 'u2_step'; 'u2_points'; ...
                   'total_combinations'; 'success_count'; 'success_rate_percent'};
    
    param_values = [params.q1_min; params.q1_max; params.q1_step; length(params.q1_range); ...
                    params.u1_min; params.u1_max; params.u1_step; length(params.u1_range); ...
                    params.q2_min; params.q2_max; params.q2_step; length(params.q2_range); ...
                    params.u2_min; params.u2_max; params.u2_step; length(params.u2_range); ...
                    params.total_combinations; size(results, 1); ...
                    100 * size(results, 1) / params.total_combinations];
    
    stats_table.Parameter = param_names;
    stats_table.Value = param_values;
    
    % 統計情報を追加
    if size(results, 1) > 0
        % 成功パラメータの統計
        success_stats_names = {'best_q1'; 'best_u1'; 'best_q2'; 'best_u2'; 'best_eigenvalue'; ...
                              'mean_q1'; 'mean_u1'; 'mean_q2'; 'mean_u2'; 'mean_eigenvalue'; ...
                              'std_q1'; 'std_u1'; 'std_q2'; 'std_u2'; 'std_eigenvalue'};
        
        success_stats_values = [results(1, 1:5)'; ...  % 最良条件
                               mean(results(:, 1:5))'; ...  % 平均
                               std(results(:, 1:5))'];      % 標準偏差
        
        success_stats_table = table(success_stats_names, success_stats_values, ...
                                   'VariableNames', {'Statistic', 'Value'});
        
        % 統計を保存
        stats_filename = sprintf('walker_gridsearch_stats_%s.xlsx', timestamp);
        try
            writetable(stats_table, stats_filename, 'Sheet', 'Search_Parameters');
            writetable(success_stats_table, stats_filename, 'Sheet', 'Success_Statistics');
            fprintf('✅ 統計情報を保存: %s\n', stats_filename);
        catch
            % CSVで保存
            stats_filename_csv = sprintf('walker_gridsearch_stats_%s.csv', timestamp);
            writetable(stats_table, stats_filename_csv);
            success_stats_filename_csv = sprintf('walker_gridsearch_success_stats_%s.csv', timestamp);
            writetable(success_stats_table, success_stats_filename_csv);
            fprintf('✅ 統計情報を保存: %s, %s\n', stats_filename_csv, success_stats_filename_csv);
        end
    end
    
    %% 4. 要約レポートを作成
    report_filename = sprintf('walker_gridsearch_report_%s.txt', timestamp);
    fid = fopen(report_filename, 'w');
    
    fprintf(fid, '=== パッシブウォーカー グリッドサーチ結果レポート ===\n');
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
    fprintf(fid, '成功条件数: %d / %d (%.1f%%)\n', ...
            size(results, 1), params.total_combinations, ...
            100 * size(results, 1) / params.total_combinations);
    
    if size(results, 1) > 0
        fprintf(fid, '\n【最良条件】\n');
        fprintf(fid, 'q1 = %.6f\n', results(1, 1));
        fprintf(fid, 'u1 = %.6f\n', results(1, 2));
        fprintf(fid, 'q2 = %.6f\n', results(1, 3));
        fprintf(fid, 'u2 = %.6f\n', results(1, 4));
        fprintf(fid, '最大固有値 = %.6f\n', results(1, 5));
        
        fprintf(fid, '\n【生成ファイル】\n');
        if exist(success_filename, 'file')
            fprintf(fid, '- %s: 成功条件の詳細\n', success_filename);
        end
        if exist(all_filename, 'file')
            fprintf(fid, '- %s: 全探索結果\n', all_filename);
        end
        if exist(stats_filename, 'file')
            fprintf(fid, '- %s: 統計情報\n', stats_filename);
        end
    end
    
    fclose(fid);
    fprintf('✅ レポートを保存: %s\n', report_filename);
    
    fprintf('\n📁 保存完了！以下のファイルが作成されました:\n');
    if size(results, 1) > 0
        fprintf('  📊 成功条件: %s\n', success_filename);
    end
    fprintf('  📋 全結果: %s\n', all_filename);
    if size(results, 1) > 0
        fprintf('  📈 統計: %s\n', stats_filename);
    end
    fprintf('  📄 レポート: %s\n', report_filename);
end