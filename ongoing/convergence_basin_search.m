function convergence_basin_search(params)
% 実際に収束する初期値を探索する関数
% params: 探索パラメータ構造体（grid_runnerと同じ形式）

    fprintf('\n=== 収束保証付き初期値探索 ===\n');
    
    % Walker設定
    walker.M = 1000; walker.m = 1.0; walker.I = 0.00; walker.l = 1.0; walker.w = 0.0; 
    walker.c = 1.0;  walker.r = 0.3; walker.g = 1.0; walker.gam = 0.009;
    
    % 探索範囲の表示
    fprintf('\n📋 探索範囲:\n');
    fprintf('q1: %.3f ～ %.3f (%.3f刻み, %d点)\n', params.q1_min, params.q1_max, params.q1_step, length(params.q1_range));
    fprintf('u1: %.3f ～ %.3f (%.3f刻み, %d点)\n', params.u1_min, params.u1_max, params.u1_step, length(params.u1_range));
    fprintf('q2: %.3f ～ %.3f (%.3f刻み, %d点)\n', params.q2_min, params.q2_max, params.q2_step, length(params.q2_range));
    fprintf('u2: %.3f ～ %.3f (%.3f刻み, %d点)\n', params.u2_min, params.u2_max, params.u2_step, length(params.u2_range));
    
    % 全組み合わせ作成
    [Q1, U1, Q2, U2] = ndgrid(params.q1_range, params.u1_range, ...
                             params.q2_range, params.u2_range);
    all_conditions = [Q1(:), U1(:), Q2(:), U2(:)];
    total = size(all_conditions, 1);
    
    fprintf('総探索数: %d\n', total);
    
    % 結果保存用
    convergent_conditions = [];
    convergence_info = [];
    
    % 並列処理の設定
    use_parallel = false;
    if total > 100
        answer = input('\n並列処理を使用しますか？ (y/n) [推奨: y]: ', 's');
        if isempty(answer) || strcmpi(answer, 'y')
            pool = gcp('nocreate');
            if isempty(pool)
                parpool;
            end
            use_parallel = true;
        end
    end
    
    fprintf('\n🚀 収束テスト開始...\n');
    tic;
    
    % 収束テストパラメータ
    test_steps = 30;  % 収束判定のための歩数
    convergence_threshold = 1e-3;  % 収束判定閾値
    
    if use_parallel
        % 並列処理
        results_cell = cell(total, 1);
        parfor idx = 1:total
            z0 = all_conditions(idx, :);
            result = test_convergence_from_initial(z0, walker, test_steps, convergence_threshold);
            if result.converged
                results_cell{idx} = result;
            end
        end
        
        % 結果集計
        for idx = 1:total
            if ~isempty(results_cell{idx})
                result = results_cell{idx};
                convergent_conditions = [convergent_conditions; result.initial_condition];
                convergence_info = [convergence_info; result];
            end
        end
    else
        % 逐次処理
        for idx = 1:total
            if mod(idx, max(1, floor(total/20))) == 0
                fprintf('進捗: %d/%d (%.1f%%)\n', idx, total, 100*idx/total);
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
    
    elapsed = toc;
    
    % 結果まとめ
    num_convergent = size(convergent_conditions, 1);
    fprintf('\n%s\n', repmat('=', 1, 60));
    fprintf('🎉 探索完了！\n');
    fprintf('実行時間: %.1f秒\n', elapsed);
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
        
        % ワークスペースに保存
        assignin('base', 'convergent_initials', sorted_conditions);
        assignin('base', 'convergence_details', sorted_info);
        assignin('base', 'best_convergent_initial', best);
        
        fprintf('\n📁 結果がワークスペースに保存されました:\n');
        fprintf('- convergent_initials: 収束する全初期値\n');
        fprintf('- convergence_details: 収束の詳細情報\n');
        fprintf('- best_convergent_initial: 最速収束する初期値\n');
        
        % 可視化
        visualize_convergence_results(sorted_info, all_conditions, params);
        
        % CSVに保存
        save_convergence_results(sorted_info, params);
        
    else
        fprintf('\n❌ 収束する初期値が見つかりませんでした。\n');
        fprintf('💡 対策案:\n');
        fprintf('  - 探索範囲を論文の値 [0.2, -0.2, 0.4, -0.3] の近くに設定\n');
        fprintf('  - 刻み幅を細かくする\n');
        fprintf('  - 収束判定の歩数を増やす\n');
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
    
    % サブプロット1: 収束速度の分布（色で表示）
    subplot(2, 2, 1);
    scatter3(conv_initials(:,1), conv_initials(:,2), conv_initials(:,3), ...
             80, conv_steps, 'filled', 'MarkerEdgeColor', 'k');
    xlabel('q1'); ylabel('u1'); zlabel('q2');
    title('収束する初期値（色：収束歩数）');
    colorbar;
    grid on;
    view(45, 30);
    
    % サブプロット2: 収束歩数のヒストグラム
    subplot(2, 2, 2);
    histogram(conv_steps, 20);
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
    for i = 1:min(5, length(convergence_info))
        history = convergence_info(i).convergence_history;
        plot(1:length(history), history, 'LineWidth', 2);
    end
    xlabel('歩数');
    ylabel('固定点からの距離');
    title('収束履歴（最速5例）');
    set(gca, 'YScale', 'log');
    grid on;
    legend('1位', '2位', '3位', '4位', '5位', 'Location', 'northeast');
    
    fprintf('\n📊 統計情報:\n');
    fprintf('平均収束歩数: %.1f歩\n', mean(conv_steps));
    fprintf('最速収束: %d歩\n', min(conv_steps));
    fprintf('最遅収束: %d歩\n', max(conv_steps));
    fprintf('標準偏差: %.1f歩\n', std(conv_steps));
end

function save_convergence_results(convergence_info, params)
% 収束結果をCSVファイルに保存
    
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    filename = sprintf('convergent_initials_%s.csv', timestamp);
    
    % テーブル作成
    T = table();
    T.No = (1:length(convergence_info))';
    
    % 初期値
    initials = vertcat(convergence_info.initial_condition);
    T.q1_initial = initials(:, 1);
    T.u1_initial = initials(:, 2);
    T.q2_initial = initials(:, 3);
    T.u2_initial = initials(:, 4);
    
    % 固定点
    fixed = vertcat(convergence_info.fixed_point);
    T.q1_fixed = fixed(:, 1);
    T.u1_fixed = fixed(:, 2);
    T.q2_fixed = fixed(:, 3);
    T.u2_fixed = fixed(:, 4);
    
    % 収束情報
    T.steps_to_converge = [convergence_info.steps_to_converge]';
    T.final_distance = [convergence_info.final_distance]';
    
    % 初期距離
    initial_distances = zeros(length(convergence_info), 1);
    for i = 1:length(convergence_info)
        initial_distances(i) = norm(initials(i,:) - fixed(i,:));
    end
    T.initial_distance = initial_distances;
    
    % CSVに保存
    writetable(T, filename);
    fprintf('\n✅ 収束する初期値を保存: %s\n', filename);
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

% single_stance, collision, heelstrike関数も必要（grid_runnerからコピー）
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