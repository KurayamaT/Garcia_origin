function run_grid_search_standalone()
% 完全スタンドアロン版グリッドサーチ - 必要な関数をすべて内蔵

    % Walker設定
    walker.M = 1000; walker.m = 1.0; walker.I = 0.00; walker.l = 1.0; walker.w = 0.0;
    walker.c = 1.0;  walker.r = 0.3; walker.g = 1.0; walker.gam = 0.009;
    
    % 探索範囲設定（元の論文条件周辺）
    q1_range = 0.15:0.05:0.25;  % 元の条件0.2を含む
    u1_range = -0.25:0.05:-0.15;  % 元の条件-0.2を含む  
    q2_range = 0.35:0.05:0.45;  % 元の条件0.4を含む
    u2_range = -0.35:0.05:-0.25;  % 元の条件-0.3を含む
    
    fprintf('\n=== 直立歩行グリッドサーチ ===\n');
    fprintf('探索範囲:\n');
    fprintf('q1: %.2f から %.2f まで (%.2f刻み, %d点)\n', min(q1_range), max(q1_range), q1_range(2)-q1_range(1), length(q1_range));
    fprintf('u1: %.2f から %.2f まで (%.2f刻み, %d点)\n', min(u1_range), max(u1_range), u1_range(2)-u1_range(1), length(u1_range));
    fprintf('q2: %.2f から %.2f まで (%.2f刻み, %d点)\n', min(q2_range), max(q2_range), q2_range(2)-q2_range(1), length(q2_range));
    fprintf('u2: %.2f から %.2f まで (%.2f刻み, %d点)\n', min(u2_range), max(u2_range), u2_range(2)-u2_range(1), length(u2_range));
    
    % 全組み合わせ作成
    [Q1, U1, Q2, U2] = ndgrid(q1_range, u1_range, q2_range, u2_range);
    all_z0 = [Q1(:), U1(:), Q2(:), U2(:)];
    total = size(all_z0, 1);
    
    fprintf('総探索数: %d\n', total);
    
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
    
    fprintf('\n探索中... (成功例は即座に表示されます)\n\n');
    
    % 探索実行
    results = cell(total, 1);
    success_count = 0;
    tic;
    
    if use_parallel
        % 並列処理
        parfor idx = 1:total
            z0 = all_z0(idx, :);
            result = evaluate_walker_condition(z0, walker);
            
            if result.success
                fprintf('【成功】 q1=%6.3f, u1=%6.3f, q2=%6.3f, u2=%6.3f | λ_max=%6.4f | θ_max=%5.1f°\n', ...
                        z0(1), z0(2), z0(3), z0(4), result.max_eig, result.max_angle*180/pi);
            end
            
            results{idx} = result;
        end
    else
        % 逐次処理
        for idx = 1:total
            if mod(idx, 50) == 0
                fprintf('進捗: %d/%d (%.1f%%)\n', idx, total, 100*idx/total);
            end
            
            z0 = all_z0(idx, :);
            result = evaluate_walker_condition(z0, walker);
            
            if result.success
                success_count = success_count + 1;
                fprintf('【成功 #%d】 q1=%6.3f, u1=%6.3f, q2=%6.3f, u2=%6.3f | λ_max=%6.4f | θ_max=%5.1f°\n', ...
                        success_count, z0(1), z0(2), z0(3), z0(4), result.max_eig, result.max_angle*180/pi);
            end
            
            results{idx} = result;
        end
    end
    
    elapsed = toc;
    
    % 結果集計と表示
    results_array = [results{:}];
    successful = results_array([results_array.success]);
    
    fprintf('\n' + repmat('=', 1, 60) + '\n');
    fprintf('🎉 グリッドサーチ完了！\n');
    fprintf('時間: %.1f秒\n', elapsed);
    fprintf('成功: %d/%d (%.1f%%)\n', length(successful), total, 100*length(successful)/total);
    fprintf(repmat('=', 1, 60) + '\n');
    
    if ~isempty(successful)
        % 安定性でソート
        [~, sort_idx] = sort([successful.max_eig]);
        sorted_successful = successful(sort_idx);
        best = sorted_successful(1);
        
        fprintf('\n🏆 【最良の初期条件】\n');
        fprintf('q1=%.3f, u1=%.3f, q2=%.3f, u2=%.3f\n', ...
                best.q1, best.u1, best.q2, best.u2);
        fprintf('最大固有値: %.6f (安定性指標)\n', best.max_eig);
        fprintf('最大角度: %.1f° (転倒回避)\n', best.max_angle*180/pi);
        fprintf('固定点: [%.4f, %.4f, %.4f, %.4f]\n', best.zstar);
        
        % 成功例のサマリー
        fprintf('\n📊 【成功例一覧】（安定性順）\n');
        fprintf('No. | q1      | u1      | q2      | u2      | λ_max    | θ_max   | 評価\n');
        fprintf('----|---------|---------|---------|---------|----------|---------|------\n');
        
        for i = 1:min(15, length(sorted_successful))
            s = sorted_successful(i);
            stability_rating = get_stability_rating(s.max_eig);
            fprintf('%3d | %7.3f | %7.3f | %7.3f | %7.3f | %8.4f | %6.1f° | %s\n', ...
                    i, s.q1, s.u1, s.q2, s.u2, s.max_eig, s.max_angle*180/pi, stability_rating);
        end
        
        if length(successful) > 15
            fprintf('... 他 %d 個の成功例\n', length(successful)-15);
        end
        
        % 統計情報
        fprintf('\n📈 【統計情報】\n');
        fprintf('固有値範囲: %.4f 〜 %.4f\n', min([successful.max_eig]), max([successful.max_eig]));
        fprintf('角度範囲: %.1f° 〜 %.1f°\n', min([successful.max_angle])*180/pi, max([successful.max_angle])*180/pi);
        fprintf('平均固有値: %.4f\n', mean([successful.max_eig]));
        
        % ワークスペースに保存
        assignin('base', 'grid_search_results', results_array);
        assignin('base', 'successful_conditions', sorted_successful);
        assignin('base', 'best_walker_condition', best);
        
        fprintf('\n変数がワークスペースに保存されました:\n');
        fprintf('- grid_search_results: 全結果\n');
        fprintf('- successful_conditions: 成功例（安定性順）\n');
        fprintf('- best_walker_condition: 最良の条件\n');
        
    else
        fprintf('\n❌ 成功例が見つかりませんでした。\n');
        fprintf('💡 対策案:\n');
        fprintf('  - 探索範囲を広げる\n');
        fprintf('  - 刻み幅を細かくする\n');
        fprintf('  - 既知の良い条件 [0.2, -0.2, 0.4, -0.3] 周辺を重点的に探索\n');
        
        assignin('base', 'grid_search_results', results_array);
    end
end

%% 個別評価関数
function result = evaluate_walker_condition(z0, walker)
    % 初期化
    result.q1 = z0(1); 
    result.u1 = z0(2);
    result.q2 = z0(3); 
    result.u2 = z0(4);
    result.success = false;
    result.max_eig = inf;
    result.max_angle = inf;
    result.zstar = [];
    
    try
        % ステップ1: 固定点探索
        options = optimset('TolFun',1e-12,'TolX',1e-12,'Display','off');
        [zstar, ~, exitflag] = fsolve(@(z) fixedpt(z, walker), z0, options);
        
        if exitflag == 1
            result.zstar = zstar;
            
            % ステップ2: 安定性チェック
            J = partialder(@(z) onestep(z, walker), zstar, walker);
            eigenvalues = eig(J);
            result.max_eig = max(abs(eigenvalues));
            
            if result.max_eig < 1
                % ステップ3: 歩行テスト（3ステップで高速化）
                [z, ~] = onestep(z0, walker, 3);
                
                % 転倒チェック（修正版 - 角度のみで判定）
                result.max_angle = max(abs(z(:,1)));
                result.success = (result.max_angle < pi/2);  % 90度未満で成功
            end
        end
    catch
        % エラーの場合はデフォルト値のまま
    end
end

%% 安定性評価関数
function rating = get_stability_rating(max_eig)
    if max_eig < 0.5
        rating = '★★★';  % 非常に安定
    elseif max_eig < 0.7
        rating = '★★☆';  % 安定
    elseif max_eig < 0.9
        rating = '★☆☆';  % やや安定
    else
        rating = '☆☆☆';  % 不安定寄り
    end
end

%% ===== PASSIVE WALKER FUNCTIONS =====
% Copy_of_pranav_passivewalker_originから必要な関数を抽出

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
        J(:,i) = (feval(FUN, ztemp1, walker) - feval(FUN, ztemp2, walker));
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
        options = odeset('abstol',1e-13,'reltol',1e-13,'events',@(t,z) collision(t,z,walker));
        tspan = linspace(t0,t0+dt,time_stamps);
        [t_temp, z_temp] = ode113(@(t,z) single_stance(t,z,walker), tspan, z0, options);
        
        zplus = heelstrike(t_temp(end), z_temp(end,:), walker); 
        
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

function zdot = single_stance(t, z, walker)  
    q1 = z(1);   u1 = z(2);                         
    q2 = z(3);   u2 = z(4);                         
    xh = z(6);  vxh = z(7);                       
    yh = z(8);  vyh = z(9);                     

    M = walker.M;  m = walker.m; I = walker.I;   
    l = walker.l;  c = walker.c; w = walker.w;   
    r = walker.r;  g = walker.g; gam = walker.gam;

    Th = 0; % external hip torque, if needed               

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

function [gstop, isterminal, direction] = collision(t, z, walker)
    q1 = z(1); q2 = z(3); 

    gstop = -q2 + 2*q1;
    if (q2 > -0.05) % allow legs to pass through for small hip angles
        isterminal = 0;
    else
        isterminal = 1; % ode should terminate is conveyed by 1
    end
    direction = -1; % The t_final can be approached by any direction
end

function zplus = heelstrike(t, z, walker)      
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