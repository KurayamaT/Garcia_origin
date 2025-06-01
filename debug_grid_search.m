function debug_grid_search_complete()
% 完全版デバッグスクリプト - 必要な関数をすべて含む
% 既知の安定な初期条件をテストして、なぜ成功例が出ないか調査

    %% Walker設定
    walker.M = 1000; walker.m = 1.0; walker.I = 0.00; walker.l = 1.0; walker.w = 0.0;
    walker.c = 1.0;  walker.r = 0.3; walker.g = 1.0; walker.gam = 0.009;
    
    %% テスト1: 元の論文の初期条件を直接テスト
    fprintf('=== テスト1: 元の論文の初期条件 ===\n');
    z0_original = [0.2, -0.2, 0.4, -0.3];
    test_single_condition(z0_original, walker, '元の論文');
    
    %% テスト2: Garcia論文の正確な固定点をテスト
    fprintf('\n=== テスト2: Garcia論文の固定点 ===\n');
    zstar_garcia = [0.200161072169750, -0.199906060087682, 0.400322144339512, -0.015805473227965];
    test_single_condition(zstar_garcia, walker, 'Garcia固定点');
    
    %% テスト3: 探索範囲に含まれているか確認
    fprintf('\n=== テスト3: 探索範囲の確認 ===\n');
    if evalin('base', 'exist(''q1_range'', ''var'')') && ...
       evalin('base', 'exist(''u1_range'', ''var'')') && ...
       evalin('base', 'exist(''q2_range'', ''var'')') && ...
       evalin('base', 'exist(''u2_range'', ''var'')')
        
        q1_range = evalin('base', 'q1_range');
        u1_range = evalin('base', 'u1_range');
        q2_range = evalin('base', 'q2_range');
        u2_range = evalin('base', 'u2_range');
        
        fprintf('現在の探索範囲:\n');
        fprintf('q1: [%.2f, %.2f] (%.3f刻み)\n', min(q1_range), max(q1_range), q1_range(2)-q1_range(1));
        fprintf('u1: [%.2f, %.2f] (%.3f刻み)\n', min(u1_range), max(u1_range), u1_range(2)-u1_range(1));
        fprintf('q2: [%.2f, %.2f] (%.3f刻み)\n', min(q2_range), max(q2_range), q2_range(2)-q2_range(1));
        fprintf('u2: [%.2f, %.2f] (%.3f刻み)\n', min(u2_range), max(u2_range), u2_range(2)-u2_range(1));
        
        % 元の条件が範囲内か確認
        q1_in = any(abs(q1_range - 0.2) < 1e-6);
        u1_in = any(abs(u1_range - (-0.2)) < 1e-6);
        q2_in = any(abs(q2_range - 0.4) < 1e-6);
        u2_in = any(abs(u2_range - (-0.3)) < 1e-6);
        
        fprintf('\n元の条件が範囲内にあるか:\n');
        fprintf('q1=0.2: %s\n', iif(q1_in, '○', '×'));
        fprintf('u1=-0.2: %s\n', iif(u1_in, '○', '×'));
        fprintf('q2=0.4: %s\n', iif(q2_in, '○', '×'));
        fprintf('u2=-0.3: %s\n', iif(u2_in, '○', '×'));
    else
        fprintf('探索範囲が設定されていません。\n');
    end
    
    %% テスト4: 評価関数の動作確認
    fprintf('\n=== テスト4: 評価関数の詳細テスト ===\n');
    result = evaluate_single_debug(z0_original, walker);
    
    %% テスト5: 範囲を含むように設定して再テスト
    fprintf('\n=== テスト5: 元の条件を含む範囲で小規模テスト ===\n');
    test_q1 = 0.1:0.1:0.3;
    test_u1 = -0.3:0.1:-0.1;
    test_q2 = 0.3:0.1:0.5;
    test_u2 = -0.4:0.1:-0.2;
    
    fprintf('テスト範囲:\n');
    fprintf('q1: '); disp(test_q1);
    fprintf('u1: '); disp(test_u1);
    fprintf('q2: '); disp(test_q2);
    fprintf('u2: '); disp(test_u2);
    
    success_count = 0;
    total_tests = length(test_q1)*length(test_u1)*length(test_q2)*length(test_u2);
    fprintf('総テスト数: %d\n', total_tests);
    
    test_count = 0;
    for q1 = test_q1
        for u1 = test_u1
            for q2 = test_q2
                for u2 = test_u2
                    test_count = test_count + 1;
                    z0 = [q1, u1, q2, u2];
                    result = evaluate_single(z0, walker);
                    if result.success
                        success_count = success_count + 1;
                        fprintf('成功 #%d: q1=%.1f, u1=%.1f, q2=%.1f, u2=%.1f (λ_max=%.4f)\n', ...
                                success_count, q1, u1, q2, u2, result.max_eig);
                    end
                    
                    % 進捗表示
                    if mod(test_count, 100) == 0
                        fprintf('進捗: %d/%d\n', test_count, total_tests);
                    end
                end
            end
        end
    end
    fprintf('\n成功数: %d / %d (%.1f%%)\n', success_count, total_tests, 100*success_count/total_tests);
    
    %% テスト6: 探索範囲の問題を特定
    fprintf('\n=== テスト6: 探索範囲の問題分析 ===\n');
    if evalin('base', 'exist(''q1_range'', ''var'')')
        analyze_search_ranges();
    end
end

%% 詳細なテスト関数
function test_single_condition(z0, walker, label)
    fprintf('\n%s: z0 = [%.3f, %.3f, %.3f, %.3f]\n', label, z0);
    
    try
        % 固定点探索
        options = optimset('TolFun',1e-12,'TolX',1e-12,'Display','off');
        [zstar, fval, exitflag] = fsolve(@(z) fixedpt(z, walker), z0, options);
        
        if exitflag == 1
            fprintf('  固定点: [%.6f, %.6f, %.6f, %.6f]\n', zstar);
            fprintf('  残差: %.2e\n', norm(fval));
            
            % 安定性解析
            J = partialder(@(z) onestep(z, walker), zstar, walker);
            eigenvalues = eig(J);
            max_eig = max(abs(eigenvalues));
            fprintf('  最大固有値: %.6f %s\n', max_eig, iif(max_eig < 1, '(安定)', '(不安定)'));
            
            if max_eig < 1
                % 歩行テスト
                try
                    [z, t] = onestep(z0, walker, 5);
                    max_angle = max(abs(z(:,1)));
                    min_height = min(z(:,8));
                    
                    fprintf('  最大角度: %.1f°\n', max_angle*180/pi);
                    fprintf('  最小高さ: %.3f\n', min_height);
                    
                    if (max_angle < pi/2) && (min_height > 0.1)
                        fprintf('  結果: ✓ 成功\n');
                    else
                        fprintf('  結果: × 失敗（転倒: 角度=%.1f°, 高さ=%.3f）\n', ...
                                max_angle*180/pi, min_height);
                    end
                catch ME
                    fprintf('  歩行テストエラー: %s\n', ME.message);
                end
            else
                fprintf('  結果: × 失敗（不安定な固定点）\n');
            end
        else
            fprintf('  結果: × 失敗（固定点が見つからない, exitflag=%d）\n', exitflag);
        end
    catch ME
        fprintf('  エラー: %s\n', ME.message);
    end
end

%% デバッグ版の評価関数
function result = evaluate_single_debug(z0, walker)
    fprintf('\n詳細評価: z0 = [%.3f, %.3f, %.3f, %.3f]\n', z0);
    
    result.q1 = z0(1); 
    result.u1 = z0(2);
    result.q2 = z0(3); 
    result.u2 = z0(4);
    result.success = false;
    result.max_eig = inf;
    result.max_angle = inf;
    
    % ステップ1: 固定点探索
    fprintf('  1. 固定点探索中...\n');
    try
        options = optimset('TolFun',1e-12,'TolX',1e-12,'Display','off');
        [zstar, fval, exitflag] = fsolve(@(z) fixedpt(z, walker), z0, options);
        
        if exitflag ~= 1
            fprintf('     → 失敗（exitflag=%d, 残差=%.2e）\n', exitflag, norm(fval));
            return;
        end
        
        result.zstar = zstar;
        fprintf('     → 成功: zstar = [%.6f, %.6f, %.6f, %.6f]\n', zstar);
        fprintf('     → 残差: %.2e\n', norm(fval));
    catch ME
        fprintf('     → エラー: %s\n', ME.message);
        return;
    end
    
    % ステップ2: 安定性チェック
    fprintf('  2. 安定性解析中...\n');
    try
        J = partialder(@(z) onestep(z, walker), zstar, walker);
        eigenvalues = eig(J);
        result.max_eig = max(abs(eigenvalues));
        fprintf('     → 最大固有値: %.6f\n', result.max_eig);
        
        if result.max_eig >= 1
            fprintf('     → 失敗（不安定）\n');
            return;
        end
    catch ME
        fprintf('     → エラー: %s\n', ME.message);
        return;
    end
    
    % ステップ3: 歩行テスト
    fprintf('  3. 歩行テスト中...\n');
    try
        [z, t] = onestep(z0, walker, 5);
        result.max_angle = max(abs(z(:,1)));
        result.min_height = min(z(:,8));
        
        fprintf('     → 最大角度: %.1f°\n', result.max_angle*180/pi);
        fprintf('     → 最小高さ: %.3f\n', result.min_height);
        
        result.success = (result.max_angle < pi/2) && (result.min_height > 0.1);
        
        if result.success
            fprintf('     → 成功！\n');
        else
            fprintf('     → 失敗（転倒）\n');
        end
    catch ME
        fprintf('     → エラー: %s\n', ME.message);
    end
end

%% 簡略版の評価関数（比較用）
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
        [zstar, ~, exitflag] = fsolve(@(z) fixedpt(z, walker), z0, options);
        
        if exitflag == 1
            result.zstar = zstar;
            
            % 安定性チェック
            J = partialder(@(z) onestep(z, walker), zstar, walker);
            eigenvalues = eig(J);
            result.max_eig = max(abs(eigenvalues));
            
            if result.max_eig < 1
                % 実際の歩行テスト
                [z, ~] = onestep(z0, walker, 5);
                
                % 転倒チェック
                result.max_angle = max(abs(z(:,1)));
                result.min_height = min(z(:,8));
                
                result.success = (result.max_angle < pi/2) && (result.min_height > 0.1);
            end
        end
    catch
        % エラーは無視
    end
end

%% 探索範囲分析
function analyze_search_ranges()
    q1_range = evalin('base', 'q1_range');
    u1_range = evalin('base', 'u1_range');
    q2_range = evalin('base', 'q2_range');
    u2_range = evalin('base', 'u2_range');
    
    fprintf('問題の可能性:\n');
    
    % 刻み幅チェック
    if length(q1_range) > 1
        q1_step = q1_range(2) - q1_range(1);
        if q1_step > 0.1
            fprintf('- q1の刻み幅が大きすぎる可能性 (%.3f)\n', q1_step);
        end
    end
    
    if length(u1_range) > 1
        u1_step = u1_range(2) - u1_range(1);
        if u1_step > 0.1
            fprintf('- u1の刻み幅が大きすぎる可能性 (%.3f)\n', u1_step);
        end
    end
    
    % 範囲チェック
    if min(u2_range) > -0.1
        fprintf('- u2の範囲が正の値に偏っている可能性\n');
    end
    
    if max(q2_range) < 0.3
        fprintf('- q2の範囲が小さすぎる可能性\n');
    end
    
    % 総数チェック
    total = length(q1_range) * length(u1_range) * length(q2_range) * length(u2_range);
    if total > 100000
        fprintf('- 探索数が多すぎる可能性 (%d)\n', total);
    end
    
    fprintf('\n推奨改善案:\n');
    fprintf('- 刻み幅を0.05以下に設定\n');
    fprintf('- u2範囲を-0.5から0に設定\n');
    fprintf('- q2範囲を0から0.6に設定\n');
end

%% Utility functions
function result = iif(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end

%% ===== PASSIVE WALKER FUNCTIONS =====
% 元のpassive walkerコードから必要な関数を抽出

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