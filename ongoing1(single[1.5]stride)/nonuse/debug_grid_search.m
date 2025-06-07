function debug_grid_search_fixed()
% 修正版デバッグスクリプト - partialder関数の呼び出しエラーを修正

    %% Walker設定
    walker.M = 1000; walker.m = 1.0; walker.I = 0.00; walker.l = 1.0; walker.w = 0.0;
    walker.c = 1.0;  walker.r = 0.3; walker.g = 1.0; walker.gam = 0.009;
    
    %% テスト1: 元の論文の初期条件を直接テスト
    fprintf('=== テスト1: 元の論文の初期条件 ===\n');
    z0_original = [0.2, -0.2, 0.4, -0.3];
    test_single_condition_fixed(z0_original, walker, '元の論文');
    
    %% テスト2: Garcia論文の正確な固定点をテスト
    fprintf('\n=== テスト2: Garcia論文の固定点 ===\n');
    zstar_garcia = [0.200161072169750, -0.199906060087682, 0.400322144339512, -0.015805473227965];
    test_single_condition_fixed(zstar_garcia, walker, 'Garcia固定点');
    
    %% テスト3: より広い範囲で小規模テスト
    fprintf('\n=== テスト3: 改善された範囲で小規模テスト ===\n');
    test_q1 = [0.15, 0.2, 0.25];
    test_u1 = [-0.25, -0.2, -0.15];
    test_q2 = [0.35, 0.4, 0.45];
    test_u2 = [-0.35, -0.3, -0.25];
    
    fprintf('テスト範囲:\n');
    fprintf('q1: '); disp(test_q1);
    fprintf('u1: '); disp(test_u1);
    fprintf('q2: '); disp(test_q2);
    fprintf('u2: '); disp(test_u2);
    
    success_count = 0;
    total_tests = length(test_q1)*length(test_u1)*length(test_q2)*length(test_u2);
    fprintf('総テスト数: %d\n', total_tests);
    
    successful_conditions = [];
    
    for q1 = test_q1
        for u1 = test_u1
            for q2 = test_q2
                for u2 = test_u2
                    z0 = [q1, u1, q2, u2];
                    result = evaluate_single_fixed(z0, walker);
                    if result.success
                        success_count = success_count + 1;
                        successful_conditions = [successful_conditions; z0, result.max_eig, result.max_angle];
                        fprintf('成功 #%d: q1=%.2f, u1=%.2f, q2=%.2f, u2=%.2f (λ_max=%.4f, θ_max=%.1f°)\n', ...
                                success_count, q1, u1, q2, u2, result.max_eig, result.max_angle*180/pi);
                    end
                end
            end
        end
    end
    
    fprintf('\n成功数: %d / %d (%.1f%%)\n', success_count, total_tests, 100*success_count/total_tests);
    
    if success_count > 0
        fprintf('\n=== 成功した初期条件一覧 ===\n');
        fprintf('No. | q1    | u1    | q2    | u2    | λ_max  | θ_max\n');
        fprintf('----|-------|-------|-------|-------|--------|-------\n');
        for i = 1:size(successful_conditions, 1)
            fprintf('%3d | %5.2f | %5.2f | %5.2f | %5.2f | %6.4f | %5.1f°\n', ...
                    i, successful_conditions(i,1:6));
        end
        
        % 最良の条件を特定
        [~, best_idx] = min(successful_conditions(:,5)); % 最小固有値
        best_condition = successful_conditions(best_idx, :);
        fprintf('\n最良の条件: q1=%.2f, u1=%.2f, q2=%.2f, u2=%.2f\n', best_condition(1:4));
        
        % ワークスペースに保存
        assignin('base', 'successful_conditions', successful_conditions);
        assignin('base', 'best_condition', best_condition);
    else
        fprintf('\n成功例が見つかりませんでした。\n');
        fprintf('推奨対策:\n');
        fprintf('1. 探索範囲をさらに広げる\n');
        fprintf('2. 刻み幅を細かくする\n');
        fprintf('3. 初期条件の設定方法を見直す\n');
    end
end

%% 修正版テスト関数
function test_single_condition_fixed(z0, walker, label)
    fprintf('\n%s: z0 = [%.3f, %.3f, %.3f, %.3f]\n', label, z0);
    
    try
        % 固定点探索
        options = optimset('TolFun',1e-12,'TolX',1e-12,'Display','off');
        [zstar, fval, exitflag] = fsolve(@(z) fixedpt(z, walker), z0, options);
        
        if exitflag == 1
            fprintf('  固定点: [%.6f, %.6f, %.6f, %.6f]\n', zstar);
            fprintf('  残差: %.2e\n', norm(fval));
            
            % 安定性解析（修正版）
            J = partialder(@onestep, zstar, walker);
            eigenvalues = eig(J);
            max_eig = max(abs(eigenvalues));
            fprintf('  最大固有値: %.6f %s\n', max_eig, iif(max_eig < 1, '(安定)', '(不安定)'));
            
            if max_eig < 1
                % 歩行テスト
                try
                    [z, t] = onestep(z0, walker, 5);
                    max_angle = max(abs(z(:,1)));
                    if size(z, 2) >= 8
                        min_height = min(z(:,8));
                    else
                        min_height = 1.0; % デフォルト値
                    end
                    
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

%% 修正版の評価関数
function result = evaluate_single_fixed(z0, walker)
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
            
            % 安定性チェック（修正版）
            J = partialder(@onestep, zstar, walker);
            eigenvalues = eig(J);
            result.max_eig = max(abs(eigenvalues));
            
            if result.max_eig < 1
                % 実際の歩行テスト
                [z, ~] = onestep(z0, walker, 3);
                
                % 転倒チェック（簡易版）
                result.max_angle = max(abs(z(:,1)));
                result.success = (result.max_angle < pi/2);
            end
        end
    catch
        % エラーは無視
    end
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
        flag = 0;
        steps = 1;
    end

    q1 = z0(1); u1 = z0(2); q2 = z0(3); u2 = z0(4);

    % Derived variables
    TE = 1/2*m*(((-l*cos(q1)-r)*u1-u1*(-c*cos(q1)+w*sin(q1)))^2+(-l*sin(q1)*u1+u1*(c*sin(q1)+w*cos(q1)))^2)+1/2*m*(((-l*cos(q1)-r)*u1-(u1-u2)*(-c*cos(q1-q2)+w*sin(q1-q2)))^2+(-l*sin(q1)*u1+(u1-u2)*(c*sin(q1-q2)+w*cos(q1-q2)))^2)+1/2*M*((-l*cos(q1)-r)^2*u1^2+l^2*sin(q1)^2*u1^2)+1/2*I*(u1^2+(u1-u2)^2)+2*m*g*cos(gam)*r+2*m*g*l*cos(gam-q1)-m*g*c*cos(gam-q1)-m*g*w*sin(gam-q1)+2*m*g*sin(gam)*r*q1-m*g*c*cos(gam-q1+q2)-m*g*w*sin(gam-q1+q2)+M*g*cos(gam)*r+M*g*l*cos(gam-q1)+M*g*sin(gam)*r*q1; 
    xp1 = 0;
    xh = -l*sin(q1) - r*q1 + xp1;
    vxh = (-l*cos(q1)-r)*u1; 
    yh =  l*cos(q1) + r;
    vyh = -l*sin(q1)*u1; 

    z0 = [q1 u1 q2 u2 TE xh vxh yh vyh];

    t0 = 0; dt = 5; time_stamps = 100;
    t_ode = t0; z_ode = z0;

    for i = 1:steps
        options = odeset('abstol',1e-13,'reltol',1e-13,'events',@(t,z) collision(t,z,walker));
        tspan = linspace(t0,t0+dt,time_stamps);
        [t_temp, z_temp] = ode113(@(t,z) single_stance(t,z,walker), tspan, z0, options);
        
        zplus = heelstrike(t_temp(end), z_temp(end,:), walker); 
        z0 = zplus; t0 = t_temp(end);
        
        t_ode = [t_ode; t_temp(2:end)];
        z_ode = [z_ode; z_temp(2:end,:)];
    end

    z = zplus(1:4);
    if flag == 1
       z = z_ode; t = t_ode;
    end
end

function zdot = single_stance(t, z, walker)  
    q1 = z(1); u1 = z(2); q2 = z(3); u2 = z(4);                         
    xh = z(6); vxh = z(7); yh = z(8); vyh = z(9);                     

    M = walker.M; m = walker.m; I = walker.I;   
    l = walker.l; c = walker.c; w = walker.w;   
    r = walker.r; g = walker.g; gam = walker.gam;

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

    ud1 = X(1); ud2 = X(2);                                       

    DTE = -ud1*I*u2+2*ud1*m*u1*r^2+m*u1*r*u2^2*c*sin(q1-q2)+m*u1*r*u2^2*w*cos(q1-q2)-m*u2^2*l*u1*c*sin(q2)+u2*m*g*c*sin(gam-q1+q2)-u2*m*g*w*cos(gam-q1+q2)+2*ud1*I*u1+ud2*I*u2+m*u2^2*l*u1*w*cos(q2)+2*ud1*m*u1*c^2+ud2*m*u2*c^2-ud2*I*u1+ud1*m*u2*c*l*cos(q2)+ud1*m*u2*w*l*sin(q2)-2*ud1*m*u1*l*c*cos(q2)-2*ud1*m*l*u1*w*sin(q2)-m*u2*u1^2*w*l*cos(q2)+m*u2*u1^2*c*l*sin(q2)+2*ud1*m*u1*w^2+ud2*m*u2*w^2+ud2*m*u1*l*c*cos(q2)+ud1*M*l^2*u1-ud2*m*u1*w^2-ud1*m*u2*c^2-ud2*m*u1*c^2+2*ud1*m*l^2*u1-ud1*m*u2*w^2-2*ud1*m*l*u1*c+2*ud1*M*l*cos(q1)*u1*r+4*ud1*m*l*cos(q1)*u1*r-2*ud1*m*u1*r*c*cos(q1)+2*ud1*m*u1*r*w*sin(q1)-2*ud1*m*u1*r*c*cos(q1-q2)-2*m*u1^3*r*l*sin(q1)+m*u1^3*r*c*sin(q1)+m*u1^3*r*w*cos(q1)+m*u1^3*r*c*sin(q1-q2)+m*u1^3*r*w*cos(q1-q2)-2*m*u1^2*r*u2*c*sin(q1-q2)-2*m*u1^2*r*u2*w*cos(q1-q2)-M*u1^3*r*l*sin(q1)+2*u1*m*g*l*sin(gam-q1)-u1*m*g*c*sin(gam-q1)+u1*m*g*w*cos(gam-q1)+2*u1*m*g*sin(gam)*r+ud2*m*l*u1*w*sin(q2)+ud1*M*u1*r^2-u1*m*g*c*sin(gam-q1+q2)+u1*m*g*w*cos(gam-q1+q2)+u1*M*g*l*sin(gam-q1)+u1*M*g*sin(gam)*r+2*ud1*m*u1*r*w*sin(q1-q2)+ud1*m*u2*c*cos(q1-q2)*r-ud1*m*u2*w*sin(q1-q2)*r+ud2*m*u1*r*c*cos(q1-q2)-ud2*m*u1*r*w*sin(q1-q2); 
    axh = l*sin(q1)*u1^2+(-l*cos(q1)-r)*ud1; 
    ayh = -l*cos(q1)*u1^2-l*sin(q1)*ud1; 

    zdot = [u1 ud1 u2 ud2 DTE vxh axh vyh ayh]';  
end

function [gstop, isterminal, direction] = collision(t, z, walker)
    q1 = z(1); q2 = z(3); 
    gstop = -q2 + 2*q1;
    if (q2 > -0.05)
        isterminal = 0;
    else
        isterminal = 1;
    end
    direction = -1;
end

function zplus = heelstrike(t, z, walker)      
    r1 = z(1); v1 = z(2); r2 = z(3); v2 = z(4);                         
    xh = z(6); yh = z(8);                       

    q1 = r1 - r2; q2 = -r2;                                       

    M = walker.M; m = walker.m; I = walker.I;   
    l = walker.l; c = walker.c; w = walker.w;   
    r = walker.r; g = walker.g; gam = walker.gam; 

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