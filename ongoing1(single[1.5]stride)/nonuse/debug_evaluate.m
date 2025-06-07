function debug_evaluate()
% 詳細デバッグ版の評価関数テスト

    fprintf('=== 詳細デバッグ評価 ===\n');
    
    % 元の論文の条件
    z0 = [0.2, -0.2, 0.4, -0.3];
    
    fprintf('テスト条件: [%.1f, %.1f, %.1f, %.1f]\n', z0);
    
    % Walker設定
    walker.M = 1000; walker.m = 1.0; walker.I = 0.00; walker.l = 1.0; walker.w = 0.0; 
    walker.c = 1.0;  walker.r = 0.3; walker.g = 1.0; walker.gam = 0.009; 
    
    % ステップ1: 固定点探索
    fprintf('\n--- ステップ1: 固定点探索 ---\n');
    options = optimset('TolFun',1e-12,'TolX',1e-12,'Display','off');
    [zstar, fval, exitflag] = fsolve(@(z) fixedpt(z, walker), z0, options);
    
    fprintf('exitflag: %d\n', exitflag);
    fprintf('残差ノルム: %.2e\n', norm(fval));
    fprintf('固定点: [%.6f, %.6f, %.6f, %.6f]\n', zstar);
    
    if exitflag ~= 1
        fprintf('❌ 固定点探索失敗\n');
        return;
    end
    
    if norm(fval) > 1e-6
        fprintf('❌ 固定点精度不足\n');
        return;
    end
    
    fprintf('✓ 固定点探索成功\n');
    
    % ステップ2: ヤコビアン計算
    fprintf('\n--- ステップ2: ヤコビアン計算 ---\n');
    
    try
        % onestep関数のテスト
        fprintf('onestep関数テスト...\n');
        z_test = onestep(zstar, walker);
        fprintf('onestep結果: [%.6f, %.6f, %.6f, %.6f]\n', z_test);
        
        % 関数ハンドルのテスト
        fprintf('関数ハンドルテスト...\n');
        fun_handle = @(z) onestep(z, walker);
        z_test2 = fun_handle(zstar);
        fprintf('関数ハンドル結果: [%.6f, %.6f, %.6f, %.6f]\n', z_test2);
        
        % partialder関数の詳細テスト
        fprintf('ヤコビアン計算中...\n');
        J = partialder_debug(@(z) onestep(z, walker), zstar, walker);
        
        fprintf('ヤコビアン サイズ: %d × %d\n', size(J));
        fprintf('ヤコビアン 範囲: %.2e ~ %.2e\n', min(J(:)), max(J(:)));
        
        % NaNやInfのチェック
        if any(any(~isfinite(J)))
            fprintf('❌ ヤコビアンに無限大またはNaNが含まれています\n');
            nan_count = sum(sum(isnan(J)));
            inf_count = sum(sum(isinf(J)));
            fprintf('NaN数: %d, Inf数: %d\n', nan_count, inf_count);
            
            % 問題のある要素を表示
            [nan_i, nan_j] = find(isnan(J));
            [inf_i, inf_j] = find(isinf(J));
            if ~isempty(nan_i)
                fprintf('NaN位置: '); 
                for k = 1:min(5, length(nan_i))
                    fprintf('(%d,%d) ', nan_i(k), nan_j(k));
                end
                fprintf('\n');
            end
            if ~isempty(inf_i)
                fprintf('Inf位置: '); 
                for k = 1:min(5, length(inf_i))
                    fprintf('(%d,%d) ', inf_i(k), inf_j(k));
                end
                fprintf('\n');
            end
            return;
        end
        
        fprintf('✓ ヤコビアン計算成功\n');
        
        % ステップ3: 固有値計算
        fprintf('\n--- ステップ3: 固有値計算 ---\n');
        eigenvalues = eig(J);
        
        fprintf('固有値数: %d\n', length(eigenvalues));
        fprintf('固有値: ');
        for i = 1:length(eigenvalues)
            if imag(eigenvalues(i)) == 0
                fprintf('%.4f ', real(eigenvalues(i)));
            else
                fprintf('%.4f%+.4fi ', real(eigenvalues(i)), imag(eigenvalues(i)));
            end
        end
        fprintf('\n');
        
        abs_eigenvalues = abs(eigenvalues);
        fprintf('固有値絶対値: ');
        for i = 1:length(abs_eigenvalues)
            fprintf('%.4f ', abs_eigenvalues(i));
        end
        fprintf('\n');
        
        max_eigenvalue = max(abs_eigenvalues);
        fprintf('最大固有値: %.6f\n', max_eigenvalue);
        
        if any(~isfinite(eigenvalues))
            fprintf('❌ 固有値に無限大またはNaNが含まれています\n');
            return;
        end
        
        if max_eigenvalue < 1
            fprintf('✓ 安定 (最大固有値 < 1)\n');
            fprintf('🎉 評価成功！この条件は安定です\n');
        else
            fprintf('❌ 不安定 (最大固有値 >= 1)\n');
        end
        
    catch ME
        fprintf('❌ エラー発生: %s\n', ME.message);
        fprintf('スタック:\n');
        for i = 1:length(ME.stack)
            fprintf('  %s (行 %d)\n', ME.stack(i).name, ME.stack(i).line);
        end
    end
end

function J = partialder_debug(FUN, z, walker)
% デバッグ版のpartialder関数

    pert = 1e-5;
    n = length(z);
    J = zeros(n, n);
    
    fprintf('  摂動による数値微分開始 (摂動: %.1e)\n', pert);
    
    for i = 1:n
        if mod(i, 2) == 0
            fprintf('  列 %d/%d 計算中...\n', i, n);
        end
        
        ztemp1 = z; ztemp2 = z;
        ztemp1(i) = ztemp1(i) + pert;
        ztemp2(i) = ztemp2(i) - pert;
        
        try
            f1 = FUN(ztemp1);
            f2 = FUN(ztemp2);
            
            % 結果のチェック
            if any(~isfinite(f1))
                fprintf('    警告: f1に無限大/NaNが含まれています (i=%d)\n', i);
            end
            if any(~isfinite(f2))
                fprintf('    警告: f2に無限大/NaNが含まれています (i=%d)\n', i);
            end
            
            J(:,i) = (f1 - f2);
            
        catch ME
            fprintf('    エラー: 列%d計算中にエラー発生: %s\n', i, ME.message);
            J(:,i) = NaN;
        end
    end
    
    J = J / (2*pert);
    fprintf('  数値微分完了\n');
end

% 必要な関数
function zdiff = fixedpt(z0, walker)
    zdiff = onestep(z0, walker) - z0; 
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
        options = odeset('abstol',1e-13,'reltol',1e-13,'events',@(t,z) collision(t,z));
        tspan = linspace(t0,t0+dt,time_stamps);
        [t_temp, z_temp] = ode113(@(t,z) single_stance(t,z,walker), tspan, z0, options);
        
        zplus = heelstrike(t_temp(end),z_temp(end,:),walker); 
        z0 = zplus; t0 = t_temp(end);
        
        t_ode = [t_ode; t_temp(2:end)];
        z_ode = [z_ode; z_temp(2:end,:)];
    end

    z = zplus(1:4);
    if flag == 1
       z = z_ode; t = t_ode;
    end
end

function zdot = single_stance(t,z,walker)  
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

function [gstop, isterminal,direction] = collision(t,z)
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