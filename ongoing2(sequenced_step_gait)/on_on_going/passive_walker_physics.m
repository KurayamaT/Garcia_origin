% PASSIVE_WALKER_PHYSICS - パッシブウォーカーの物理計算関数群
% このファイルには、single_stride_passivewalker2.mと共通の物理計算関数を含みます

%===================================================================
function [z_traj, t_traj, z_final, z_midpoint] = one_and_half_stride_detailed(z0, walker)
%===================================================================
% 1.5ストライド（スイング脚がスタンス脚に追いつくまで）の詳細な軌道を返す関数

M = walker.M;  m = walker.m; I = walker.I;   
l = walker.l;  c = walker.c; w = walker.w;   
r = walker.r;  g = walker.g; gam = walker.gam;

q1 = z0(1); u1 = z0(2); q2 = z0(3); u2 = z0(4);

% エネルギーと位置の計算
TE = calculate_energy(z0, walker);
xp1 = 0;
xh = -l*sin(q1) - r*q1 + xp1;
vxh = (-l*cos(q1)-r)*u1; 
yh =  l*cos(q1) + r;
vyh = -l*sin(q1)*u1; 

z0_extended = [q1 u1 q2 u2 TE xh vxh yh vyh];

% === 1ステップ目 ===
fprintf('1ステップ目（最初の脚）を実行中...\n');
t0 = 0; 
dt = 5; % 1ステップの最大時間
time_stamps = 100; % 軌道の解像度

% Single stance phase
options = odeset('abstol',1e-13,'reltol',1e-13,'events',@collision);
tspan = linspace(t0, t0+dt, time_stamps);
[t_step1, z_step1] = ode113(@single_stance, tspan, z0_extended, options, walker);

% 1ステップ目のHeel strike
if length(t_step1) < time_stamps
    fprintf('1ステップ目完了: t = %.4f秒で次の脚が地面に接触\n', t_step1(end));
    z_after_collision1 = heelstrike(t_step1(end), z_step1(end,:), walker);
    z_midpoint = z_after_collision1(1:4);
else
    error('1ステップ目: 指定時間内に衝突が発生しませんでした');
end

% === 0.5ステップ目（スイング脚がスタンス脚に追いつくまで） ===
fprintf('0.5ステップ目（スイング脚がスタンス脚に追いつくまで）を実行中...\n');
z1_extended = z_after_collision1;
t1 = t_step1(end);

% 新しい衝突条件：スイング脚がスタンス脚に追いつく（地面接触なし）
options_catchup = odeset('abstol',1e-13,'reltol',1e-13,'events',@swing_catches_stance);
tspan = linspace(t1, t1+dt, time_stamps);
[t_step1_5, z_step1_5] = ode113(@single_stance, tspan, z1_extended, options_catchup, walker);

% 時間を調整（連続的にする）
t_step1_5 = t_step1_5 - t1 + t_step1(end);

% 0.5ステップ目の終了（スイング脚追いつき、Heel Strikeなし）
if length(t_step1_5) < time_stamps
    fprintf('0.5ステップ目完了: t = %.4f秒でスイング脚がスタンス脚に追いつきました\n', t_step1_5(end));
    fprintf('最終角度: q1=%.4f, q2=%.4f, 差=%.6f\n', z_step1_5(end,1), z_step1_5(end,3), z_step1_5(end,3)-z_step1_5(end,1));
    z_final = z_step1_5(end, 1:4);
else
    fprintf('警告: 指定時間内にスイング脚が追いつきませんでした\n');
    fprintf('最終角度: q1=%.4f, q2=%.4f, 差=%.6f\n', z_step1_5(end,1), z_step1_5(end,3), z_step1_5(end,3)-z_step1_5(end,1));
    z_final = z_step1_5(end, 1:4);
end

% === 軌道データの結合 ===
% 重複する時間点を除去して連結
t_traj = [t_step1; t_step1_5(2:end)];
z_traj = [z_step1; z_step1_5(2:end,:)];

fprintf('1.5ストライド（スイング脚追いつき）完了\n');

end

%===================================================================
function TE = calculate_energy(z, walker)
%===================================================================
% 総エネルギーを計算

q1 = z(1); u1 = z(2); q2 = z(3); u2 = z(4);
M = walker.M;  m = walker.m; I = walker.I;   
l = walker.l;  c = walker.c; w = walker.w;   
r = walker.r;  g = walker.g; gam = walker.gam;

TE = 1/2*m*(((-l*cos(q1)-r)*u1-u1*(-c*cos(q1)+w*sin(q1)))^2+(-l*sin(q1)*u1+u1*(c*sin(q1)+w*cos(q1)))^2)+1/2*m*(((-l*cos(q1)-r)*u1-(u1-u2)*(-c*cos(q1-q2)+w*sin(q1-q2)))^2+(-l*sin(q1)*u1+(u1-u2)*(c*sin(q1-q2)+w*cos(q1-q2)))^2)+1/2*M*((-l*cos(q1)-r)^2*u1^2+l^2*sin(q1)^2*u1^2)+1/2*I*(u1^2+(u1-u2)^2)+2*m*g*cos(gam)*r+2*m*g*l*cos(gam-q1)-m*g*c*cos(gam-q1)-m*g*w*sin(gam-q1)+2*m*g*sin(gam)*r*q1-m*g*c*cos(gam-q1+q2)-m*g*w*sin(gam-q1+q2)+M*g*cos(gam)*r+M*g*l*cos(gam-q1)+M*g*sin(gam)*r*q1;

end

%===================================================================
function [gstop, isterminal,direction]=swing_catches_stance(t,z,walker)
%===================================================================
% スイング脚がスタンス脚に追いつく条件

q1 = z(1); q2 = z(3); 

% スイング脚角度がスタンス脚角度に追いつく条件
gstop = q2 - q1;
isterminal = 1; % 条件に達したら積分を停止
direction = 0;  % どちらの方向からでも検出（0 = 両方向）

% デバッグ用
fprintf('Debug: t=%.4f, q1=%.4f, q2=%.4f, diff=%.6f\n', t, q1, q2, gstop);

end

%===================================================================
function zdot=single_stance(t,z,walker)  
%===================================================================

q1 = z(1);   u1 = z(2);                         
q2 = z(3);   u2 = z(4);                         
xh = z(6);  vxh = z(7);                       
yh = z(8);  vyh = z(9);                     

M = walker.M;  m = walker.m; I = walker.I;   
l = walker.l;  c = walker.c; w = walker.w;   
r = walker.r;  g = walker.g; gam = walker.gam;

Th=0;

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

%===================================================================
function [gstop, isterminal,direction]=collision(t,z,walker)
%===================================================================

q1 = z(1); q2 = z(3); 

gstop = -q2 + 2*q1;
if (q2>-0.05)
    isterminal = 0;
else
    isterminal=1;
end
direction=-1;

end

%===================================================================
function zplus=heelstrike(t,z,walker)      
%===================================================================

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

u1 = X(1); u2 = X(2);                                      

TE = 1/2*m*(((-l*cos(q1)-r)*u1-u1*(-c*cos(q1)+w*sin(q1)))^2+(-l*sin(q1)*u1+u1*(c*sin(q1)+w*cos(q1)))^2)+1/2*m*(((-l*cos(q1)-r)*u1-(u1-u2)*(-c*cos(q1-q2)+w*sin(q1-q2)))^2+(-l*sin(q1)*u1+(u1-u2)*(c*sin(q1-q2)+w*cos(q1-q2)))^2)+1/2*M*((-l*cos(q1)-r)^2*u1^2+l^2*sin(q1)^2*u1^2)+1/2*I*(u1^2+(u1-u2)^2)+2*g*m*cos(gam)*r+2*m*g*l*cos(gam-q1)-m*g*c*cos(gam-q1)-m*g*w*sin(gam-q1)+2*g*m*sin(gam)*r*q1-m*g*c*cos(gam-q1+q2)-m*g*w*sin(gam-q1+q2)+g*M*cos(gam)*r+M*g*l*cos(gam-q1)+g*M*sin(gam)*r*q1; 
vxh = (-l*cos(q1)-r)*u1; 
vyh = -l*sin(q1)*u1; 

zplus = [q1 u1 q2 u2 TE xh vxh yh vyh];                     

end