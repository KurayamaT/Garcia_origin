% SINGLE_STRIDE_WALKER - 1.5ストライド（スイング脚がスタンス脚に追いつくまで）を実行するパッシブウォーカー
% 使い方:
% single_stride_walker(1, [0.2, -0.25, 0.35, -0.4])  % 指定初期値から1.5歩
% single_stride_walker(1)                              % デフォルト初期値から1.5歩

function single_stride_walker(flag, z0_input)

clc
close all
format long

% 引数の処理
if nargin == 0
    flag = 1; % デフォルトは最も単純なウォーカー
    z0_input = []; % 初期値は後で入力
elseif nargin == 1
    z0_input = []; % 初期値は後で入力
end

% インタラクティブな初期値入力
if isempty(z0_input)
    fprintf('=== パッシブウォーカー 1.5ストライド（スイング脚追いつきまで）シミュレーション ===\n\n');
    
    if flag == 1
        fprintf('ウォーカータイプ: Garcia''s Simplest Walker\n');
        fprintf('推奨固定点: [0.2002, -0.1999, 0.4003, -0.0158]\n');
    else
        fprintf('ウォーカータイプ: General Round Feet Walker\n');
        fprintf('推奨固定点: [0.1895, -0.2391, 0.3789, -0.0537]\n');
    end
    
    fprintf('\n初期値の入力方法を選択してください:\n');
    fprintf('1. デフォルト初期値を使用\n');
    fprintf('2. 固定点から開始（安定歩行）\n');
    fprintf('3. 手動で初期値を入力\n');
    fprintf('4. プリセット初期値を選択\n');
    
    choice = input('選択 (1-4): ');
    
    switch choice
        case 1
            % デフォルト初期値
            if flag == 1
                z0_input = [0.2, -0.2, 0.4, -0.3];
            else
                z0_input = [0.2, -0.3, 0.4, -0.3];
            end
            fprintf('デフォルト初期値を使用: [%.3f, %.3f, %.3f, %.3f]\n', z0_input);
            
        case 2
            % 固定点から開始
            if flag == 1
                z0_input = [0.200161072169750, -0.199906060087682, 0.400322144339512, -0.015805473227965];
            else
                z0_input = [0.189472782205104, -0.239124222551699, 0.378945564410209, -0.053691703909393];
            end
            fprintf('固定点から開始: [%.6f, %.6f, %.6f, %.6f]\n', z0_input);
            
        case 3
            % 手動入力
            fprintf('\n初期状態を入力してください:\n');
            fprintf('q1: スタンス脚角度 (rad) - 推奨範囲: [-0.5, 0.5]\n');
            fprintf('u1: スタンス脚角速度 (rad/s) - 推奨範囲: [-1.0, 1.0]\n');
            fprintf('q2: スイング脚角度 (rad) - 推奨範囲: [-0.5, 0.8]\n');
            fprintf('u2: スイング脚角速度 (rad/s) - 推奨範囲: [-1.0, 1.0]\n\n');
            
            q1 = input('q1 (スタンス脚角度): ');
            u1 = input('u1 (スタンス脚角速度): ');
            q2 = input('q2 (スイング脚角度): ');
            u2 = input('u2 (スイング脚角速度): ');
            
            z0_input = [q1, u1, q2, u2];
            fprintf('入力された初期値: [%.6f, %.6f, %.6f, %.6f]\n', z0_input);
            
        case 4
            % プリセット初期値
            fprintf('\nプリセット初期値を選択してください:\n');
            if flag == 1
                presets = {
                    [0.1, -0.1, 0.2, -0.2], '安定な初期値1';
                    [0.3, -0.3, 0.5, -0.1], '安定な初期値2';
                    [0.0, -0.15, 0.0, -0.35], '揃え足歩行1';
                    [-.01, -0.1, -.1, -0.3], '揃え足歩行2';
                    [0.4, -0.1, 0.6, -0.4], '不安定な初期値1';
                    [0.1, -0.5, 0.3, -0.6], '不安定な初期値2'
                };
            else
                presets = {
                    [0.15, -0.2, 0.3, -0.3], '安定な初期値1';
                    [0.25, -0.3, 0.4, -0.2], '安定な初期値2';
                    [0.18, -0.24, 0.36, -0.06], '固定点近傍1';
                    [0.21, -0.22, 0.40, -0.04], '固定点近傍2';
                    [0.35, -0.1, 0.5, -0.5], '不安定な初期値1';
                    [0.1, -0.4, 0.25, -0.7], '不安定な初期値2'
                };
            end
            
            for i = 1:size(presets, 1)
                fprintf('%d. %s: [%.3f, %.3f, %.3f, %.3f]\n', i, presets{i,2}, presets{i,1});
            end
            
            preset_choice = input(sprintf('プリセット選択 (1-%d): ', size(presets, 1)));
            
            if preset_choice >= 1 && preset_choice <= size(presets, 1)
                z0_input = presets{preset_choice, 1};
                fprintf('選択されたプリセット: %s\n', presets{preset_choice, 2});
                fprintf('初期値: [%.6f, %.6f, %.6f, %.6f]\n', z0_input);
            else
                fprintf('無効な選択です。デフォルト初期値を使用します。\n');
                if flag == 1
                    z0_input = [0.2, -0.2, 0.4, -0.3];
                else
                    z0_input = [0.2, -0.3, 0.4, -0.3];
                end
            end
            
        otherwise
            fprintf('無効な選択です。デフォルト初期値を使用します。\n');
            if flag == 1
                z0_input = [0.2, -0.2, 0.4, -0.3];
            else
                z0_input = [0.2, -0.3, 0.4, -0.3];
            end
    end
    
    % 初期値の妥当性チェック
    if any(abs(z0_input(1:2:3)) > 1.0)  % 角度チェック
        fprintf('\n⚠️  警告: 角度が大きすぎます（|q1|, |q2| > 1.0 rad）\n');
        fprintf('   シミュレーションが不安定になる可能性があります。\n');
    end
    
    if any(abs(z0_input(2:2:4)) > 2.0)  % 角速度チェック
        fprintf('\n⚠️  警告: 角速度が大きすぎます（|u1|, |u2| > 2.0 rad/s）\n');
        fprintf('   シミュレーションが不安定になる可能性があります。\n');
    end
    
    fprintf('\n続行しますか？ (y/n): ');
    continue_sim = input('', 's');
    if ~strcmpi(continue_sim, 'y') && ~strcmpi(continue_sim, 'yes') && ~isempty(continue_sim)
        fprintf('シミュレーションを中止しました。\n');
        return;
    end
    
    fprintf('\n');
end

% アニメーション保存オプション（初期値設定後、シミュレーション前）
fprintf('\nアニメーションを保存しますか？\n');
fprintf('1. 保存しない（表示のみ）\n');
fprintf('2. GIFファイルとして保存\n');
fprintf('3. AVIファイルとして保存\n');
fprintf('4. 両方の形式で保存\n');
save_choice = input('選択 (1-4): ');

% ファイル名の設定
if save_choice > 1
    default_filename = sprintf('passive_walker_%s', datestr(now, 'yyyymmdd_HHMMSS'));
    fprintf('\n保存ファイル名を入力してください（拡張子なし）\n');
    fprintf('デフォルト: %s\n', default_filename);
    filename_input = input('ファイル名: ', 's');
    if isempty(filename_input)
        save_filename = default_filename;
    else
        save_filename = filename_input;
    end
else
    save_filename = '';
end

if flag == 1
    %% Garcia's simplest walker
    walker.M = 1000; walker.m = 1.0; walker.I = 0.00; walker.l = 1.0; walker.w = 0.0; 
    walker.c = 1.0;  walker.r = 0.0; walker.g = 1.0; walker.gam = 0.009; 
    
    %%%% Initial State %%%%%
    % ユーザー指定の初期値
    if length(z0_input) ~= 4
        error('初期値は4要素のベクトル [q1, u1, q2, u2] である必要があります');
    end
    z0 = z0_input;
    
    % 固定点（参考）
    zstar = [0.200161072169750  -0.199906060087682   0.400322144339512  -0.015805473227965];
    
else 
    %%  More General round feet walker
    walker.M = 1.0; walker.m = 0.5; walker.I = 0.02; walker.l = 1.0; walker.w = 0.0; 
    walker.c = 0.5; walker.r = 0.2; walker.g = 1.0; walker.gam = 0.01; 
    
    %%%% Initial State %%%%%
    % ユーザー指定の初期値
    if length(z0_input) ~= 4
        error('初期値は4要素のベクトル [q1, u1, q2, u2] である必要があります');
    end
    z0 = z0_input;
    
    % 固定点（参考）
    zstar = [0.189472782205104  -0.239124222551699   0.378945564410209  -0.053691703909393];
end

%=========================================================================
% 1.5ストライド（スイング脚がスタンス脚に追いつくまで）のシミュレーション
fps = 10;

fprintf('\n=== 1.5ストライド（スイング脚追いつきまで）シミュレーション開始 ===\n');

%%%% 1.5ストライドの実行 %%%%
[z_trajectory, t_trajectory, z_final, z_midpoint] = one_and_half_stride_detailed(z0, walker);

fprintf('\n=== 1.5ストライド（スイング脚追いつき）完了 ===\n');
fprintf('初期状態      : [%.6f, %.6f, %.6f, %.6f]\n', z0);
fprintf('1ステップ後   : [%.6f, %.6f, %.6f, %.6f]\n', z_midpoint);
fprintf('最終状態(1.5歩): [%.6f, %.6f, %.6f, %.6f]\n', z_final);

% 固定点との比較
fprintf('\n=== 固定点との比較 ===\n');
fprintf('固定点     : [%.6f, %.6f, %.6f, %.6f]\n', zstar);
diff_from_fixed = z_final - zstar;
fprintf('最終差分   : [%.6f, %.6f, %.6f, %.6f]\n', diff_from_fixed);
fprintf('最終距離   : %.6f\n', norm(diff_from_fixed));

% 周期性チェック（1.5ストライド後の状態確認）
fprintf('\n=== 歩行状態チェック ===\n');
fprintf('スタンス脚角度 (q1): %.6f rad\n', z_final(1));
fprintf('スイング脚角度 (q2): %.6f rad\n', z_final(3));
fprintf('脚間角度差: %.6f rad\n', abs(z_final(1) - z_final(3)));

if abs(z_final(1) - z_final(3)) < 0.05
    fprintf('✅ スイング脚がスタンス脚に正常に追いつきました\n');
elseif abs(z_final(1) - z_final(3)) < 0.2
    fprintf('⚠️ スイング脚がスタンス脚にほぼ追いついています\n');
else
    fprintf('❌ スイング脚がスタンス脚に追いついていません\n');
end

%%%% エネルギー解析 %%%%
fprintf('\n=== エネルギー解析 ===\n');
TE_initial = calculate_energy(z0, walker);
TE_midpoint = calculate_energy(z_midpoint, walker);
TE_final = calculate_energy(z_final, walker);
fprintf('初期エネルギー: %.6f\n', TE_initial);
fprintf('中間エネルギー: %.6f\n', TE_midpoint);
fprintf('最終エネルギー: %.6f\n', TE_final);
fprintf('1ステップ目の損失: %.6f\n', TE_initial - TE_midpoint);
fprintf('0.5ステップ目の損失: %.6f\n', TE_midpoint - TE_final);
fprintf('総エネルギー損失: %.6f\n', TE_initial - TE_final);

%%%% 歩行成功判定 %%%%
fprintf('\n=== 歩行成功判定 ===\n');
if abs(z_final(1)) < 1.0 && abs(z_final(3)) < 1.0  % 角度が妥当な範囲内
    fprintf('✅ 1.5ストライド（スイング脚追いつき）が正常に完了しました\n');
else
    fprintf('❌ 歩行が不安定です（角度が大きすぎます）\n');
end

%%%% アニメーション（保存機能付き） %%%%
fprintf('\n1.5ストライド（スイング脚追いつきまで）のアニメーションを表示します...\n');
animate_single_stride(t_trajectory, z_trajectory, walker, save_choice, save_filename);

%%%% プロット %%%%
fprintf('\n軌道をプロットします...\n');
plot_single_stride(t_trajectory, z_trajectory, z0, z_final, zstar, z_midpoint);

end  % メイン関数の終了

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FUNCTIONS START HERE %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    % Heel Strikeは発生させない（スイング脚が追いついただけ）
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
function plot_single_stride(t, z, z0, z_final, zstar, z_midpoint)
%===================================================================

% 1ステップ目と0.5ステップ目の境界を見つける
% 時間が大きく跳ぶ点や、x座標のヒップ位置の不連続点を探す
step_boundary = length(t) / 2; % 近似的に中点を境界とする
for i = 2:length(t)-1
    if i < length(t) && (t(i+1) - t(i) > mean(diff(t))*2 || abs(z(i+1,6) - z(i,6)) > 0.1)
        step_boundary = i;
        break;
    end
end

figure('Name', '1.5ストライド（スイング脚追いつき）の軌道', 'Position', [100, 100, 1400, 900]);

% 状態変数の時間変化
subplot(3,3,1)
plot(t(1:step_boundary), z(1:step_boundary,1), 'r-', 'LineWidth', 2)
hold on
plot(t(step_boundary+1:end), z(step_boundary+1:end,1), 'r--', 'LineWidth', 2)
plot([t(1), t(end)], [z0(1), z0(1)], 'g:', 'LineWidth', 1)
plot([t(1), t(end)], [z_final(1), z_final(1)], 'b:', 'LineWidth', 1)
plot([t(1), t(end)], [zstar(1), zstar(1)], 'k--', 'LineWidth', 1)
xlabel('時間 (s)')
ylabel('q1 (rad)')
title('スタンス脚角度')
legend('1ステップ目', '0.5ステップ目', '初期値', '最終値', '固定点', 'Location', 'best')
grid on

subplot(3,3,2)
plot(t(1:step_boundary), z(1:step_boundary,2), 'b-', 'LineWidth', 2)
hold on
plot(t(step_boundary+1:end), z(step_boundary+1:end,2), 'b--', 'LineWidth', 2)
plot([t(1), t(end)], [z0(2), z0(2)], 'g:', 'LineWidth', 1)
plot([t(1), t(end)], [z_final(2), z_final(2)], 'b:', 'LineWidth', 1)
plot([t(1), t(end)], [zstar(2), zstar(2)], 'k--', 'LineWidth', 1)
xlabel('時間 (s)')
ylabel('u1 (rad/s)')
title('スタンス脚角速度')
legend('1ステップ目', '0.5ステップ目', '初期値', '最終値', '固定点', 'Location', 'best')
grid on

subplot(3,3,3)
plot(t(1:step_boundary), z(1:step_boundary,3), 'g-', 'LineWidth', 2)
hold on
plot(t(step_boundary+1:end), z(step_boundary+1:end,3), 'g--', 'LineWidth', 2)
plot([t(1), t(end)], [z0(3), z0(3)], 'g:', 'LineWidth', 1)
plot([t(1), t(end)], [z_final(3), z_final(3)], 'b:', 'LineWidth', 1)
plot([t(1), t(end)], [zstar(3), zstar(3)], 'k--', 'LineWidth', 1)
xlabel('時間 (s)')
ylabel('q2 (rad)')
title('スイング脚角度')
legend('1ステップ目', '0.5ステップ目', '初期値', '最終値', '固定点', 'Location', 'best')
grid on

subplot(3,3,4)
plot(t(1:step_boundary), z(1:step_boundary,4), 'm-', 'LineWidth', 2)
hold on
plot(t(step_boundary+1:end), z(step_boundary+1:end,4), 'm--', 'LineWidth', 2)
plot([t(1), t(end)], [z0(4), z0(4)], 'g:', 'LineWidth', 1)
plot([t(1), t(end)], [z_final(4), z_final(4)], 'b:', 'LineWidth', 1)
plot([t(1), t(end)], [zstar(4), zstar(4)], 'k--', 'LineWidth', 1)
xlabel('時間 (s)')
ylabel('u2 (rad/s)')
title('スイング脚角速度')
legend('1ステップ目', '0.5ステップ目', '初期値', '最終値', '固定点', 'Location', 'best')
grid on

% エネルギーの変化
subplot(3,3,5)
plot(t(1:step_boundary), z(1:step_boundary,5), 'k-', 'LineWidth', 2)
hold on
plot(t(step_boundary+1:end), z(step_boundary+1:end,5), 'k--', 'LineWidth', 2)
xlabel('時間 (s)')
ylabel('総エネルギー')
title('エネルギー変化')
legend('1ステップ目', '0.5ステップ目', 'Location', 'best')
grid on

% 位相図 (q1-u1)
subplot(3,3,6)
plot(z(1:step_boundary,1), z(1:step_boundary,2), 'r-', 'LineWidth', 2)
hold on
plot(z(step_boundary+1:end,1), z(step_boundary+1:end,2), 'r--', 'LineWidth', 2)
plot(z0(1), z0(2), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g')
plot(z_midpoint(1), z_midpoint(2), 'yo', 'MarkerSize', 8, 'MarkerFaceColor', 'y')
plot(z_final(1), z_final(2), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b')
plot(zstar(1), zstar(2), 'k*', 'MarkerSize', 12, 'MarkerFaceColor', 'k')
xlabel('q1 (rad)')
ylabel('u1 (rad/s)')
title('位相図 (q1-u1)')
legend('1ステップ目', '0.5ステップ目', '初期点', '中間点', '最終点', '固定点', 'Location', 'best')
grid on

% 位相図 (q2-u2)
subplot(3,3,7)
plot(z(1:step_boundary,3), z(1:step_boundary,4), 'g-', 'LineWidth', 2)
hold on
plot(z(step_boundary+1:end,3), z(step_boundary+1:end,4), 'g--', 'LineWidth', 2)
plot(z0(3), z0(4), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g')
plot(z_midpoint(3), z_midpoint(4), 'yo', 'MarkerSize', 8, 'MarkerFaceColor', 'y')
plot(z_final(3), z_final(4), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b')
plot(zstar(3), zstar(4), 'k*', 'MarkerSize', 12, 'MarkerFaceColor', 'k')
xlabel('q2 (rad)')
ylabel('u2 (rad/s)')
title('位相図 (q2-u2)')
legend('1ステップ目', '0.5ステップ目', '初期点', '中間点', '最終点', '固定点', 'Location', 'best')
grid on

% ヒップ軌道
subplot(3,3,8)
plot(z(1:step_boundary,6), z(1:step_boundary,8), 'c-', 'LineWidth', 2)
hold on
plot(z(step_boundary+1:end,6), z(step_boundary+1:end,8), 'c--', 'LineWidth', 2)
plot(z(1,6), z(1,8), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g')
plot(z(step_boundary,6), z(step_boundary,8), 'yo', 'MarkerSize', 8, 'MarkerFaceColor', 'y')
plot(z(end,6), z(end,8), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b')
xlabel('x (m)')
ylabel('y (m)')
title('ヒップ軌道')
legend('1ステップ目', '0.5ステップ目', '開始', '中間', '終了', 'Location', 'best')
grid on
axis equal

% 状態空間での軌道比較
subplot(3,3,9)
% 全状態の正規化距離を計算
state_norm = sqrt(sum(z(:,1:4).^2, 2));
plot(t(1:step_boundary), state_norm(1:step_boundary), 'k-', 'LineWidth', 2)
hold on
plot(t(step_boundary+1:end), state_norm(step_boundary+1:end), 'k--', 'LineWidth', 2)
plot([t(1), t(end)], [norm(z0), norm(z0)], 'g:', 'LineWidth', 1)
plot([t(1), t(end)], [norm(z_final), norm(z_final)], 'b:', 'LineWidth', 1)
plot([t(1), t(end)], [norm(zstar), norm(zstar)], 'k:', 'LineWidth', 1)
xlabel('時間 (s)')
ylabel('状態ベクトルのノルム')
title('状態空間での軌道')
legend('1ステップ目', '0.5ステップ目', '初期ノルム', '最終ノルム', '固定点ノルム', 'Location', 'best')
grid on

sgtitle('1.5ストライド（スイング脚追いつき）の詳細解析')

end

%===================================================================
function animate_single_stride(t_all, z_all, walker, save_choice, save_filename)
%===================================================================
% アニメーション関数（保存機能付き）

% 引数チェック（後方互換性のため）
if nargin < 4
    save_choice = 1;  % デフォルトは保存しない
    save_filename = '';
end

fps = 10;
z_all_plot = [z_all(:,1) z_all(:,3) z_all(:,6) z_all(:,8)];
nn = size(z_all_plot,2);
total_frames = round(t_all(end)*fps);
t = linspace(0,t_all(end),total_frames);
z = zeros(total_frames,nn);
for i=1:nn
    z(:,i) = interp1(t_all,z_all_plot(:,i),t);
end

figure('Name', '1.5ストライド（スイング脚追いつき）アニメーション', 'Position', [200, 200, 800, 600]);
clf

M = walker.M;  m = walker.m; I = walker.I;   
l = walker.l;  c = walker.c; w = walker.w;   
r = walker.r;  g = walker.g; gam = walker.gam; 

mm = size(z,1);

min_xh = min(z(:,3)); max_xh = max(z(:,3)); 
window_xmin = min_xh - 1*l; window_xmax = max_xh + 1*l;
window_ymin = -0.1; window_ymax = 1.1*(l+r);

axis('equal')
axis([window_xmin window_xmax window_ymin window_ymax])
axis off
set(gcf,'Color',[1,1,1])

% 軌跡を描画するための準備
trajectory_x = z(:,3);
trajectory_y = z(:,4);

lines_for_feet = 4;
counter = 2 + 2*lines_for_feet;
th = 0.25;

%%% create object for hinge %%%%%
hingepic=line('xdata',0,'ydata',0, 'marker','.','markersize',[20], ...
          'color','black');

%%% create trajectory line %%%
trajpic=line('xdata',[],'ydata',[], 'linewidth', 1, 'color','green', 'linestyle', '--');
   
%%%% create object for legs and feet %%%%
barref = [0 0; 0 -1];
y = [0;-1];
O = [0; 0];

%%%% legs in red %%%      
for p = 1:2
    barpic(p)=line('xdata',barref(1,:),'ydata',barref(2,:),'linewidth', 3, 'color','red');
end

%%% feet in blue %%%
for p = 3:counter
    barpic(p)=line('xdata',barref(1,:),'ydata',barref(2,:),'linewidth', 2, 'color','blue');
end

%%%% create ramp %%%%
rampref=[window_xmin window_xmax ; 0 0];
line('xdata',rampref(1,:),'ydata',rampref(2,:), 'linewidth', 2,'color','black');

title('1.5ストライド（スイング脚追いつき）アニメーション', 'FontSize', 14)

% GIF保存の準備
if save_choice == 2 || save_choice == 4
    gif_filename = [save_filename '.gif'];
    fprintf('GIFファイルを保存中: %s\n', gif_filename);
end

% AVI保存の準備
if save_choice == 3 || save_choice == 4
    avi_filename = [save_filename '.avi'];
    fprintf('AVIファイルを保存中: %s\n', avi_filename);
    v = VideoWriter(avi_filename);
    v.FrameRate = fps;
    open(v);
end

% アニメーションループ
for i=1:mm
    
   q1 = z(i,1); q2 = z(i,2); 
   xh = z(i,3); yh = z(i,4);  

   %%% hinge coordinates
   hinge=[xh; yh];   
   
   %%% leg coordinates
   A = [q1 -(q2-q1)];
   
   for p = 1:2
       bar(:,:,p) = [hinge, hinge] + (l+r)*R(A(p))*barref;
       center(:,:,p) = hinge + l*R(A(p))*y;
   end
   
   %%%% feet coordinates
   B = [-th -2*th; 0 -th; th 0; 2*th th];

   incr = 3;
   for p=1:2
       for q=1:4
              C = A(p) + B(q,:);
              bar(:,:,incr) = [center(:,:,p), center(:,:,p)] + r*R(C(1))*[O,y] + r*R(C(2))*[y, O]; 
              incr = incr + 1;
       end
   end
          
   %%% animate now    
   set(hingepic,'xdata',hinge(1),'ydata',hinge(2));
   set(trajpic,'xdata',trajectory_x(1:i),'ydata',trajectory_y(1:i));

   for p=1:counter
       set(barpic(p),'xdata',bar(1,:,p),'ydata',bar(2,:,p));
   end
   
   pause(0.1)
   drawnow
   
   % フレームを保存
   if save_choice > 1
       frame = getframe(gcf);
       
       % GIF保存
       if save_choice == 2 || save_choice == 4
           im = frame2im(frame);
           [imind,cm] = rgb2ind(im,256);
           if i == 1
               imwrite(imind,cm,gif_filename,'gif','Loopcount',inf,'DelayTime',1/fps);
           else
               imwrite(imind,cm,gif_filename,'gif','WriteMode','append','DelayTime',1/fps);
           end
       end
       
       % AVI保存
       if save_choice == 3 || save_choice == 4
           writeVideo(v,frame);
       end
   end
  
end

% AVIファイルを閉じる
if save_choice == 3 || save_choice == 4
    close(v);
    fprintf('AVIファイル保存完了: %s\n', avi_filename);
end

% GIF保存完了メッセージ
if save_choice == 2 || save_choice == 4
    fprintf('GIFファイル保存完了: %s\n', gif_filename);
end

end

%===================================================================
function rotation = R(A)
%===================================================================
rotation = [cos(A) -sin(A); sin(A) cos(A)];
end

%===================================================================
function [gstop, isterminal,direction]=swing_catches_stance(t,z,walker)
%===================================================================
% スイング脚がスタンス脚に追いつく条件

q1 = z(1); q2 = z(3); 

% スイング脚角度がスタンス脚角度に追いつく条件
% q2 - q1 = 0 になったときに停止
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