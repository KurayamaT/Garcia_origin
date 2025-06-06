% SINGLE_STRIDE_WALKER - 1.5ストライド（スイング脚がスタンス脚に追いつくまで）を実行するパッシブウォーカー
% 使い方:
% single_stride_walker(1, [0.2, -0.25, 0.35, -0.4])  % 指定初期値から1.5歩
% single_stride_walker(1)                              % デフォルト初期値から1.5歩
% 成功例　single_stride_passivewalker(1, [0, -0.1, 0, -0.35])

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
                    [0.05, -0.15, 0.35, -0.25], '固定点近傍1';
                    [0.25, -0.25, 0.45, -0.05], '固定点近傍2';
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

if flag == 1
    %% Garcia's simplest walker
    walker.M = 1000; walker.m = 1.0; walker.I = 0.00; walker.l = 1.0; walker.w = 0.0; 
    walker.c = 1.0;  walker.r = 0.3; walker.g = 1.0; walker.gam = 0.009; 
    
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
% 1.5ストライド × 5回の連続歩行シミュレーション
fps = 10;
num_cycles = 5; % 歩行サイクル数
pause_duration = 0.5; % 各サイクル間の静止時間（秒）

fprintf('\n=== 1.5ストライド × %d回 連続歩行シミュレーション開始 ===\n', num_cycles);

% 全サイクルのデータを保存する変数
all_trajectories = [];
all_times = [];
all_final_states = [];
all_midpoint_states = [];
current_z0 = z0; % 現在の初期値

for cycle = 1:num_cycles
    fprintf('\n--- サイクル %d/%d ---\n', cycle, num_cycles);
    fprintf('初期状態: [%.6f, %.6f, %.6f, %.6f]\n', current_z0);
    
    %%%% 1.5ストライドの実行 %%%%
    [z_trajectory, t_trajectory, z_final, z_midpoint] = one_and_half_stride_detailed(current_z0, walker);
    
    % 時間オフセットを追加（前のサイクルからの連続性）
    if cycle > 1
        time_offset = all_times(end) + pause_duration; % 静止時間を追加
        t_trajectory = t_trajectory + time_offset;
    end
    
    % データを保存
    if cycle == 1
        all_trajectories = z_trajectory;
        all_times = t_trajectory;
    else
        % 静止期間を追加
        pause_start = all_times(end);
        pause_end = pause_start + pause_duration;
        pause_time = linspace(pause_start, pause_end, 50)'; % 静止期間のタイムステップ
        
        % 静止期間中は最後の状態を維持
        last_state = all_trajectories(end, :);
        pause_states = repmat(last_state, length(pause_time), 1);
        
        % データを連結
        all_trajectories = [all_trajectories; pause_states(2:end,:); z_trajectory(2:end,:)];
        all_times = [all_times; pause_time(2:end); t_trajectory(2:end)];
    end
    
    all_final_states = [all_final_states; z_final];
    all_midpoint_states = [all_midpoint_states; z_midpoint];
    
    fprintf('最終状態: [%.6f, %.6f, %.6f, %.6f]\n', z_final);
    fprintf('脚間角度差: %.6f rad\n', abs(z_final(1) - z_final(3)));
    
    % 次のサイクルの初期条件を設定
    % 位置角度は終了時の角度、速度は初期条件と同じに設定
    next_q1 = z_final(1);
    next_q2 = z_final(3);
    next_u1 = z0(2); % 元の初期角速度を使用
    next_u2 = z0(4); % 元の初期角速度を使用
    
    current_z0 = [next_q1, next_u1, next_q2, next_u2];
    
    if cycle < num_cycles
        fprintf('0.5秒静止後、次のサイクルへ...\n');
    end
end

fprintf('\n=== 1.5ストライド × %d回 連続歩行完了 ===\n', num_cycles);
fprintf('全体のシミュレーション時間: %.3f秒\n', all_times(end));

% 各サイクルの結果を表示
fprintf('\n=== 各サイクルの結果 ===\n');
for cycle = 1:num_cycles
    final_state = all_final_states(cycle, :);
    midpoint_state = all_midpoint_states(cycle, :);
    
    fprintf('サイクル %d:\n', cycle);
    fprintf('  1ステップ後: [%.4f, %.4f, %.4f, %.4f]\n', midpoint_state);
    fprintf('  最終状態    : [%.4f, %.4f, %.4f, %.4f]\n', final_state);
    fprintf('  脚間角度差  : %.6f rad\n', abs(final_state(1) - final_state(3)));
end

% 最終サイクルの詳細解析
z_final = all_final_states(end, :);
z_midpoint = all_midpoint_states(end, :);

% 固定点との比較
fprintf('\n=== 固定点との比較 ===\n');
fprintf('固定点     : [%.6f, %.6f, %.6f, %.6f]\n', zstar);
diff_from_fixed = z_final - zstar;
fprintf('最終差分   : [%.6f, %.6f, %.6f, %.6f]\n', diff_from_fixed);
fprintf('最終距離   : %.6f\n', norm(diff_from_fixed));

% 全体の歩行状態チェック
fprintf('\n=== 全体の歩行状態チェック ===\n');
fprintf('最終サイクルの脚間角度差: %.6f rad\n', abs(z_final(1) - z_final(3)));

success_count = 0;
for cycle = 1:num_cycles
    final_state = all_final_states(cycle, :);
    if abs(final_state(1) - final_state(3)) < 0.05
        success_count = success_count + 1;
    end
end

fprintf('成功したサイクル: %d/%d\n', success_count, num_cycles);
if success_count == num_cycles
    fprintf('✅ 全サイクルでスイング脚がスタンス脚に正常に追いつきました\n');
elseif success_count > num_cycles/2
    fprintf('⚠️ 大部分のサイクルでスイング脚が追いつきました\n');
else
    fprintf('❌ 多くのサイクルでスイング脚が追いついていません\n');
end

%%%% エネルギー解析 %%%%
fprintf('\n=== エネルギー解析 ===\n');
TE_initial = calculate_energy(z0, walker);
TE_final = calculate_energy(z_final, walker);
fprintf('初期エネルギー（1サイクル目開始時）: %.6f\n', TE_initial);
fprintf('最終エネルギー（%dサイクル目終了時）: %.6f\n', num_cycles, TE_final);
fprintf('全体のエネルギー損失: %.6f\n', TE_initial - TE_final);
fprintf('サイクルあたり平均エネルギー損失: %.6f\n', (TE_initial - TE_final)/num_cycles);

% 各サイクルのエネルギー変化
fprintf('\n各サイクルのエネルギー変化:\n');
prev_energy = TE_initial;
for cycle = 1:num_cycles
    cycle_final_energy = calculate_energy(all_final_states(cycle, :), walker);
    energy_loss = prev_energy - cycle_final_energy;
    fprintf('サイクル %d: %.6f → %.6f (損失: %.6f)\n', cycle, prev_energy, cycle_final_energy, energy_loss);
    prev_energy = cycle_final_energy;
end

%%%% 歩行成功判定 %%%%
fprintf('\n=== 全体の歩行成功判定 ===\n');
if abs(z_final(1)) < 1.0 && abs(z_final(3)) < 1.0  % 角度が妥当な範囲内
    fprintf('✅ %d回の1.5ストライド連続歩行が正常に完了しました\n', num_cycles);
else
    fprintf('❌ 歩行が不安定です（最終角度が大きすぎます）\n');
end

%%%% アニメーション %%%%
fprintf('\n%d回の1.5ストライド連続歩行のアニメーションを表示します...\n', num_cycles);
animate_continuous_walking(all_times, all_trajectories, walker, num_cycles);

%%%% プロット %%%%
fprintf('\n連続歩行の軌道をプロットします...\n');
plot_continuous_walking(all_times, all_trajectories, z0, z_final, zstar, all_final_states, all_midpoint_states, num_cycles);

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
function animate_continuous_walking(t_all, z_all, walker, num_cycles)
%===================================================================

fps = 10;
z_all_plot = [z_all(:,1) z_all(:,3) z_all(:,6) z_all(:,8)];
nn = size(z_all_plot,2);
total_frames = round(t_all(end)*fps);
t = linspace(0,t_all(end),total_frames);
z = zeros(total_frames,nn);
for i=1:nn
    z(:,i) = interp1(t_all,z_all_plot(:,i),t);
end

figure('Name', sprintf('1.5ストライド × %d回 連続歩行アニメーション', num_cycles), 'Position', [200, 200, 1000, 600]);
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

title(sprintf('1.5ストライド × %d回 連続歩行アニメーション', num_cycles), 'FontSize', 14)

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
   
   pause(0.05)  % 少し早くして全体を見やすく
   drawnow  
  
end

end

%===================================================================
function plot_continuous_walking(t, z, z0, z_final, zstar, all_final_states, all_midpoint_states, num_cycles)
%===================================================================

figure('Name', sprintf('1.5ストライド × %d回 連続歩行の軌道', num_cycles), 'Position', [100, 100, 1600, 1000]);

% 状態変数の時間変化
subplot(3,4,1)
plot(t, z(:,1), 'r-', 'LineWidth', 2)
hold on
plot([t(1), t(end)], [z0(1), z0(1)], 'g:', 'LineWidth', 1)
plot([t(1), t(end)], [z_final(1), z_final(1)], 'b:', 'LineWidth', 1)
plot([t(1), t(end)], [zstar(1), zstar(1)], 'k--', 'LineWidth', 1)
xlabel('時間 (s)')
ylabel('q1 (rad)')
title('スタンス脚角度')
legend('軌道', '初期値', '最終値', '固定点', 'Location', 'best')
grid on

subplot(3,4,2)
plot(t, z(:,2), 'b-', 'LineWidth', 2)
hold on
plot([t(1), t(end)], [z0(2), z0(2)], 'g:', 'LineWidth', 1)
plot([t(1), t(end)], [z_final(2), z_final(2)], 'b:', 'LineWidth', 1)
plot([t(1), t(end)], [zstar(2), zstar(2)], 'k--', 'LineWidth', 1)
xlabel('時間 (s)')
ylabel('u1 (rad/s)')
title('スタンス脚角速度')
legend('軌道', '初期値', '最終値', '固定点', 'Location', 'best')
grid on

subplot(3,4,3)
plot(t, z(:,3), 'g-', 'LineWidth', 2)
hold on
plot([t(1), t(end)], [z0(3), z0(3)], 'g:', 'LineWidth', 1)
plot([t(1), t(end)], [z_final(3), z_final(3)], 'b:', 'LineWidth', 1)
plot([t(1), t(end)], [zstar(3), zstar(3)], 'k--', 'LineWidth', 1)
xlabel('時間 (s)')
ylabel('q2 (rad)')
title('スイング脚角度')
legend('軌道', '初期値', '最終値', '固定点', 'Location', 'best')
grid on

subplot(3,4,4)
plot(t, z(:,4), 'm-', 'LineWidth', 2)
hold on
plot([t(1), t(end)], [z0(4), z0(4)], 'g:', 'LineWidth', 1)
plot([t(1), t(end)], [z_final(4), z_final(4)], 'b:', 'LineWidth', 1)
plot([t(1), t(end)], [zstar(4), zstar(4)], 'k--', 'LineWidth', 1)
xlabel('時間 (s)')
ylabel('u2 (rad/s)')
title('スイング脚角速度')
legend('軌道', '初期値', '最終値', '固定点', 'Location', 'best')
grid on

% エネルギーの変化
subplot(3,4,5)
plot(t, z(:,5), 'k-', 'LineWidth', 2)
xlabel('時間 (s)')
ylabel('総エネルギー')
title('エネルギー変化')
grid on

% 各サイクル終了時の状態変化
subplot(3,4,6)
cycle_numbers = 1:num_cycles;
plot(cycle_numbers, all_final_states(:,1), 'ro-', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(cycle_numbers, all_final_states(:,3), 'bo-', 'LineWidth', 2, 'MarkerSize', 8)
plot([1, num_cycles], [zstar(1), zstar(1)], 'r--', 'LineWidth', 1)
plot([1, num_cycles], [zstar(3), zstar(3)], 'b--', 'LineWidth', 1)
xlabel('サイクル番号')
ylabel('角度 (rad)')
title('各サイクル終了時の角度')
legend('q1', 'q2', 'q1固定点', 'q2固定点', 'Location', 'best')
grid on

% 脚間角度差の変化
subplot(3,4,7)
angle_differences = abs(all_final_states(:,1) - all_final_states(:,3));
plot(cycle_numbers, angle_differences, 'ko-', 'LineWidth', 2, 'MarkerSize', 8)
xlabel('サイクル番号')
ylabel('脚間角度差 (rad)')
title('各サイクル終了時の脚間角度差')
grid on

% ヒップ軌道
subplot(3,4,8)
plot(z(:,6), z(:,8), 'c-', 'LineWidth', 2)
xlabel('x (m)')
ylabel('y (m)')
title('ヒップ軌道（全体）')
grid on
axis equal

% 位相図 (q1-u1)
subplot(3,4,9)
plot(z(:,1), z(:,2), 'r-', 'LineWidth', 2)
hold on
plot(z0(1), z0(2), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g')
plot(z_final(1), z_final(2), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b')
plot(zstar(1), zstar(2), 'k*', 'MarkerSize', 12, 'MarkerFaceColor', 'k')
xlabel('q1 (rad)')
ylabel('u1 (rad/s)')
title('位相図 (q1-u1)')
legend('軌道', '初期点', '最終点', '固定点', 'Location', 'best')
grid on

% 位相図 (q2-u2)
subplot(3,4,10)
plot(z(:,3), z(:,4), 'g-', 'LineWidth', 2)
hold on
plot(z0(3), z0(4), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g')
plot(z_final(3), z_final(4), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b')
plot(zstar(3), zstar(4), 'k*', 'MarkerSize', 12, 'MarkerFaceColor', 'k')
xlabel('q2 (rad)')
ylabel('u2 (rad/s)')
title('位相図 (q2-u2)')
legend('軌道', '初期点', '最終点', '固定点', 'Location', 'best')
grid on

% エネルギー変化（サイクル別）
subplot(3,4,11)
walker_temp.M = 1000; walker_temp.m = 1.0; walker_temp.I = 0.00; walker_temp.l = 1.0; walker_temp.w = 0.0; 
walker_temp.c = 1.0;  walker_temp.r = 0.3; walker_temp.g = 1.0; walker_temp.gam = 0.009;
cycle_energies = zeros(num_cycles, 1);
for i = 1:num_cycles
    cycle_energies(i) = calculate_energy(all_final_states(i,:), walker_temp);
end
plot(cycle_numbers, cycle_energies, 'ko-', 'LineWidth', 2, 'MarkerSize', 8)
xlabel('サイクル番号')
ylabel('エネルギー')
title('各サイクル終了時のエネルギー')
grid on

% 状態空間での軌道
subplot(3,4,12)
state_norm = sqrt(sum(z(:,1:4).^2, 2));
plot(t, state_norm, 'k-', 'LineWidth', 2)
xlabel('時間 (s)')
ylabel('状態ベクトルのノルム')
title('状態空間での軌道')
grid on

sgtitle(sprintf('1.5ストライド × %d回 連続歩行の詳細解析', num_cycles))

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

% 以下、元のコードから必要な関数をコピー
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

end