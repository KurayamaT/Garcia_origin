function grid_walker()
    % パッシブウォーカーのグリッドサーチGUIプログラム
    % 初期条件と walker パラメータを GUI で設定し、グリッドサーチを実行
    
    clear; close all;
    
    fprintf('=== パッシブウォーカー グリッドサーチ GUI ===\n');
    
    % GUI データ構造体の初期化
    gui_data = struct();
    
    % デフォルトのwalkerパラメータ設定
    gui_data.walker = struct();
    gui_data.walker.M = 1000;     % Hip mass [kg]
    gui_data.walker.m = 1.0;      % Leg mass [kg]
    gui_data.walker.I = 0.00;     % Leg inertia [kg*m^2]
    gui_data.walker.l = 1.0;      % Leg length [m]
    gui_data.walker.w = 0.0;      % Leg width [m]
    gui_data.walker.c = 1.0;      % Leg CoM distance [m]
    gui_data.walker.r = 0.0;      % Foot radius [m]
    gui_data.walker.g = 1.0;      % Gravity [m/s^2]
    gui_data.walker.gam = 0.009;  % Slope angle [rad]
    
    % デフォルトのグリッドサーチパラメータ
    gui_data.grid = struct();
    gui_data.grid.q1_min = -0.5;
    gui_data.grid.q1_max = 0.5;
    gui_data.grid.q1_step = 0.1;
    gui_data.grid.u1_min = -1.0;
    gui_data.grid.u1_max = 1.0;
    gui_data.grid.u1_step = 0.2;
    gui_data.grid.q2_min = -0.5;
    gui_data.grid.q2_max = 0.5;
    gui_data.grid.q2_step = 0.1;
    gui_data.grid.u2_min = -1.0;
    gui_data.grid.u2_max = 1.0;
    gui_data.grid.u2_step = 0.2;
    gui_data.grid.thin_factor = 1;  % 間引き係数（1=間引きなし）
    
    % デフォルトの解析パラメータ
    gui_data.analysis = struct();
    gui_data.analysis.max_time = 100000000;      % 最大シミュレーション時間
    gui_data.analysis.num_cycles = 3;      % 解析サイクル数
    gui_data.analysis.success_threshold = 0.1;  % 成功判定閾値
    gui_data.analysis.target_q2 = 0.05;     % 目標スイング脚角度 [rad]
    
    % メインGUIの作成
    create_main_gui(gui_data);
end

function create_main_gui(gui_data)
    % メインGUIウィンドウの作成
    
    fig = figure('Name', 'パッシブウォーカー グリッドサーチ設定', ...
                 'Position', [100 100 800 600], ...
                 'MenuBar', 'none', ...
                 'NumberTitle', 'off', ...
                 'Resize', 'off');
    
    % タブの代わりにボタンで切り替える方式に変更
    % タブボタンの作成
    tab_buttons = [];
    tab_names = {'Walker パラメータ', 'グリッドサーチ設定', '解析設定', '実行・結果'};
    for i = 1:4
        tab_buttons(i) = uicontrol(fig, 'Style', 'pushbutton', ...
                                  'String', tab_names{i}, ...
                                  'Position', [20 + (i-1)*180, 570, 160, 25], ...
                                  'Tag', sprintf('tab_button_%d', i), ...
                                  'Callback', @(src,evt) switch_tab(fig, i));
    end
    
    % タブコンテンツ用のパネル
    for i = 1:4
        tab_panels(i) = uipanel(fig, 'Position', [0.02 0.02 0.96 0.9], ...
                               'Visible', 'off', ...
                               'Tag', sprintf('tab_panel_%d', i));
    end
    
    % 各タブの内容を作成
    create_walker_tab(tab_panels(1), gui_data);
    create_grid_tab(tab_panels(2), gui_data);
    create_analysis_tab(tab_panels(3), gui_data);
    create_execution_tab(tab_panels(4), gui_data);
    
    % 最初のタブを表示
    set(tab_panels(1), 'Visible', 'on');
    set(tab_buttons(1), 'BackgroundColor', [0.8 0.8 1.0]);
    
    % データをfigureに保存
    setappdata(fig, 'gui_data', gui_data);
    setappdata(fig, 'tab_panels', tab_panels);
    setappdata(fig, 'tab_buttons', tab_buttons);
    setappdata(fig, 'current_tab', 1);
end

function switch_tab(fig, tab_num)
    % タブ切り替え関数
    
    tab_panels = getappdata(fig, 'tab_panels');
    tab_buttons = getappdata(fig, 'tab_buttons');
    current_tab = getappdata(fig, 'current_tab');
    
    % 現在のタブを非表示
    if current_tab <= length(tab_panels)
        set(tab_panels(current_tab), 'Visible', 'off');
        set(tab_buttons(current_tab), 'BackgroundColor', get(0,'defaultuicontrolbackgroundcolor'));
    end
    
    % 新しいタブを表示
    if tab_num <= length(tab_panels)
        set(tab_panels(tab_num), 'Visible', 'on');
        set(tab_buttons(tab_num), 'BackgroundColor', [0.8 0.8 1.0]);
        setappdata(fig, 'current_tab', tab_num);
    end
end

function create_walker_tab(parent, gui_data)
    % Walker パラメータ設定タブ
    
    % ウォーカータイプ選択
    uicontrol(parent, 'Style', 'text', 'String', 'ウォーカータイプ:', ...
              'Position', [20 520 100 20], 'HorizontalAlignment', 'left');
    
    type_popup = uicontrol(parent, 'Style', 'popupmenu', ...
                          'String', {'Garcia''s Simplest Walker', 'General Round Feet Walker'}, ...
                          'Position', [130 520 200 25], ...
                          'Callback', @(src,evt) update_walker_preset(src, evt, parent));
    
    % パラメータ入力フィールド
    param_names = {'M', 'm', 'I', 'l', 'w', 'c', 'r', 'g', 'gam'};
    param_labels = {'Hip mass [kg]:', 'Leg mass [kg]:', 'Leg inertia [kg*m^2]:', ...
                   'Leg length [m]:', 'Leg width [m]:', 'Leg CoM distance [m]:', ...
                   'Foot radius [m]:', 'Gravity [m/s^2]:', 'Slope angle [rad]:'};
    param_values = [gui_data.walker.M, gui_data.walker.m, gui_data.walker.I, ...
                   gui_data.walker.l, gui_data.walker.w, gui_data.walker.c, ...
                   gui_data.walker.r, gui_data.walker.g, gui_data.walker.gam];
    
    y_pos = 480;
    for i = 1:length(param_names)
        uicontrol(parent, 'Style', 'text', 'String', param_labels{i}, ...
                  'Position', [20 y_pos 150 20], 'HorizontalAlignment', 'left');
        
        edit_handle = uicontrol(parent, 'Style', 'edit', ...
                               'String', num2str(param_values(i)), ...
                               'Position', [180 y_pos 100 25], ...
                               'Tag', ['walker_' param_names{i}]);
        
        y_pos = y_pos - 35;
    end
    
    % プリセット適用ボタン
    uicontrol(parent, 'Style', 'pushbutton', 'String', 'プリセット適用', ...
              'Position', [350 520 100 25], ...
              'Callback', @(src,evt) apply_walker_preset(src, evt, parent, type_popup));
    
    % 現在の設定表示
    uicontrol(parent, 'Style', 'text', 'String', '現在の設定:', ...
              'Position', [400 480 100 20], 'HorizontalAlignment', 'left', ...
              'FontWeight', 'bold');
    
    create_walker_display(parent, gui_data);
end

function create_grid_tab(parent, gui_data)
    % グリッドサーチ設定タブ
    
    % タイトル
    uicontrol(parent, 'Style', 'text', 'String', '初期条件グリッドサーチ設定', ...
              'Position', [20 540 200 20], 'HorizontalAlignment', 'left', ...
              'FontWeight', 'bold', 'FontSize', 12);
    
    % q1 設定（スタンス脚角度）
    y_start = 480;
    uicontrol(parent, 'Style', 'text', 'String', 'q1 (スタンス脚角度) [rad]:', ...
              'Position', [20 y_start+40 180 20], 'HorizontalAlignment', 'left', ...
              'FontWeight', 'bold', 'FontSize', 10);
    create_variable_controls(parent, 'q1', gui_data.grid.q1_min, gui_data.grid.q1_max, ...
                           gui_data.grid.q1_step, y_start);
    
    % u1 設定（スタンス脚角速度）
    y_start = 380;
    uicontrol(parent, 'Style', 'text', 'String', 'u1 (スタンス脚角速度) [rad/s]:', ...
              'Position', [20 y_start+40 180 20], 'HorizontalAlignment', 'left', ...
              'FontWeight', 'bold', 'FontSize', 10);
    create_variable_controls(parent, 'u1', gui_data.grid.u1_min, gui_data.grid.u1_max, ...
                           gui_data.grid.u1_step, y_start);
    
    % q2 設定（スイング脚角度）
    y_start = 280;
    uicontrol(parent, 'Style', 'text', 'String', 'q2 (スイング脚角度) [rad]:', ...
              'Position', [20 y_start+40 180 20], 'HorizontalAlignment', 'left', ...
              'FontWeight', 'bold', 'FontSize', 10);
    create_variable_controls(parent, 'q2', gui_data.grid.q2_min, gui_data.grid.q2_max, ...
                           gui_data.grid.q2_step, y_start);
    
    % u2 設定（スイング脚角速度）
    y_start = 180;
    uicontrol(parent, 'Style', 'text', 'String', 'u2 (スイング脚角速度) [rad/s]:', ...
              'Position', [20 y_start+40 180 20], 'HorizontalAlignment', 'left', ...
              'FontWeight', 'bold', 'FontSize', 10);
    create_variable_controls(parent, 'u2', gui_data.grid.u2_min, gui_data.grid.u2_max, ...
                           gui_data.grid.u2_step, y_start);
    
    % 間引き設定
    uicontrol(parent, 'Style', 'text', 'String', '間引き係数:', ...
              'Position', [20 120 80 20], 'HorizontalAlignment', 'left', ...
              'FontWeight', 'bold');
    uicontrol(parent, 'Style', 'edit', 'String', num2str(gui_data.grid.thin_factor), ...
              'Position', [110 120 60 25], 'Tag', 'thin_factor');
    uicontrol(parent, 'Style', 'text', 'String', '(1=間引きなし)', ...
              'Position', [180 120 80 20], 'HorizontalAlignment', 'left', ...
              'FontSize', 8);
    
    % グリッドサイズ表示
    uicontrol(parent, 'Style', 'text', 'String', 'グリッドサイズ:', ...
              'Position', [280 120 80 20], 'HorizontalAlignment', 'left', ...
              'FontWeight', 'bold');
    uicontrol(parent, 'Style', 'text', 'String', '計算中...', ...
              'Position', [370 120 100 20], 'HorizontalAlignment', 'left', ...
              'Tag', 'grid_size_display', 'FontWeight', 'bold', ...
              'ForegroundColor', [0 0 0.8]);
    
    % サイズ更新ボタン
    uicontrol(parent, 'Style', 'pushbutton', 'String', 'サイズ更新', ...
              'Position', [480 120 80 25], ...
              'Callback', @(src,evt) update_grid_size_display(parent));
    
    % プレビューボタン
    uicontrol(parent, 'Style', 'pushbutton', 'String', 'グリッドプレビュー', ...
              'Position', [570 120 100 25], ...
              'Callback', @(src,evt) preview_grid(parent));
    
    % 推奨値設定ボタン
    uicontrol(parent, 'Style', 'pushbutton', 'String', '推奨値設定', ...
              'Position', [480 80 80 25], 'BackgroundColor', [0.9 0.9 0.6], ...
              'Callback', @(src,evt) set_recommended_grid(parent));
    
    % 全リセットボタン
    uicontrol(parent, 'Style', 'pushbutton', 'String', '全リセット', ...
              'Position', [570 80 80 25], 'BackgroundColor', [1.0 0.8 0.8], ...
              'Callback', @(src,evt) reset_grid_settings(parent));
end

function create_analysis_tab(parent, gui_data)
    % 解析設定タブ
    
    % タイトル
    uicontrol(parent, 'Style', 'text', 'String', '解析・計算設定', ...
              'Position', [20 540 150 20], 'HorizontalAlignment', 'left', ...
              'FontWeight', 'bold', 'FontSize', 12);
    
    % 最大シミュレーション時間
    uicontrol(parent, 'Style', 'text', 'String', '最大シミュレーション時間 [s]:', ...
              'Position', [20 500 180 20], 'HorizontalAlignment', 'left');
    uicontrol(parent, 'Style', 'edit', 'String', num2str(gui_data.analysis.max_time), ...
              'Position', [210 500 80 25], 'Tag', 'max_time');
    
    % 解析サイクル数
    uicontrol(parent, 'Style', 'text', 'String', '解析サイクル数:', ...
              'Position', [20 460 120 20], 'HorizontalAlignment', 'left');
    uicontrol(parent, 'Style', 'edit', 'String', num2str(gui_data.analysis.num_cycles), ...
              'Position', [150 460 80 25], 'Tag', 'num_cycles');
    
    % 成功判定閾値
    uicontrol(parent, 'Style', 'text', 'String', '成功判定閾値:', ...
              'Position', [20 420 120 20], 'HorizontalAlignment', 'left');
    uicontrol(parent, 'Style', 'edit', 'String', num2str(gui_data.analysis.success_threshold), ...
              'Position', [150 420 80 25], 'Tag', 'success_threshold');
    
    % 目標スイング脚角度
    uicontrol(parent, 'Style', 'text', 'String', '目標スイング脚角度 [rad]:', ...
              'Position', [20 380 150 20], 'HorizontalAlignment', 'left');
    uicontrol(parent, 'Style', 'edit', 'String', num2str(gui_data.analysis.target_q2), ...
              'Position', [180 380 80 25], 'Tag', 'target_q2');
    uicontrol(parent, 'Style', 'text', 'String', '(0.5ストライド終了条件)', ...
              'Position', [270 380 150 20], 'HorizontalAlignment', 'left', ...
              'FontSize', 8, 'ForegroundColor', [0.5 0.5 0.5]);
    
    % 並列処理設定
    uicontrol(parent, 'Style', 'text', 'String', '並列計算設定', ...
              'Position', [20 340 120 20], 'HorizontalAlignment', 'left', ...
              'FontWeight', 'bold', 'FontSize', 11);
    
    uicontrol(parent, 'Style', 'checkbox', 'String', '並列処理を使用（parfor）', ...
              'Position', [20 310 180 25], 'Tag', 'use_parallel', 'Value', 1);
    
    % 並列プール情報表示
    uicontrol(parent, 'Style', 'text', 'String', '現在の並列プール状況:', ...
              'Position', [20 280 150 20], 'HorizontalAlignment', 'left');
    
    pool_status_text = uicontrol(parent, 'Style', 'text', 'String', get_pool_status(), ...
                                'Position', [170 280 300 20], 'HorizontalAlignment', 'left', ...
                                'Tag', 'pool_status', 'FontSize', 9);
    
    % プール管理ボタン
    uicontrol(parent, 'Style', 'pushbutton', 'String', 'プール開始', ...
              'Position', [20 250 80 25], ...
              'Callback', @(src,evt) start_pool_callback(parent));
    
    uicontrol(parent, 'Style', 'pushbutton', 'String', 'プール停止', ...
              'Position', [110 250 80 25], ...
              'Callback', @(src,evt) stop_pool_callback(parent));
    
    uicontrol(parent, 'Style', 'pushbutton', 'String', '状況更新', ...
              'Position', [200 250 80 25], ...
              'Callback', @(src,evt) update_pool_status(parent));
    
    % 出力設定
    uicontrol(parent, 'Style', 'text', 'String', '出力・ログ設定', ...
              'Position', [20 210 120 20], 'HorizontalAlignment', 'left', ...
              'FontWeight', 'bold', 'FontSize', 11);
    
    uicontrol(parent, 'Style', 'checkbox', 'String', '詳細ログ出力', ...
              'Position', [20 180 120 25], 'Tag', 'verbose_output');
    
    uicontrol(parent, 'Style', 'checkbox', 'String', '結果をCSVで保存', ...
              'Position', [150 180 140 25], 'Tag', 'save_csv', 'Value', 1);
    
    uicontrol(parent, 'Style', 'checkbox', 'String', '進捗バー表示', ...
              'Position', [20 150 120 25], 'Tag', 'show_progress', 'Value', 1);
    
    % パフォーマンス設定
    uicontrol(parent, 'Style', 'text', 'String', 'パフォーマンス設定', ...
              'Position', [20 110 120 20], 'HorizontalAlignment', 'left', ...
              'FontWeight', 'bold', 'FontSize', 11);
    
    uicontrol(parent, 'Style', 'text', 'String', 'ODE許容誤差 (abstol):', ...
              'Position', [20 80 140 20], 'HorizontalAlignment', 'left');
    uicontrol(parent, 'Style', 'edit', 'String', '1e-10', ...
              'Position', [170 80 80 25], 'Tag', 'ode_abstol');
    
    uicontrol(parent, 'Style', 'text', 'String', '時間解像度:', ...
              'Position', [20 50 80 20], 'HorizontalAlignment', 'left');
    uicontrol(parent, 'Style', 'edit', 'String', '50', ...
              'Position', [110 50 60 25], 'Tag', 'time_resolution');
    
    % 推定計算時間表示
    uicontrol(parent, 'Style', 'text', 'String', '推定計算時間:', ...
              'Position', [400 200 100 20], 'HorizontalAlignment', 'left', ...
              'FontWeight', 'bold');
    uicontrol(parent, 'Style', 'text', 'String', '設定後に計算されます', ...
              'Position', [400 180 200 20], 'HorizontalAlignment', 'left', ...
              'Tag', 'estimated_time', 'FontSize', 9);
    
    % 計算時間見積もりボタン
    uicontrol(parent, 'Style', 'pushbutton', 'String', '計算時間見積もり', ...
              'Position', [400 150 120 25], ...
              'Callback', @(src,evt) estimate_computation_time(parent));
end

function pool_status = get_pool_status()
    % 現在の並列プール状況を取得
    
    try
        current_pool = gcp('nocreate');
        if isempty(current_pool)
            pool_status = 'プールなし（並列処理無効）';
        else
            pool_status = sprintf('%d workers (アクティブ)', current_pool.NumWorkers);
        end
    catch
        pool_status = 'Parallel Computing Toolbox 未利用可能';
    end
end

function start_pool_callback(parent)
    % 並列プール開始コールバック
    
    try
        current_pool = gcp('nocreate');
        if isempty(current_pool)
            fprintf('並列計算プールを開始中...\n');
            parpool('local');
            fprintf('並列計算プール開始完了\n');
        else
            fprintf('並列計算プールは既に動作中です\n');
        end
        update_pool_status(parent);
    catch me
        fprintf('並列プール開始エラー: %s\n', me.message);
    end
end

function stop_pool_callback(parent)
    % 並列プール停止コールバック
    
    try
        current_pool = gcp('nocreate');
        if ~isempty(current_pool)
            fprintf('並列計算プールを停止中...\n');
            delete(current_pool);
            fprintf('並列計算プール停止完了\n');
        else
            fprintf('停止する並列計算プールがありません\n');
        end
        update_pool_status(parent);
    catch me
        fprintf('並列プール停止エラー: %s\n', me.message);
    end
end

function update_pool_status(parent)
    % プール状況表示の更新
    
    pool_status_text = findobj(parent, 'Tag', 'pool_status');
    if ~isempty(pool_status_text)
        set(pool_status_text, 'String', get_pool_status());
    end
end

function estimate_computation_time(parent)
    % 計算時間の見積もり
    
    % 現在の設定を読み込み
    fig = ancestor(parent, 'figure');
    gui_data = getappdata(fig, 'gui_data');
    gui_data = read_current_settings(fig, gui_data);
    
    % グリッドサイズを計算
    q1_range = gui_data.grid.q1_min:gui_data.grid.q1_step:gui_data.grid.q1_max;
    u1_range = gui_data.grid.u1_min:gui_data.grid.u1_step:gui_data.grid.u1_max;
    q2_range = gui_data.grid.q2_min:gui_data.grid.q2_step:gui_data.grid.q2_max;
    u2_range = gui_data.grid.u2_min:gui_data.grid.u2_step:gui_data.grid.u2_max;
    
    total_points = length(q1_range) * length(u1_range) * length(q2_range) * length(u2_range);
    
    % 間引きを考慮
    if gui_data.grid.thin_factor > 1
        total_points = ceil(total_points / gui_data.grid.thin_factor);
    end
    
    % 計算時間の見積もり（経験的な値）
    time_per_point_sequential = 0.1; % 秒/点（逐次計算）
    time_per_point_parallel = 0.03;  % 秒/点（並列計算、4コア想定）
    
    if gui_data.analysis.use_parallel
        try
            current_pool = gcp('nocreate');
            if ~isempty(current_pool)
                estimated_time = total_points * time_per_point_parallel / current_pool.NumWorkers;
                mode_str = sprintf('並列処理 (%d workers)', current_pool.NumWorkers);
            else
                estimated_time = total_points * time_per_point_sequential;
                mode_str = '逐次処理（プールなし）';
            end
        catch
            estimated_time = total_points * time_per_point_sequential;
            mode_str = '逐次処理';
        end
    else
        estimated_time = total_points * time_per_point_sequential;
        mode_str = '逐次処理';
    end
    
    % 時間の表示形式を調整
    if estimated_time < 60
        time_str = sprintf('約%.1f秒', estimated_time);
    elseif estimated_time < 3600
        time_str = sprintf('約%.1f分', estimated_time/60);
    else
        time_str = sprintf('約%.1f時間', estimated_time/3600);
    end
    
    % 結果表示
    estimated_time_text = findobj(parent, 'Tag', 'estimated_time');
    if ~isempty(estimated_time_text)
        result_str = sprintf('%s (%s)', time_str, mode_str);
        set(estimated_time_text, 'String', result_str);
    end
    
    fprintf('計算時間見積もり: %s (総点数: %d, %s)\n', time_str, total_points, mode_str);
end

function create_execution_tab(parent, gui_data)
    % 実行・結果タブ
    
    % 設定確認表示エリア
    uicontrol(parent, 'Style', 'text', 'String', '設定確認:', ...
              'Position', [20 520 80 20], 'HorizontalAlignment', 'left', ...
              'FontWeight', 'bold');
    
    summary_text = uicontrol(parent, 'Style', 'text', 'String', '設定を読み込み中...', ...
                            'Position', [20 450 350 60], 'HorizontalAlignment', 'left', ...
                            'Tag', 'settings_summary');
    
    % 実行ボタン
    uicontrol(parent, 'Style', 'pushbutton', 'String', 'グリッドサーチ実行', ...
              'Position', [400 500 120 40], 'FontSize', 12, 'FontWeight', 'bold', ...
              'BackgroundColor', [0.2 0.8 0.2], ...
              'Callback', @(src,evt) execute_grid_search(src, evt, parent));
    
    % 進捗表示エリア
    uicontrol(parent, 'Style', 'text', 'String', '進捗:', ...
              'Position', [20 400 50 20], 'HorizontalAlignment', 'left');
    
    % 進捗バー用のaxes（簡略化）
    progress_axes = axes('Parent', parent, 'Position', [0.1 0.65 0.8 0.05], ...
                        'XLim', [0 1], 'YLim', [0 1], 'Box', 'on', ...
                        'XTick', [], 'YTick', [], 'Tag', 'progress_axes');
    
    % ログ表示エリア
    uicontrol(parent, 'Style', 'text', 'String', 'ログ:', ...
              'Position', [20 350 50 20], 'HorizontalAlignment', 'left');
    
    log_listbox = uicontrol(parent, 'Style', 'listbox', ...
                           'Position', [20 50 700 290], ...
                           'Tag', 'log_listbox');
    
    % 結果表示ボタン
    uicontrol(parent, 'Style', 'pushbutton', 'String', '結果可視化', ...
              'Position', [550 500 80 25], ...
              'Callback', @(src,evt) visualize_results(parent));
    
    uicontrol(parent, 'Style', 'pushbutton', 'String', '結果エクスポート', ...
              'Position', [640 500 80 25], ...
              'Callback', @(src,evt) export_results(parent));
end

function create_variable_controls(parent, var_name, min_val, max_val, step_val, y_pos)
    % 変数設定コントロールの作成
    
    % 最小値
    uicontrol(parent, 'Style', 'text', 'String', '最小値:', ...
              'Position', [60 y_pos 50 20], 'HorizontalAlignment', 'left');
    uicontrol(parent, 'Style', 'edit', 'String', num2str(min_val), ...
              'Position', [110 y_pos 70 25], 'Tag', [var_name '_min'], ...
              'Callback', @(src,evt) update_variable_count(parent, var_name));
    
    % 最大値
    uicontrol(parent, 'Style', 'text', 'String', '最大値:', ...
              'Position', [190 y_pos 50 20], 'HorizontalAlignment', 'left');
    uicontrol(parent, 'Style', 'edit', 'String', num2str(max_val), ...
              'Position', [240 y_pos 70 25], 'Tag', [var_name '_max'], ...
              'Callback', @(src,evt) update_variable_count(parent, var_name));
    
    % 刻み幅
    uicontrol(parent, 'Style', 'text', 'String', '刻み幅:', ...
              'Position', [320 y_pos 50 20], 'HorizontalAlignment', 'left');
    uicontrol(parent, 'Style', 'edit', 'String', num2str(step_val), ...
              'Position', [370 y_pos 70 25], 'Tag', [var_name '_step'], ...
              'Callback', @(src,evt) update_variable_count(parent, var_name));
    
    % サイズ表示
    num_points = calculate_points(min_val, max_val, step_val);
    uicontrol(parent, 'Style', 'text', 'String', sprintf('点数: %d', num_points), ...
              'Position', [450 y_pos 80 20], 'HorizontalAlignment', 'left', ...
              'Tag', [var_name '_count'], 'ForegroundColor', [0 0.6 0]);
    
    % 範囲表示
    range_str = sprintf('[%.3f, %.3f]', min_val, max_val);
    uicontrol(parent, 'Style', 'text', 'String', range_str, ...
              'Position', [540 y_pos 120 20], 'HorizontalAlignment', 'left', ...
              'Tag', [var_name '_range'], 'FontSize', 8, 'ForegroundColor', [0.5 0.5 0.5]);
end

function num_points = calculate_points(min_val, max_val, step_val)
    % 点数計算
    if step_val > 0 && max_val > min_val
        num_points = length(min_val:step_val:max_val);
    else
        num_points = 0;
    end
end

function update_variable_count(parent, var_name)
    % 個別変数の点数更新
    
    min_edit = findobj(parent, 'Tag', [var_name '_min']);
    max_edit = findobj(parent, 'Tag', [var_name '_max']);
    step_edit = findobj(parent, 'Tag', [var_name '_step']);
    count_text = findobj(parent, 'Tag', [var_name '_count']);
    range_text = findobj(parent, 'Tag', [var_name '_range']);
    
    if ~isempty(min_edit) && ~isempty(max_edit) && ~isempty(step_edit)
        min_val = str2double(get(min_edit, 'String'));
        max_val = str2double(get(max_edit, 'String'));
        step_val = str2double(get(step_edit, 'String'));
        
        num_points = calculate_points(min_val, max_val, step_val);
        
        if ~isempty(count_text)
            set(count_text, 'String', sprintf('点数: %d', num_points));
        end
        
        if ~isempty(range_text)
            range_str = sprintf('[%.3f, %.3f]', min_val, max_val);
            set(range_text, 'String', range_str);
        end
    end
    
    % 全体のグリッドサイズも更新
    update_grid_size_display(parent);
end

function create_walker_display(parent, gui_data)
    % Walker パラメータ表示エリア
    
    y_pos = 450;
    param_names = fieldnames(gui_data.walker);
    
    for i = 1:length(param_names)
        param_name = param_names{i};
        param_value = gui_data.walker.(param_name);
        
        text_str = sprintf('%s: %.6g', param_name, param_value);
        uicontrol(parent, 'Style', 'text', 'String', text_str, ...
                  'Position', [400 y_pos 150 20], 'HorizontalAlignment', 'left', ...
                  'Tag', ['display_' param_name]);
        
        y_pos = y_pos - 25;
    end
end

function update_walker_preset(src, evt, parent)
    % ウォーカープリセットの更新
    % この関数は実際のプリセット適用のトリガーとして機能
end

function apply_walker_preset(src, evt, parent, type_popup)
    % ウォーカープリセットの適用
    
    walker_type = get(type_popup, 'Value');
    
    if walker_type == 1
        % Garcia's simplest walker
        preset_values = [1000, 1.0, 0.00, 1.0, 0.0, 1.0, 0.0, 1.0, 0.009];
    else
        % General round feet walker
        preset_values = [1.0, 0.5, 0.02, 1.0, 0.0, 0.5, 0.2, 1.0, 0.01];
    end
    
    param_names = {'M', 'm', 'I', 'l', 'w', 'c', 'r', 'g', 'gam'};
    
    for i = 1:length(param_names)
        edit_handle = findobj(parent, 'Tag', ['walker_' param_names{i}]);
        if ~isempty(edit_handle)
            set(edit_handle, 'String', num2str(preset_values(i)));
        end
    end
    
    popup_strings = get(type_popup, 'String');
    fprintf('プリセット適用完了: %s\n', popup_strings{walker_type});
end

function update_grid_size_display(parent)
    % グリッドサイズ表示の更新
    
    vars = {'q1', 'u1', 'q2', 'u2'};
    total_points = 1;
    
    for i = 1:length(vars)
        var_name = vars{i};
        
        min_edit = findobj(parent, 'Tag', [var_name '_min']);
        max_edit = findobj(parent, 'Tag', [var_name '_max']);
        step_edit = findobj(parent, 'Tag', [var_name '_step']);
        count_text = findobj(parent, 'Tag', [var_name '_count']);
        
        if ~isempty(min_edit) && ~isempty(max_edit) && ~isempty(step_edit)
            min_val = str2double(get(min_edit, 'String'));
            max_val = str2double(get(max_edit, 'String'));
            step_val = str2double(get(step_edit, 'String'));
            
            if step_val > 0 && max_val > min_val
                num_points = length(min_val:step_val:max_val);
                total_points = total_points * num_points;
                
                if ~isempty(count_text)
                    set(count_text, 'String', sprintf('点数: %d', num_points));
                end
            end
        end
    end
    
    % 間引き係数を考慮
    thin_edit = findobj(parent, 'Tag', 'thin_factor');
    if ~isempty(thin_edit)
        thin_factor = str2double(get(thin_edit, 'String'));
        if thin_factor > 1
            total_points = ceil(total_points / thin_factor);
        end
    end
    
    size_display = findobj(parent, 'Tag', 'grid_size_display');
    if ~isempty(size_display)
        set(size_display, 'String', sprintf('総計: %d点', total_points));
    end
end

function set_recommended_grid(parent)
    % 推奨グリッド設定の適用
    
    % パッシブウォーカーでよく使われる範囲の推奨値
    recommended_settings = struct();
    recommended_settings.q1_min = -0.3;
    recommended_settings.q1_max = 0.3;
    recommended_settings.q1_step = 0.05;
    
    recommended_settings.u1_min = -0.8;
    recommended_settings.u1_max = 0.8;
    recommended_settings.u1_step = 0.1;
    
    recommended_settings.q2_min = -0.4;
    recommended_settings.q2_max = 0.4;
    recommended_settings.q2_step = 0.05;
    
    recommended_settings.u2_min = -0.6;
    recommended_settings.u2_max = 0.6;
    recommended_settings.u2_step = 0.1;
    
    % 設定を適用
    vars = {'q1', 'u1', 'q2', 'u2'};
    for i = 1:length(vars)
        var_name = vars{i};
        
        min_edit = findobj(parent, 'Tag', [var_name '_min']);
        max_edit = findobj(parent, 'Tag', [var_name '_max']);
        step_edit = findobj(parent, 'Tag', [var_name '_step']);
        
        if ~isempty(min_edit)
            set(min_edit, 'String', num2str(recommended_settings.([var_name '_min'])));
        end
        if ~isempty(max_edit)
            set(max_edit, 'String', num2str(recommended_settings.([var_name '_max'])));
        end
        if ~isempty(step_edit)
            set(step_edit, 'String', num2str(recommended_settings.([var_name '_step'])));
        end
        
        % 個別の点数も更新
        update_variable_count(parent, var_name);
    end
    
    % 間引き係数も推奨値に設定
    thin_edit = findobj(parent, 'Tag', 'thin_factor');
    if ~isempty(thin_edit)
        set(thin_edit, 'String', '2');  % 推奨間引き係数
    end
    
    update_grid_size_display(parent);
    fprintf('推奨グリッド設定を適用しました\n');
end

function reset_grid_settings(parent)
    % グリッド設定のリセット
    
    % デフォルト値
    default_settings = struct();
    default_settings.q1_min = -0.5;
    default_settings.q1_max = 0.5;
    default_settings.q1_step = 0.1;
    
    default_settings.u1_min = -1.0;
    default_settings.u1_max = 1.0;
    default_settings.u1_step = 0.2;
    
    default_settings.q2_min = -0.5;
    default_settings.q2_max = 0.5;
    default_settings.q2_step = 0.1;
    
    default_settings.u2_min = -1.0;
    default_settings.u2_max = 1.0;
    default_settings.u2_step = 0.2;
    
    % 設定を適用
    vars = {'q1', 'u1', 'q2', 'u2'};
    for i = 1:length(vars)
        var_name = vars{i};
        
        min_edit = findobj(parent, 'Tag', [var_name '_min']);
        max_edit = findobj(parent, 'Tag', [var_name '_max']);
        step_edit = findobj(parent, 'Tag', [var_name '_step']);
        
        if ~isempty(min_edit)
            set(min_edit, 'String', num2str(default_settings.([var_name '_min'])));
        end
        if ~isempty(max_edit)
            set(max_edit, 'String', num2str(default_settings.([var_name '_max'])));
        end
        if ~isempty(step_edit)
            set(step_edit, 'String', num2str(default_settings.([var_name '_step'])));
        end
        
        % 個別の点数も更新
        update_variable_count(parent, var_name);
    end
    
    % 間引き係数もリセット
    thin_edit = findobj(parent, 'Tag', 'thin_factor');
    if ~isempty(thin_edit)
        set(thin_edit, 'String', '1');  % デフォルト間引き係数
    end
    
    update_grid_size_display(parent);
    fprintf('グリッド設定をデフォルト値にリセットしました\n');
end

function preview_grid(parent)
    % グリッドプレビュー表示
    
    vars = {'q1', 'u1', 'q2', 'u2'};
    grid_data = struct();
    
    % グリッドデータの読み込み
    for i = 1:length(vars)
        var_name = vars{i};
        
        min_edit = findobj(parent, 'Tag', [var_name '_min']);
        max_edit = findobj(parent, 'Tag', [var_name '_max']);
        step_edit = findobj(parent, 'Tag', [var_name '_step']);
        
        min_val = str2double(get(min_edit, 'String'));
        max_val = str2double(get(max_edit, 'String'));
        step_val = str2double(get(step_edit, 'String'));
        
        grid_data.(var_name) = min_val:step_val:max_val;
    end
    
    % プレビューウィンドウの作成
    preview_fig = figure('Name', 'グリッドプレビュー', 'Position', [200 200 800 600]);
    
    % 2Dプロット（q1-q2, u1-u2）
    subplot(2,2,1);
    [Q1, Q2] = meshgrid(grid_data.q1, grid_data.q2);
    plot(Q1(:), Q2(:), 'bo', 'MarkerSize', 3);
    xlabel('q1 [rad]'); ylabel('q2 [rad]');
    title('角度グリッド (q1-q2)');
    grid on; axis equal;
    
    subplot(2,2,2);
    [U1, U2] = meshgrid(grid_data.u1, grid_data.u2);
    plot(U1(:), U2(:), 'ro', 'MarkerSize', 3);
    xlabel('u1 [rad/s]'); ylabel('u2 [rad/s]');
    title('角速度グリッド (u1-u2)');
    grid on; axis equal;
    
    % 統計情報
    subplot(2,2,3);
    text(0.1, 0.8, sprintf('q1 点数: %d', length(grid_data.q1)), 'FontSize', 12);
    text(0.1, 0.7, sprintf('u1 点数: %d', length(grid_data.u1)), 'FontSize', 12);
    text(0.1, 0.6, sprintf('q2 点数: %d', length(grid_data.q2)), 'FontSize', 12);
    text(0.1, 0.5, sprintf('u2 点数: %d', length(grid_data.u2)), 'FontSize', 12);
    text(0.1, 0.3, sprintf('総格子点数: %d', length(grid_data.q1) * length(grid_data.u1) * ...
                           length(grid_data.q2) * length(grid_data.u2)), 'FontSize', 12, 'FontWeight', 'bold');
    xlim([0 1]); ylim([0 1]); axis off;
    title('グリッド統計');
    
    % 範囲情報
    subplot(2,2,4);
    text(0.1, 0.8, sprintf('q1: [%.3f, %.3f]', min(grid_data.q1), max(grid_data.q1)), 'FontSize', 10);
    text(0.1, 0.7, sprintf('u1: [%.3f, %.3f]', min(grid_data.u1), max(grid_data.u1)), 'FontSize', 10);
    text(0.1, 0.6, sprintf('q2: [%.3f, %.3f]', min(grid_data.q2), max(grid_data.q2)), 'FontSize', 10);
    text(0.1, 0.5, sprintf('u2: [%.3f, %.3f]', min(grid_data.u2), max(grid_data.u2)), 'FontSize', 10);
    xlim([0 1]); ylim([0 1]); axis off;
    title('範囲情報');
end

function execute_grid_search(src, evt, parent)
    % グリッドサーチの実行
    
    % GUIからデータを読み込み
    fig = ancestor(parent, 'figure');
    gui_data = getappdata(fig, 'gui_data');
    
    % 最新の設定を読み込み
    gui_data = read_current_settings(fig, gui_data);
    
    % ログ初期化
    log_listbox = findobj(parent, 'Tag', 'log_listbox');
    set(log_listbox, 'String', {});
    add_log(log_listbox, 'グリッドサーチを開始します...');
    
    try
        % グリッドサーチ実行
        results = perform_grid_search(gui_data, log_listbox, parent);
        
        % 結果をGUIデータに保存
        setappdata(fig, 'search_results', results);
        
        add_log(log_listbox, 'グリッドサーチが完了しました！');
        
        % 結果の簡単な統計を表示
        success_rate = sum([results.success]) / length(results) * 100;
        add_log(log_listbox, sprintf('成功率: %.1f%% (%d/%d)', success_rate, sum([results.success]), length(results)));
        
    catch me
        add_log(log_listbox, ['エラー: ' me.message]);
        fprintf('グリッドサーチエラー: %s\n', me.message);
    end
end

function gui_data = read_current_settings(fig, gui_data)
    % 現在のGUI設定を読み込み
    
    % Walker パラメータ
    walker_params = {'M', 'm', 'I', 'l', 'w', 'c', 'r', 'g', 'gam'};
    for i = 1:length(walker_params)
        param_name = walker_params{i};
        edit_handle = findobj(fig, 'Tag', ['walker_' param_name]);
        if ~isempty(edit_handle)
            gui_data.walker.(param_name) = str2double(get(edit_handle, 'String'));
        end
    end
    
    % グリッドパラメータ
    grid_vars = {'q1', 'u1', 'q2', 'u2'};
    for i = 1:length(grid_vars)
        var_name = grid_vars{i};
        
        min_edit = findobj(fig, 'Tag', [var_name '_min']);
        max_edit = findobj(fig, 'Tag', [var_name '_max']);
        step_edit = findobj(fig, 'Tag', [var_name '_step']);
        
        if ~isempty(min_edit)
            gui_data.grid.([var_name '_min']) = str2double(get(min_edit, 'String'));
        end
        if ~isempty(max_edit)
            gui_data.grid.([var_name '_max']) = str2double(get(max_edit, 'String'));
        end
        if ~isempty(step_edit)
            gui_data.grid.([var_name '_step']) = str2double(get(step_edit, 'String'));
        end
    end
    
    % 間引き係数
    thin_edit = findobj(fig, 'Tag', 'thin_factor');
    if ~isempty(thin_edit)
        gui_data.grid.thin_factor = str2double(get(thin_edit, 'String'));
    end
    
    % 解析パラメータ
    analysis_params = {'max_time', 'num_cycles', 'success_threshold', 'target_q2'};
    for i = 1:length(analysis_params)
        param_name = analysis_params{i};
        edit_handle = findobj(fig, 'Tag', param_name);
        if ~isempty(edit_handle)
            gui_data.analysis.(param_name) = str2double(get(edit_handle, 'String'));
        end
    end
    
    % チェックボックス設定
    checkbox_params = {'use_parallel', 'verbose_output', 'save_csv', 'show_progress'};
    for i = 1:length(checkbox_params)
        param_name = checkbox_params{i};
        checkbox_handle = findobj(fig, 'Tag', param_name);
        if ~isempty(checkbox_handle)
            gui_data.analysis.(param_name) = get(checkbox_handle, 'Value');
        end
    end
end

function results = perform_grid_search(gui_data, log_listbox, parent)
    % 実際のグリッドサーチ実行（並列計算対応）
    
    % グリッドの生成
    q1_range = gui_data.grid.q1_min:gui_data.grid.q1_step:gui_data.grid.q1_max;
    u1_range = gui_data.grid.u1_min:gui_data.grid.u1_step:gui_data.grid.u1_max;
    q2_range = gui_data.grid.q2_min:gui_data.grid.q2_step:gui_data.grid.q2_max;
    u2_range = gui_data.grid.u2_min:gui_data.grid.u2_step:gui_data.grid.u2_max;
    
    % 全組み合わせ生成
    [Q1, U1, Q2, U2] = ndgrid(q1_range, u1_range, q2_range, u2_range);
    
    total_points = numel(Q1);
    
    % 間引き処理
    if gui_data.grid.thin_factor > 1
        indices = 1:gui_data.grid.thin_factor:total_points;
        Q1 = Q1(indices);
        U1 = U1(indices);
        Q2 = Q2(indices);
        U2 = U2(indices);
        total_points = length(indices);
    end
    
    add_log(log_listbox, sprintf('総計算点数: %d', total_points));
    
    % 並列処理の設定
    if gui_data.analysis.use_parallel
        setup_parallel_pool(log_listbox);
    end
    
    % 結果格納用配列の事前割り当て
    q1_vals = zeros(total_points, 1);
    u1_vals = zeros(total_points, 1);
    q2_vals = zeros(total_points, 1);
    u2_vals = zeros(total_points, 1);
    success_vals = false(total_points, 1);
    final_energy_vals = NaN(total_points, 1);
    max_cycle_vals = zeros(total_points, 1);
    error_messages = cell(total_points, 1);
    
    % 進捗バー設定
    progress_axes = findobj(parent, 'Tag', 'progress_axes');
    if ~isempty(progress_axes)
        cla(progress_axes);
        progress_bar = rectangle(progress_axes, 'Position', [0, 0.2, 0, 0.6], ...
                               'FaceColor', 'blue', 'EdgeColor', 'none');
    end
    
    % 進捗管理用変数
    progress_update_interval = max(1, floor(total_points/100)); % 1%ごとに更新
    
    add_log(log_listbox, sprintf('並列計算を開始します（並列処理: %s）', ...
                                 mat2str(gui_data.analysis.use_parallel)));
    
    % グリッドサーチのメインループ
    if gui_data.analysis.use_parallel
        % === 並列計算版 ===
        add_log(log_listbox, '並列計算プールで実行中...');
        
        % walkerにtarget_q2を追加
        walker_with_target = gui_data.walker;
        walker_with_target.target_q2 = gui_data.analysis.target_q2;
        
        % 並列ループ
        parfor i = 1:total_points
            % 初期条件設定
            z0 = [Q1(i); U1(i); Q2(i); U2(i)];
            
            % シミュレーション実行
            try
                [success, final_energy, max_cycle, error_msg] = simulate_walker_single_parallel(z0, walker_with_target, gui_data.analysis);
                
                % 結果保存（parforでは直接配列に代入）
                q1_vals(i) = Q1(i);
                u1_vals(i) = U1(i);
                q2_vals(i) = Q2(i);
                u2_vals(i) = U2(i);
                success_vals(i) = success;
                final_energy_vals(i) = final_energy;
                max_cycle_vals(i) = max_cycle;
                error_messages{i} = error_msg;
                
            catch me
                % エラー時のデフォルト値
                q1_vals(i) = Q1(i);
                u1_vals(i) = U1(i);
                q2_vals(i) = Q2(i);
                u2_vals(i) = U2(i);
                success_vals(i) = false;
                final_energy_vals(i) = NaN;
                max_cycle_vals(i) = 0;
                error_messages{i} = me.message;
            end
        end
        
        add_log(log_listbox, '並列計算完了！');
        
    else
        % === 逐次計算版（進捗表示付き） ===
        add_log(log_listbox, '逐次計算で実行中...');
        
        % walkerにtarget_q2を追加
        walker_with_target = gui_data.walker;
        walker_with_target.target_q2 = gui_data.analysis.target_q2;
        
        for i = 1:total_points
            % 初期条件設定
            z0 = [Q1(i); U1(i); Q2(i); U2(i)];
            
            % シミュレーション実行
            try
                [success, final_energy, max_cycle, error_msg] = simulate_walker_single_parallel(z0, walker_with_target, gui_data.analysis);
                
                % 結果保存
                q1_vals(i) = Q1(i);
                u1_vals(i) = U1(i);
                q2_vals(i) = Q2(i);
                u2_vals(i) = U2(i);
                success_vals(i) = success;
                final_energy_vals(i) = final_energy;
                max_cycle_vals(i) = max_cycle;
                error_messages{i} = error_msg;
                
            catch me
                % エラー時のデフォルト値
                q1_vals(i) = Q1(i);
                u1_vals(i) = U1(i);
                q2_vals(i) = Q2(i);
                u2_vals(i) = U2(i);
                success_vals(i) = false;
                final_energy_vals(i) = NaN;
                max_cycle_vals(i) = 0;
                error_messages{i} = me.message;
            end
            
            % 進捗更新
            if gui_data.analysis.show_progress && mod(i, progress_update_interval) == 0
                progress = i / total_points;
                if ~isempty(progress_axes) && isvalid(progress_bar)
                    set(progress_bar, 'Position', [0, 0.2, progress, 0.6]);
                    drawnow;
                end
                
                % ログ更新（5%ごと）
                if mod(i, max(1, floor(total_points/20))) == 0
                    progress_percent = i / total_points * 100;
                    success_so_far = sum(success_vals(1:i));
                    add_log(log_listbox, sprintf('進捗: %.1f%% (%d/%d), 成功: %d', ...
                           progress_percent, i, total_points, success_so_far));
                end
            end
        end
        
        add_log(log_listbox, '逐次計算完了！');
    end
    
    % 構造体配列に変換
    results = struct('q1', num2cell(q1_vals), 'u1', num2cell(u1_vals), ...
                     'q2', num2cell(q2_vals), 'u2', num2cell(u2_vals), ...
                     'success', num2cell(success_vals), 'final_energy', num2cell(final_energy_vals), ...
                     'max_cycle', num2cell(max_cycle_vals), 'error_message', error_messages);
    
    % CSV保存
    if gui_data.analysis.save_csv
        save_results_to_csv(results, gui_data);
        add_log(log_listbox, 'Results saved to CSV file.');
    end
end

function setup_parallel_pool(log_listbox)
    % 並列計算プールの設定
    
    try
        % 現在のプールの確認
        current_pool = gcp('nocreate');
        
        if isempty(current_pool)
            add_log(log_listbox, '並列計算プールを開始中...');
            
            % プールを開始（自動的に最適なワーカー数を決定）
            pool = parpool('local');
            add_log(log_listbox, sprintf('並列計算プール開始完了: %d workers', pool.NumWorkers));
            
        else
            add_log(log_listbox, sprintf('既存の並列計算プールを使用: %d workers', current_pool.NumWorkers));
        end
        
    catch me
        add_log(log_listbox, ['並列計算プールエラー: ' me.message]);
        add_log(log_listbox, '逐次計算に切り替えます');
    end
end

function [success, final_energy, max_cycle, error_msg] = simulate_walker_single_parallel(z0, walker, analysis_params)
    % 単一初期条件でのウォーカーシミュレーション（並列計算用）
    % この関数は並列実行に最適化されており、外部変数への依存を最小化
    
    success = false;
    final_energy = NaN;
    max_cycle = 0;
    error_msg = '';
    
    try
        for cycle = 1:analysis_params.num_cycles
            % 1.5ストライドシミュレーション
            [z_traj, t_traj, z_final, z_midpoint] = one_and_half_stride_detailed_parallel(z0, walker);
            
            if isempty(z_traj) || size(z_traj, 1) < 2
                error_msg = sprintf('Failed at cycle %d: trajectory empty', cycle);
                max_cycle = cycle - 1;
                return;
            end
            
            % エネルギー計算
            final_energy = calculate_energy_parallel(z_final, walker);
            
            % 成功判定（最終状態が初期状態に近い）
            deviation = norm(z_final - z0);
            if deviation < analysis_params.success_threshold
                success = true;
                max_cycle = cycle;
                return;
            end
            
            % 次のサイクルのために角速度を軽く調整
            z0 = apply_simple_correction_parallel(z_final, z0);
            max_cycle = cycle;
            
            % 最大時間チェック
            if t_traj(end) > analysis_params.max_time
                error_msg = 'Exceeded maximum time';
                return;
            end
        end
        
        % 全サイクル完了
        success = true;
        
    catch me
        error_msg = me.message;
    end
end

function [z_traj, t_traj, z_final, z_midpoint] = one_and_half_stride_detailed_parallel(z0, walker)
    % 並列計算用の1.5ストライド関数（最適化版）
    
    M = walker.M;  m = walker.m; I = walker.I;   
    l = walker.l;  c = walker.c; w = walker.w;   
    r = walker.r;  g = walker.g; gam = walker.gam;
    
    % target_q2をwalkerから取得
    if isfield(walker, 'target_q2')
        target_q2 = walker.target_q2;
    else
        target_q2 = 0.5;  % デフォルト値
    end
    
    q1 = z0(1); u1 = z0(2); q2 = z0(3); u2 = z0(4);
    
    % エネルギーと位置の計算
    TE = calculate_energy_parallel(z0, walker);
    xp1 = 0;
    xh = -l*sin(q1) - r*q1 + xp1;
    vxh = (-l*cos(q1)-r)*u1; 
    yh =  l*cos(q1) + r;
    vyh = -l*sin(q1)*u1; 
    
    z0_extended = [q1 u1 q2 u2 TE xh vxh yh vyh];
    
    % === 1ステップ目 ===
    t0 = 0; 
    dt = 3; % 最大時間を短縮（並列処理では高速化優先）
    time_stamps = 50; % 解像度を下げて高速化
    
    % Single stance phase
    options = odeset('abstol',1e-10,'reltol',1e-10,'events',@collision);
    tspan = linspace(t0, t0+dt, time_stamps);
    [t_step1, z_step1] = ode113(@(t,z) single_stance(t,z,walker), tspan, z0_extended, options);
    
    % 1ステップ目のHeel strike
    if length(t_step1) < time_stamps
        z_after_collision1 = heelstrike_parallel(t_step1(end), z_step1(end,:), walker);
        z_midpoint = z_after_collision1(1:4);
    else
        error('1ステップ目: 指定時間内に衝突が発生しませんでした');
    end
    
    % === 0.5ステップ目（スイング脚がスタンス脚に追いつくまで） ===
    z1_extended = z_after_collision1;
    t1 = t_step1(end);
    
    % 新しい衝突条件：スイング脚がスタンス脚に追いつく
    options_catchup = odeset('abstol',1e-10,'reltol',1e-10,'events',@swing_catches_stance);
    tspan = linspace(t1, t1+dt, time_stamps);
    [t_step1_5, z_step1_5] = ode113(@(t,z) single_stance(t,z,walker), tspan, z1_extended, options_catchup);
    
    % 時間を調整
    t_step1_5 = t_step1_5 - t1 + t_step1(end);
    
    % 0.5ステップ目の終了
    if length(t_step1_5) < time_stamps
        z_final = z_step1_5(end, 1:4);
    else
        z_final = z_step1_5(end, 1:4);
    end
    
    % === 軌道データの結合 ===
    t_traj = [t_step1; t_step1_5(2:end)];
    z_traj = [z_step1; z_step1_5(2:end,:)];
end

function TE = calculate_energy_parallel(z, walker)
    % 並列計算用のエネルギー計算（最適化版）
    
    q1 = z(1); u1 = z(2); q2 = z(3); u2 = z(4);
    M = walker.M;  m = walker.m; I = walker.I;   
    l = walker.l;  c = walker.c; w = walker.w;   
    r = walker.r;  g = walker.g; gam = walker.gam;
    
    % 運動エネルギー項
    KE1 = 1/2*m*(((-l*cos(q1)-r)*u1-u1*(-c*cos(q1)+w*sin(q1)))^2+(-l*sin(q1)*u1+u1*(c*sin(q1)+w*cos(q1)))^2);
    KE2 = 1/2*m*(((-l*cos(q1)-r)*u1-(u1-u2)*(-c*cos(q1-q2)+w*sin(q1-q2)))^2+(-l*sin(q1)*u1+(u1-u2)*(c*sin(q1-q2)+w*cos(q1-q2)))^2);
    KE3 = 1/2*M*((-l*cos(q1)-r)^2*u1^2+l^2*sin(q1)^2*u1^2);
    KE4 = 1/2*I*(u1^2+(u1-u2)^2);
    
    % 位置エネルギー項
    PE1 = 2*m*g*cos(gam)*r+2*m*g*l*cos(gam-q1)-m*g*c*cos(gam-q1)-m*g*w*sin(gam-q1);
    PE2 = 2*m*g*sin(gam)*r*q1-m*g*c*cos(gam-q1+q2)-m*g*w*sin(gam-q1+q2);
    PE3 = M*g*cos(gam)*r+M*g*l*cos(gam-q1)+M*g*sin(gam)*r*q1;
    
    TE = KE1 + KE2 + KE3 + KE4 + PE1 + PE2 + PE3;
end

function z_corrected = apply_simple_correction_parallel(z_final, z_original)
    % 並列計算用の簡単な角速度補正
    
    energy_factor = 1.02;  % より控えめな補正
    angle_adjustment = 0.005;  % より小さな角度調整
    
    z_corrected = [z_final(1) + angle_adjustment; 
                   z_original(2) * energy_factor;
                   z_final(3) - angle_adjustment; 
                   z_original(4) * energy_factor];
end

function zplus=heelstrike_parallel(t,z,walker)      
    % 並列計算用のヒールストライク（最適化版）
    
    r1 = z(1);   v1 = z(2);                         
    r2 = z(3);   v2 = z(4);                         
    xh = z(6);   yh = z(8);                       
    
    q1 = r1 - r2;                         
    q2 = -r2;                                       
    
    M = walker.M;  m = walker.m; I = walker.I;   
    l = walker.l;  c = walker.c; w = walker.w;   
    r = walker.r;  g = walker.g; gam = walker.gam; 
    
    % 質量行列
    M11 = 2*m*l^2-2*m*l*c+2*m*c^2+2*m*w^2+2*m*r^2+4*m*r*l*cos(q1)-2*m*r*c*cos(q1)+2*m*w*sin(q1)*r-2*m*l*c*cos(q2)-2*m*l*w*sin(q2)-2*m*r*c*cos(q1-q2)+2*m*sin(q1-q2)*w*r+M*l^2+2*M*r*l*cos(q1)+M*r^2+2*I; 
    M12 = m*l*c*cos(q2)+m*l*w*sin(q2)-m*c^2-m*w^2+m*r*c*cos(q1-q2)-m*sin(q1-q2)*w*r-I; 
    M21 = -m*l*c*cos(q2)-m*l*w*sin(q2)+m*c^2+m*w^2-m*r*c*cos(q1-q2)+m*sin(q1-q2)*w*r+I; 
    M22 = -m*w^2-m*c^2-I; 
    
    % 右辺項
    RHS1 = m*l*v2*c-m*c^2*v2+M*v1*r^2-2*m*c*l*v1+M*v1*l^2*cos(r2)+2*m*l^2*v1*cos(r2)+2*I*v1-I*v2+2*m*v1*r^2-2*m*l*v1*c*cos(r2)+2*m*c^2*v1+2*m*w^2*v1-m*w^2*v2+2*m*r*v1*w*sin(r1)+M*r*v1*l*cos(r1)+M*l*cos(-r1+r2)*v1*r-2*m*r*v1*c*cos(r1)+2*m*l*cos(-r1+r2)*v1*r+2*m*r*v1*l*cos(r1)-2*m*r*v1*c*cos(-r1+r2)-2*m*r*v1*w*sin(-r1+r2)+m*r*v2*c*cos(-r1+r2)+m*r*v2*w*sin(-r1+r2); 
    RHS2 = m*r*v1*w*sin(r1)-m*r*v1*c*cos(r1)+I*v1-I*v2+m*w^2*v1-m*c*l*v1+m*c^2*v1; 
    
    MM = [M11 M12; M21 M22];    
    RHS = [RHS1; RHS2];                      
    X = MM \ RHS;                                    
    
    u1 = X(1); u2 = X(2);                                      
    
    % エネルギー計算（簡略化）
    TE = calculate_energy_parallel([q1; u1; q2; u2], walker);
    vxh = (-l*cos(q1)-r)*u1; 
    vyh = -l*sin(q1)*u1; 
    
    zplus = [q1 u1 q2 u2 TE xh vxh yh vyh];                     
end

function save_results_to_csv(results, gui_data)
    % 結果をCSVファイルに保存
    
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    filename = sprintf('grid_search_results_%s.csv', timestamp);
    
    % CSVファイル作成
    fid = fopen(filename, 'w');
    
    % ヘッダー行
    fprintf(fid, 'q1,u1,q2,u2,success,final_energy,max_cycle,error_message\n');
    
    % データ行
    for i = 1:length(results)
        if iscell(results(i).error_message)
            error_msg = results(i).error_message{1};
        else
            error_msg = results(i).error_message;
        end
        
        fprintf(fid, '%.6f,%.6f,%.6f,%.6f,%d,%.6f,%d,"%s"\n', ...
                results(i).q1, results(i).u1, results(i).q2, results(i).u2, ...
                results(i).success, results(i).final_energy, results(i).max_cycle, ...
                error_msg);
    end
    
    fclose(fid);
    
    fprintf('Results saved to: %s\n', filename);
end

function add_log(log_listbox, message)
    % ログメッセージの追加
    
    if isvalid(log_listbox)
        current_logs = get(log_listbox, 'String');
        timestamp = datestr(now, 'HH:MM:SS');
        new_message = sprintf('[%s] %s', timestamp, message);
        
        if ischar(current_logs)
            current_logs = {current_logs};
        elseif isempty(current_logs)
            current_logs = {};
        end
        
        updated_logs = [current_logs; {new_message}];
        set(log_listbox, 'String', updated_logs);
        set(log_listbox, 'Value', length(updated_logs));  % 最新行を表示
        drawnow;
    end
end

function visualize_results(parent)
    % 結果の可視化
    
    fig = ancestor(parent, 'figure');
    results = getappdata(fig, 'search_results');
    
    if isempty(results)
        msgbox('実行結果がありません。まずグリッドサーチを実行してください。', 'Warning');
        return;
    end
    
    % 結果可視化ウィンドウ
    vis_fig = figure('Name', 'グリッドサーチ結果', 'Position', [300 300 1000 700]);
    
    % 成功/失敗の分布
    subplot(2,3,1);
    success_idx = [results.success];
    q1_vals = [results.q1];
    q2_vals = [results.q2];
    
    if any(success_idx)
        scatter(q1_vals(success_idx), q2_vals(success_idx), 'go', 'filled');
        hold on;
    end
    if any(~success_idx)
        scatter(q1_vals(~success_idx), q2_vals(~success_idx), 'rx');
    end
    xlabel('q1 [rad]'); ylabel('q2 [rad]');
    title('成功/失敗分布 (q1-q2)');
    legend('Success', 'Failure', 'Location', 'best');
    grid on;
    
    % エネルギー分布
    subplot(2,3,2);
    final_energies = [results.final_energy];
    valid_energy_idx = ~isnan(final_energies);
    
    if any(valid_energy_idx)
        scatter(q1_vals(valid_energy_idx), final_energies(valid_energy_idx), 'filled');
        xlabel('q1 [rad]'); ylabel('Final Energy');
        title('最終エネルギー分布');
        grid on;
    end
    
    % 成功率統計
    subplot(2,3,3);
    success_rate = sum(success_idx) / length(results) * 100;
    success_count = sum(success_idx);
    failure_count = sum(~success_idx);
    
    if success_count > 0 || failure_count > 0
        pie([success_count, failure_count], {'Success', 'Failure'});
    end
    title(sprintf('成功率: %.1f%%', success_rate));
    
    % その他の統計や詳細プロットを追加可能
    
    % 統計情報テキスト
    subplot(2,3,4);
    stats_text = sprintf('総計算点数: %d\n成功: %d\n失敗: %d\n成功率: %.1f%%', ...
                        length(results), success_count, failure_count, success_rate);
    text(0.1, 0.5, stats_text, 'FontSize', 12);
    xlim([0 1]); ylim([0 1]); axis off;
    title('統計情報');
end

function export_results(parent)
    % 結果のエクスポート
    
    fig = ancestor(parent, 'figure');
    results = getappdata(fig, 'search_results');
    
    if isempty(results)
        msgbox('エクスポートする結果がありません。', 'Warning');
        return;
    end
    
    % ファイル保存ダイアログ
    [filename, pathname] = uiputfile('*.csv', 'Save Results As');
    
    if filename ~= 0
        full_filename = fullfile(pathname, filename);
        
        % CSVエクスポート
        fid = fopen(full_filename, 'w');
        fprintf(fid, 'q1,u1,q2,u2,success,final_energy,max_cycle,error_message\n');
        
        for i = 1:length(results)
            if iscell(results(i).error_message)
                error_msg = results(i).error_message{1};
            else
                error_msg = results(i).error_message;
            end
            
            fprintf(fid, '%.6f,%.6f,%.6f,%.6f,%d,%.6f,%d,"%s"\n', ...
                    results(i).q1, results(i).u1, results(i).q2, results(i).u2, ...
                    results(i).success, results(i).final_energy, results(i).max_cycle, ...
                    error_msg);
        end
        
        fclose(fid);
        msgbox(sprintf('Results exported to: %s', full_filename), 'Export Complete');
    end
end

% ========== 以下、元のシミュレーション関数群 ==========
% multi_cycle_passivewalker_3.m から必要な関数をコピー

function [z_traj, t_traj, z_final, z_midpoint] = one_and_half_stride_detailed(z0, walker)
    % 1.5ストライド（スイング脚が目標角度に到達するまで）の詳細な軌道を返す関数
    
    M = walker.M;  m = walker.m; I = walker.I;   
    l = walker.l;  c = walker.c; w = walker.w;   
    r = walker.r;  g = walker.g; gam = walker.gam;
    
    % target_q2をwalkerから取得
    if isfield(walker, 'target_q2')
        target_q2 = walker.target_q2;
    else
        target_q2 = 0.5;  % デフォルト値
    end
    
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
    t0 = 0; 
    dt = 5; % 1ステップの最大時間
    time_stamps = 100; % 軌道の解像度
    
    % Single stance phase
    options = odeset('abstol',1e-13,'reltol',1e-13,'events',@collision);
    tspan = linspace(t0, t0+dt, time_stamps);
    [t_step1, z_step1] = ode113(@(t,z) single_stance(t,z,walker), tspan, z0_extended, options);
    
    % 1ステップ目のHeel strike
    if length(t_step1) < time_stamps
        z_after_collision1 = heelstrike(t_step1(end), z_step1(end,:), walker);
        z_midpoint = z_after_collision1(1:4);
    else
        error('1ステップ目: 指定時間内に衝突が発生しませんでした');
    end
    
    % === 0.5ステップ目（スイング脚がスタンス脚に追いつくまで） ===
    z1_extended = z_after_collision1;
    t1 = t_step1(end);
    
    % 新しい衝突条件：スイング脚がスタンス脚に追いつく（地面接触なし）
    options_catchup = odeset('abstol',1e-13,'reltol',1e-13,'events',@swing_catches_stance);
    tspan = linspace(t1, t1+dt, time_stamps);
    [t_step1_5, z_step1_5] = ode113(@(t,z) single_stance(t,z,walker), tspan, z1_extended, options_catchup);
    
    % 時間を調整（連続的にする）
    t_step1_5 = t_step1_5 - t1 + t_step1(end);
    
    % 0.5ステップ目の終了（スイング脚追いつき、Heel Strikeなし）
    if length(t_step1_5) < time_stamps
        z_final = z_step1_5(end, 1:4);
    else
        z_final = z_step1_5(end, 1:4);
    end
    
    % === 軌道データの結合 ===
    % 重複する時間点を除去して連結
    t_traj = [t_step1; t_step1_5(2:end)];
    z_traj = [z_step1; z_step1_5(2:end,:)];
end

function TE = calculate_energy(z, walker)
    % 総エネルギーを計算
    
    q1 = z(1); u1 = z(2); q2 = z(3); u2 = z(4);
    M = walker.M;  m = walker.m; I = walker.I;   
    l = walker.l;  c = walker.c; w = walker.w;   
    r = walker.r;  g = walker.g; gam = walker.gam;
    
    TE = 1/2*m*(((-l*cos(q1)-r)*u1-u1*(-c*cos(q1)+w*sin(q1)))^2+(-l*sin(q1)*u1+u1*(c*sin(q1)+w*cos(q1)))^2)+1/2*m*(((-l*cos(q1)-r)*u1-(u1-u2)*(-c*cos(q1-q2)+w*sin(q1-q2)))^2+(-l*sin(q1)*u1+(u1-u2)*(c*sin(q1-q2)+w*cos(q1-q2)))^2)+1/2*M*((-l*cos(q1)-r)^2*u1^2+l^2*sin(q1)^2*u1^2)+1/2*I*(u1^2+(u1-u2)^2)+2*m*g*cos(gam)*r+2*m*g*l*cos(gam-q1)-m*g*c*cos(gam-q1)-m*g*w*sin(gam-q1)+2*m*g*sin(gam)*r*q1-m*g*c*cos(gam-q1+q2)-m*g*w*sin(gam-q1+q2)+M*g*cos(gam)*r+M*g*l*cos(gam-q1)+M*g*sin(gam)*r*q1;
end

function [gstop, isterminal,direction]=swing_catches_stance(t,z,walker)
    % スイング脚が目標角度に到達する条件
    
    q1 = z(1); q2 = z(3); 
    
    % 目標角度を取得（walkerに含まれている場合）
    if isfield(walker, 'target_q2')
        target_q2 = walker.target_q2;
    else
        target_q2 = 0.5;  % デフォルト値
    end
    
    % スイング脚角度が目標値に到達する条件
    gstop = q2 - target_q2;
    isterminal = 1; % 条件に達したら積分を停止
    direction = 0;  % どちらの方向からでも検出
end

function zdot=single_stance(t,z,walker)  
    % 単一スタンス期の微分方程式
    
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

function [gstop, isterminal,direction]=collision(t,z,walker)
    % 衝突検知
    
    q1 = z(1); q2 = z(3); 
    
    gstop = -q2 + 2*q1;
    if (q2>-0.05)
        isterminal = 0;
    else
        isterminal=1;
    end
    direction=-1;
end

function zplus=heelstrike(t,z,walker)      
    % ヒールストライク
    
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

function rotation = R(A)
    % 回転行列
    rotation = [cos(A) -sin(A); sin(A) cos(A)];
end