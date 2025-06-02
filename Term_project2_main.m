%1st commit
clc;
clear;
close all;

[r_ECI, v_ECI,main_h_direction] = read_r_v('FENGYUN 1C.txt');

trash_location = [r_ECI(:,1), r_ECI(:,2), r_ECI(:,3)];
trash_velocity = [v_ECI(:,1), v_ECI(:,2), v_ECI(:,3)];
% most, h_unit = [0.3048 0.9402 -0.1522]

main_h_direction = [0.3048, 0.9402, -0.1522];

[n, i, RAAN] = h_to_n_i_RAAN(main_h_direction);

%%
y=kepler_simulation(r_ECI,v_ECI);

function y=kepler_simulation(r_list, v_list)
    % --- 상수 ---
    mu = 398600.4418;
    Re = 6378.137;

    % --- 위성 초기 궤도 ---
    rp = Re + 250;
    ra = Re + 800;
    a = (rp + ra)/2;
    e = (ra - rp)/(ra + rp);
    i_deg = 98.7542;
    RAAN_deg = 162.0380;  
    omega_deg = 14;
    nu0_deg = 0;
    elements = [a, e, i_deg, RAAN_deg, omega_deg, nu0_deg];

    % --- 위성 초기 상태 계산 ---
    [r0_sat, v0_sat] = kepler_to_rv(elements, mu);
    y0 = [r0_sat; v0_sat];
    

    k=3;
    tspan = linspace(0, 86400*k, k*8640);
    dt = tspan(2) - tspan(1);   % ← 여기에 정확한 시간 간격 계산
    t0 = tspan(1);
    t_array = tspan(:);         % [N x 1]
    N = length(t_array);                     % 총 시점 수
    K_array = (0:N-1)';                      % 인덱스 배열

    
    % --- 위성 궤도 전파 ---
    opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
    [t, y] = ode45(@(t,y) two_body_j2_ode(t, y, mu), tspan, y0, opts);
    r_sat_all = y(:,1:3);
    v_sat_all = y(:,4:6);

    % --- 객체 궤도 생성 ---
    object_positions_list = generate_object_positions(r_list, v_list, mu, t_array, t0);

    % --- 접근 이벤트 분석 ---
    M = size(r_list, 1);
    threshold = 50;  % [km]

    a=0;

    for i = 1:M
        r0_obj = r_list(i,:)';
        v0_obj = v_list(i,:)';

        % 위치 및 속도 배열
        r_obj_all = get_position_on_circular_orbit_vec(r0_obj, v0_obj, K_array, dt);
        v_obj_all = get_velocity_on_circular_orbit_vec(r0_obj, v0_obj, mu, K_array, dt);

        % 이벤트 분석
        [match_idx, rel_v, total_area, event_count, event_area_array, ...
         start_times, end_times, durations] = ...
            find_proximity_and_area_with_events( ...
                r_sat_all, v_sat_all, r_obj_all, v_obj_all, t_array, threshold);
        % 이벤트 마커 저장
        event_markers_list{i} = r_obj_all(match_idx, :);

        if event_count >= 1
            a=a+1;
        end
       
        
        % 출력
        fprintf('▶ 객체 %d\n', i);
        fprintf('  접근 이벤트: %d회\n', event_count);
        fprintf('  총 면적: %.2f km^2\n', total_area);
        fprintf('  각 이벤트 면적: %s\n', mat2str(event_area_array', 4));
    end
    
    disp(a);

    % --- 시각화 ---
    plot_orbit_3D_all(r_sat_all, object_positions_list,event_markers_list);
end

% ------------------------ 서브 함수들 ------------------------

function dydt = two_body_j2_ode(~, y, mu)
    r_vec = y(1:3);
    v_vec = y(4:6);
    x = r_vec(1); y1 = r_vec(2); z = r_vec(3);
    r = norm(r_vec);
    
    a_grav = -mu / r^3 * r_vec;
    a_total = a_grav;
    dydt = [v_vec; a_total];
end

function [r_eci, v_eci] = kepler_to_rv(elements, mu)
    a = elements(1);
    e = elements(2);
    i = deg2rad(elements(3));
    RAAN = deg2rad(elements(4));
    omega = deg2rad(elements(5));
    nu = deg2rad(elements(6));

    p = a * (1 - e^2);
    r_pqw = (p / (1 + e*cos(nu))) * [cos(nu); sin(nu); 0];
    v_pqw = sqrt(mu/p) * [-sin(nu); e + cos(nu); 0];
    

    R3_Omega = [cos(-RAAN), -sin(-RAAN), 0;
                sin(-RAAN),  cos(-RAAN), 0;
                0,           0,          1];

    R1_i = [1, 0, 0;
            0, cos(-i), -sin(-i);
            0, sin(-i),  cos(-i)];

    R3_omega = [cos(-omega), -sin(-omega), 0;
                sin(-omega),  cos(-omega), 0;
                0,            0,           1];

    Q = R3_Omega * R1_i * R3_omega;
   
    r_eci = Q * r_pqw;
    v_eci = Q * v_pqw;
end

function [a, e, i, RAAN, omega, nu] = rv_to_elements(r, v, mu)
    h = cross(r, v);
    n = cross([0;0;1], h);
    e_vec = (1/mu) * (cross(v, h) - mu*r/norm(r));
    e = norm(e_vec);
    i = acos(h(3)/norm(h));
    RAAN = atan2(n(2), n(1));
    omega = atan2(dot(cross(n, e_vec), h)/norm(h), dot(n, e_vec));
    nu = atan2(dot(cross(e_vec, r), h)/norm(h), dot(e_vec, r));
    a = 1 / (2/norm(r) - norm(v)^2/mu);
end

function plot_orbit_3D(t, r)
    figure;
    plot3(r(:,1), r(:,2), r(:,3), 'b', 'LineWidth', 1.2);
    hold on;

    % 지구 그리기 (선택)
    [X,Y,Z] = sphere(50);
    Re = 6378.137; % 지구 반지름 [km]
    surf(Re*X, Re*Y, Re*Z, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', [0.4 0.7 1]);
    
    % 시작점 표시
    plot3(r(1,1), r(1,2), r(1,3), 'go', 'MarkerFaceColor', 'g', 'DisplayName','Start');
    % 종료점 표시
    plot3(r(end,1), r(end,2), r(end,3), 'ro', 'MarkerFaceColor', 'r', 'DisplayName','End');

    % 설정
    grid on;
    axis equal;
    xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
    title('3D Orbit in ECI Frame');
    legend show;
    view(30, 30); % 시각적 각도 조절
end

function positions = get_position_on_circular_orbit_vec(r0, v0, t_array, t0)
    % get_position_on_circular_orbit_vec
    % ----------------------------------
    % 입력:
    %   r0      : 초기 위치 벡터 [3x1] (km)
    %   v0      : 초기 속도 벡터 [3x1] (km/s)
    %   K_array : 시점 인덱스 배열 [1xN] 또는 [Nx1]
    %   dt      : 시간 간격 (초)
    %
    % 출력:
    %   positions : 각 시점의 위치 [Nx3] 배열 (km)

    % 1. 초기 단위 벡터들
    r_mag = norm(r0);
    r_hat = r0 / r_mag;

    h_vec = cross(r0, v0);
    h_hat = h_vec / norm(h_vec);
    t_hat = cross(h_hat, r_hat);



    % 2. 각속도
    v_mag = norm(v0);
    omega = v_mag / r_mag;

    % 3. 회전각 배열
    theta = omega * (t_array(:) - t0);  % Nx1

    % 4. 위치 계산 (벡터화)
    cos_theta = cos(theta);  % Nx1
    sin_theta = sin(theta);  % Nx1

    positions = r_mag * (cos_theta .* r_hat' + sin_theta .* t_hat');  % Nx3
end

function velocities = get_velocity_on_circular_orbit_vec(r0, v0, mu, K_array, dt)
    % get_velocity_on_circular_orbit_vec
    % ----------------------------------
    % 입력:
    %   r0      : 초기 위치 벡터 [3x1] (km)
    %   v0      : 초기 속도 벡터 [3x1] (km/s)
    %   mu      : 중심체 중력 상수 (km^3/s^2)
    %   K_array : 시점 인덱스 배열 [1xN] 또는 [Nx1]
    %   dt      : 시간 간격 (초)
    %
    % 출력:
    %   velocities : 각 시점의 속도 [Nx3] 배열 (km/s)

    % 1. 초기 단위 벡터들
    r_mag = norm(r0);
    r_hat = r0 / r_mag;

    h_vec = cross(r0, v0);
    h_hat = h_vec / norm(h_vec);
    t_hat = cross(h_hat, r_hat);  % 초기 속도 방향

    % 2. 속도 크기 (원운동)
    v_mag = sqrt(mu / r_mag);

    % 3. 회전각
    theta = v_mag / r_mag * (K_array(:) * dt);  % Nx1

    % 4. 속도 계산 (회전 후 접선 방향)
    % v(t) = -sin(theta)*r̂ + cos(theta)*t̂
    sin_theta = sin(theta);
    cos_theta = cos(theta);

    velocities = v_mag * (-sin_theta .* r_hat' + cos_theta .* t_hat');  % Nx3
end


function find_proximity_and_area(r_sat_all, v_sat_all, ...
                             r_obj_all, v_obj_all, ...
                             threshold, dt)
    % 입력:
    %   r_sat_all : 위성 위치 [N x 3]
    %   v_sat_all : 위성 속도 [N x 3]
    %   r_obj_all : 객체 위치 [N x 3]
    %   v_obj_all : 객체 속도 [N x 3]
    %   threshold : 거리 임계값 [km]
    %   dt        : 시간 간격 [s]

    % 거리 계산 (벡터화)
    d = vecnorm(r_sat_all - r_obj_all, 2, 2);  % N x 1

    % 조건 만족하는 인덱스
    matching_indices = find(d < threshold);

    % 상대 속도 계산 (Nx3 → Mx3)
    v_rel = v_sat_all(matching_indices,:) - v_obj_all(matching_indices,:);
    relative_velocities = vecnorm(v_rel, 2, 2);  % M x 1

    % 면적 겹침 구간 찾기
    inside = d < threshold;  % N x 1 logical
    shift_inside = [false; inside(1:end-1)];
    event_start = find(~shift_inside & inside);  % 0→1
    event_end   = find(shift_inside & ~inside);  % 1→0

    % 예외 처리: 열려 있는 끝단 처리
    if ~isempty(event_start) && (isempty(event_end) || event_end(1) < event_start(1))
        event_end = [event_end; length(inside)];  % 마지막까지 붙어있을 경우
    end
    if length(event_end) > length(event_start)
        event_end(end) = [];  % 이상한 끼어듦 방지
    end

    % 각 이벤트별 면적 계산
    event_count = length(event_start);
    event_area_array = zeros(event_count, 1);

    for i = 1:event_count
        k1 = event_start(i);
        k2 = event_end(i);
        v_rel_seg = v_sat_all(k1:k2,:) - v_obj_all(k1:k2,:);
        rel_speed = vecnorm(v_rel_seg, 2, 2);
        seg_dist  = d(k1:k2);
        event_area_array(i) = sum(rel_speed .* seg_dist * dt);
    end

    return;
end

function [r_list, v_list] = split_state_matrix(data)
    % split_state_matrix
    % -------------------
    % 입력:
    %   data: M x 6 행렬 [r_x, r_y, r_z, v_x, v_y, v_z]
    %
    % 출력:
    %   r_list: M x 3 위치 행렬
    %   v_list: M x 3 속도 행렬

    if size(data,2) ~= 6
        error('입력 행렬은 반드시 6열이어야 합니다. [r_x, r_y, r_z, v_x, v_y, v_z]');
    end

    r_list = data(:, 1:3);
    v_list = data(:, 4:6);
    return;
end

function object_positions_list = generate_object_positions(r_list, v_list, mu, t_array, t0)
    M = size(r_list, 1);
    object_positions_list = cell(M, 1);
    for i = 1:M
        r0 = r_list(i,:)';
        v0 = v_list(i,:)';
        object_positions_list{i} = get_position_on_circular_orbit_vec(r0, v0, t_array, t0);
    end
end

function [match_idx, rel_v, total_area, event_count, event_area_array, ...
          start_times, end_times, durations] = ...
    find_proximity_and_area_with_events(r_sat, v_sat, r_obj, v_obj, t_array, threshold)

    d = vecnorm(r_sat - r_obj, 2, 2);
    inside = d < threshold;
    shifted = [false; inside(1:end-1)];
    event_start = find(~shifted & inside);
    event_end   = find(shifted & ~inside);

    if ~isempty(event_start) && (isempty(event_end) || event_end(1) < event_start(1))
        event_end = [event_end; length(inside)];
    end
    if length(event_end) > length(event_start)
        event_end(end) = [];
    end

    % 면적 계산
    event_count = length(event_start);
    event_area_array = zeros(event_count,1);
    for i = 1:event_count
        k1 = event_start(i); k2 = event_end(i);
        v_rel = v_sat(k1:k2,:) - v_obj(k1:k2,:);
        rel_speed = vecnorm(v_rel, 2, 2);
        event_area_array(i) = sum(rel_speed .* d(k1:k2) * (t_array(2)-t_array(1)));
    end

    total_area = sum(event_area_array);
    match_idx = find(inside);
    rel_v = vecnorm(v_sat(inside,:) - v_obj(inside,:), 2, 2);
    start_times = t_array(event_start);
    end_times = t_array(event_end);
    durations = end_times - start_times;
end

function plot_orbit_3D_all(r_sat_all, object_positions_list, event_markers_list)
    % plot_orbit_3D_all
    % ------------------
    % 위성 + 객체 궤도 + 지구를 함께 3D 시각화
    % 사용자 입력으로 표시 여부 제어

    if nargin < 3
        event_markers_list = {};  % 이벤트 마커가 없다면 빈 셀
    end

    % 사용자 입력
    show_objects = input('객체 궤도를 표시할까요? (1 = 예, 0 = 아니오): ');
    show_events  = input('접근 이벤트 지점을 표시할까요? (1 = 예, 0 = 아니오): ');

    % 시각화 준비
    figure;
    hold on;
    grid on;
    axis equal;

    % 1. 지구 표시
    Re = 6378.137;
    [X, Y, Z] = sphere(50);
    surf(Re*X, Re*Y, Re*Z, ...
         'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', [0.4 0.7 1]);

    % 2. 위성 궤도
    plot3(r_sat_all(:,1), r_sat_all(:,2), r_sat_all(:,3), ...
          'b-', 'LineWidth', 2.0, 'DisplayName', 'Satellite');

    % 위성 시작점
    plot3(r_sat_all(1,1), r_sat_all(1,2), r_sat_all(1,3), ...
          'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'Start (Satellite)');

    % 3. 객체 궤도 (조건부)
    if show_objects
        M = length(object_positions_list);
        cmap = lines(M);
        for i = 1:M
            r_obj = object_positions_list{i};

            % 궤도
            plot3(r_obj(:,1), r_obj(:,2), r_obj(:,3), ...
                  'Color', cmap(i,:), 'LineWidth', 1.2);

            % 시작점 마커
            %plot3(r_obj(1,1), r_obj(1,2), r_obj(1,3), ...
             %     'o', 'MarkerEdgeColor', cmap(i,:), ...
              %    'MarkerFaceColor', cmap(i,:), ...
               %   'HandleVisibility', 'off');
        end
    end

    % 4. 접근 이벤트 마커 (조건부)
    if show_events && ~isempty(event_markers_list)
        for i = 1:length(event_markers_list)
            event_pts = event_markers_list{i};
            if ~isempty(event_pts)
                scatter3(event_pts(:,1), event_pts(:,2), event_pts(:,3), ...
                         25, 'r', 'filled');
            end
        end
    end
   

    % 5. 라벨 및 보기
    xlabel('X [km]');
    ylabel('Y [km]');
    zlabel('Z [km]');
    title('3D Orbits of Satellite and Objects');
    view(30, 30);
    
end

Re = 6378.137;
mu = 398600.4418;
% Debris 새로운 타원궤도의 목표 근지점 Re+120 [km]
r_p = Re + 120;
r_a = Re + 250;
a = (r_a + r_p)/2;
e = (r_a - r_p)/(r_a + r_p);
i = rad2deg(i);
RAAN = rad2deg(RAAN);
nu = 180;
omega = 14;
p = a * (1 - e^2);
v_sat = y(1,4:6);
v_sat_mag = norm(y(1,4:6));
v_new_scalar = sqrt(mu*(2/r_a-1/a));
del_v_mag = v_new_scalar - v_sat_mag;
v_unit = v_sat/v_sat_mag;
% Debris 새로운 타원궤도의 원지점에서 위치, 속도
v_debris = v_sat + del_v_mag .* v_unit;
r_debris = y(1,1:3);

tspan = linspace(0, 86400*3, 3*8640);
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);
f = @(t, X) [X(4); X(5); X(6); ...
 -mu*X(1)/norm(X(1:3))^3; ...
 -mu*X(2)/norm(X(1:3))^3; ...
 -mu*X(3)/norm(X(1:3))^3];
% f = [v_x, v_y, v_z, a_x, a_y, a_z]'
X0 = [r_debris, v_debris];
[t, X] = ode45(f, tspan, X0, options);

% 근지점 거리 계산 및 찾기
distances = sqrt(X(:,1).^2 + X(:,2).^2 + X(:,3).^2);
perigee_target = Re + 120;
tolerance = 5; % 5km 허용 오차

% 근지점에 도달하는 첫 번째 인덱스 찾기 (시작점 이후)
perigee_indices = find(abs(distances - perigee_target) < tolerance);
% 시작점(첫 번째 점) 제외하고 두 번째 근지점 찾기
if length(perigee_indices) > 1
    stop_index = perigee_indices(2);
else
    % 근지점을 찾지 못한 경우 최소 거리 지점 사용
    [~, stop_index] = min(distances(100:end)); % 처음 100개 점 제외
    stop_index = stop_index + 99;
end

%% 애니메이션 설정
figure;
hold on;
grid on;
axis equal;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('Orbital Animation: Debris and Satellite');

% 지구 그리기
[l, Y, Z] = sphere(50);
surf(Re*l, Re*Y, Re*Z, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', [0.4 0.7 1]);

% 위성 전체 궤도 미리 그리기 (얇은 선)
plot3(y(:,1), y(:,2), y(:,3), 'b-', 'LineWidth', 1.0);

% 초기 점들 설정
debris_line = plot3(X(1,1), X(1,2), X(1,3), 'g-', 'LineWidth', 2.0);
debris_point = plot3(X(1,1), X(1,2), X(1,3), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
sat_point = plot3(y(1,1), y(1,2), y(1,3), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');

% 시작점 표시
plot3(y(1,1), y(1,2), y(1,3), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'Start Point');


view(3);

% 축 범위 설정
max_range = max([max(abs(X(:,1:3)), [], 'all'), max(abs(y(:,1:3)), [], 'all')]) * 1.1;
xlim([-max_range, max_range]);
ylim([-max_range, max_range]);
zlim([-max_range, max_range]);

%% 동시 애니메이션 설정
fprintf('Starting simultaneous animation of satellite and debris...\n');

% 시간 동기화를 위한 설정
debris_x = [];
debris_y = [];
debris_z = [];
sat_x = [];
sat_y = [];
sat_z = [];

% 디브리스는 stop_index까지만 진행
debris_frames = stop_index;

% 위성은 더 빠른 속도로 여러 바퀴 돌게 설정(위성의 근지점, debris의 원지점 속도 비율로 배속)
sat_speed_multiplier = norm(y(1,4:6))/v_sat_mag;
max_frames = debris_frames * 3;

% 디브리스 궤적 선을 위한 추가 변수
debris_trail = plot3(X(1,1), X(1,2), X(1,3), 'g-', 'LineWidth', 1.5);
sat_trail = plot3(y(1,1), y(1,2), y(1,3), 'b-', 'LineWidth', 1.5);

% 디브리스 상태 플래그
debris_stopped = false;
debris_stop_frame = 0;

for frame = 1:max_frames
    % 디브리스 인덱스 계산 (stop_index까지만)
    if ~debris_stopped
        debris_idx = min(frame, stop_index);
        if debris_idx >= stop_index
            debris_stopped = true;
            debris_stop_frame = frame;           
        end
    else
        debris_idx = stop_index; % 멈춘 상태 유지
    end
    
    % 위성 인덱스 계산 (빠른 속도로 순환)
    sat_position = mod((frame * sat_speed_multiplier - 1), length(y)) + 1;
    sat_idx = round(sat_position);
    
    % 디브리스 궤적 누적 (멈추지 않았을 때만)
    if ~debris_stopped || frame == debris_stop_frame
        debris_x = [debris_x, X(debris_idx,1)];
        debris_y = [debris_y, X(debris_idx,2)];
        debris_z = [debris_z, X(debris_idx,3)];
        set(debris_trail, 'XData', debris_x, 'YData', debris_y, 'ZData', debris_z);
    end
    
    % 위성 궤적 누적 (계속 돌면서 궤적 갱신)
    sat_x = [sat_x, y(sat_idx,1)];
    sat_y = [sat_y, y(sat_idx,2)];
    sat_z = [sat_z, y(sat_idx,3)];
    
    % 위성 궤적이 너무 길어지면 일부 제거 (최근 궤도만 표시)
    trail_length = 200; % 표시할 궤적 길이
    if length(sat_x) > trail_length
        sat_x = sat_x(end-trail_length+1:end);
        sat_y = sat_y(end-trail_length+1:end);
        sat_z = sat_z(end-trail_length+1:end);
    end
    
    % 궤적 선 업데이트
    set(sat_trail, 'XData', sat_x, 'YData', sat_y, 'ZData', sat_z);
    
    % 현재 위치 점 업데이트
    set(debris_point, 'XData', X(debris_idx,1), 'YData', X(debris_idx,2), 'ZData', X(debris_idx,3));
    set(sat_point, 'XData', y(sat_idx,1), 'YData', y(sat_idx,2), 'ZData', y(sat_idx,3));
    
    drawnow;
    pause(0.05); % 애니메이션 속도 조절
end
