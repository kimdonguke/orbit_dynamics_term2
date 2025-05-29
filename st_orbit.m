%1st commit
clear;
clc;

data = readmatrix('orbit_data.xlsx');
% 1. 엑셀 파일에서 데이터 읽기

% 2. 위치/속도 분리
[r_list, v_list] = split_state_matrix(data);  % 각각 [M x 3]

% 3. 시뮬레이션 실행
kepler_simulation(r_list, v_list);

function kepler_simulation(r_list, v_list)
    % --- 상수 ---
    mu = 398600.4418;
    Re = 6378.137;

    % --- 위성 초기 궤도 ---
    rp = Re + 250;
    ra = Re + 800;
    a = (rp + ra)/2;
    e = (ra - rp)/(ra + rp);
    i_deg = 53;
    RAAN_deg = 40;
    omega_deg = 60;
    nu0_deg = 0;
    elements = [a, e, i_deg, RAAN_deg, omega_deg, nu0_deg];

    % --- 위성 초기 상태 계산 ---
    [r0_sat, v0_sat] = kepler_to_rv(elements, mu);
    y0 = [r0_sat; v0_sat];

    % --- 시간 설정 ---
    tspan = linspace(0, 86400*10, 30000);
    dt = tspan(2) - tspan(1);   % ← 여기에 정확한 시간 간격 계산
    t0 = tspan(1);
    t_array = tspan(:);         % [N x 1]

    % --- 위성 궤도 전파 ---
    opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
    [t, y] = ode45(@(t,y) two_body_j2_ode(t, y, mu), tspan, y0, opts);
    r_sat_all = y(:,1:3);
    v_sat_all = y(:,4:6);

    % --- 객체 궤도 생성 ---
    object_positions_list = generate_object_positions(r_list, v_list, mu, t_array, t0);

    % --- 접근 이벤트 분석 ---
    M = size(r_list, 1);
    threshold = 100;  % [km]

    for i = 1:M
        r0_obj = r_list(i,:)';
        v0_obj = v_list(i,:)';

        % 위치 및 속도 배열
        r_obj_all = get_position_on_circular_orbit_vec(r0_obj, v0_obj, mu, t_array, t0);
        v_obj_all = get_velocity_on_circular_orbit_vec(r0_obj, v0_obj, mu, t_array, t0);

        % 이벤트 분석
        [match_idx, rel_v, total_area, event_count, event_area_array, ...
         start_times, end_times, durations] = ...
            find_proximity_and_area_with_events( ...
                r_sat_all, v_sat_all, r_obj_all, v_obj_all, t_array, threshold);

        % 출력
        fprintf('▶ 객체 %d\n', i);
        fprintf('  접근 이벤트: %d회\n', event_count);
        fprintf('  총 면적: %.2f km^2\n', total_area);
        fprintf('  각 이벤트 면적: %s\n', mat2str(event_area_array', 4));
    end

    % --- 시각화 ---
    plot_orbit_3D_all(r_sat_all, object_positions_list);
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

function positions = get_position_on_circular_orbit_vec(r0, v0, mu, K_array, dt)
    % get_position_on_circular_orbit_vec
    % ----------------------------------
    % 입력:
    %   r0      : 초기 위치 벡터 [3x1] (km)
    %   v0      : 초기 속도 벡터 [3x1] (km/s)
    %   mu      : 중심체 중력 상수 (km^3/s^2)
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
    theta = omega * (K_array(:) * dt);  % Nx1

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
        object_positions_list{i} = get_position_on_circular_orbit_vec(r0, v0, mu, t_array, t0);
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
