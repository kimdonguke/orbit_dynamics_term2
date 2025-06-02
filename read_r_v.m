function [r_ECI, v_ECI,main_h_direction] = read_r_v(inFile)
% parse_TLE_to_ECI - TLE 파일에서 ECI 위치·속도 벡터 계산 및 시각화
%
% 입력:
%   inFile - TLE 텍스트 파일 이름 (예: 'FENGYUN 1C.txt')
%
% 출력:
%   r_ECI  - 위치 벡터 (ECI, km)
%   v_ECI  - 속도 벡터 (ECI, km/s)

    %% 상수
    mu  = 398600.4418; % [km^3/s^2] 지구 중력상수
    Re  = 6371;         % [km] 지구 반경

    %% 1) TLE 파일 읽기 및 PQW 계산
    fid = fopen(inFile, 'r');
    assert(fid ~= -1, 'TLE 파일 열기 실패');

    Incl = []; RAAN = []; Ecc = []; Aop = [];
    Mdeg = []; n_rev = []; r_PQW = []; v_PQW = [];

    while ~feof(fid)
        name = strtrim(fgetl(fid));  % 이름 (무시)
        if isempty(name), continue, end
        fgetl(fid);                  % Line 1 (무시)
        L2 = fgetl(fid);             % Line 2 (데이터)

        % 고정 칼럼 파싱
        incl = str2double(L2(9:16));        
        raan = str2double(L2(18:25));       
        ecc  = str2double(['0.' L2(27:33)]);
        aop  = str2double(L2(35:42));       
        M0   = str2double(L2(44:51));       
        n0   = str2double(L2(53:63));       

        M     = deg2rad(M0);                
        n_rad = n0 * 2*pi/86400;            

        % Newton-Raphson으로 E 계산
        E_prev = M;
        dE = E_prev;
        while abs(dE) > 0.001
            E = M + ecc * sin(E_prev);
            dE = E - E_prev;
            E_prev = E;
        end

        % PQW 요소 계산
        nu   = atan2(sqrt(1 - ecc^2) * sin(E), cos(E) - ecc);
        a    = (mu / n_rad^2)^(1/3);
        p    = a * (1 - ecc^2);
        rmag = p / (1 + ecc * cos(nu));
        rPQW = [rmag * cos(nu), rmag * sin(nu), 0];
        vPQW = [-sqrt(mu/p) * sin(nu), sqrt(mu/p) * (ecc + cos(nu)), 0];

        % 누적
        Incl(end+1) = incl; RAAN(end+1) = raan;
        Ecc(end+1)  = ecc;  Aop(end+1)  = aop;
        Mdeg(end+1) = M0;   n_rev(end+1) = n0;
        r_PQW(end+1, :) = rPQW; v_PQW(end+1, :) = vPQW;
    end
    fclose(fid);

    %% 2) PQW → ECI 변환
    r_ECI = zeros(size(r_PQW));
    v_ECI = zeros(size(v_PQW));
    for i = 1:numel(Incl)
        R = DCM(deg2rad(RAAN(i)), deg2rad(Aop(i)), deg2rad(Incl(i)));
        r_ECI(i,:) = (R * r_PQW(i,:).').';
        v_ECI(i,:) = (R * v_PQW(i,:).').';
    end
    h = cross(r_ECI, v_ECI, 2);

% 3. 단위 벡터화 (방향성만 추출)
h_unit = h ./ vecnorm(h, 2, 2);

% 4. 최적의 클러스터 수 k 결정 (실루엣 점수 기반)
k = 20;

% [idx, ~] = kmeans(h_unit, k, 'Replicates', 10, 'Display', 'off');
% [s, ~] = silhouette(h_unit, idx);
    

% 5. 최적 k로 클러스터링 수행
[idx, C] = kmeans(h_unit, k, 'Replicates', 10, 'Display', 'off');

% 6. 가장 빈도가 많은 클러스터의 중심 벡터 추출
counts = histcounts(idx, 1:21);
[~, main_cluster_idx] = max(counts);
main_h_direction = C(main_cluster_idx, :);
main_h_direction = main_h_direction / norm(main_h_direction);


end
