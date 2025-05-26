%% 1) 파일에서 모든 줄 읽기
fid = fopen('C:\Users\wsdvy\OneDrive\바탕 화면\tlefile.txt','r');
C = textscan(fid,'%s','Delimiter','\n');
fclose(fid);
lines = C{1};

% 위성 당 3줄씩(이름, line1, line2) 이라고 가정
N = numel(lines)/3;

%% 2) 결과 저장용 배열 미리 할당
trash_e    = nan(N,1);
trash_E    = nan(N,1);
trash_i    = nan(N,1);
trash_RAAN = nan(N,1);
trash_w    = nan(N,1);
trash_M    = nan(N,1);
trash_n    = nan(N,1);
trash_nu   = nan(N,1);
trash_a    = nan(N,1);

mu = 398600.4418;  % km^3/s^2

for k = 1:N
    line2 = lines{3*(k-1)+3};

    %--- 파싱 (고정폭 포맷) ---
    i_deg      = str2double( line2(9:16) );          % [deg]
    Omega_deg  = str2double( line2(18:25) );         % RAAN [deg]
    e_str      = line2(27:33);                       % 소수점 생략
    e          = str2double(['0.' e_str]);           % 이심률
    omega_deg  = str2double( line2(35:42) );         % arg. of perigee [deg]
    M_deg      = str2double( line2(44:51) );         % mean anomaly [deg]
    n_rev_day  = str2double( line2(53:63) );         % mean motion [rev/day]

    %--- 라디안 변환 ---
    i_rad     = deg2rad(i_deg);
    Omega_rad = deg2rad(Omega_deg);
    omega_rad = deg2rad(omega_deg);
    M_rad     = deg2rad(M_deg);

    %--- 반장반경 a 계산 ---
    n   = n_rev_day * 2*pi/86400;        % [rad/s]
    a   = (mu / n^2)^(1/3);              % [km]

    %--- Newton–Raphson 으로 Eccentric Anomaly E 구하기 ---
    E = M_rad;
    for iter = 1:50
        f  = E - e*sin(E) - M_rad;
        fp = 1 - e*cos(E);
        dE = -f/fp;
        E  = E + dE;
        if abs(dE) < 1e-12, break; end
    end

    %--- True Anomaly ν ---
    nu = 2*atan2( sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2) );

    %--- 결과 저장 ---
    trash_e(k)    = e;
    trash_E(k)    = E;          % [rad]
    trash_i(k)    = i_rad;      % [rad]
    trash_RAAN(k) = Omega_rad;  % [rad]
    trash_w(k)    = omega_rad;  % [rad]
    trash_M(k)    = M_rad;      % [rad]
    trash_n(k)    = n;          % [rad/s]
    trash_nu(k)   = nu;         % [rad]
    trash_a(k)    = a;          % [km]
end

%% 3) 테이블로 묶어서 Excel로 저장
T = table( trash_e, trash_E, trash_i, trash_RAAN, trash_w, trash_M, trash_n, trash_nu, trash_a, ...
           'VariableNames', {'trash_e','trash_E','trash_i','trash_RAAN','trash_w','trash_M','trash_n','trash_nu','trash_a'} );

writetable(T, 'trash_elements.xlsx', 'Sheet', 'Sheet1');
disp('trash_elements.xlsx 에 궤도 요소가 저장되었습니다.');
