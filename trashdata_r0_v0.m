
T = readtable('trash_elements.xlsx','Sheet','Sheet1');
N= height(T);

mu = 398600.4418;   % gravitational parameter 

%결과 저장용
r_eci = zeros(N,3);
v_eci = zeros(N,3);

for k= 1:N
    %PQW
    e     = T.trash_e(k);      % 이심률
    a     = T.trash_a(k);      % 반장반경 [km]
    nu    = T.trash_nu(k);     % 진근점이각 True Anomaly [rad]
    i     = T.trash_i(k);      % 경사각 [rad]
    RAAN  = T.trash_RAAN(k);    % RAAN [rad]
    w     = T.trash_w(k);      % w[rad]

    p = a * (1 - e^2);
    r = p / (1 + e*cos(nu));

    r_pqw = [r*cos(nu); r*sin(nu); 0];
    v_pqw = [-sqrt(mu/p)*sin(nu); sqrt(mu/p)*(e + cos(nu)); 0];

% PQW -> ECI
R3_RAAN = [cos(RAAN), -sin(RAAN), 0;
        sin(RAAN),  cos(RAAN), 0;
        0,           0,          1];

R1_i = [1, 0, 0;
        0, cos(i), -sin(i);
        0, sin(i),  cos(i)];

R3_w = [cos(w), -sin(w), 0;
        sin(w),  cos(w), 0;
        0,         0,      1];

Q_pqw2eci = R3_RAAN * R1_i * R3_w;

    state_r     = Q_pqw2eci * r_pqw;    % 3×1
    state_v     = Q_pqw2eci * v_pqw;    % 3×1
    r_eci(k,:)  = state_r.';    % 1×3 로 저장
    v_eci(k,:)  = state_v.'; 
end

%-------------------------------------------------------------------

% 3) 원본 테이블에 r,v 컬럼 붙이고 Excel 로 내보내기

T_out = [ T, ...
          array2table(r_eci, 'VariableNames', {'r_x','r_y','r_z'}), ...
          array2table(v_eci, 'VariableNames', {'v_x','v_y','v_z'}) ];

writetable(T_out, 'trash_elements_with_state.xlsx', 'Sheet','Sheet1');

