%% ──────────────────────────────────────────────────────────────
%  0) 기본 상수·파일 이름
% ───────────────────────────────────────────────────────────────
clear; clc;
mu      = 398600.4418;         % [km^3/s^2] 지구 중력상수
inFile  = 'FENGYUN 1C.txt';    % TLE 원본
outFile = 'TLE_ECI_2.xlsx';      % 결과 엑셀
Re      = 6371;                % [km] 지구 반경

%% ──────────────────────────────────────────────────────────────
%  1) TLE 읽어 오면서 궤도요소, r_PQW, v_PQW 계산
% ───────────────────────────────────────────────────────────────
fid = fopen(inFile,'r');  
assert(fid~=-1,'TLE 열기 실패');

Names={}; Incl=[]; RAAN=[]; Ecc=[]; Aop=[];
Mdeg=[]; n_rev=[]; nu_deg=[]; a_km=[]; p_km=[];
r_PQW=[]; v_PQW=[];

while ~feof(fid)
    name = strtrim(fgetl(fid));        % Line-0 (이름)
    if isempty(name),   continue, end
    fgetl(fid);                         % Line-1 (unused)
    L2 = fgetl(fid);                    % Line-2 (데이터)

    %── 1-1) 고정 칼럼 파싱 ──────────────────────────────
    incl = str2double(L2(9:16));          % [deg]
    raan = str2double(L2(18:25));         % [deg]
    ecc  = str2double(['0.' L2(27:33)]);
    aop  = str2double(L2(35:42));         % [deg]
    M0   = str2double(L2(44:51));         % [deg]
    n0   = str2double(L2(53:63));         % [rev/day]

    %── 1-2) 단위 변환 및 뉴턴-랩슨으로 E 찾기 ───────────
    M     = deg2rad(M0);                 % [rad]
    n_rad = n0 * 2*pi/86400;             % [rad/s]


    %여기야!
    E_prev = M;
    dE = E_prev;
    while (dE > 0.001)
        E = M + ecc*sin(E_prev);
        E_prev = E;
        dE = (E - E_prev);
    end

    %── 1-3) True anomaly, a, p, r, v (PQW) ──────────────
    nu   = atan2( sqrt(1-ecc^2)*sin(E), cos(E)-ecc );   % [rad]
    a    = (mu/n_rad^2)^(1/3);                          % [km]
    p    = a*(1 - ecc^2);                               % [km]
    rmag = p/(1 + ecc*cos(nu));                         % [km]
    rPQW = [ rmag*cos(nu)  ,  rmag*sin(nu)  , 0 ];
    vPQW = [-sqrt(mu/p)*sin(nu) ,  sqrt(mu/p)*(ecc+cos(nu)) , 0];

    %── 1-4) 결과 누적 ───────────────────────────────────
    Names{end+1,1} = name;
    Incl(end+1)    = incl;   RAAN(end+1) = raan;
    Ecc(end+1)     = ecc;    Aop(end+1)  = aop;
    Mdeg(end+1)    = M0;     n_rev(end+1)= n0;
    nu_deg(end+1)  = rad2deg(nu);
    a_km(end+1)    = a;      p_km(end+1) = p;
    r_PQW(end+1,:) = rPQW;   v_PQW(end+1,:)=vPQW;
end
fclose(fid);

%% ──────────────────────────────────────────────────────────────
%  2) PQW → ECI 변환 (DCM 함수 사용)
% ───────────────────────────────────────────────────────────────
r_ECI = zeros(size(r_PQW));
v_ECI = zeros(size(v_PQW));
for i = 1:numel(Incl)
    R   = DCM( deg2rad(RAAN(i)), deg2rad(Aop(i)), deg2rad(Incl(i)) );
    r_ECI(i,:) = (R * r_PQW(i,:).').';
    v_ECI(i,:) = (R * v_PQW(i,:).').';
end

%% ──────────────────────────────────────────────────────────────
%  3) 엑셀 저장 (6-parameter + r,v 모두 포함)
% ───────────────────────────────────────────────────────────────
T = table(Names, Incl.', RAAN.', Ecc.', Aop.', Mdeg.', n_rev.',...
          nu_deg.', a_km.', p_km.', ...
          r_ECI(:,1), r_ECI(:,2), r_ECI(:,3), ...
          v_ECI(:,1), v_ECI(:,2), v_ECI(:,3), ...
      'VariableNames',{'Name','Incl_deg','RAAN_deg','Ecc',...
      'ArgPer_deg','MeanAno_deg','MeanMot_revday','TrueAno_deg',...
      'a_km','p_km','x_ECI_km','y_ECI_km','z_ECI_km',...
      'vx_ECI_kms','vy_ECI_kms','vz_ECI_kms'});

writetable(T,outFile,'FileType','spreadsheet');
fprintf("✔ %d개 객체의 궤도요소·ECI r,v를 '%s'에 저장!\n",height(T),outFile);

%% ──────────────────────────────────────────────────────────────
%  4) plotting
% ───────────────────────────────────────────────────────────────
figure('Color','w'); hold on; grid on; axis equal
title('Earth + 10 cm Space-Debris Distribution (ECI)')
xlabel('X  [km]'); ylabel('Y  [km]'); zlabel('Z  [km]')

% 4-1) 지구 
[xe,ye,ze] = sphere(72);
surf(Re*xe, Re*ye, Re*ze, ...
    'EdgeColor','none','FaceAlpha',0.25,'FaceColor',[0.2 0.6 1]);

% 4-2) 우주쓰레기 plotting
scatter3(r_ECI(:,1), r_ECI(:,2), r_ECI(:,3), ...
         8, 'filled', 'MarkerFaceColor',[1 0 0]);

view(35,22);   
