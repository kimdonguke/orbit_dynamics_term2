clear;

% 쓰레기 개수
N = 10000;

% 쓰레기 위치할 반지름 [km]
r1 = 6778;
r2 = 7178;
r = rand(N,1) * (r2 - r1) + r1;

% 지구  반지름 [km]
earth_radius = 6378 ;  

% 무작위 각도
theta = 2 * pi * rand(N, 1);      % azimuth (0 to 2pi)
phi = pi/2 + randn(N, 1) * (pi/8);           % 


% 지구 
[xe, ye, ze] = sphere(100);  % 100은 해상도
xe = xe * earth_radius;
ye = ye * earth_radius;
ze = ze * earth_radius;

% 쓰레기 위치
x = r .* sin(phi) .* cos(theta);
y = r .* sin(phi) .* sin(theta);
z = r .* cos(phi);

% K-Means 기법으로 쓰레기 위치 clustering
k = 5;
trash_location = [x y z];
opts = statset('MaxIter', 1000);
[idx,centroid] = kmeans(trash_location, k,'Options',opts);
colormap(hsv(k));


% 시각화
scatter3(x, y, z, 10, idx , 'filled');
caxis([1 k]); % 색상 범위 조정
axis equal;
hold on;
earth_surf = surf(xe, ye, ze);
set(earth_surf, 'FaceColor', [0 0.3 0.3], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
xlabel('X'); ylabel('Y'); zlabel('Z');

plot3(centroid(:,1), centroid(:,2), centroid(:,3), 'kx', 'MarkerSize', 15, 'LineWidth', 5);

[most, frequency] = mode(idx);

perigee = centroid(most,:);