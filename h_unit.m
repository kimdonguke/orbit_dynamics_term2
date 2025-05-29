
% 1. 데이터 불러오기
data = readtable('/content/orbit_states_only.csv');
r = [data.r_x, data.r_y, data.r_z];
v = [data.v_x, data.v_y, data.v_z];

% 2. 각운동량 벡터 h 계산
h = cross(r, v, 2);

% 3. 단위 벡터화 (방향성만 추출)
h_unit = h ./ vecnorm(h, 2, 2);

% 4. 최적의 클러스터 수 k 결정 (실루엣 점수 기반)
silhouette_scores = zeros(1, 19);
k_range = 2:20;
for k = k_range
    [idx, ~] = kmeans(h_unit, k, 'Replicates', 10, 'Display', 'off');
    [s, ~] = silhouette(h_unit, idx);
    silhouette_scores(k-1) = mean(s);
end

% 최적의 k 선택
[~, optimal_k_idx] = max(silhouette_scores);
optimal_k = k_range(optimal_k_idx);
disp(['최적의 클러스터 수 (k): ', num2str(optimal_k)]);

% 5. 최적 k로 클러스터링 수행
[idx, C] = kmeans(h_unit, optimal_k, 'Replicates', 10, 'Display', 'off');

% 6. 가장 빈도가 많은 클러스터의 중심 벡터 추출
counts = histcounts(idx, 1:optimal_k+1);
[~, main_cluster_idx] = max(counts);
main_h_direction = C(main_cluster_idx, :);
main_h_direction = main_h_direction / norm(main_h_direction);

disp('대표 궤도면 방향 벡터 (법선):');
disp(main_h_direction);