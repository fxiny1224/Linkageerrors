clear; clc;

mindist = ...; nsim = 30000; target_dist = ...; r0 = [x0; 0; z0];
nseg = 3;

param_lists = generate_random_params(nsim, nseg);  
noise_sigmas = [sigma_fab, sigma_fc, ...];  

param_names = {'L1', 'L2', 'L3', 'phi1', 'phi2', 'phi3', 'theta1', 'theta2', 'theta3'};
mi_matrix = load_mi_matrix(9, 1);  
error_mi = mi_matrix(1:9, end);
param_ranges = [7.0, 7.0, 8.5, 45, 60, 112.6, 45, 360, 360];

max_probs = compute_max_prob_values(param_lists, param_names);

candidate_pairs = [4,5; 5,6; 4,6; 4,7];

analysis_results = struct();

for p = 1:size(candidate_pairs, 1)
    idx1 = candidate_pairs(p, 1); idx2 = candidate_pairs(p, 2);
    
    [range1, is_angle1] = get_param_range_static(idx1, param_ranges);
    [range2, is_angle2] = get_param_range_static(idx2, param_ranges);
    
    n_points = ; sampled_range1 = linspace(min(range1), max(range1), n_points);
    sampled_range2 = linspace(min(range2), max(range2), n_points);
    
    mean_errors = zeros(length(sampled_range1), length(sampled_range2));
    
    for i = 1:length(sampled_range1)
        for j = 1:length(sampled_range2)
            val1 = sampled_range1(i); val2 = sampled_range2(j);
            if is_angle1, val1 = deg2rad(val1); end; if is_angle2, val2 = deg2rad(val2); end;
            
            current_params = create_fixed_params(idx1, val1, idx2, val2, max_probs);
            
            temp_errors = zeros(1, 1000);
            for k = 1:1000
                chain_params = assemble_params(current_params);
                XYZ_rel = kinkedchain(chain_params);
                epitope_pos = select_random_epitope(nsim);
                XYZ = XYZ_rel + repmat(epitope_pos, 1, size(XYZ_rel, 2));
                dist = compute_dist(XYZ(:, end), r0);
                offset = (target_dist - dist) / 2;
                XYZ_adj = XYZ + offset;
                dist_adj = compute_dist(XYZ_adj(:, end), r0);
                linkage_error = (dist_adj - ref_dist) / 2;
                temp_errors(k) = linkage_error;
            end
            mean_errors(i, j) = mean(temp_errors);
        end
    end
    
    [X, Y] = meshgrid(sampled_range1, sampled_range2); Z = mean_errors';
    x_data = X(:); y_data = Y(:); z_data = Z(:);
    valid_idx = ~isnan(z_data); x_data = x_data(valid_idx); y_data = y_data(valid_idx); z_data = z_data(valid_idx);
    
    X_design = [ones(size(x_data)), x_data, y_data, x_data.*y_data, x_data.^2, y_data.^2];
    [beta, ~, ~, ~, stats] = regress(z_data, X_design);
    R2 = stats(1);
    abs_beta = abs(beta(2:end)); total_effect = sum(abs_beta); contributions = abs_beta / total_effect * 100;
    
    analysis_results(p).regression.R2 = R2; analysis_results(p).regression.contributions = contributions;
    analysis_results(p).regression.param_names = {};
    analysis_results(p).regression.beta = beta;
    
    error_range = max(Z(:)) - min(Z(:)); mean_val = mean(Z(:)); min_val = min(Z(:));
    [dx, dy] = gradient(Z); gradient_norm = sqrt(dx.^2 + dy.^2);
    [dxx, dxy] = gradient(dx); [dyx, dyy] = gradient(dy);
    curvature_matrix = zeros(size(Z));
    for m = 1:size(Z, 1)
        for n = 1:size(Z, 2)
            H = [dxx(m,n), dxy(m,n); dyx(m,n), dyy(m,n)];
            curvature_matrix(m,n) = norm(H, 'fro');
        end
    end
    mean_curvature = mean(curvature_matrix(:));
    CSI = (mean_val - min_val) / abs(mean_val); ISR = sum(gradient_norm(:) > quantile(gradient_norm(:), 0.9)) / numel(Z);
    
    analysis_results(p).pair = [param_names{idx1} '-' param_names{idx2}]; analysis_results(p).error_range = error_range;
    analysis_results(p).mean_error = mean_val; analysis_results(p).min_error = min_val; analysis_results(p).mean_curvature = mean_curvature;
    analysis_results(p).CSI = CSI; analysis_results(p).ISR = ISR;
    
    plot_response_surface(X, Y, Z, param_names{idx1}, param_names{idx2}, error_range, CSI);
    plot_regression_fit(x_data, y_data, z_data, beta, param_names{idx1}, param_names{idx2}, R2);
end

summarize_results(analysis_results, param_names, candidate_pairs);
