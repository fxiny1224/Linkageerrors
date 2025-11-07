clear; clc;

mindist = ...; nsim = 30000; target_dist = ...; r0 = [x0; 0; z0];
nseg = 1;

param_lists = generate_random_params(nsim, nseg);  

param_names = {'L1', 'phi1', 'theta1'};
param_ranges = {linspace(minL, maxL, 15), linspace(min_phi, max_phi, 15), linspace(min_theta, max_theta, 15)};

max_probs = compute_max_prob_values(param_lists, param_names);

param_pairs = nchoosek(1:3, 2);

analysis_results = struct();

for p = 1:size(param_pairs, 1)
    idx1 = param_pairs(p, 1); idx2 = param_pairs(p, 2);
    fixed_idx = setdiff(1:3, [idx1, idx2]);
    fixed_val = max_probs(fixed_idx);
    n_points = ; sampled_range1 = linspace(min(range1), max(range1), n_points);
    range1 = param_ranges{idx1}; range2 = param_ranges{idx2};
    mean_errors = zeros(length(range1), length(range2));
    
    for i = 1:length(range1)
        for j = 1:length(range2)
            params = zeros(3, 1); params(fixed_idx) = fixed_val;
            params(idx1) = convert_to_rad_if_angle(range1(i), idx1);
            params(idx2) = convert_to_rad_if_angle(range2(j), idx2);
            
            temp_errors = zeros(1, 500);
            for k = 1:500
                XYZ_rel = kinkedchain(params);
                epitope_pos = select_random_epitope(nsim);
                XYZ = XYZ_rel + repmat(epitope_pos, 1, size(XYZ_rel, 2));
                offset = compute_offset(XYZ, target_dist, r0);
                XYZ_adj = XYZ + offset;
                dist_adj = compute_dist(XYZ_adj(:, end), r0);
                linkage_error = (dist_adj - ref_dist) / 2;
                temp_errors(k) = linkage_error;
            end
            mean_errors(i, j) = mean(temp_errors);
        end
    end
    
    [X, Y] = meshgrid(range1, range2); Z = mean_errors';
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

summarize_results(analysis_results, param_names, param_pairs);
