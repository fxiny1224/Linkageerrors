clear; clc;

mindist = ...;  
nsim = 30000;   
target_dist = ...;  
r0 = [x0; 0; z0];  
nseg = 1;

param_lists = generate_random_params(nsim, nseg);  

endpoints = zeros(3, nsim);
distances = zeros(1, nsim);
validconf = true(1, nsim);
for i = 1:nsim
    params_i = param_lists(:, i);  
    XYZ = kinkedchain(params_i) + r0_offset;  
    validconf(i) = ~isclashing(XYZ, mindist);
    if validconf(i)
        distances(i) = compute_dist_to_mt(XYZ(:,end), r0);  
        endpoints(:,i) = adjust_offset(XYZ(:,end), target_dist, distances(1:i));  
    end
end

valid_idx = find(validconf);
linkage_errors = compute_errors(endpoints(:,valid_idx), r0, target_dist);  
mean_err = mean(linkage_errors); std_err = std(linkage_errors);
disp(['mean_std: ' num2str(mean_err) ' ± ' num2str(std_err)]);
plot_distributions(param_lists, linkage_errors);  

param_types = {'l1', 'psi1', 'theta1'};  
for p_type = param_types
    fixed_vals = linspace(min_p, max_p, 5);  
    err_means = zeros(size(fixed_vals)); err_stds = zeros(size(fixed_vals));
    for val = fixed_vals
        param_fixed = generate_fixed_params(nsim, p_type, val, nseg);
        [endpoints_fixed, validconf_fixed] = run_simulation(param_fixed, r0, mindist, target_dist);
        errors_fixed = compute_errors(endpoints_fixed(:,find(validconf_fixed)), r0, target_dist);
        err_means(k) = mean(errors_fixed); err_stds(k) = std(errors_fixed);
    end
    plot_sensitivity(fixed_vals, err_means, err_stds, p_type);  
end
