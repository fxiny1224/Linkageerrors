clear; clc;

mindist = ...; nsim = 30000; target_dist = ...; r0 = [x0; 0; z0];
nseg = 6;  

param_lists = generate_random_params(nsim, nseg);  
noise_sigmas = [sigma_fab, sigma_fc, ...];  

endpoints = zeros(3, nsim); validconf = true(1, nsim);
for i = 1:nsim
    params_i = param_lists(:, i);
    XYZ = kinkedchain(params_i) + r0_offset;  
    XYZ = add_hinge_noise(XYZ, noise_sigmas);  
    validconf(i) = ~isclashing(XYZ, mindist);
    if validconf(i)
        endpoints(:,i) = adjust_radial_offset(XYZ(:,end), target_dist, r0);  
    end
end

valid_idx = find(validconf);
linkage_errors = compute_errors(endpoints(:,valid_idx), r0, target_dist);
mean_err = mean(linkage_errors); std_err = std(linkage_errors);
disp(['mean_std: ' num2str(mean_err) ' ± ' num2str(std_err)]);
plot_distributions(param_lists, linkage_errors);  

param_types = {'l1','l2','l3','l4','ldye', 'psi1','psi2',...};  
for p_type = param_types
    fixed_vals = linspace(min_p, max_p, 5);  
    err_means = zeros(size(fixed_vals)); err_stds = zeros(size(fixed_vals));
    for val = fixed_vals
        param_fixed = generate_fixed_params(nsim, p_type, val, nseg);  
        [endpoints_fixed, validconf_fixed] = run_simulation_with_noise(param_fixed, r0, mindist, target_dist);
        errors_fixed = compute_errors(endpoints_fixed(:,find(validconf_fixed)), r0, target_dist);
        err_means(k) = mean(errors_fixed); err_stds(k) = std(errors_fixed);
    end
    plot_sensitivity(fixed_vals, err_means, err_stds, p_type);  
end
