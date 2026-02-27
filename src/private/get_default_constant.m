function constant_value = get_default_constant(constant_name)
switch constant_name
    case {"MaxFunctionEvaluations_dim_factor"}
        constant_value = 500;
    case {"ftarget"}
        constant_value = -inf;
    case {"StepTolerance"}
        constant_value = 1e-6;
    case {"use_function_value_stop"}
        constant_value = false;
    case {"func_window_size"}
        constant_value = 20;
    case {"func_tol"}
        constant_value = 1e-6;
    case {"use_estimated_gradient_stop"}
        constant_value = false;
    case {"grad_window_size"}
        constant_value = 1;
    case {"grad_tol"}
        constant_value = 1e-6;
    case {"lipschitz_constant"}
        constant_value = 1e3;
    case {"block_visiting_pattern"}
        constant_value = "sorted";
    case {"alpha_init"}
        constant_value = 1;
    case {"ds_expand_small"}
        constant_value = 1.25;
    case {"ds_shrink_small"}
        constant_value = 0.4;
    case {"ds_expand_small_noisy"}
        constant_value = 2;
    case {"ds_shrink_small_noisy"}
        constant_value = 0.5;
    case {"ds_expand_big"}
        constant_value = 1.25;
    case {"ds_shrink_big"}
        constant_value = 0.4;
    case {"ds_expand_big_noisy"}
        constant_value = 1.25;
    case {"ds_shrink_big_noisy"}
        constant_value = 0.4;
    case {"expand_small"}
        constant_value = 2;
    case {"shrink_small"}
        constant_value = 0.5;
    case {"expand_big"}
        constant_value = 1.25;
    case {"shrink_big"}
        constant_value = 0.65;
    case {"expand_big_noisy"}
        constant_value = 1.25;
    case {"shrink_big_noisy"}
        constant_value = 0.85;
    case {"is_noisy"}
        constant_value = false;
    case {"forcing_function"}
        constant_value = @(alpha) alpha^2;
    case {"reduction_factor"}
        constant_value = [0, eps, eps];
    case {"polling_inner"}
        constant_value = "opportunistic";
    case {"cycling_inner"}
        constant_value = 1;
    case {"seed"}
        constant_value = "shuffle";
    case {"output_xhist"}
        constant_value = false;
    case {"output_alpha_hist"}
        constant_value = false;
    case {"output_block_hist"}
        constant_value = false;
    case {"output_grad_hist"}
        constant_value = false;
    case {"iprint"}
        constant_value = 0;
    case {"debug_flag"}
        constant_value = false;
    case {"gradient_estimation_complete"}
        constant_value = false;
    otherwise
        error("Unknown constant name '%s'.", constant_name);
end
end