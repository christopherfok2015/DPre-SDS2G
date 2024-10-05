% num_of_BN_matrices is the number of columns in
% list_of_BN_matrices_in_vec_form.
function [list_of_BN_coefficients, objective_func_val, exit_flag] = momp_quadprog(list_of_BN_matrices_in_vec_form, ...
                                                                                  input_PBN_matrix_flattened, ...
                                                                                  algo_choice, ...
                                                                                  quadprog_initial_point_argument, ...
                                                                                  num_of_BN_matrices)
    H = list_of_BN_matrices_in_vec_form' * list_of_BN_matrices_in_vec_form;
    f = - list_of_BN_matrices_in_vec_form' * input_PBN_matrix_flattened;
    Aeq = ones(1, num_of_BN_matrices);
    beq = 1;
    lb = zeros(num_of_BN_matrices, 1);
    options = optimoptions('quadprog', 'Algorithm', algo_choice);

    if quadprog_initial_point_argument == "uniform"
        x0 = ones(num_of_BN_matrices, 1) / num_of_BN_matrices;
    elseif quadprog_initial_point_argument == "concentrated at first position"
        x0 = zeros(num_of_BN_matrices, 1);
        x0(1) = 1;
    end

    [list_of_BN_coefficients, objective_func_val, exit_flag, ~] = quadprog(H, f, [], [], Aeq, beq, lb, [], x0, options);
end