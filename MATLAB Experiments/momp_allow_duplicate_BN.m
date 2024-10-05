function [list_of_BN_coefficients, list_of_BN_matrices_in_terms_of_pos, list_of_exit_flags, ...
          duplicate_BN_matrix_itf_position, presence_of_duplicate_BN] = momp_allow_duplicate_BN(input_PBN_matrix, ...
                                            input_matrix_row_num, input_matrix_col_num, algo_choice, ...
                                            initial_guess_is_uniform, initial_guess_coefficients, ...
                                            initial_guess_BN_matrices_in_terms_of_nonzero_pos, ...
                                            quadprog_initial_point_argument, ...
                                            stopping_threshold, stopping_criteria_type)
    % S^k will be implemented as the variable list_of_BN_matrices_in_vec_form.
    % There will be two variables representing the list of BN matrices: 
    % (1) list_of_BN_matrices_in_vec_form; (2) list_of_BN_matrices_in_terms_of_pos

    % Step 0 (modified): create and initialize some variables, and compute b – Ax0 (in unflattened form)
    list_of_BN_matrices_in_vec_form = [];
    list_of_BN_matrices_in_terms_of_pos = [];
    list_of_BN_coefficients = [];  % This output will be a column vector.

    list_of_exit_flags = [];

    duplicate_BN_matrix_itf_position = zeros(input_matrix_col_num, 1);
    presence_of_duplicate_BN = 0;

    input_PBN_matrix_flattened = flatten_matrix_to_col_vec(input_PBN_matrix, ...
                                                           input_matrix_row_num, input_matrix_col_num);

    k = 0; % count the number of MOMP main loops executed.
    matrix_of_nonzero_positions = form_matrix_of_nonzero_positions(input_PBN_matrix, ...
                                                                   input_matrix_row_num, input_matrix_col_num);
    num_of_nonzero_entries_each_col = count_num_of_nonzero_entries_each_col(input_PBN_matrix, ...
                                                                            input_matrix_row_num, input_matrix_col_num);

    current_PBN_matrix_found = zeros(input_matrix_row_num, input_matrix_col_num);

    if initial_guess_is_uniform == true
        for col_count = 1 : input_matrix_col_num
            for row_count = 1 : num_of_nonzero_entries_each_col(col_count)
                current_row_index = matrix_of_nonzero_positions(row_count, col_count);
                current_PBN_matrix_found(current_row_index, col_count) = 1 / num_of_nonzero_entries_each_col(col_count);
            end
        end
    else  % when initial_guess_is_uniform == false
        current_PBN_matrix_found = sum_up_several_BN_matrices(initial_guess_coefficients, ...
            initial_guess_BN_matrices_in_terms_of_nonzero_pos, input_matrix_row_num, input_matrix_col_num);
    end
    residue_matrix_b_minus_currentPBN = input_PBN_matrix - current_PBN_matrix_found;
    current_error = 1 + stopping_threshold;  % any value greater than stopping_theshold will do – I need to ensure that the main loop of MOMP is entered.

    % Step 1-3 of MOMP
    while (current_error > stopping_threshold) && (presence_of_duplicate_BN == 0)
        % Step 1
        chosen_BN_matrix_in_terms_of_positions = choose_a_new_BN_matrix_allow_duplicate_BN(residue_matrix_b_minus_currentPBN, ...
                                                                                           input_matrix_col_num, ...
                                                                                           matrix_of_nonzero_positions, ...
                                                                                           num_of_nonzero_entries_each_col);
        
        if ~isempty(list_of_BN_matrices_in_terms_of_pos) && ismember(chosen_BN_matrix_in_terms_of_positions', ...
                                                                     list_of_BN_matrices_in_terms_of_pos', 'rows')
            duplicate_BN_matrix_itf_position = chosen_BN_matrix_in_terms_of_positions;
            presence_of_duplicate_BN = 1;

        else % if the chosen BN matrix is not a previously chosen BN matrix
            list_of_BN_matrices_in_terms_of_pos(:, k + 1) = chosen_BN_matrix_in_terms_of_positions;
            chosen_BN_matrix_in_matrix_form = sum_up_several_BN_matrices(1, chosen_BN_matrix_in_terms_of_positions, ...
                                                                         input_matrix_row_num, input_matrix_col_num);
            chosen_BN_matrix_in_vec_form = flatten_matrix_to_col_vec(chosen_BN_matrix_in_matrix_form, ...
                                                                     input_matrix_row_num, input_matrix_col_num);
            list_of_BN_matrices_in_vec_form(:, k + 1) = chosen_BN_matrix_in_vec_form;
    
            % Step 2
            [list_of_BN_coefficients, objective_func_val, current_exit_flag] = momp_quadprog(list_of_BN_matrices_in_vec_form, ...
                input_PBN_matrix_flattened, algo_choice, quadprog_initial_point_argument, k + 1);
    
            % Step 3
            k = k + 1;
            current_PBN_matrix_found = sum_up_several_BN_matrices(list_of_BN_coefficients, ...
                                                                  list_of_BN_matrices_in_terms_of_pos, ...
                                                                  input_matrix_row_num, input_matrix_col_num);
            residue_matrix_b_minus_currentPBN = input_PBN_matrix - current_PBN_matrix_found;
            list_of_exit_flags(k) = current_exit_flag;
            current_error = compute_current_error(residue_matrix_b_minus_currentPBN, ...
                                                  list_of_BN_matrices_in_terms_of_pos, ...
                                                  current_PBN_matrix_found, ...
                                                  matrix_of_nonzero_positions, ...
                                                  num_of_nonzero_entries_each_col, ...
                                                  stopping_criteria_type, ...
                                                  input_matrix_row_num, input_matrix_col_num);
        end  % end if
    end  % end while

end % end function