function the_matrix_sum = sum_up_several_BN_matrices(list_of_BN_coefficients, ...
                                                     list_of_BN_matrices_in_terms_of_pos, ...
                                                     matrix_row_num, matrix_col_num)
    
    num_of_BN_matrices_involved = length(list_of_BN_coefficients);
    the_matrix_sum = zeros(matrix_row_num, matrix_col_num);
    
    for BN_count = 1 : num_of_BN_matrices_involved
        current_coefficient = list_of_BN_coefficients(BN_count);

        for col_count = 1 : matrix_col_num
            current_row_pos = list_of_BN_matrices_in_terms_of_pos(col_count, BN_count);
            the_matrix_sum(current_row_pos, col_count) = the_matrix_sum(current_row_pos, col_count) + current_coefficient;
        end
    end

end