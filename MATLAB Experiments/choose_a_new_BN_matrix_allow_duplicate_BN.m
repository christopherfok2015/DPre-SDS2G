function chosen_BN_matrix_in_terms_of_positions = choose_a_new_BN_matrix_allow_duplicate_BN(residue_matrix_b_minus_currentPBN, ...
                                                                                            input_matrix_col_num, ...
                                                                                            matrix_of_nonzero_positions, ...
                                                                                            num_of_nonzero_entries_each_col)
    chosen_BN_matrix_in_terms_of_positions = zeros(input_matrix_col_num, 1);

    for col_count = 1 : input_matrix_col_num
        row_index_of_largest_element_so_far = matrix_of_nonzero_positions(1, col_count);
        largest_element_so_far = residue_matrix_b_minus_currentPBN(row_index_of_largest_element_so_far, col_count);

        for row_count = 2 : num_of_nonzero_entries_each_col(col_count)
            current_row_index = matrix_of_nonzero_positions(row_count, col_count);
            current_element = residue_matrix_b_minus_currentPBN(current_row_index, col_count);
            
            if largest_element_so_far < current_element
                largest_element_so_far = current_element;
                row_index_of_largest_element_so_far = current_row_index;
            end
        end  % end inner for
        chosen_BN_matrix_in_terms_of_positions(col_count) = row_index_of_largest_element_so_far;
    end  % end outer for
end