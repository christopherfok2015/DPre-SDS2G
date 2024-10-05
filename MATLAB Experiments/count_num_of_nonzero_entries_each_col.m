function num_of_nonzero_entries_each_col = count_num_of_nonzero_entries_each_col(input_PBN_matrix, ...
                                                                                 input_matrix_row_num, ...
                                                                                 input_matrix_col_num)
    
    % count the number of nonzero entries in each column of input PBN matrix.
    num_of_nonzero_entries_each_col = zeros(1, input_matrix_col_num);
    for col_count = 1 : input_matrix_col_num
        for row_count = 1 : input_matrix_row_num
            if input_PBN_matrix(row_count, col_count) > 0
                num_of_nonzero_entries_each_col(col_count) = num_of_nonzero_entries_each_col(col_count) + 1;
            end
        end
    end

end