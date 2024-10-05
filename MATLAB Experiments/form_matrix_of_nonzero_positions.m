function matrix_of_nonzero_positions = form_matrix_of_nonzero_positions(input_PBN_matrix, ...
                                                                        input_matrix_row_num, ...
                                                                        input_matrix_col_num)
    
    matrix_of_nonzero_positions = zeros(input_matrix_row_num, input_matrix_col_num);
    
    for col_count = 1 : input_matrix_col_num
        pointer = 1;
        for row_count = 1 : input_matrix_row_num
            if input_PBN_matrix(row_count, col_count) > 0
                matrix_of_nonzero_positions(pointer, col_count) = row_count;
		        pointer = pointer + 1;
            end
        end
    end
end