function output_col_vec = flatten_matrix_to_col_vec(input_matrix, matrix_row_num, matrix_col_num)
    output_col_vec = zeros(matrix_row_num * matrix_col_num, 1);
    
    for count = 1 : matrix_col_num
        output_col_vec(matrix_row_num * (count - 1) + 1 : matrix_row_num * count) = input_matrix(:, count);
    end
end