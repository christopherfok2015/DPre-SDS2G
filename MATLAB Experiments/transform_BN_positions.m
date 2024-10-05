function BN_itf_positions_transformed = transform_BN_positions(BN_itf_positions_untransformed, input_matrix_col_num, ...
                                                               matrix_of_nonzero_positions)
    
    BN_itf_positions_transformed = zeros(input_matrix_col_num, 1);
    for col_count = 1 : input_matrix_col_num
        BN_itf_positions_transformed(col_count) = ...
            matrix_of_nonzero_positions(BN_itf_positions_untransformed(col_count), col_count);
    end
end