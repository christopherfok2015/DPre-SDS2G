function number_in_special_num_system_form = convert_from_ordinary_to_special(number_in_ordinary_form, ...
                                                                              input_matrix_col_num, ...
                                                                              num_of_nonzero_entries_each_col)
    number_in_special_num_system_form = zeros(input_matrix_col_num, 1);
    current_remaining_value = number_in_ordinary_form;

    for col_count = input_matrix_col_num : -1 : 2
        current_divisor = num_of_nonzero_entries_each_col(col_count);
        digit = rem(current_remaining_value, current_divisor);
        number_in_special_num_system_form(col_count) = digit;
        current_remaining_value = (current_remaining_value - digit) / current_divisor;
    end

    number_in_special_num_system_form(1) = current_remaining_value;
end