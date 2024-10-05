function number_in_ordinary_form = convert_from_special_to_ordinary(number_in_special_num_system_form, ...
                                                                    input_matrix_col_num, ...
                                                                    num_of_nonzero_entries_each_col)
    number_in_ordinary_form = number_in_special_num_system_form(1);

    for count = 2 : input_matrix_col_num
        number_in_ordinary_form = num_of_nonzero_entries_each_col(count) * number_in_ordinary_form + ...
                                  number_in_special_num_system_form(count);
    end
end