% specify the input PBN matrix and matrix size (column number or row number).
input_PBN_matrix = 0.0001 * [1599, 2501, 1271, 1964, 2501,    0,    0,  242;
                              184, 1618, 1906,    0,    0, 1906, 1308, 2501;
                             1308,    0,  224,  198,  173, 2494, 2469, 1895;
                                0, 1308,    0, 1599, 1906, 1271, 2494, 2520;
                                0,  198, 2494, 2480, 1290, 1567,    0, 1578;
                             2494, 2480,    0,    0, 1592, 2538, 1927, 1264;
                             1895,    0, 2538, 2469, 2538,  224, 1578,    0;
                             2520, 1895, 1567, 1290,    0,    0,  224,    0];
input_matrix_size = size(input_PBN_matrix);
input_matrix_row_num = input_matrix_size(1);
input_matrix_col_num = input_matrix_size(2);
algo_choice = 'active-set'; % either interior-point-convex or active-set.
matrix_of_nonzero_positions = form_matrix_of_nonzero_positions(input_PBN_matrix, ...
                                                               input_matrix_row_num, ...
                                                               input_matrix_col_num);

% Set the initial guess x^0 for the MOMP algorithm. 
% Note that x^0 here cannot involve too many nonzero positive entries; 
% otherwise, memory error will be resulted. If needed at a later time, 
% I will modify the MOMP function so that we can set an initial guess x^0 
% involving many nonzero positive entries (for example, setting x^0 to be 
% the uniform distribution).
initial_guess_is_uniform = true;
initial_guess_coefficients = 1;
initial_guess_BN_matrices_in_terms_of_nonzero_pos = zeros(input_matrix_col_num, 1);

for col_count = 1 : input_matrix_col_num
    initial_guess_BN_matrices_in_terms_of_nonzero_pos(col_count) = matrix_of_nonzero_positions(1, col_count);
end

% set the initial point argument of the quadprog function.
% quadprog_initial_point_argument can be set to 2 values: either "uniform"
% or "concentrated at first position".
quadprog_initial_point_argument = "uniform";

% set the stopping_threshold argument of the momp function.
% Please see page 9 of the MOMP paper.
% stopping_criteria_type can be either "the_obvious_difference_Ax_minus_b"
% or "according_to_page_9_of_the_paper".
stopping_threshold = 10^(-7);
stopping_criteria_type = "according_to_page_9_of_the_paper";

% Perform MOMP PBN construction.
%
% If step 1 of MOMP produces a BN vector that is already present in S^k 
% (list_of_BN_matrices_in_terms_of_pos), the value of presence_of_duplicate_BN
% will be set to 1. Otherwise, the value of presence_of_duplicate_BN will
% be set to 0.
%
% If no duplicate BN vector is generated throughout the execution of 
% momp_allow_duplicate_BN, duplicate_BN_matrix_itf_position will be a zero
% column vector of length input_matrix_col_num.
tic
[list_of_BN_coefficients, list_of_BN_matrices_in_terms_of_pos, list_of_exit_flags, ...
 duplicate_BN_matrix_itf_position, presence_of_duplicate_BN] = momp_allow_duplicate_BN(input_PBN_matrix, ...
                                   input_matrix_row_num, input_matrix_col_num, algo_choice, ...
                                   initial_guess_is_uniform, initial_guess_coefficients, ...
                                   initial_guess_BN_matrices_in_terms_of_nonzero_pos, quadprog_initial_point_argument, ...
                                   stopping_threshold, stopping_criteria_type);
toc

% Display presence_of_duplicate_BN.
presence_of_duplicate_BN
% Display duplicate_BN_matrix_itf_position.
duplicate_BN_matrix_itf_position

% Find out the number of positive entries, the number of zero entries,
% and the number of negative entries in list_of_BN_coefficients.
num_of_positive_entries_in_decomposition_found = sum(list_of_BN_coefficients > 0);
num_of_zero_entries_in_decomposition_found = sum(list_of_BN_coefficients == 0);
num_of_negative_entries_in_decomposition_found = sum(list_of_BN_coefficients < 0);

% display the value of these variables.
num_of_positive_entries_in_decomposition_found
num_of_zero_entries_in_decomposition_found
num_of_negative_entries_in_decomposition_found

% check that the sum of entries of list_of_BN_coefficients equals 1.
sum_of_coefficients_in_decomposition_found = sum(list_of_BN_coefficients);
% Display the value of the variable.
sum_of_coefficients_in_decomposition_found

% check that Ax (in unflattened form) looks very similar to (or even the same
% as) the original PBN matrix b (in unflattened form).
sum_of_decomposition_found = sum_up_several_BN_matrices(list_of_BN_coefficients, ...
                             list_of_BN_matrices_in_terms_of_pos, input_matrix_row_num, ...
                             input_matrix_col_num);
% display Ax (in unflattened form).
sum_of_decomposition_found

% Check that ||Ax - b||_{2} is small.
Ax_minus_b_vector_form = flatten_matrix_to_col_vec(sum_of_decomposition_found - input_PBN_matrix, ...
                                                   input_matrix_row_num, ...
                                                   input_matrix_col_num);
norm_of_the_vec_Ax_minus_b = norm(Ax_minus_b_vector_form);
% display the value of ||Ax - b||_{2}.
norm_of_the_vec_Ax_minus_b

% Display the exit flags' values for each quadratic optimization performed 
% (i.e., each iteration of step 2 of the MOMP algorithm).
list_of_exit_flags

% Display the BN matrices found and their corresponding coefficients, in
% the order of the iterations of steps 1 to 3 of MOMP.
% Note that list_of_BN_coefficients may contain zero and/or negative
% entries.
list_of_BN_coefficients
list_of_BN_matrices_in_terms_of_pos

% Sort the positive coefficients in x in descending order. 
% Sort the list of BN matrices accordingly.
[list_of_BN_coefficients_sorted, indices_of_sorted_coefficients] = sort(list_of_BN_coefficients, ...
                                                                        'descend');
list_of_BN_matrices_itf_pos_sorted = list_of_BN_matrices_in_terms_of_pos(:, indices_of_sorted_coefficients);

% Display the sorted list of coefficients and the sorted list of BN
% matrices.
list_of_BN_coefficients_sorted
list_of_BN_matrices_itf_pos_sorted