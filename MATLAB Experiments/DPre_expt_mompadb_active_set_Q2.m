% specify the input PBN matrix and matrix size (column number or row number).
input_PBN_matrix = (1 / 86605) * [    0,  5750, 13587,  3407, 13616, 18379,     0, 12038, 13574, 13566,  3217,  6668,     0, 18387, 13577,  4819;
                                   4813,  9470,  2153,  3260,  1660, 13992,  3352, 18385,  2137,   829,  3751,  1617,  7117,   806,  3325,   143;
                                    143,  4590,  3334,  3769,  9511, 15837, 18533,  2269,  2292, 11971, 13597,  2003,     0,  3733,  6668,     0;
                                  15173,  4851,   825,  9453,  4833,  2311,  4838,  6683,  2279,  3734,  3343,  2271,  9601,  4789,  2262,  6658;
                                    844,     0,  2279, 14398,  2272,  4797,  4555,  2271, 15305,  2059,  4788,  2128, 13593,  2272,  2004,  2178;
                                   3792, 11989, 11978,   154,  2270,  9461,  8877,   158,  4799,  6658,  6695,  2325,  2055,  1646, 12038,  3788;
                                   2282, 13567,  6695,     0,   816,     0,  3193,  1607, 18422,  4856, 18362,  3755,  4852,  2149,  3755,  2306;
                                   3334,   197,  1996,  2012, 12003,   134,     0,  3213,  3213,  2288,  2045, 18387,  3263, 11995,  9453,  2028;
                                  13574,  3202,  1664, 14256,  1995,  3370,     0,   843,  1673,   158,  1605,  9468, 13576,  2288,  4821, 11995;
                                   1627, 20507, 18404, 18381,  3326,  1664, 17301,  3725,  9476,  1646,  9453,  3201,  2139,  3237,  3237,  3352;
                                   6683,     0,  3256,  3907,  3213,  3732,  9488, 13630,   178,  3203, 11981, 11971,     0,   201,  2282,     0;
                                  11731,  3396,  2294,     0, 18381,   862,   816,  2045,  7512,  2279,  2294,  3336,  6651, 13630,   862,  4017;
                                   2146,  1629,  9474,     0,  6681,  3213, 12043,  4815,  2020, 18399,   799,   160,  3067,  2005,  1646, 11040;
                                   2045,  7457,  4799,  6661,   143,     0,  2005,  9461,     0,  9476,  2325,  4856,     0,  3356, 18381,  2272;
                                  18418,     0,   143,  6947,  2153,  2192,     0,  2128,     0,  2149,  2196,   848,  2306,  9443,   141, 13587;
                                      0,     0,  3724,     0,  3732,  6661,  1604,  3334,  3725,  3334,   154, 13611, 18385,  6668,  2153, 18422];

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
[list_of_BN_coefficients, list_of_BN_matrices_in_terms_of_pos, list_of_exit_flags, ...
 duplicate_BN_matrix_itf_position, presence_of_duplicate_BN] = momp_allow_duplicate_BN(input_PBN_matrix, ...
                                   input_matrix_row_num, input_matrix_col_num, algo_choice, ...
                                   initial_guess_is_uniform, initial_guess_coefficients, ...
                                   initial_guess_BN_matrices_in_terms_of_nonzero_pos, quadprog_initial_point_argument, ...
                                   stopping_threshold, stopping_criteria_type);

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