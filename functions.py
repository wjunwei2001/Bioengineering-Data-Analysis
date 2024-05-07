import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import math
from SurvivalFunctions import compute_survival , compare_survivals
from scipy.optimize import minimize, curve_fit

def confidence_interval(arr, confidence_level):
    avg = np.mean(arr) # calculate mean
    sd = np.std(arr, ddof=1) # calculate standard deviation
    dof = len(arr) - 1
    crit_val = stats.t.ppf((1+confidence_level)/2, dof)
    print("Critical value: ", crit_val)
    lower_bound = avg - crit_val * sd / np.sqrt(len(arr))
    upper_bound = avg + crit_val * sd / np.sqrt(len(arr))
    print ("[" + str(lower_bound) + ", " + str(upper_bound) + "]")

# confidence_interval([20.1, 25.2, 22.4, 23.1, 24.8, 25.6, 22.0, 26.1], 0.95)

def t_test(arr, confidence_level, mil, two_sided=True):
    avg = np.mean(arr)
    sd = np.std(arr, ddof=1)
    dof = len(arr) - 1
    print("Mean: " + str(avg))
    print("Sample SD: " + str(sd))
    if two_sided == True:
        t_stats = (avg - mil) / (sd/np.sqrt(len(arr)))
        crit_val = stats.t.ppf((1 + confidence_level)/2, dof)
        print("T statistic: " + str(t_stats))
        print("T critical: " + str(crit_val))
        p_value = 2 * (1 - stats.t.cdf(abs(t_stats), dof))
        print("P-value: " + str(p_value))
    else:
        t_stats = (avg - mil) / (sd/np.sqrt(len(arr)))
        crit_val = stats.t.ppf(confidence_level, dof)
        print("T statistic: " + str(t_stats))
        print("T critical: " + str(crit_val))
        p_value = (1 - stats.t.cdf(abs(t_stats), dof))
        print("P-value: " + str(p_value))
    return

# t_test([20.1, 25.2, 22.4, 23.1, 24.8, 25.6, 22.0, 26.1], 0.95, 20)

'''####################### Hypo Testing on 2 means ########################'''
############# Case 1: population variance unknown, equal variance

# arr1 = []
# arr2 = []

# arr1_mean = np.mean(arr1)
# arr1_sd = np.std(arr1, ddof=1)
# arr2_mean = np.mean(arr2)
# arr2_sd = np.std(arr2, ddof=1)

# var1 = np.var(arr1, ddof=1)
# var2 = np.var(arr2, ddof=1)
# l1 = len(arr1)
# l2 = len(arr2)
# combined_var = ((l1-1)*var1 + (l2-1)*var2) / (l1+l2-2)
# # calc_comb_var(arr1,arr2)
# dof = l1+l2-2
# t_stats = ((arr1_mean - arr2_mean) - 0) / (np.sqrt(combined_var) * np.sqrt((1/l1) + (1/l2)))
# crit_val = stats.t.ppf(0.975, dof)

# p_value = 2*(1 - stats.t.cdf(abs(t_stats), dof))

# print("T statistics: ", t_stats)
# print("T Critical value: ", crit_val)
# print("P value: ", p_value)

########### Case 2: population variance unknown and unequal variance
# arr1 = [1,2,3,4]
# arr2 = [1,2,3,4]
# arr1_mean = np.mean(arr1)
# arr1_sd = np.std(arr1, ddof=1)

# arr2_mean = np.mean(arr2)
# arr2_sd = np.std(arr2, ddof=1)

# def calc_dof(arr1, arr2):
#     var1 = np.var(arr1, ddof=1)
#     var2 = np.var(arr2, ddof=1)
#     num1 = len(arr1)
#     num2 = len(arr2)
#     result = (var1/num1 + var2/num2)**2 / (((var1/num1)**2/(num1 - 1)) + ((var2/num2)**2/(num2 - 1)))
#     return math.ceil(result)

# dof = calc_dof(arr1,arr2)
# t_stats = ((arr1_mean - arr2_mean) - 0) / np.sqrt(arr1_sd**2/len(arr1) + arr2_sd**2/len(arr2))
# print("T statistics: ", t_stats)
# crit_val = stats.t.ppf(confidence_level, dof)
# p_value = 2*(1 - stats.t.cdf(abs(t_stats), dof))

# print("Critical value: ", crit_val)
# print("P value: ", p_value)

'''######################################## Power ###################################################'''
# delta = 
# sigma = 
# n1 = 
# n2 = 
# D = delta/(sigma*np.sqrt(1/n1 + 1/n2))
# t_star = 1.98 - D
# power = 1 - stats.t.cdf(t_star, n1+n1-2)
# print("T star: ", t_star)
# print("Power: ", power)


################ Finding min sample size for power of at least p ####################
# delta = 
# std = np.std(, ddof=1)
# n1 = 
# n2 = 0 
# power = 0
# desired_power = 0.5
# alpha = 0.05

# while (power < desired_power):  #iterate from 1 onwards until power >= 0.5
#     n2 = n2 + 1
#     D = delta/(std*np.sqrt(1/n1 + 1/n2))
#     t_crit = stats.t.ppf(1 - (alpha / 2), n1 + n2 - 2)
#     t_star = t_crit - D
#     power = 1 - stats.t.cdf(t_star, n1+n1-2)

# print(f"Min sample size: {n2}")


def paired_t_distr(arr1, arr2, confidence_level):
    diff_arr = arr1 - arr2
    length = len(arr1)
    avg = np.mean(diff_arr)  # calculate the mean of the difference
    sd = np.std(diff_arr, ddof=1)  # calculate the standard deviation of the difference
    t_stats = avg/(sd/np.sqrt(length))
    t_crit = stats.t.ppf((1+confidence_level)/2, length-1) 
    print("T statistics: " + str(t_stats))
    print("T critical: " + str(t_crit))
    p_value_paired = 2*(1 - stats.t.cdf(np.abs(t_stats), length-1))
    print("P value: " + str(p_value_paired))
    return


'''############## ANOVA #####################'''
# s1 = [4.4764, 4.7195, 3.6412, 4.8303, 4.5301, 4.1661, 4.1024, 4.0912]
# s2 = [4.2484, 4.4115, 4.4256, 4.8130, 5.1635, 4.6335, 4.2951, 4.7126]
# s3 = [4.2916, 3.6851, 4.6746, 4.3535, 4.2676, 4.4576, 4.3307, 3.7293]
# s4 = [4.1188, 4.1270, 3.9340, 5.0411, 3.7621, 3.9581, 3.8440, 3.6129]
# m = ??

# s1_mean = np.mean(s1)
# s1_variance = np.std(s1, ddof=1)**2

# s2_mean = np.mean(s2)
# s2_variance = np.std(s2, ddof=1)**2

# s3_mean = np.mean(s3)
# s3_variance = np.std(s3, ddof=1)**2

# s4_mean = np.mean(s4)
# s4_variance = np.std(s4, ddof=1)**2

# s_wit = (s1_variance + s2_variance + s3_variance + s4_variance) / m
# cluster = np.array([s1_mean, s2_mean, s3_mean, s4_mean])
# cluster_variance = np.std(cluster, ddof=1)**2
# s_bet = len(s1)*cluster_variance
# f_value = s_bet/s_wit

# confidence_level = 0.95
# f_crit = stats.f.ppf(confidence_level, m-1, m*(len(s1)-1))
# p_value = 1 - stats.f.cdf(f_value, m-1, m*(len(s1)-1))
# print("S bet: ", s_bet)
# print("S wit: ", s_wit)
# print("F statistics: ", f_value)
# print("F critical value: ", f_crit)
# print("P value: ", p_value)

############## Bonferroni t test #######################
# denominator = np.sqrt(s_wit/len(s1) + s_wit/len(s1))
# t_stat_1 = (s1_mean - s2_mean)/denominator
# t_stat_2 = (s1_mean - s3_mean)/denominator
# t_stat_3 = (s1_mean - s4_mean)/denominator
# t_stat_4 = (s2_mean - s3_mean)/denominator
# t_stat_5 = (s2_mean - s4_mean)/denominator
# t_stat_6 = (s3_mean - s4_mean)/denominator
# t_crit_BONF = stats.t.ppf(?, m*(n-1))
# print("T critical (Bonferroni): ", t_crit_BONF)

################# ANOVA with unequal sample size ##################
# d1 = 
# d2 = 
# d3 = 

# d1_mean = np.mean(d1)
# d1_variance = np.std(d1, ddof=1)**2

# d2_mean = np.mean(d2)
# d2_variance = np.std(d2, ddof=1)**2

# d3_mean = np.mean(d3)
# d3_variance = np.std(d3, ddof=1)**2

# N = len(d1) + len(d2) + len(d3)
# m = 

# s_wit = ((len(d1)-1)*d1_variance + (len(d2)-1)*d2_variance + (len(d3)-1)*d3_variance) / (N-m)
# s_bet = ((len(d1)*d1_mean**2 + len(d2)*d2_mean**2 + len(d3)*d3_mean**2) - ((len(d1)*d1_mean + len(d2)*d2_mean + len(d3)*d3_mean)**2 / N)) / (m-1)
# f_value = s_bet/s_wit

# confidence_level = 0.90
# f_crit = stats.f.ppf(confidence_level, m-1, N-m)
# p_value = 1 - stats.f.cdf(f_value, m-1, N-m)
# print("P value: ", p_value)

# t_crit_BONF = stats.t.ppf(1 - 0.1/3/2, N-m)
# t_stat_1 = (d1_mean - d2_mean)/np.sqrt(s_wit/len(d1) + s_wit/len(d2))
# t_stat_2 = (d1_mean - d3_mean)/np.sqrt(s_wit/len(d1) + s_wit/len(d3))
# t_stat_3 = (d2_mean - d3_mean)/np.sqrt(s_wit/len(d2) + s_wit/len(d3))


'''################ Chi Square test ######################'''
def chi_square_test_manual(observed, confidence_level):
    row_sum = [sum(i) for i in observed]
    col_sum = [sum(i) for i in zip(*observed)]
    total = sum(col_sum)

    # 1. Building the "Expected Table"
    expected_matrix = []
    for i in range(len(observed)):
        expected_matrix.append([])
        for j in range(len(observed[0])):
            expected_matrix[i].append(row_sum[i] * col_sum[j]/total)
    
    # 2. Calculate Chi square value for this problem
    chi_stats = 0
    for i in range(len(observed)):
        for j in range(len(observed[0])):
            chi_stats += (observed[i][j] - expected_matrix[i][j])**2 / expected_matrix[i][j]
    
    # 3. Calculate Chi Square critical value and p-value
    dof = (len(observed)-1) * (len(observed[1])-1)
    chi_crit = stats.chi2.ppf(confidence_level, dof)
    print("Chi Square statistics: ", chi_stats)
    print("Chi Square critical value: ", chi_crit)
    p_value = 1 - stats.chi2.cdf(chi_stats, dof)
    print("P value: ", p_value)
    
    return p_value

# observed = []
# chi_square_test_manual(observed, ?)

##### Alternative #####
# print(stats.chi2_contingency(observed))

######### Special Chi Square stat ##############
# Nominal data chapter slide 40

'''If null hypothesis is rejected, proceed to pairwise chi square test '''
def pairwise_chi_square_stats(array): 
    # Calculate Chi Square statistics for the pairwise comparison
    row_sum = [sum(i) for i in array]
    col_sum = [sum(i) for i in zip(*array)]
    total = sum(row_sum)
    result = 0
    for i in range(len(array)):
        for j in range(len(array[0])):
            result += (np.abs(array[i][j] - row_sum[i]*col_sum[j]/total)-0.5)**2 / (row_sum[i]*col_sum[j]/total)
    return result

# array = []
# # array.append(observed[0])
# # array.append(observed[1])

# confidence_level = 
# alpha = 1 - confidence_level
# dof = (len(array)-1) * (len(array[1])-1)
# chi_crit_BONF = stats.chi2.ppf(1 - alpha/???, dof)  # chisquare is ALWAYS one-sided
# print("Chi Critical value (Bonferroni): ", chi_crit_BONF)
# print(pairwise_chi_square_stats(array))

###### Alternative ####
# print(stats.chi2_contingency(array))


'''################## Survival Analysis (Tutorial 6) #######################'''
# Refer to SurvivalFunctions.py for working

# existing_results = compute_survival(unsorted_time_of_events, unsorted_type_of_events)
# new_results = compute_survival(unsorted_time_of_events, unsorted_type_of_events)

# result = compare_survivals(existing_results, new_results)
# test_stats = result[0]/np.sqrt(result[1])
# crit_value = stats.norm.ppf(0.975)
# print(f"Test statistics is {test_stats}")
# print(f"Critical value is {crit_value}")
# p_value = 2 * (1 - stats.norm.cdf(np.abs(test_stats)))

'''################## Linear Regression #####################'''
# # Find x bar and y bar
# y = 
# x = 
# mean_y = np.mean(y)
# mean_x = np.mean(x)

# # Compute slope B1 and intercept B0
# B1 = sum((y[i]-mean_y)*(x[i]-mean_x) for i in range(len(y))) / sum((i-mean_x)**2 for i in x)
# B0 = B0 = mean_y - B1*mean_x
# print("B1 (slope): ", B1)
# print("B0 (intercept): ", B0)

# #Use for checking (m=B1), (b=B0)
# m, b = np.polyfit(x, y, 1)

# sse = sum(( (B1*x[i]+B0) - y[i])**2 for i in range(len(x)))
# sst = sum((i-mean_y)**2 for i in y)
# r_squared = 1 - sse/sst

# dof = len(x) - 2
# Sy = np.sqrt(sse/dof)
# print("Residual Variance: ", Sy**2)
# SB1 = Sy*np.sqrt(1/sum((i-mean_x)**2 for i in x))
# SB0 = Sy*np.sqrt(1/len(x) + mean_x**2 / sum((i-mean_x)**2 for i in x))

# ########### Hypo testing on the slope ##########
# t_stat = B1/SB1
# t_crit = stats.t.ppf(0.975, dof)
# p_value = 2 * (1 - stats.t.cdf(abs(t_stats), dof))
# print("T statistics: ", t_stat)
# print("T critical value: ", t_crit)
# print("P value: ", p_value)
# # Data supports an increase/decrease in ___ as ___ increases. The trend was significant
# # Data fails to support any difference in ___ with ___

######## Confidence interval ###########
# t_crit = stats.t.ppf(0.975, dof)
# confidence_interval_intercept = [B0-t_crit*SB0, B0+t_crit*SB0]
# confidence_interval_slope= [B1-t_crit*SB1, B1+t_crit*SB1]
# print("Confidence interval for intercept (B0): ", confidence_interval_intercept)
# print("Confidence interval for slope (B1): ", confidence_interval_slope)

# ######## Prediction Interval ##########
# Xh = 
# Yh = B0 + B1*Xh
# Syh = Sy*np.sqrt(1/len(x) + (Xh-mean_x)**2/sum((j-mean_x)**2 for j in x))
# print(Yh - t_crit*Syh)
# print(Yh + t_crit*Syh)

########## Estimation of Ypred ###########
# Xh = 
# Yh = B0 + B1*Xh
# Syh = Sy*np.sqrt(1 + 1/len(x) + (Xh-mean_x)**2/sum((j-mean_x)**2 for j in x))
# print(Yh - t_crit*Syh)
# print(Yh + t_crit*Syh)

'''###################### General case (matrix-based) linear regression #####################'''
###### Obtaining coefficient values
# x = 
# y = 

# # Create matrix A
# order = ??
# matrix_A = []
# for i in range(len(x)):
#     lst = []
#     for j in range(order+1):
#         lst.append(np.power(x[i], j))
#     matrix_A.append(lst)
# matrix_A

######### Backup ##########
# matrix_A = np.zeros((len(x),2))
# for i in range(0,len(x)):
#     matrix_A[i,0] = 1 #coefficient of beta_0
#     matrix_A[i,1] = np.exp(x[i]) #coefficient of beta_1
############################

# K = np.dot(np.linalg.inv(np.dot(np.transpose(matrix_A), matrix_A)),np.transpose(matrix_A))
# B = np.dot(K,y) # array of B0, B1 ...
# print(K)
# print(B)

# dof = len(x) - len(B)    # number of data points - number of parameters in B
# sse = sum((y[i] - (Y_MODEL))**2 for i in range(len(y))) 
# Sy_squared = sse / (dof)
# Sy = np.sqrt(Sy_squared)

# SB0 = np.sqrt(sum(K[0][i]**2 * Sy_squared for i in range(len(K[0]))))
# SB1 = np.sqrt(sum(K[1][i]**2 * Sy_squared for i in range(len(K[1]))))
# SB2 = np.sqrt(sum(K[2][i]**2 * Sy_squared for i in range(len(K[2]))))

# ###### Hypothesis testing on regression parameters
# t_stat = B[1]/SB1
# t_crit = stats.t.ppf(0.975, dof)
# print("T statistics: ", t_stat)
# print("T critical value: ", t_crit)

# ###### Confidence intervals for regression parameters
# t_crit = stats.t.ppf(0.975, dof)
# upperSB0 = B[0] + t_crit*SB0
# lowerSB0 = B[0] - t_crit*SB0

# upperSB1 = B[1] + t_crit*SB1
# lowerSB1 = B[1] - t_crit*SB1

# upperSB2 = B[2] + t_crit*SB2
# lowerSB2 = B[2] - t_crit*SB2
# print(f"BO confidence interval: [{lowerSB0}, {upperSB0}]")
# print(f"B1 confidence interval: [{lowerSB1}, {upperSB1}]")
# print(f"B2 confidence interval: [{lowerSB2}, {upperSB2}]")

'''##################### Single variable non-linear regression ##################'''
################## Manual calculation ##################################
# x = []
# y = []
# a = 
# sse = sum((y[i] - (Y_MODEL) )**2 for i in range(len(y)))
# print(sse) 
#########################################################################

def objective_equation(ec50, concentration): # TO BE REPLACED WITH EQUATION IN THE CONTEXT
    return 100 / (1+(ec50/concentration)**2)

def calculate_sse(equation, beta, x, y):
    return sum((y[i] - equation(beta, x[i]))**2 for i in range(len(y)))

def golden_section(lower_bound, upper_bound, tol, equation, input_x, input_y):
    r = (np.sqrt(5)-1)/2
    
    # while upper_bound-lower_bound >= tol:
    #     x1 = lower_bound + (upper_bound-lower_bound)*(1-r)
    #     x2 = upper_bound - (upper_bound-lower_bound)*(1-r)
    #     sse_x1 = calculate_sse(equation, x1, input_x, input_y)
    #     sse_x2 = calculate_sse(equation, x2, input_x, input_y)
    #     if sse_x1 > sse_x2:
    #         lower_bound = x1
    #     else:
    #         upper_bound = x2
    
    max_iter = 10000000

    x1 = lower_bound + (upper_bound-lower_bound)*(1-r)
    x2 = upper_bound - (upper_bound-lower_bound)*(1-r)
    sse_x1 = calculate_sse(equation, x1, input_x, input_y)
    sse_x2 = calculate_sse(equation, x2, input_x, input_y)
    for i in range(0, max_iter):
        if sse_x1 > sse_x2:
            lower_bound = x1
            x1 = x2
            x2 = r * (upper_bound-lower_bound) + lower_bound
            sse_x1 = sse_x2
            sse_x2 = calculate_sse(equation, x2, input_x, input_y)
        else:
            upper_bound = x2
            x2 = x1
            x1 = lower_bound + (upper_bound-lower_bound)*(1-r)
            sse_x2 = sse_x1
            sse_x1 = calculate_sse(equation, x1, input_x, input_y)
        if (np.abs(upper_bound-lower_bound) < tol):
            break
            
    return (lower_bound+upper_bound)/2


####### Root finding
def obj_fun(x): # REPLACE WITH THE ACTUAL EQUATION
    return np.exp(x) - 25*x

def root_finding(lower_bound, upper_bound, tol, equation, input_x, input_y):
    r = (np.sqrt(5)-1)/2

    max_iter = 10000000

    x1 = lower_bound + (upper_bound-lower_bound)*(1-r)
    x2 = upper_bound - (upper_bound-lower_bound)*(1-r)
    f_x1 = equation(x1)
    f_x2 = equation(x2)
    for i in range(0, max_iter):
        if np.abs(f_x1) > np.abs(f_x2):
            lower_bound = x1
            x1 = x2
            x2 = r * (upper_bound-lower_bound) + lower_bound
            f_x1 = f_x2
            f_x2 = equation(x2)
        else:
            upper_bound = x2
            x2 = x1
            x1 = lower_bound + (upper_bound-lower_bound)*(1-r)
            f_x2 = f_x1
            f_x1 = equation(x1)
        if (np.abs(upper_bound-lower_bound) < tol):
            break
            
    return (lower_bound+upper_bound)/2

'''####################### Multiple variables non linear regression ##################'''
# x = ??
# y = ??
# def equation1(x, params):
#     a1 = params[0]
#     a2 = params[1]
#     mu = params[2]
#     sigma = params[3]
#     return (a1/np.sqrt(2*np.pi*sigma**2)) * np.exp(-((x-mu)**2 / sigma**2)) + a2

# def objective_func1(params):
#     fitted_y = equation1(x, params)
#     sse = np.sum((y-fitted_y)**2)
#     return sse

# initial_guess = [10, 10, 15, 1]
# result = minimize(objective_func1, initial_guess, method="Nelder-Mead", options = {'xatol':1e-8, 'disp': True})
# optimal_params = result.x
# a1_optimal, a2_optimal, mu_optimal, sigma_optimal = optimal_params

# # Print optimal parameters
# print(f'Optimal parameters:')
# print(f'a1: {a1_optimal}')
# print(f'a2: {a2_optimal}')
# print(f'mu: {mu_optimal}')
# print(f'sigma: {sigma_optimal}')