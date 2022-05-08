The following two Matlab script files are examples used in the paper "Using connectome-based predictive modeling to predict individual behavior from brain connectivity" by Xilin Shen, et al.

codeshare_behavioralprediction.m implements the connectome based model to predict an external measure, e.g. behavioral scores in a leave-one-subject-out scheme. 
permutation_test_example.m estimates the distribution of the test statistics. 1000 CPM models are built based on randomly shuffled pairing between the connectome and the dependent variable. The statistical significance of the true model is calculated as the percentage of cases that exceeding the performance of the true model.





