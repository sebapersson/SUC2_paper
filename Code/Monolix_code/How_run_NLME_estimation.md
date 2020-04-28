This file contains a description on how to run the NLME parameter estimation the same way as in the paper for the *Simple_feedback* and *Snf1_feedback* models. 

## Simple feedback model (feedback cascade in the paper) 

1. Use the **Data_monolix_SUC2.csv** as data-file. 
2. Load the Simple_feedback_model.txt files from the **Monolix_code** directory. 
3. Use the initial values in the table below. 
4. Choose a constant error model. 
5. Use a full correlation model, except for k5 and k8 that should be fixed. 
6. In population parameter settings: 
    * Set max number of iterations to 1500 
    * For *Methodology for parameters without variability* choose *Variability at the first stage*. 
7. Besides point 4, 5 and 6 use default settings and run the full analysis. 
8. Save the mlxtran file as **Simple_feedback.mlxtran** in the **Code/Monolix_code/Simple_feedback** directory. 
9. Export charts data. 

## Snf1 feedback model (feedback mediated in the paper)
 
1. Use the **Data_monolix_SUC2.csv** as data-file.  
2. Load the Snf1_feedback_model.txt files from the **Monolix_code** directory. 
3. Use the initial values in the table below. 
4. Choose a constant error model. 
5. Use a full correlation model. 
6. In population parameter settings: 
    * Set max number of iterations to 1500 
    * Disable simulated annealing. Recommended by Monolix when the start guess should be close to an optimum. 
7. Besides point 4, 5 and 6 use default settings and run the full analysis. 
8. Save the mlxtran file as **Snf1_feedback.mlxtran** in the **Code/Monolix_code/Snf1_feedback** directory. 
9. Export charts data. 


| Parameter | Simple feedbac | Snf1 feedback |
| -------- | -------- | -------- |
| k1     | 0.027     | 0.47     |
| k2     | -     | -     |
| k3     | 0.22     | 0.048     |
| k4     | 10.66     | 1.19     |
| k5     | 1.35    | 0.17     |
| k6     | -     | 0.071     |
| k7     | 6.29     | -     |
| k8     | 5.35     | 1.16     |
| k9     | 10.93     | 0.0074     |
| tau_X     | 183.8     | -     |
| SUC20     | 4.05     | 4.1     |
