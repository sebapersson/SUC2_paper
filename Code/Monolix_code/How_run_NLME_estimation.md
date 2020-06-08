This file contains a description on how to run the NLME parameter estimation the same way as in the paper for the *Simple_feedback* and *Snf1_feedback* models. 

## Simple feedback model (feedback cascade in the paper) 

1. Use the **Data_monolix_SUC2.csv** as data-file. 
2. Load the Simple_feedback_model.txt files from the **Monolix_code** directory. 
3. Use the initial values in the table below. 
4. Choose a constant error model. 
5. Use a full correlation model, except for k5 and k8 that should be fixed. 
6. In population parameter settings: 
    * Set max number of iterations to 1000 
    * For *Methodology for parameters without variability* choose *No variability*
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
    * Set max number of iterations to 1000 
    * Disable simulated annealing. Recommended by Monolix when the start guess should be close to an optimum. 
7. Besides point 4, 5 and 6 use default settings and run the full analysis. 
8. Save the mlxtran file as **Snf1_feedback.mlxtran** in the **Code/Monolix_code/Snf1_feedback** directory. 
9. Export charts data. 


| Parameter | Simple feedback | Snf1 feedback |
| -------- | -------- | -------- |
| k1     | 0.087     | 0.25     |
| k2     | 0.055     | 1.33     |
| k3     | 4.83     | 0.09     |
| k4     | 12.73     | 1.69     |
| k5     | 2.09    | 0.15     |
| k6     | 0.87     | 0.81     |
| k7     | 0.62     | 0.02     |
| k8     | 9.07     | 0.17     |
| k9     | 5.17     | 0.00064     |
| k10     | -     | 0.017     |
| tau_X     | 218.8     | -     |

