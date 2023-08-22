import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf


# Sample data
data = {
    'Category': catagory_sm,
    'Y': F_data 
}

df = pd.DataFrame(data)

# Convert the 'Category' column to dummy variables
df = pd.get_dummies(df, columns=['Category'], drop_first=True)

# Since 'Category_A' is our reference category (dropped), we include 'Category_B' and 'Category_C' in the model
model = smf.quantreg('Y ~ Category_y', df)

# Fit the model at the desired quantile (e.g., 0.25 for the lower quartile)
res = model.fit(q=0.99)

print(res.summary())


# zero filtered out

data = {
    'sm': sm_valid,
    'Y': F_w_sm_stress 
}

df = pd.DataFrame(data)

model = smf.quantreg('Y ~ sm', df)

# Fit the model at the desired quantile (e.g., 0.25 for the lower quartile)
res = model.fit(q=0.99)

print(res.summary())

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

# Assuming you've already defined your data and created the df DataFrame
# as well as fitted your quantile regression model

# Plotting the data points
plt.scatter(df['sm'], df['Y'], label='Data', alpha=0.05)

# Generating values for the predictor to visualize the model
x_vals = np.linspace(df['sm'].min(), df['sm'].max(), 100)

# Using the model to predict y values
y_vals = res.predict(pd.DataFrame({'sm': x_vals}))

# Plotting the model
plt.plot(x_vals, y_vals, color='red', label='0.95 Quantile Regression')

# Additional plot settings
plt.title('Quantile Regression of Y on sm')
plt.xlabel('sm')
plt.ylabel('Y')
plt.legend()
plt.grid(True)
plt.show()




# Sample data
data = {
    'Category': [int(i) for i in biome_list],
    'Y': F_by_biome 
}

df = pd.DataFrame(data)

# Convert the 'Category' column to dummy variables
df = pd.get_dummies(df, columns=['Category'], drop_first=True)

# Since 'Category_A' is our reference category (dropped), we include 'Category_B' and 'Category_C' in the model
model = smf.quantreg('Y ~ Category_2 + Category_3 + Category_4 + Category_5+Category_6 +Category_7 + Category_8', df)

# Fit the model at the desired quantile (e.g., 0.25 for the lower quartile)
res = model.fit(q=0.95)

print(res.summary())
