""" least squares model """
from sklearn import linear_model

reg = linear_model.LinearRegression()
x = [[0,0],[1,1],[2,2]]
y = [0,1,2]
reg.fit(x,y)
print(reg.coef_)

# If X has a linear denpendence ->
# least square estimate is highly sensitive to random errors in the observed responce.


""" ridge regression = least squares + regularization """
reg = linear_model.Ridge(alpha = 0.5)
x = [[0,0],[0,0],[1,1]]
y = [0,0.1,1]
reg.fit(x,y)
print('w:', reg.coef_)
print('w0:', reg.intercept_)

# set regularization strength with Cross-Validation: alpha
reg = linear_model.RidgeCV(alphas = [0.1, 1.0, 10.0])
reg.fit(x,y)
print('best alpha', reg.alpha_)



