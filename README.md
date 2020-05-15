COVID-19 Case and Death in USA Projection
================
Tahmidul Islam

## Application of Gaussian Process Regression with Richard’s Growth Curve Prior

This project aims at building a forecast model for cumulative COVID-19
cases in USA. To view projection on cumulative cases in US, go to:
<https://tahmid-usc.github.io/covidGP/> For projection in CSV format, go
to <https://github.com/tahmid-usc/covidGP/tree/master/docs/prediction>

### Model

Parametric nonlinear growth curve models have been used to model
COVID-19 cumulative confirmed cases and deaths with data obtained from
all over the world. Logistic curve (Verhulst, 1838), generalized
logistic curve, Richard’s curve(Richards (1959)), Gompertz curve
(Gompertz, 1825; Winsor, 1932), cumulative distribution function method
etc. were applied by various researchers to fit the COVID-19 data and
forecast the trend in the future. The forecast is tremendously important
to assess and configure emergency life saving resources like PPE,
ventilator, medical supply and medicine.

The parametric models are somewhat rigid in terms of deviation of
observed data from hypothetical model. If parametric assumption is
violated then the usefulness of the model declines specially for
forecasting purpose. But completely nonparametric model also suffers
from the inability to make valid prediction outside the support of the
data.

Here we propose a Bayesian nonlinear regression model using Gaussian
process regression where we have incorporated a modified Richard’s curve
as the prior mean function propagating our belief about the shape of the
growth curve. But this model allows deviation from the parametric model
by updating the prior with observed data. The used modified Richad’s
curve is:

\[ \mu(t) = \frac{A}{(1 + B_0 e^{(- B_1 t)})^v}.\]

Here, \(\mu(t)\) is the cumulative number of cases/deaths at time \(t\).
\(A, B_0, B_1,\) and \(v\) are the parameters of this nonlinear growth
function. We use this as the prior mean function in GP regression
framework.

Here \(K\) is the kernel function, usually chosen to be Gaussian and
used in this work. We have used numerical optimization with multiple
starts to estimate the marginal likelihood function to estimate all the
parameters in the model. The posterior predictive distribution was
derived and the mean of this distribution is used as the point
estimate/forecast of the response. The diagonal of the posterior
covariance was used to construct the credible error band.

A sample projection of cumulative number of confirmed cases for South
Carolina is shown.

![](docs/plot/South%20Carolina.png)
