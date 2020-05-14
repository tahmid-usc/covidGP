# COVID-19 Confirmed Cumulative Case Projection
## Application of Gaussian Process Regression with Richard's Growth Curve Prior

This project aims at building a forecast model for cumulative COVID-19 cases in USA.
To view projection on cumulative cases in US, go to: https://tahmid-usc.github.io/covidGP/
For projection in CSV format, go to https://github.com/tahmid-usc/covidGP/tree/master/docs/prediction

### Model

Parametric nonlinear growth curve models have been used to model COVID-19 cumulative confirmed cases and deaths with data obtained from all over the world. Logistic curve (Verhulst, 1838), generalized logistic curve, Richard's curve(Richards (1959)), Gompertz curve (Gompertz, 1825; Winsor, 1932), cumulative distribution function method etc. were applied by various researchers to fit the COVID-19 data and forecast the trend in the future. The forecast is tremendously important to assess and configure emergency life saving resources like PPE, ventilator, medical supply and medicine 
