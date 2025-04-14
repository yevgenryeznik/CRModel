# CRModel

## `Julia` routine to calculate Locally Optimal Restricted Designs for Continuation-Ratio Model. 

Here, we have a collection of `Julia` scripts to calculate locally optmal designs for the _**continuation-ratio model**_.
The model is a four-parameter model that might be used to model _**multinomial probabilites**_ in the _**dose-efficacy-toxiciity**_ studies.

$$\boldsymbol{\eta}(x; \boldsymbol{\theta}) = \left\\{
\begin{aligned}
\pi_0(x; \boldsymbol{\theta}) &= \frac{1}{1+\exp{(\theta_1+ + \theta_2 x)}}\frac{1}{1+\exp{(\theta_3+ + \theta_3 x)}} \\
\pi_1(x; \boldsymbol{\theta}) &= \frac{\exp{(\theta_1+ + \theta_2 x)}}{1+\exp{(\theta_1+ + \theta_2 x)}}\frac{1}{1+\exp{(\theta_3+ + \theta_3 x)}} \\
\pi_2(x; \boldsymbol{\theta}) &= \frac{\exp{(\theta_3+ + \theta_4 x)}}{1+\exp{(\theta_3+ + \theta_3 x)}}.
\end{aligned}\right.$$
