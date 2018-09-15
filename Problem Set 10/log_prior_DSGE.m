function LP=log_prior_DSGE(theta)

LP=0;

% %(1) sigma  - risk aversion
pmean=2; pstdd=.1;
LP = LP + lpdfNormal(theta(1),pmean,pstdd); 

%(2) beta - discount factor
pmean=0.99; pstdd=0.05;
a = (1-pmean)*pmean^2/pstdd^2 - pmean;
b = a*(1/pmean - 1);
LP = LP + lpdfBeta(theta(2),a,b); 

% %(3) phi - inverse of elasticity of labor supply
pmean=4; pstdd=.1;
LP = LP + lpdfNormal(theta(3),pmean,pstdd); 

%(4) phi_pi - coefficient on inflation in Taylor rule
pmean=1.5; pstdd=0.1;
LP = LP + lpdfNormal(theta(5),pmean,pstdd) ;

%(5) phi_y - coefficient on output in Taylor rule
pmean=.125; pstdd=0.1;
LP = LP + lpdfNormal(theta(6),pmean,pstdd) ;

%(6) theta - price stickiness
pmean=0.75; pstdd=0.1;
a = (1-pmean)*pmean^2/pstdd^2 - pmean;
b = a*(1/pmean - 1);
LP = LP + lpdfBeta(theta(7),a,b); 

%(7) a - labor share in GDP
pmean=1/3; pstdd=0.1;
a = (1-pmean)*pmean^2/pstdd^2 - pmean;
b = a*(1/pmean - 1);
LP = LP + lpdfBeta(theta(8),a,b); 
 
