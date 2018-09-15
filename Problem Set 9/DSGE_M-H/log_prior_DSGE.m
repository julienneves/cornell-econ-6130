function LP=log_prior_DSGE(theta)


LP=0;

% %(1) tau - consumption curvature
pmean=4; pstdd=1;
LP = LP + lpdfNormal(theta(1),4,1); 

%(2) beta - discount factor
pmean=2; pstdd=0.05;
LP = LP + lpdfNormal(theta(2),pmean,pstdd) ;

%(3) d - Calvo price setting parameter
pmean=0.75; pstdd=0.05;
a = (1-pmean)*pmean^2/pstdd^2 - pmean;
b = a*(1/pmean - 1);
LP = LP + lpdfBeta(theta(3),a,b); 

%(4) f - coefficient on inflation in Taylor rule
pmean=1.5; pstdd=0.1;
LP = LP + lpdfNormal(theta(4),pmean,pstdd) ;

%(5) f - coefficient on output in Taylor rule
pmean=.125; pstdd=0.1;
LP = LP + lpdfNormal(theta(5),pmean,pstdd) ;

% %(6) varphi - disutility of labor curvature
pmean=4; pstdd=1;
LP = LP + lpdfNormal(theta(6),4,1); 

%(7) d - Calvo price setting parameter
pmean=1/3; pstdd=0.1;
a = (1-pmean)*pmean^2/pstdd^2 - pmean;
b = a*(1/pmean - 1);
LP = LP + lpdfBeta(theta(7),a,b); 
 
