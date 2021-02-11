function TimeConstant_us = FindTimeConstant_us(ratio,ExpTime_us)
t_idf= 0.18; % interframe deadtime in microseconds (for pco.1600)
TimeConstant_us= fzero(@(tau)lifetime(tau,ratio,ExpTime_us,ExpTime_us+t_idf),5);
end
function err=lifetime(tau,ratio,delta,Delta)
err= (1-exp(-delta/tau))*exp(Delta/tau) - ratio;
end

