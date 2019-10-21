function [newParam]=chooseStep(param,min,max)
%A function that chooses a parameter randomly shifted by a step within the
%the parameter range.

%define step range
n_min=min-param;
n_max=max-param;



%pick step size
step=(n_max-n_min)*rand()+n_min;

newParam=param+step;

end