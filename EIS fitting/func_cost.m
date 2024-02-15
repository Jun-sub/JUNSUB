function [cost] = func_cost(x_data,y_data,para)

y_model = func_model(x_data,para);

residual = y_model - y_data;

cost = sqrt(sum(residual.^2));



end

