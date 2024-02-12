function [derivative] = calc_derivative(data,time_difference)
%CALC_DERIVATIVE Calculates derivative of the data
derivative = (data(1:end-time_difference, :) - ...
              data(1+time_difference:end, :))/time_difference;
end
