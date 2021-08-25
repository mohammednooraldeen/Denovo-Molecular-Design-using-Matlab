function [c, ceq] = simple_constraint(x)
%SIMPLE_CONSTRAINT Summary of this function goes here
%   Detailed explanation goes here
c = [1.5 + x(1)*x(2) + x(1) - x(2);   -x(1)*x(2) + 10];
ceq = [];

end

