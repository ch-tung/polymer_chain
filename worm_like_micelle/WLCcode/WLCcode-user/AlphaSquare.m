function y = AlphaSquare(x)
if x > 10
y = (1 + (x/3.12).^2 + (x/8.67).^3).^(0.176/3);
else
y = (1 + (x/3.12).^2 + (x/8.67).^3).^(0.170/3);
end

% Modified by Guan-Rong Huang Oct. 16 2019