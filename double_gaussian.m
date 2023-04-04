function z = gaussian(x,y,sigma)

     z=exp(-.5 .* (x.^2+y.^2)/sigma^2)/(sigma*(sqrt(2*pi)));
