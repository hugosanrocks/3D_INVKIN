       function res = interpol(x_in,y_in,x,y,val)

       %estimate weights
       w1 = ((y(2) - y(3))*(x_in-x(3))+(x(3)-x(2))*(y_in-y(3))) / ((y(2) - y(3))*(x(1)-x(3))+(x(3)-x(2))*(y(1)-y(3)));
       w2 = ((y(3) - y(1))*(x_in-x(3))+(x(1)-x(3))*(y_in-y(3))) / ((y(2) - y(3))*(x(1)-x(3))+(x(3)-x(2))*(y(1)-y(3)));
       w3 = 1. - w1 - w2;
       res = w1*val(1) + w2*val(2) + w3*val(3);
      
       return
