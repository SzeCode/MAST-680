%RK4 method
function [t,w] = RK4_Sys_Inputs(f, to, yo, n, tf, inputs)
  
  t = [to];
  w = [yo];
  tol = 1.0e-10;
  
  i = 1;
 
  h = (tf - to)/n;
  Time_Increment = h
  
  tic;
  while(tf - t(i) > tol)
  
    k1 = f(t(i),w(i,:), inputs);
    k2 = f(t(i) + h/2.0, w(i,:) + h/2.0*k1, inputs);
    k3 = f(t(i) + h/2.0, w(i,:) + h/2.0*k2, inputs);
    k4 = f(t(i) + h, w(i,:) + h*k3, inputs);
  
    t = [t ; t(i) + h];
    w = [w ; w(i,:) + h/6.0 * (k1 + 2*k2 + 2*k3 + k4)];
 
    i++;
  
  endwhile
  
  disp("");
  disp("Computation Time:");
  toc;

endfunction
