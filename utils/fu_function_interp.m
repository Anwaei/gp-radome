function out1 = fu_function_interp(out, delta_t_disc)
E = out.yout.numElements;
t = out.tout;
tt = 0:delta_t_disc:t(end);
out1.tout = tt;
for i = 1:E
    y = out.yout{i}.Values.Data;
    y1 = interp1(t,y,tt,'spline');
    out1.yout{i}.Values.Data = y1;
end
end