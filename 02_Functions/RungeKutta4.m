function Y_next = RungeKutta4(model, t, Y, dt)
    k1 = model(t, Y);
    k2 = model(t + 0.5*dt, Y + 0.5*dt*k1);
    k3 = model(t + 0.5*dt, Y + 0.5*dt*k2);
    k4 = model(t + dt, Y + dt*k3);
    Y_next = Y + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4);
end