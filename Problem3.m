u10 = 1;
u20 = 2;
u0 = [u10; u20];
f = @(t, u) [-1000 1; 0 (-1/10)] * u;

T = 1;
N = 400;

h = T / N;

[t, u] = fourthOrderRK(f, u0, h, N);

[T, U] = ode23s(f, [0, 1], u0);
size(T)

u1true = @(t) (9979/9999)*exp(-1000*t) + (20/9999)*exp(-0.1*t);
u2true = @(t) 2*exp(-0.1*t);

%plot(t, u(2,:), '-r', "LineWidth", 2)
hold on;
plot(t, u1true(t), '-b', "LineWidth", 2)
plot(T, U(:, 1), '-r', "LineWidth", 2)
title("ode23s Solution vs True Solution")
legend("ode23s Solution", "True Solution")
xlabel("t")
ylabel("u_1(t)")

function [x, u] = fourthOrderRK(f, u0, h, N)
    x = (h .* (0:N)).';
    u1 = [u0(1) zeros(1, N)];
    u2 = [u0(2) zeros(1, N)];
    u = [u1; u2];
    for j = 2:(N+1)
        k1 = f(x(j-1), u(:, j-1));
        k2 = f(x(j-1) + (h/2), u(:, j-1) + (h/2).*k1);
        k3 = f(x(j-1) + (h/2), u(:, j-1) + (h/2).*k2);
        k4 = f(x(j-1) + h, u(:, j-1) + h.*k3);
        u(:, j) = u(:, j-1) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
end