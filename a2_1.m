close all;
clear all;

% !!! Set y_bound to 0 for part 1a) 1 for part 1b)

%set simulation parameters
L = 20; % x-axis
W = 20; % y-axis
v0x = 1; % BC of x axis (2 sides)
v0y = 0;
y_bound = 0; %y boundary status: 1=bounded, 0 = unbound ( for part 1_a): set to 0, part 1_b): set to 1 )
R1 = 1; %resistance outside the box (= 1/ conductivity)
R2 = 1/(10^-2); %reisitance inside the box, 10^-2 = counductivity
%inserted boxes define by interval of x and y axis [x_L x_R y_down y_up]
box_num = 2;
box = [8 12 12 W; 8 12 0 8];

%param prepare
nx = 50;
ny = 50; %num of data points 
%eig_num = 20;
x = linspace(0, L, nx);
y = linspace(0, W, ny);
[xx yy] = meshgrid(x, y);
F = zeros(ny*nx, 1);
G = zeros(nx*ny, nx*ny);
V_fd = zeros(ny, nx);


%V equation solution 
V_eq = zeros(ny, nx);
a = W;
b = L/2;
v0=v0x;
iter = 100; %num of iteration
n = 1;
for i=1:iter
    n = [n (2*i+1)];
end
V_eq = zeros(ny, nx);
for i=1:2:iter
    for j=1:nx
        for k=1:ny
            temp1 = cosh(i .* pi .* (x(j)-(L/2)) ./ a) / cosh(i .* pi .* b ./ a);
            temp2 = sin(i .* pi .* y(k)/a) ./ i;
            V_eq(k, j) = V_eq(k, j)  + (4*v0 / pi) * temp1.*temp2;
        end
    end
    figure(2);
    subplot(1, 2, 1);
    surf(xx, yy, V_eq);
    title('V from equations (100 iteration');
    xlabel('x-axis (L)');
    ylabel('y-axis (W)');
    zlabel('Potential V');
end

% FD method to find potential
%G matrix for potential
for i=1:(nx)
    for j=1:(ny)
        index = j + (i-1)*ny;
        L = j + (i-2)*ny;
        R = j + (i)*ny;
        up = j+1 +(i-1)*ny;
        down = j-1 + (i-1)*ny;      
        %let dx^2 = 1
        if(i==1 || i==nx)
            G(index, index) = 1;
            if(y_bound)
                F(index) = 1;
            else
                if(i==1)
                    F(index) = 1;
                end
            end
        elseif(j==1)
            G(index, index) = 1;
            if (y_bound == 1)
                G(index, up) = 1;
            else
                G(index, up) = -1;
            end
        elseif(j==ny)
            G(index, index) = 1;
            if (y_bound == 1)
                G(index, down) = 1;
            else
                G(index, down) = -1;
            end
        else
            G(index, index) = -4;
            G(index, up) = 1;
            G(index, down) = 1;
            G(index, R) = 1;
            G(index, L) = 1;
        end         
    end
end
% GV  = F  V = G\F  V = inv(G)*F;
V_fd = solV(G, F, nx, ny);

%comment this out for part 1_b):
V_x = V_fd(ny/2, :);
figure(1);
plot(x, V_x);
title('2D plot of Vx');
xlabel('x-axis');
ylabel('V');

figure(2);
subplot(1, 2, 2);
surf(xx, yy, V_fd-V_eq);
title('Difference between FD and Eq');
xlabel('x-axis (L)');
ylabel('y-axis (W)');
zlabel('Potential V');
figure(3);
subplot(1, 2, 1);
surf(xx, yy, V_fd);
title('Surf plot of V');
xlabel('x-axis (L)');
ylabel('y-axis (W)');
zlabel('Potential V');
[dvx dvy] = gradient(V_fd);
subplot(1, 2, 2);
quiver(xx, yy, dvx, dvy);
title('E vector');
xlabel('x-axis (L)');
ylabel('y-axis (W)');




% solve potential from G and F
function V = solV(G, F, nx, ny)
    V = zeros(ny, nx);
    vtemp = G\F;
    for i=1:nx
        for j=1:ny
            ind = j + (i-1)*ny;
            V(j, i) = vtemp(ind);
        end
    end
end

