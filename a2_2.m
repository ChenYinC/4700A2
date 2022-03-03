close all;
clear all;

%set simulation parameters
L = 20; % x-axis
W = 20; % y-axis
Wb = 8; % hight of the box
Lb = 4; % length of the box
y_bound = 0; %y boundary status: 1=bounded, 0 = unbound
c1 = 1; %conductivity outside the box 
c2 = (10^-2); %conductivity inside the box
%inserted boxes define by interval of x and y axis [x_L x_R y_down y_up]
box = [8 12 12 W; 8 12 0 8]; %incerted boxes: [upper_box; lower_box]

%param prepare
nx = 50;
ny = 50; %num of data points 
%eig_num = 20;
x = linspace(0, L, nx);
y = linspace(0, W, ny);
[xx yy] = meshgrid(x, y);

%main calculation part
cmap = makeCmap(x, y, c1, c2, L, W, nx, ny, Wb, Lb);
[G F] = createG(nx, ny, cmap);
V = solV(G, F, nx, ny);
[dvx dvy] =  gradient(V);
Ex = dvx .* (-1);
Ey = dvy .* (-1);
jx = cmap .* Ex;
jy = cmap .* Ey;
jxC = jx ./ cmap;
jyC = jy ./ cmap; % current / conductivity

%plot the result
figure(1);
subplot(2, 3, 1);
surf(xx, yy, cmap);
title('Conductivity Map');
xlabel('x-axis (L)');
ylabel('y-axis (W)');
zlabel('Condictivity');
subplot(2, 3, 2);
surf(xx, yy, V);
title('Potential Map');
xlabel('x-axis (L)');
ylabel('y-axis (W)');
zlabel('Potential V');
subplot(2, 3, 3);
quiver(Ex, Ey, 3);
title('E Vector');
xlabel('x-axis (L)');
ylabel('y-axis (W)');
subplot(2, 3, 4);
quiver(jx, jy, 3);
title('Current Vector');
xlabel('x-axis (L)');
ylabel('y-axis (W)');
subplot(2, 3, 5);
surf(xx, yy, jxC);
title('Current (Jx) / Conductivity');
xlabel('x-axis (L)');
ylabel('y-axis (W)');
zlabel('Jx / Conductivity');
subplot(2, 3, 6);
surf(xx, yy, jyC);
title('Current (Jy) / Conductivity');
xlabel('x-axis (L)');
ylabel('y-axis (W)');
zlabel('Jy / Conductivity');

%current / mesh size
meshS = [];
jSum = [];
for i=5:1:50
    meshS = [meshS i];
    cmap = makeCmap(x, y, c1, c2, L, W, i, i, Wb, Lb);
    [G F] = createG(i, i, cmap);
    V = solV(G, F, i, i);
    [dvx dvy] =  gradient(V);
    Ex = dvx .* (-1);
    jx = cmap .* Ex;
    jxC = jx ./ cmap;
    jSum = [jSum sum(jxC(:, i))];
end
figure(2);
plot(meshS, jSum);
title('Current(Jx) / mesh size');
xlabel('Mesh Size');
ylabel('Current (Jx)');

%current / bottleneck(height)
width = [];
jSumB = [];
for i=1:9
    width = [width i];
    cmap = makeCmap(x, y, c1, c2, L, W, nx, ny, i, Lb);
    [G F] = createG(nx, ny, cmap);
    V = solV(G, F, nx, ny);
    [dvx dvy] =  gradient(V);
    Ex = dvx .* (-1);
    jx = cmap .* Ex;
    jxC = jx ./ cmap;
    jSumB = [jSumB sum(jxC(:, i))];
end
figure(3);
plot(width, jSumB);
title('Current(Jx) / Bottleneck (height)');
xlabel('Bottlenexh Width');
ylabel('Current (Jx)');







%function that get V from G and F then convert it to 2D matrix
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

function cmap = makeCmap(x, y, c1, c2, L, W, nx, ny, Wb, Lb)
    box = [(L/2 - Lb/2) (L/2 + Lb/2) (W - Wb) W; (L/2 - Lb/2) (L/2 + Lb/2) 0 Wb]; %incerted boxes: [upper_box; lower_box]
    cmap = zeros(ny, nx);
    for i=1:nx
        for j=1:ny
                if (x(i) >= box(1, 1)) && (x(i) < box(1, 2))
                    if (y(j) > box(2, 4)) && (y(j) < box(1, 3))
                        cmap(j, i) = c1;
                    else
                        cmap(j, i) = c2;
                    end
                else
                    cmap(j, i) = c1;
                end
        end
    end
end

function [G F] = createG(nx, ny, cmap)
    F = zeros(ny*nx, 1);
    G = zeros(nx*ny, nx*ny);
    %G matrix formation
    for i=1:nx
        for j=1:ny
            ct = j + (i-1)*ny;
            L = j + (i-2)*ny;
            R = j + (i)*ny;
            up = j+1 +(i-1)*ny;
            down = j-1 + (i-1)*ny;
    
            if (i==1) %lefttest
                G(ct, ct) = 1;
                F(ct) = 1;
            elseif (i==nx) %Rightest
                G(ct, ct) = 1;
            elseif (j==1) %bottom
                % 3 resistors: up , left , right
                r_up = (cmap(j, i) + cmap(j+1, i)) / 2;
                r_r = (cmap(j, i) + cmap(j, i+1)) / 2;
                r_l = (cmap(j, i) + cmap(j, i-1)) / 2;
    
                G(ct, ct) = -(r_up + r_r + r_l);
                G(ct, up) = r_up;
                G(ct, R) = r_r;
                G(ct, L) = r_l;
            elseif (j==ny) %top
                % 3 resistors: down, left, right
                r_down = (cmap(j, i) + cmap(j-1, i)) / 2;
                r_r = (cmap(j, i) + cmap(j, i+1)) / 2;
                r_l = (cmap(j, i) + cmap(j, i-1)) / 2;
    
                G(ct, ct) = -(r_down + r_r + r_l);
                G(ct, down) = r_down;
                G(ct, R) = r_r;
                G(ct, L) = r_l;
            else
                % 4 resistors: up, down, left, right
                r_up = (cmap(j, i) + cmap(j+1, i)) / 2;
                r_down = (cmap(j, i) + cmap(j-1, i)) / 2;
                r_r = (cmap(j, i) + cmap(j, i+1)) / 2;
                r_l = (cmap(j, i) + cmap(j, i-1)) / 2;
    
                G(ct, ct) = -(r_up + r_down + r_r + r_l);
                G(ct, up) = r_up;
                G(ct, down) = r_down;
                G(ct, R) = r_r;
                G(ct, L) = r_l;
            end
        end
    end
end