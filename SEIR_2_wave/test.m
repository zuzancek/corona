figure;
drawArrow = @(x,y) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0, 'HeadWidth', 10 );   

x1 = [10 30];
y1 = [10 30];

drawArrow(x1,y1); hold on

x2 = [25 15];
y2 = [15 25];

drawArrow(x2,y2); hold on;

x3 = [20 20];
y3 = [26 16];

drawArrow(x3,y3)

f=quiver( x2(1),y2(1),x2(2)-x2(1),y2(2)-y2(1),0 );