function [] = drawPipe(pos,r,color)
x = linspace(-1,1,100);
k = 0.22;
%function describing shape of pipe
y = r + (k-r)*x.^4;
x = x/2;

X = [x,x;y,-y];
theta = atan((pos(4)-pos(2))/(pos(3)-pos(1)));


V = [(pos(1)+pos(3))/2,(pos(2)+pos(4))/2]';

% generate face from function
Xf=[X(1,1:100),fliplr(X(1,1:100))];                
Yf=[X(2,1:100),fliplr(X(2,101:200))]; 
%rotate face
Xff = Xf*cos(theta)-Yf*sin(theta);
Yff = Xf*sin(theta)+Yf*cos(theta);
%move face
Xff = Xff+V(1);
Yff = Yff+V(2);
%fill face
fill(Xff,Yff,color,'EdgeColor','none');
hold on;
end

