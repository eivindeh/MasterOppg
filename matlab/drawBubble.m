function [patch] = drawBubble(pos,r,color,start,stop)
x = linspace(2*start-1,2*stop-1,20);
k = 0.22;
%function describing shape of pipe
y = r + (k-r)*x.^4;
x = x/2;

X = [x,x;y,-y];
theta = atan2((pos(4)-pos(2)),(pos(3)-pos(1)));


V = [(pos(1)+pos(3))/2,(pos(2)+pos(4))/2]';

% generate face from function
Xf=[X(1,1:100/5),fliplr(X(1,1:100/5))];                
Yf=[X(2,1:100/5),fliplr(X(2,100/5+1:200/5))]; 
%rotate face
Xff = Xf*cos(theta)-Yf*sin(theta);
Yff = Xf*sin(theta)+Yf*cos(theta);
%move face
Xff = Xff+V(1);
Yff = Yff+V(2);
%fill face
patch = fill(Xff,Yff,color,'EdgeColor','none','Visible','off');
end

