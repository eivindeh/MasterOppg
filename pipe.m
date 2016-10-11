x = linspace(-1,1,100);
r = 0.08;
k = 0.17;

pos = [1,1,2,2];

y = r + (k-r)*x.^4;
plot(x,y);
hold on;
plot(x,-y);

X = [x,x;y,-y];
theta = atan((pos(4)-pos(2))/(pos(3)-pos(1)));
V = [(pos(1)+pos(3))/2,(pos(2)+pos(4))/2]';
X = (X'*[cos(theta),sin(theta);-sin(theta),cos(theta)])';
for i = 1:200
    X(:,i) = X(:,i)+V;
end
plot(X(1,1:100),X(2,1:100),'b',X(1,101:200),X(2,101:200),'b');