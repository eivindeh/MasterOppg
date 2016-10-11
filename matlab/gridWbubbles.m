clf;
N = 6*3;
M = 4*3;
K = 104;
data = importdata('Output.txt');
bubbleStartData = importdata('bStart.txt');
bubbleStopData  = importdata('bStop.txt');
bubbleStart = zeros(N*M*3/2,15,K);
bubbleStop = zeros(N*M*3/2,15,K);
for i = 1:N*M*3/2
    for j = 1:K
    bubbleStart(:,:,j) = bubbleStartData((j-1)*N*M*3/2+1:j*N*M*3/2,:);
    bubbleStop(:,:,j) = bubbleStopData((j-1)*N*M*3/2+1:j*N*M*3/2,:);
    end
end

Q = data(1,:);
X = data(2:5,:);

Qmax = max(abs(Q));
figure(1)
for i = 1:N*M/2*3
    color = [0.5,0.8,1];
    if(abs(X(1,i)-X(3,i))>1)
         drawPipe(X(:,i)'+[N*sqrt(3)/2,0,0,0],0.1,color);
    elseif(abs(X(2,i)-X(4,i))>1)
         drawPipe(X(:,i)'+[0,0,0,M*1.5],0.1,color);
    else
        drawPipe(X(:,i)',0.1,color);
    end     
end

for i = 1:N*M*3/2
    if(X(2,i)-X(4,i) == 0.5)
        temp = X(1:2,i);
        X(1:2,i) = X(3:4,i);
        X(3:4,i) = temp;
    end
end
l = 1;


for k = 1:K
    l
if(k ~= 1)
    for i = 1:l-1
        set(P(i),'Visible','off');
    end
end
l = 1;
for i = 1:N*M*3/2
    color = [1,0,0];
    for j = 1:15
        if(bubbleStart(i,j,k) ~= bubbleStop(i,j,k))
            if(abs(X(1,i)-X(3,i))>1)
                    if(X(1,i)-X(3,i)<0)
                        P(l)=drawBubble(X(:,i)'+[N*sqrt(3)/2,0,0,0],0.1,color,bubbleStart(i,j,k),bubbleStop(i,j,k)); 
                    else
                        P(l)=drawBubble(X(:,i)'+[0,0,N*sqrt(3)/2,0],0.1,color,bubbleStart(i,j,k),bubbleStop(i,j,k));
                    end
            elseif(abs(X(2,i)-X(4,i))>1)
                    P(l)=drawBubble(X(:,i)'+[0,0,0,M*1.5],0.1,color,bubbleStart(i,j,k),bubbleStop(i,j,k));
            else
                    P(l)=drawBubble(X(:,i)',0.1,color,bubbleStart(i,j,k),bubbleStop(i,j,k));
            end
            l= l+1;
        end
    end
end
for m = 1:l-1
    set(P(m),'Visible','on');
end
drawnow;
pause(0.0001);
end