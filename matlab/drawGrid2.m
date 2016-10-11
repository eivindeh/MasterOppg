clf;
N = 6*6;
M = 4*6;
data = importdata('Output.txt');
Q = data(1,:);
X = data(2:5,:);
Qmax = max(abs(Q));

for i = 1:N*M/2*3
    color = (1-abs(Q(i))/(Qmax))*[0.5,0.8,1+(abs(Q(i))/(Qmax))];
    if(abs(X(1,i)-X(3,i))>1)
         drawPipe(X(:,i)'+[N*sqrt(3)/2,0,0,0],0.1,color);
    elseif(abs(X(2,i)-X(4,i))>1)
         drawPipe(X(:,i)'+[0,0,0,M*1.5],0.1,color);
    else
        drawPipe(X(:,i)',0.1,color);
    end     
end

% 
% %Draq pipes
% for i = 1:N*M
%     for j = 1:N*M
%         if (R(j,i) ~= 0 && i>j)
%             if(X(1,i)-X(1,j)>1)
%                 drawPipe([X(1,i),X(2,i),X(1,j)+N*sqrt(3)/2,X(2,j)],0.1,(1-abs(Q(j,i))/(Qmax))*[0.5,0.8,1+(abs(Q(j,i))/(Qmax))]);
%             elseif(X(2,i)-X(2,j)>1)
%                 %Horizontal boundary
%                 drawPipe([X(1,i),X(2,i),X(1,j),X(2,j)+M*1.5],0.1,(1-abs(Q(j,i))/(Qmax))*[0.5,0.8,1+(abs(Q(j,i))/(Qmax))]);
%             else
%                 %vertical boundary
%                 drawPipe([X(1,i),X(2,i),X(1,j),X(2,j)],0.1,(1-abs(Q(j,i))/(Qmax))*[0.5,0.8,1+(abs(Q(j,i))/(Qmax))]);
%             end
%          end
%     end
% end