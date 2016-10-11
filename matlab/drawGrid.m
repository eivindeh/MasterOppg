clf;
X = 0;
N = 6;
M = 4;
R = importdata('R.txt');
Q = importdata('Q.txt');
Qmax = max(max(abs(Q)));
%for loop generating coordinates for pipes
%X(1,:) coresponds to X coordinates and X(2,:) are Y.
for i = 1:(M*N)
    x = mod(i-1,N)+1;
    y = ceil(i/N);
    X(1,i) = x*sqrt(3)/2;
    if (~mod(M,2) && mod(y,2))
        X(2,i) = 1.5*y +(mod(i,2))/2;
    else
        X(2,i) = 1.5*y+(1-mod(i,2))/2;
    end
end
%Draq pipes
for i = 1:N*M
    for j = 1:N*M
        if (R(j,i) ~= 0 && i>j)
            if(X(1,i)-X(1,j)>1)
                drawPipe([X(1,i),X(2,i),X(1,j)+N*sqrt(3)/2,X(2,j)],0.1,(1-abs(Q(j,i))/(Qmax))*[0.5,0.8,1+(abs(Q(j,i))/(Qmax))]);
            elseif(X(2,i)-X(2,j)>1)
                %Horizontal boundary
                drawPipe([X(1,i),X(2,i),X(1,j),X(2,j)+M*1.5],0.1,(1-abs(Q(j,i))/(Qmax))*[0.5,0.8,1+(abs(Q(j,i))/(Qmax))]);
            else
                %vertical boundary
                drawPipe([X(1,i),X(2,i),X(1,j),X(2,j)],0.1,(1-abs(Q(j,i))/(Qmax))*[0.5,0.8,1+(abs(Q(j,i))/(Qmax))]);
            end
         end
    end
end