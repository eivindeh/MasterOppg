clf;
R = importdata('R.txt');
N = sqrt(length(R));
Q = importdata('Q.txt');
R = R.^(1/4);
Qmax = max(max(abs(Q)));
%R(1:N,N^2-N+1:N^2)=zeros(N);
%R(N^2-N+1:N^2,1:N) = zeros(N);

width = 0.06;
for j= 1:N^2
    for i = 1: N^2
        xi = mod(i-1,N)+1;
        yi = ceil(i/N);
        xj = mod(j-1,N)+1;
        yj = ceil(j/N);
        if(R(j,i) ~= 0)
            width = R(j,i)*0.08;
            if(yi == yj && (xj -xi == 1 || (xi == N && xj == 1)))
                rectangle('Position',[xi+width,yi-width,1-2*width,2*width],'Curvature',[1,0.5],'FaceColor',[1-abs(Q(j,i))/(Qmax),0,abs(Q(j,i))/(Qmax)])
            end
            if(xi == xj && ((yj-yi)==1 || (yi == N && yj == 1)))
                rectangle('Position',[xi-width,yi+width,2*width,1-2*width],'Curvature',[0.5,1],'FaceColor',[1-abs(Q(j,i))/(Qmax),0,abs(Q(j,i))/(Qmax)]);
            end
        end
    end
    
end