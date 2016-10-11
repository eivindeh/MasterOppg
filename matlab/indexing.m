i = 17;

n = mod(i-1,N);
m = floor((i-1)/N);

left = (m*((N+M-1))+n)+(n==0)*(N)
right = m*((N+M-1))+n+1
above = m*((N+M-1))+N+(n+1+mod(m+1,2))/2
below = m*((N+M-1))-(M-n+mod(m+1,2))/2+(m==0)*(N+M-1)*M
