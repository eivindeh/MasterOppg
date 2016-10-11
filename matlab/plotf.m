f = importdata('f.txt');
x = linspace(0,1,1000);
plot(x,f);
hold on;
plot(x,-sin(x)+sin(1)*x);