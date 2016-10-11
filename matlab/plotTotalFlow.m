TotalFlow = importdata('TotalFlow.txt');
figure(1)
hold on;
plot(TotalFlow);
avgFlow = importdata('AvgFlow.txt');
figure(2)
hold on;
plot(avgFlow);
Flow = importdata('Flow.txt');
figure(3)
hold on;
plot(Flow)