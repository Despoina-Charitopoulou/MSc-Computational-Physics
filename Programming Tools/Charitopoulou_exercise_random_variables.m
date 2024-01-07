%% Random Variables
% Charitopoulou Despoina 4405

clear;
clf;
clc;
m=1000;
mu=0; % μέση τιμή 
sigma=1; % τυπική απόκλιση
x = normrnd(mu,sigma,m,1); % κανονική κατανομή

figure(1)
plot(x,'b.');
xlabel('m')
ylabel('x')
title('Plot 1 - 1000 random numbers from normal distribution')

figure(2)
histogram(x);
title('Histogram 1 - normal distribution with mu=0 and sigma=1')

m1=sum(x)/m; % μέση τιμή
m2=mean(x); % μέση τιμή με την εντολή mean
s1=sqrt((sum((x-m1).^2))/(m-1)); % τυπική απόκλιση 
s2=(std(x)); % τυπική απόκλιση με την εντολή std

fprintf('Mean value of x is : %g \n',m2)
fprintf('Standard deviation of x is : %g \n',s2)


% πειραματικός προσδιορισμός μέσης τιμής
y=sum(x,m)/m;
figure(3)
plot(y,'b.');
xlabel('m')
ylabel('y')
title('Plot 2 - distribution of mean value')

figure(4)
histogram(y);
title('Histogram 2 - distribution of mean value')

m3=mean(y);
s3=std(y);
sx=s3/sqrt(m); % σφάλμα

figure(5)
plot(x,'b.','DisplayName','Normal Distribution');
hold on
plot(y,'r.','DisplayName','Distribution of Mean Value');
hold off
xlabel('m')
ylabel('x,y')
title('Plot 1 and Plot 2')

fprintf('Mean value of y is : %g \n',m3)
fprintf('Standard deviation of y is : %g \n',s3)
fprintf('Error of y is : %g \n',sx)
