%% Monte Carlo υπολογισμός του π
% Charitopoulou Despoina 4405

clear;
clc;
clf;
tic
% για n=100 -> 2 δεκαδικά
% για n=1000 -> 3 δεκαδικά
% για n=10000 έως 10^8 -> 4 δεκαδικά 
% για n πάνω από 10^8 ο υπολογιστής δεν αποκρίνεται
n=10000; 
c=0;
a = -1; % lower limit
b = 1; % upper limit
x=a+(b-a).*rand(1,n); % uniform distribution
y=a+(b-a).*rand(1,n);
z=a+(b-a).*rand(1,n);

plot3(x,y,z,'.','DisplayName','n=10000')
axis equal
xlabel('x(t)')
ylabel('y(t)')
zlabel('z(t)')
hold on 
sphere
axis equal
title('uniform distribution-sphere')

for i=1:n
    if x(1,i)^2+y(1,i)^2+z(1,i)^2<=1
        c(i)=1;
    else
        c(i)=0;
    end

end
m1=sum(c)/n;% μέση τιμή του αριθμού των ανυσμάτων που βρίσκονται μέσα στη σφαίρα 
m2=mean(c);

p=6*m2 % πρόβλεψη του π
fprintf('Η πρόβλεψη του π είναι π = %g\n',p)

err=6*((std(c))/(sqrt(n))); % σφάλμα πρόβλεψης
fprintf('Το σφάλμα του π είναι error = %g\n',err)

fprintf('Η ακριβής τιμή του π αναμένεται να βρίσκεται στη περιοχή π = [%.5f,%.5f] \n',p-err,p+err)
fprintf('Η ακριβής τιμή του π είναι π=%g \n',pi)

toc

