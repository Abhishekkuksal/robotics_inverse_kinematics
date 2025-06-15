clc; clf;

a1 = 7; a2 = 7;
n = 6; %no of intermediate steps

%-----------------------------Co-Ordinates of letter 'A'--------------%

% ##x1=0; y1=0;
% ##x2=5; y2=10;
% ##x3=10; y3=0;
% ##x4=7.5; y4=5;
% ##x5=2.5; y5=5;

x1 = [0 5]; y1 = [0 10];
x2 = [5 10]; y2 = [10 0];
x3 = [10 7.5]; y3 = [0 5];
x4 = [7.5 2.5]; y4 = [5 5];
line(x1,y1,'Color','red','LineStyle','--');
line(x2,y2,'Color','red','LineStyle','--');
line(x4,y4,'Color','red','LineStyle','--');
% ##pause(10);

%-----------------------------Plotting My Alphabet------------------%
% ##plot(x1,y1,x2,y2,'--',x3,y3,':',x4,y4,'k',x5,y5,'-k')
% ##hold on;
% ##pause(5);
%------------------------------Zeroing Everything---------------------%
pXArray = zeros (n,1); pYArray = zeros(n,1); aXArray = zeros (n,1); aYArray = zeros(n,1);
theta1Array = zeros(n,1); theta2Array = zeros (n,1);
timeStart = 0; timeEnd = 1;
timeArray=linspace (timeStart, timeEnd,n);
%------------------------------Zeroing Everything---------------------%



%-------1st Animation Starts-----------------------------------------------%
p1=[0 0]; p2 = [5 10]; %#co- ordinate
deltaP=p2 - p1;
tArray=linspace(0,1,n);
for i=1:n
t = tArray(i);
p = p1+ t *deltaP;
pXArray(i) = p(1);
pYArray(i) = p(2);
end

check='1st Animation'
[~, ~] = IKin2Rfunction (a1, a2, pXArray(1), pYArray(1));
for i=1:n
px = pXArray(i); py=pYArray(i);
[theta1_1, theta2_1] = IKin2Rfunction (a1, a2, px, py);
aXArray(i)=a1*cos(theta1_1); aYArray(i)=a1*sin(theta1_1);
theta1Array(i)= theta1_1; theta2Array(i)= theta2_1;
end
xAxisArrayXcoords = [10 -10];
xAxisArrayYcoords = [0 0];
yAxixArrayXcoords = [0 0];
yAxisArrayYcoords = [-10 10];
%##plot(xAxisArrayXcoords, xAxisArrayYcoords,'r', yAxixArrayXcoords, yAxisArrayYcoords,'g')
%##hold on
for i=1:n
   ax= aXArray(i); ay=aYArray(i); bx=pXArray(i); by=pYArray(i);
   traceXArray(i) = bx; traceYArray(i) = by;
   link1Xcoords= [0 ax]; link1Ycoords = [0 ay];
   link2Xcoords= [ax bx]; link2Ycoords = [ay by];

line(x1,y1,'Color','red','LineStyle','--');
line(x2,y2,'Color','red','LineStyle','--');
line(x4,y4,'Color','red','LineStyle','--');
hold on
plot(xAxisArrayXcoords, xAxisArrayYcoords,'r', yAxixArrayXcoords, yAxisArrayYcoords,'g')
hold on
plot(link1Xcoords, link1Ycoords, 'b', link2Xcoords, link2Ycoords, 'c')
hold on
%##plot(x1,y1,x2,y2,'--',x3,y3,':',x4,y4,'k',x5,y5,'-k')
%##plot(traceXArray, traceYArray,'-.k')
hold on
pause(0.05)
hold off
end
%---------------------1st Animation Ends-----------------------------%

%--------------------2nd Animation Starts-------------------------------%
p1=[5 10]; p2 = [10 0]; %#co- ordinate
deltaP=p2 - p1;
%##tArray=linspace(0,1,n);
for i=1:n
t = tArray(i);
p = p1+ t *deltaP;
pXArray(i) = p(1);
pYArray(i) = p(2);
end
[~, ~] = IKin2Rfunction (a1, a2, pXArray(1), pYArray(1));
for i=1:n
px = pXArray(i); py=pYArray(i);
[theta1_1, theta2_1] = IKin2Rfunction (a1, a2, px, py);
aXArray(i)=a1*cos(theta1_1); aYArray(i)=a1*sin(theta1_1);
theta1Array(i)= theta1_1; theta2Array(i)= theta2_1;
end


for i=1:n
   ax= aXArray(i); ay=aYArray(i); bx=pXArray(i); by=pYArray(i);
   traceXArray(i) = bx; traceYArray(i) = by;
   link1Xcoords= [0 ax]; link1Ycoords = [0 ay];
   link2Xcoords= [ax bx]; link2Ycoords = [ay by];
   line(x1,y1,'Color','red','LineStyle','--');
line(x2,y2,'Color','red','LineStyle','--');
line(x4,y4,'Color','red','LineStyle','--');
hold on
plot(xAxisArrayXcoords, xAxisArrayYcoords,'r', yAxixArrayXcoords, yAxisArrayYcoords,'g')
hold on
plot(link1Xcoords, link1Ycoords, 'b', link2Xcoords, link2Ycoords, 'c')
hold on
%##plot(traceXArray, traceYArray,'-.k')
hold on
pause(0.05)
hold off
end
%------------------------2nd Animation Ends---------------------------%

%--------------------------3rd Animation Starts------------------------------%
p1=[10 0]; p2 = [7.5 5]; %#co- ordinate
deltaP=p2 - p1;
%##tArray=linspace(0,1,n);
for i=1:n
t = tArray(i);
p = p1+ t *deltaP;
pXArray(i) = p(1);
pYArray(i) = p(2);
end
[theta1_1, theta2_1] = IKin2Rfunction (a1, a2, pXArray(1), pYArray(1));
for i=1:n
px = pXArray(i); py=pYArray(i);
[theta1_1, theta2_1] = IKin2Rfunction (a1, a2, px, py);
aXArray(i)=a1*cos(theta1_1); aYArray(i)=a1*sin(theta1_1);
theta1Array(i)= theta1_1; theta2Array(i)= theta2_1;
end

for i=1:n
   ax= aXArray(i); ay=aYArray(i); bx=pXArray(i); by=pYArray(i);
   traceXArray(i) = bx; traceYArray(i) = by;
   link1Xcoords= [0 ax]; link1Ycoords = [0 ay];
   link2Xcoords= [ax bx]; link2Ycoords = [ay by];
   line(x1,y1,'Color','red','LineStyle','--');
line(x2,y2,'Color','red','LineStyle','--');
line(x4,y4,'Color','red','LineStyle','--');
hold on
plot(xAxisArrayXcoords, xAxisArrayYcoords,'r', yAxixArrayXcoords, yAxisArrayYcoords,'g')
hold on
plot(link1Xcoords, link1Ycoords, 'b', link2Xcoords, link2Ycoords, 'c')
hold on
%##plot(traceXArray, traceYArray,'-.k')
hold on
pause(0.05)
hold off
end
%----------------------------------3rd Animation Ends-----------------------%

%----------------------------------4th Animation Starts----------------------%
p1=[7.5 5]; p2 = [2.5 5]; %#co- ordinate
deltaP=p2 - p1;
% ##tArray=linspace(0,1,n);
for i=1:n
t = tArray(i);
p = p1+ t *deltaP;
pXArray(i) = p(1);
pYArray(i) = p(2);
end
[theta1_1, theta2_1] = IKin2Rfunction (a1, a2, pXArray(1), pYArray(1));
for i=1:n
px = pXArray(i); py=pYArray(i);
[theta1_1, theta2_1] = IKin2Rfunction (a1, a2, px, py);
aXArray(i)=a1*cos(theta1_1); aYArray(i)=a1*sin(theta1_1);
theta1Array(i)= theta1_1; theta2Array(i)= theta2_1;
end
for i=1:n
   ax= aXArray(i); ay=aYArray(i); bx=pXArray(i); by=pYArray(i);
   traceXArray(i) = bx; traceYArray(i) = by;
   link1Xcoords= [0 ax]; link1Ycoords = [0 ay];
   link2Xcoords= [ax bx]; link2Ycoords = [ay by];
   line(x1,y1,'Color','red','LineStyle','--');
line(x2,y2,'Color','red','LineStyle','--');
line(x4,y4,'Color','red','LineStyle','--');
hold on
plot(xAxisArrayXcoords, xAxisArrayYcoords,'r', yAxixArrayXcoords, yAxisArrayYcoords,'g')
hold on
plot(link1Xcoords, link1Ycoords, 'b', link2Xcoords, link2Ycoords, 'c')
hold on
% ##plot(traceXArray, traceYArray,'-.k')
hold on
pause(0.05)
hold off
end
%------------------------4th Animation Ends----------------------------------%




function [theta1_1, theta2_1] = IKin2Rfunction(a1, a2, px, py)
c2=(px*px + py*py - a1*a1 - a2*a2)/(2*a1*a2);
%if c2>1, px and py lies outside the workspace of this robot.
  if c2<= 1

    s2_1=sqrt(1-c2*c2);
    theta2_1=atan2(s2_1,c2); %first solution
    angle_1=a1*a1 + a2*a2 + 2*a1*a2*cos(theta2_1);

    s1_1= py*(a1+a2*cos(theta2_1)) - px*a2*sin(theta2_1);
    c1_1= px*(a1+a2*cos(theta2_1)) - py*a2*sin(theta2_1);

    theta1_1 = atan2(s1_1, c1_1); %second solution

  end
end
