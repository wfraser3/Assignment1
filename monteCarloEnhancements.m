%{
ELEC 4700 Assignment 1
William Fraser
101001393
%}
colours = ['b','g','r','c','m','y','k'];
maxTraces = length(colours);
numparticles = input("Please input the number of particles: ");
plotting = input('Type 1 to do live plotting, 0 for no live plotting: ');

if(numparticles < maxTraces)
    numTraces = numparticles;
else
    numTraces = maxTraces;
end

kb = 1.38064852e-23;
m0 = 9.11e-31;
m = 0.26*m0;
T = 300;
vth = sqrt((kb*T)/m);
vth = vth*1e9;
vth = vth*1e-13; %nm/0.1ps
vth2 = vth^2;

velocities = linspace(0,50,102);
v2 = velocities.^2;
coef1 = 4*pi;
kbMax = kb*1e18*1e-26; %Fixing units
coef2 = (m0/(2*pi*kbMax*T));
exponent = exp(-(m0*v2)/(2*kbMax*T));
maxBol = coef1*coef2*v2.*exponent*numparticles;
velAssigner = zeros(1,round(sum(maxBol)));
currentIndex = 1;

for i = 1:length(maxBol)
    if(round(maxBol(i))~=0)
        stopPoint = currentIndex + round(maxBol(i))-1;
        velAssigner(currentIndex:stopPoint) = velocities(i);
    else
        stopPoint = currentIndex;
    end
    currentIndex = stopPoint + 1;
end

velAssigner = velAssigner(randperm(length(velAssigner)));

grid1 = zeros(100,200);
particles = zeros(numparticles,4); %Linear Index; row; column; X velocity; Y velocity
if(numparticles <= 100)
    particles(:,2) = randperm(100,numparticles);
    particles(:,3) = randperm(200,numparticles);
else
    start = 1;
    stop = 100;
    endReached = 0;
    while(endReached==0)
        if(stop<=numparticles)
            particles(start:stop,2) = randperm(100,100);
            particles(start:stop,3) = randperm(200,100);
            start = start + 100;
            stop = stop + 100;
        elseif(stop>=numparticles)
            particles(start:numparticles,2) = randperm(100,(numparticles-start+1));
            particles(start:numparticles,3) = randperm(100,(numparticles-start+1));
            endReached = 1;
        end 
    end
end
inBoxX1 = particles(:,3) >= 75;
inBoxX2 = particles(:,3) <= 125;
inBoxX = inBoxX1.*inBoxX2;
inBoxY1 = particles(:,2) <= 40;
inBoxY2 = particles(:,2) >= 60;
inBox1 = inBoxY1.*inBoxX;
inBox2 = inBoxY2.*inBoxX;
inBox = inBox1 + inBox2;
notInBox = inBox==0;
particles(:,2) = (particles(:,2).*notInBox) + (inBox*50);
particles(:,3) = (particles(:,3).*notInBox) + (inBox*100);
particles(:,1) = sub2ind(size(grid1),particles(:,2),particles(:,3));
traces = zeros(numTraces,2);
traces(:,1) = particles(1:numTraces,1);
maxIndex = length(velAssigner);

for i = 1:length(particles(:,1))
    choice = randperm(maxIndex,1);
    tempVel = velAssigner(choice);
    tempVel2 = tempVel^2;
    xRat = rand;
    xDir = (-1)^(round(rand));
    yDir = (-1)^(round(rand));
    yRat = 1 - xRat;
    xVel = (sqrt(xRat*tempVel2));
    yVel = (sqrt(yRat*tempVel2));
    particles(i,4) = xVel*xDir;
    particles(i,5) = yVel*yDir;
end

grid1(particles(:,1)) = 1;
gridSize = size(grid1);

pScat = 1 - exp(-0.1/0.2);

clear choice coef1 coef2 currentIndex exponent i kb m0 maxBol normalize
clear stopPoint T tempVel tempVel2 v2 vth vth2 velocities histVel maxBolNorm

scatterMatrix = zeros(numparticles,1000);
averageVel = zeros(1,1000);

boxX = [75 75 125 125];
box1y = [0 40 40 0];
box2y = [100 60 60 100];

if(plotting)
    figure('units','normalized','outerposition',[0 0 1 1])
    f(1) = subplot(1,2,1);
    set(f(1),'position',[0.05 0.4 0.4 0.35])
    f(2) = subplot(1,2,2);
    set(f(2),'position',[0.55,0.4,0.4,0.35])
    handles = findobj(figure(1),'Type','axes');
    axes(handles(2))
    xlim([0 200])
    ylim([0 100])
    set(gca, 'YDir','reverse')
end

for time = 1:1000
    squaredVel = (particles(:,4).^2) + (particles(:,5).^2);
    meanVel = mean(squaredVel);
    averageVel(time) = sqrt(meanVel);
    temperature = (m*meanVel)/kbMax;
    display('Time = ',num2str(time));
    if(time~=1)
        traces(:,1) = traces(:,2);
    end
    grid2 = zeros(size(grid1));
    scatVector = rand(1,numparticles);
    willScatter = scatVector <= pScat;
    scatterMatrix(:,time) = willScatter;
    noScatter = willScatter == 0;
    reTherm = randperm(maxIndex,numparticles);
    reTherm = (velAssigner(reTherm).^2).*willScatter;
    xRatio = rand(1,numparticles);
    yRatio = abs(xRatio-1);
    xDir = (-1).^(round(rand(1,numparticles)));
    yDir = (-1).^(round(rand(1,numparticles)));
    xDir = transpose(xDir);
    yDir = transpose(yDir);
    noScatter = transpose(noScatter);
    willScatter = transpose(willScatter);
    reTherm = transpose(reTherm);
    xRatio = transpose(xRatio);
    yRatio = transpose(yRatio);
    particles(:,4) = (particles(:,4).*noScatter) + (willScatter.*sqrt(xRatio.*reTherm).*xDir);
    particles(:,5) = (particles(:,5).*noScatter) + (willScatter.*sqrt(yRatio.*reTherm).*yDir);
    inBoxX1 = particles(:,3) >= 75;
    inBoxX2 = particles(:,3) <= 125;
    inBoxX = inBoxX1.*inBoxX2;
    inBoxY1 = particles(:,2) <= 40;
    inBoxY2 = particles(:,2) >= 60;
    inBoxY = inBoxY1.*inBoxY2;
    enter1Left = particles(:,3)+particles(:,4)>=75;
    enter1Left = enter1Left.*inBoxX2.*inBoxY1;
    enter1Right = particles(:,3)+particles(:,4)<=125;
    enter1Right = enter1Right.*inBoxX1.*inBoxY1;
    enter1Sides = enter1Right + enter1Left;
    not1Sides = enter1Sides==0;
    particles(:,4) = (particles(:,4).*not1Sides) - (particles(:,4).*enter1Sides);
    enter2Left = particles(:,3)+particles(:,4)>=75;
    enter2Left = enter2Left.*inBoxX2.*inBoxY2;
    enter2Right = particles(:,3)+particles(:,4)<=125;
    enter2Right = enter2Right.*inBoxX1.*inBoxY2;
    enter2Sides = enter2Right + enter2Left;
    not2Sides = enter2Sides==0;
    particles(:,4) = (particles(:,4).*not2Sides) - (particles(:,4).*enter2Sides);
    enter1Bottom = particles(:,2)+particles(:,5)<=40;
    enter1Bottom = enter1Bottom.*inBoxX1.*inBoxX2;
    enter2Top = particles(:,2)+particles(:,5)>=60;
    enter2Top = enter2Top.*inBoxX1.*inBoxX2;
    enterHorizontal = enter1Bottom+enter2Top;
    noEnterHorizontal = enterHorizontal==0;
    particles(:,5) = (particles(:,5).*noEnterHorizontal) - (particles(:,5).*enterHorizontal);
    particles(:,3) = round(particles(:,3) + particles(:,4));
    particles(:,2) = round(particles(:,2) + particles(:,5));
    yBoundMaxed = particles(:,2) > gridSize(1);
    yBoundMined = particles(:,2) < 0;
    yIsZero = particles(:,2)==0;
    yBounds = yBoundMaxed + yBoundMined + yIsZero;
    keepY = yBounds==0;
    particles(:,2) = (particles(:,2).*keepY) - (yBoundMined.*particles(:,2)) + (yBoundMaxed.*(gridSize(1)-(particles(:,2)-gridSize(1)))) - round(yIsZero.*particles(:,5));
    particles(:,5) = (particles(:,5).*keepY) - (yBoundMined.*particles(:,5)) - (yBoundMaxed.*particles(:,5)) - (particles(:,5).*yIsZero);
    xBoundMaxed = particles(:,3) > gridSize(2);
    xBoundMined = particles(:,3) < 0;
    xIsZero = particles(:,3)==0;
    xBounds = xBoundMaxed + xBoundMined + xIsZero;
    keepX = xBounds==0;
    particles(:,3) = (particles(:,3).*keepX) + (xBoundMaxed.*(particles(:,3)-gridSize(2))) + (xBoundMined.*(particles(:,3)+gridSize(2))) + (xIsZero.*round(gridSize(2)+particles(:,4)));
    particles(:,1) = sub2ind(gridSize,particles(:,2),particles(:,3));
    traces(:,2) = particles(1:numTraces,1);
    holder = traces(:,2);
    jump = abs(traces(:,1) - traces(:,2)) > 10000;
    noJump = jump == 0;
    temp = traces(:,1).*noJump;
    traces(:,1) = temp + (jump.*traces(:,2));
    if(plotting)
        if(time~=1)
            [y1,x1] = ind2sub(gridSize,traces(:,1));
            [y2,x2] = ind2sub(gridSize,traces(:,2));
            x(1,:) = x1;
            x(2,:) = x2;
            y(1,:) = y1;
            y(2,:) = y2;

            axes(handles(1))
            xlim([0 200])
            ylim([0 100])
            title('Electron Path Tracing')
            ylabel('Y Position (nm)')
            xlabel('X Position (nm)')
            set(gca, 'YDir','reverse')
            plot(boxX,box1y,'k',boxX,box2y,'k')
            hold on
            if(numTraces==7)
                plot(x(:,1),y(:,1),colours(1),x(:,2),y(:,2),colours(2),x(:,3),y(:,3),colours(3),x(:,4),y(:,4),colours(4),x(:,5),y(:,5),colours(5),x(:,6),y(:,6),colours(6),x(:,7),y(:,7),colours(7))
            elseif(numTraces==6)
                plot(x(:,1),y(:,1),colours(1),x(:,2),y(:,2),colours(2),x(:,3),y(:,3),colours(3),x(:,4),y(:,4),colours(4),x(:,5),y(:,5),colours(5),x(:,6),y(:,6),colours(6))
            elseif(numTraces==5)
                plot(x(:,1),y(:,1),colours(1),x(:,2),y(:,2),colours(2),x(:,3),y(:,3),colours(3),x(:,4),y(:,4),colours(4),x(:,5),y(:,5),colours(5))
            elseif(numTraces==4)
                plot(x(:,1),y(:,1),colours(1),x(:,2),y(:,2),colours(2),x(:,3),y(:,3),colours(3),x(:,4),y(:,4),colours(4))
            elseif(numTraces==3)
                plot(x(:,1),y(:,1),colours(1),x(:,2),y(:,2),colours(2),x(:,3),y(:,3),colours(3))
            elseif(numTraces==2)
                plot(x(:,1),y(:,1),colours(1),x(:,2),y(:,2),colours(2))
            else
                plot(x(:,1),y(:,1),colours(1))
            end
        end
        currentTime = time/10;
        grid2(particles(:,1)) = 1;
        grid1 = grid2;
        axes(handles(2))
        plot(boxX,box1y,'k',boxX,box2y,'k')
        hold on
        spy(grid1,5)
        hold off
        title(['Temperature = ',num2str(temperature),'K at Time = ',num2str(currentTime),'ps'])
        xlabel('X Position (nm)')
        ylabel('Y Position (nm)')
        pause(0.0001)
        clear grid2
    end
end

h = ones(3,3);
density = imfilter(grid1,h);
figure(2)
imagesc(density)
colorbar
title('Electron Position Density (Electrons/nm^2)')

electronTemperatures = (squaredVel.*m)./kbMax;
tempGrid = zeros(100,200);
tempGrid(particles(:,1)) = electronTemperatures;
figure(3)
imagesc(tempGrid,[0 500])
colorbar
title('Temperature Map (K)')