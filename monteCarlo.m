%{
ELEC 4700 Assignment 1
William Fraser
101001393
%}
colours = ['b','g','r','c','m','y','k'];
maxTraces = length(colours);
numparticles = input("Please input the number of particles: ");

if(numparticles < maxTraces)
    numTraces = numparticles;
else
    numTraces = maxTraces;
end

setSpeed = input('Would you like the animation to be fast or slow? (Please only type "f" or "s"): ','s');

if(setSpeed=='f')
    speedFactor = 1e4;
elseif(setSpeed=='s')
    speedFactor = 1e5;
end

kb = 1.38064852e-23;
m0 = 9.11e-31;
m = 0.26*m0;
T = 300;
vth = sqrt((kb*T)/m);
vth2 = vth^2;

grid1 = zeros(100+numparticles,200+numparticles);
particles = zeros(numparticles,4); %Linear Index; row; column; X velocity; Y velocity
particles(:,2) = randperm(100+numparticles,numparticles);
particles(:,3) = randperm(200+numparticles,numparticles);
particles(:,1) = sub2ind(size(grid1),particles(:,2),particles(:,3));
traces = zeros(numTraces,2);
traces(:,1) = particles(1:numTraces,1);

for i = 1:length(particles(:,1))
    xRat = rand;
    xDir = (-1)^(round(rand));
    yDir = (-1)^(round(rand));
    yRat = 1 - xRat;
    xVel = (sqrt(xRat*vth2))/speedFactor;
    yVel = (sqrt(yRat*vth2))/speedFactor;
    particles(i,4) = xVel*xDir;
    particles(i,5) = yVel*yDir;
end

grid1(particles(:,1)) = 1;
gridSize = size(grid1);
f1 = figure(1);
movegui(f1,'west')
spy(grid1,10)

for time = 1:1000
    if(time~=1)
        traces(:,1) = traces(:,2);
    end
    grid2 = zeros(size(grid1));
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
    if(time~=1)
        if(time==2)
            f2 = figure(2);
            movegui(f2,'east')
            hold on
        else
            figure(2)
            hold on
        end
        for loop = 1:numTraces
            difference = abs(traces(loop,1) - traces(loop,2));
            if(difference > 5000)
                continue
            end
            [y1,x1] = ind2sub(gridSize,traces(loop,1));
            [y2,x2] = ind2sub(gridSize,traces(loop,2));
            plot([x1,x2],[y1,y2],colours(loop))
        end
        hold off
    end
    grid2(particles(:,1)) = 1;
    grid1 = grid2;
    figure(1)
    spy(grid1,10)
    pause(0.0001)
    clear grid2
end