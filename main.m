goodInput = 0;

while(goodInput == 0)
    part = input('Please select the part of the assignment you wish to run (1, 2, or 3): ');
    if(part==1)
        goodInput = 1;
    elseif(part==2)
        goodInput = 1;
    elseif(part==3)
        goodInput = 1;
    else
        fprintf('You did not enter a valid input. Please read the options and try again. \n')
    end
end

if(part==1)
    fprintf('Electon Modeling Selected! \n')
    pause(1)
    monteCarlo
elseif(part==2)
    fprintf('Collisions with Mean Free Path Selected! \n')
    pause(1)
    monteCarloScattering
elseif(part==3)
    fprintf('Enhancements Selected! \n')
    pause(1)
    monteCarloEnhancements
end
    