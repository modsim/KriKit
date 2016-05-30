function [output] = tutorialFunction(input)
    part1 = 10*input(:,1)./(5 + input(:,1));
    part2 = sin(input(:,2));
    part3 = input(:,3);
    
    output = part1.*part2 + part3;
end

