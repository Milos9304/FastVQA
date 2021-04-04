clc
clear
close all

lattices = dir('../*.l');

for i=1:1%length(lattices)
    
   path = [lattices(i).folder '/' lattices(i).name];
   B = readmatrix(path... % filename
    ,'LineEnding',{']'}... % defines ']' as the end of each row instead of a carriage return
    ,'Delimiter',{'[',' ','\r','\n'}... % define the remaining non-numeric characters as delimiters
    ,'ConsecutiveDelimitersRule','join'... % treat consecutive delimiters as one
    ,'LeadingDelimitersRule','ignore'... % ignore delimiters that start a line
    ,'FileType','text'...
    );

    G = B * transpose(B);
    dim = size(G,1);    
   
    X=sym('X', [1 dim]);
    Z=sym('Z', [1 dim]);
    
    L = 1000;
    
    sumXiZi = 0;
    sumXiZj = 0;
    for i=1:dim
       sumXiZi=sumXiZi -L * X(i) * Z(i);
       for j=i+1:dim
           sumXiZj=sumXiZj + L * Z(i) * X(j);
       end
    end
    
    sumXiZj
    
    X*G*transpose(X) + L + sumXiZi;
    
    %https://uk.mathworks.com/help/symbolic/symbolic-summation.html
   
end