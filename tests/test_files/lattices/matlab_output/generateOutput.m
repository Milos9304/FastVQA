clc
clear
close all

lattices = dir('../*.l');

for k=1:1%length(lattices)
    
   path = [lattices(k).folder '/' lattices(k).name];
   B = readmatrix(path... % filename
    ,'LineEnding',{']'}... % defines ']' as the end of each row instead of a carriage return
    ,'Delimiter',{'[',' ','\r','\n'}... % define the remaining non-numeric characters as delimiters
    ,'ConsecutiveDelimitersRule','join'... % treat consecutive delimiters as one
    ,'LeadingDelimitersRule','ignore'... % ignore delimiters that start a line
    ,'FileType','text'...
    );

    G = B * transpose(B);
    dim = size(G,1);    
   
    X=sym('x', [1 dim]);
    Z=sym('z', [1 dim]);
    
    L = 1000;
    
    sumXiZi = 0;
    sumXiZj = 0;
    for i=1:dim
       sumXiZi=sumXiZi -L * X(i) * Z(i);
       for j=i+1:dim
           sumXiZj=sumXiZj + L * Z(i) * X(j);
       end
    end
    
    f=X*G*transpose(X);    
    penalized_f = f + L + sumXiZi;
    
    subs1 = subs(penalized_f, Z(1), 1);
    subs2 = subs(penalized_f, Z(2), X(2));
    
    expand(subs2)
    fpath = [lattices(k).folder '/matlab_output/' lattices(k).name(1:end-2) '.txt'];
    fid = fopen(fpath,'wt');
    fprintf(fid,'%s', expand(subs2));
    fclose(fid);    
    
    %https://uk.mathworks.com/help/symbolic/symbolic-summation.html
   
end