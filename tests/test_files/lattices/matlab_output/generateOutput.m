clc
clear
close all

lattices = dir('../*.l');

num_qubits_per_x = 2;
L = 1000; %penalty

for k=1:length(lattices)
    
   path = [lattices(k).folder '/' lattices(k).name];
   B = readmatrix(path... % filename
    ,'LineEnding',{']'}... % defines ']' as the end of each row instead of a carriage return
    ,'Delimiter',{'[',' ','\r','\n'}... % define the remaining non-numeric characters as delimiters
    ,'ConsecutiveDelimitersRule','join'... % treat consecutive delimiters as one
    ,'LeadingDelimitersRule','ignore'... % ignore delimiters that start a line
    ,'FileType','text'...
    );

    G = B * transpose(B);
    dim = size(G,1) * num_qubits_per_x;    
   
    X=sym('x', [1 dim]);
    Z=sym('z', [1 dim]);
    
    sumXiZi = 0;
    sumXiZj = 0;
    for i=1:dim
       sumXiZi=sumXiZi -L * X(i) * Z(i);
       for j=i+1:dim
           sumXiZj=sumXiZj + L * Z(i) * X(j);
       end
    end
    
    f=X(1:size(G,1))*G*transpose(X(1:size(G,1)));

    penalized_f = f + L + sumXiZi + sumXiZj;
    
    subs1 = subs(penalized_f, Z(1), 1);
    subs2 = subs(subs1, Z(2), X(2));    
   
    subs_sq = subs2;

    if num_qubits_per_x > 1 % non binary x's
        for i=1:size(G,1)
            subs_sq = subs(expand(subs_sq), X(i), )
        end        
    end
    
    %expand(subs_sq)
    
    for i=1:dim % get rid of squares
        subs_sq=subs(expand(subs_sq), X(i)^2, X(i));
    end
    
    final_exp = subs_sq;
    vars=sort(symvar(final_exp));
    
    num_qubits = size(vars, 2);
    Q=sym('Z', [1 num_qubits]);
    
    for i=1:num_qubits
       final_exp=subs(final_exp, vars(i), (1-Q(i))/2); 
    end
    
    fpath = [lattices(k).folder '/matlab_output/' lattices(k).name(1:end-2) '_' num_qubits_per_x '.txt'];
    fid = fopen(fpath,'wt');
    fprintf(fid,'%s', vpa(expand(final_exp)));
    fclose(fid);    
     
end