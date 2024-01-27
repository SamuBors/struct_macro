function [s]=to_latex_sym(M,P)

% Get matrix dimensions
[m,n]=size(M);
c=cell(m,n);
for k=1:m
    for l=1:n
        if M(k,l)==0
            if k>l
                c{k,l}='';
            else
                c{k,l}=[num2str(M(k,l))];
            end
        elseif round(M(k,l),3)~=0
            if k>l
                c{k,l}='';
            elseif P(k,l) < 0.01
                c{k,l}=[num2str(round(M(k,l),3)) '^{***}'];
            elseif P(k,l) < 0.05
                c{k,l}=[num2str(round(M(k,l),3)) '^{**}'];
            elseif P(k,l) < 0.1
                c{k,l}=[num2str(round(M(k,l),3)) '^{*}'];
            else
                c{k,l}=[num2str(round(M(k,l),3))];
            end
        else
            if k>l
                c{k,l}='';
            elseif P(k,l) < 0.01
                c{k,l}=[num2str(M(k,l),'%.2e') '^{***}'];
            elseif P(k,l) < 0.05
                c{k,l}=[num2str(M(k,l),'%.2e') '^{**}'];
            elseif P(k,l) < 0.1
                c{k,l}=[num2str(M(k,l),'%.2e') '^{*}'];
            else
                c{k,l}=[num2str(M(k,l),'%.2e')];
            end
        end
    end
end

s = sprintf('  \\begin{bmatrix}\n  ');
for k = 1:m
    for l = 1:n
        s = sprintf('%s %s', s,c{k, l}); % print 3 decimal places, align to 6 characters
        if l < n
            s = sprintf('%s &', s);
        end
    end
    if k < m
        s = sprintf('%s \\cr', s);
    end
    s = sprintf('%s\n  ', s);
end
s = sprintf('%s\\end{bmatrix}\n', s);

end