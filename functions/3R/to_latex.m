function [s]=to_latex(M,P)

% Get matrix dimensions
[m,n]=size(M);
c=cell(m,n);
for k=1:m
    for l=1:n
        if M(k,l)==0
            c{k,l}=[num2str(M(k,l))];
        elseif round(M(k,l),3)~=0
            if P(k,l) < 0.01
                c{k,l}=[num2str(round(M(k,l),4)) '^{***}'];
            elseif P(k,l) < 0.05
                c{k,l}=[num2str(round(M(k,l),4)) '^{**}'];
            elseif P(k,l) < 0.1
                c{k,l}=[num2str(round(M(k,l),4)) '^{*}'];
            else
                c{k,l}=[num2str(round(M(k,l),4))];
            end
        else
            if P(k,l) < 0.01
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