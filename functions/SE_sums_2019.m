% standard errors of the sums
function [SE_sums]=SE_sums(SVAR)
    SE_sums={};
    c=0;
    for j=1:3
        for i=1:3
            if SVAR.B(i,j)==0
                SE_sums{1}(i,j)=NaN;
            else
                c=c+1;
                SE_sums{1}(i,j)=SVAR.SE(c);
            end
        end
    end % loop for the first matrix
    
    c_B=0;
    c_Q2=sum(sum(SVAR.B~=0));
    for j=1:3
        for i=1:3
    
            if SVAR.Q2(i,j)~=0
                c_Q2=c_Q2+1;
            end
            if SVAR.B(i,j)~=0
                c_B=c_B+1;
            end
            
            if SVAR.Q2(i,j)==0
                SE_sums{2}(i,j)=SE_sums{1}(i,j);
            elseif SVAR.B(i,j)==0
                SE_sums{2}(i,j)=SVAR.SE(c_Q2);
            else
                SE_sums{2}(i,j)=sqrt(sum(sum(SVAR.V([c_B,c_Q2],[c_B,c_Q2]))));
            end
        end
    end % loop for the second matrix
    
    c_B=0;
    c_Q2=sum(sum(SVAR.B~=0));
    c_Q3=sum(sum(SVAR.B~=0))+sum(sum(SVAR.Q2~=0));
    for j=1:3
        for i=1:3
    
            if SVAR.Q2(i,j)~=0
                c_Q2=c_Q2+1;
            end
            if SVAR.B(i,j)~=0
                c_B=c_B+1;
            end
            if SVAR.Q3(i,j)~=0
                c_Q3=c_Q3+1;
            end
    
            if SVAR.Q3(i,j)==0
                SE_sums{3}(i,j)=SE_sums{2}(i,j);
            elseif (SVAR.B(i,j)==0) && (SVAR.Q2(i,j)==0)
                SE_sums{3}(i,j)=SVAR.SE(c_Q3);
            elseif (SVAR.B(i,j)==0)
                SE_sums{3}(i,j)=sqrt(sum(sum(SVAR.V([c_Q2,c_Q3],[c_Q2,c_Q3]))));
            elseif (SVAR.Q2(i,j)==0)
                SE_sums{3}(i,j)=sqrt(sum(sum(SVAR.V([c_B,c_Q3],[c_B,c_Q3]))));
            else
                SE_sums{3}(i,j)=sqrt(sum(sum(SVAR.V([c_B,c_Q2,c_Q3],[c_B,c_Q2,c_Q3]))));
            end
        end
    end % loop for the third matrix
    
end