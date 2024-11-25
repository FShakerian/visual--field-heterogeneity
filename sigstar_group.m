function [g,p]=sigstar_group(Xvar,PP)

 g=[];p=[];
z=0;
    for i=1:length(Xvar)
        for j= i+1:length(Xvar)
            if PP(i,j)<=(0.05/1); z = z+1;
                g{z} = [Xvar(i),Xvar(j)]; p(z) = PP(i,j);
            end;
        end;
    end;
