function [comarray,widtharray] = stats(indist)

a= size(indist);
comarray=zeros(1,a(2));
widtharray = zeros(1,a(2));
for i = 1: a(2);
    tot= a(1);
    comarray(i) = sum(indist(1:a(1),i))/tot;
    rel =   indist(1:a(1),i) - comarray(i);
    widtharray(i) = sqrt(sum((rel.*rel))/tot);
end
