function x = CalcularPunto(vector,yc)
data =size(vector)
for i = 1:data(1)
    if  vector(i,2) >= yc
        x = interpolar(yc,vector(i,1),vector(i-1,1),vector(i,2),vector(i-1,2));
        break;
    end
end