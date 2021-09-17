function x = interpolar(y,x1,x2,y1,y2)
x = x1 + ((x2-x1)/(y2-y1))*(y-y1)
end