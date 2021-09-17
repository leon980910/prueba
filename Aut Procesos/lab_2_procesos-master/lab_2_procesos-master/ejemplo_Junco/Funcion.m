function res = Funcion(v1, v2, valor)
    for k=1:length(v1)
        Aux1=v1(k);
        if Aux1>valor
            res=v2(k-1)+((v2(k)-v2(k-1))/(v1(k)-v1(k-1)))*(valor-v1(k-1));
            break
        end
    end
end