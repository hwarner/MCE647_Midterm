function m=unpackSymbolics(symVals,listOfSymVarNames)


for i=1:length(listOfSymVarNames)
    if isvector(symVals)  %If it is a vector
        m.(char(listOfSymVarNames(i)))=symVals(i);
    else  %It's a mtrix
        m.(char(listOfSymVarNames(i)))=symVals(:,i);
    end
end