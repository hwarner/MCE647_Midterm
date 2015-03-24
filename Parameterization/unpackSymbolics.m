function m=unpackSymbolics(symVals,listOfSymVarNames)

%unpackSymbolics - Match a set of values returned from a symbolic function
%to list of symbolic variables and place them in a structure.  
%
%
%m=unpackSymbolics(symVals,listOfSymVarNames)
%
%       Inputs:
%               symVals - values returned from a symbolic equation.
%               listOfSymVarNames - cell array of symbolic variables.
%                   Typically this comes from: varList=symvar(equation);
%       Outputs:
%               m - structure containing an field per variable using it's
%                   name and value.  


for i=1:length(listOfSymVarNames)
    if isvector(symVals)  %If it is a vector
        m.(char(listOfSymVarNames(i)))=symVals(i);
    else  %It's a mtrix
        m.(char(listOfSymVarNames(i)))=symVals(:,i);
    end
end