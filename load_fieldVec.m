function [fieldVec] = load_fieldVec(structName,fieldName,dim)


if isfield(structName,fieldName)
    fieldVec = arrayfun(@(x) x.(fieldName), structName, 'UniformOutput', false); fieldVec = [fieldVec{:}]';
    fieldVec = reshape(fieldVec,dim,length(structName))';
else
    display(['No fieldVec named ' fieldName])
    fieldVec = [];
end