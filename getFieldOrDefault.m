function value = getFieldOrDefault(struct, fieldName, defaultValue)
    if isfield(struct, fieldName)
        value = struct.(fieldName);
    else
        value = defaultValue;
    end
end