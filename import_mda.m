function data = import_mda(filename)

fi = fopen(filename, 'r');
if fi == 0
    disp('Couldn''t open file');
    return;
end

fstart = fread(fi, 3, 'int32');

dtypes = ["complex float32",
    "byte",
    "float32",
    "int16",
    "int32",
    "uint16",
    "double",
    "uint32"];

dtype = dtypes(-fstart(1));
num_dims = fstart(3);
dim_sz = fread(fi,num_dims, 'int32');
if (sum(dim_sz == 0) > 1)
    disp "Too many unspecified dimensions"
    return;
end

overwrite_file = false;
if (prod(dim_sz) == 0)
    data = fread(fi, Inf, dtype);

    ip = 'k';
    while ip ~= 'y' && ip ~= 'n'

        ip = input('One dimension is 0. Would you like to overwrite this dimension to be the correct size (y/n)?', 's') ;
    end
    if ip == 'y'
        overwrite_file = true;
        overwrite_dim = find(dim_sz == 0);
    end
else
    data = fread(fi, prod(dim_sz), dtype);
end

%disp "dim_sz"
%disp(dim_sz);
%disp "fstart"
%disp(fstart);
%disp "data(1:10)"
%disp(data(1:10));
%disp "numel(data)"
%disp(numel(data));


function v = f(x)
    if x == 0
        v = [];
    else
        v = x;
    end
end

if length(dim_sz) > 1
    sz_arg = mat2cell(dim_sz, ones(1,numel(dim_sz)), 1);
    sz_arg = cellfun(@f, sz_arg, 'UniformOutput', false);
    data = reshape(data, sz_arg{:});
end

fclose(fi);

if overwrite_file
    fi = fopen(filename, 'w');
    if fi == 0
        disp('Couldn''t open file for overwriting');
        return;
    end
    dim_sz(overwrite_dim) = size(data, overwrite_dim);
    fwrite(fi, fstart, 'int32');
    fwrite(fi, dim_sz, 'int32');
    fwrite(fi, data, dtype);
    fclose(fi);
end

end