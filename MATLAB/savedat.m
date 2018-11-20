function savedat(SR, filename)
%SAVEDAT Summary of this function goes here
%   Detailed explanation goes here

fid = fopen(filename, 'w');
headers = {'i_Model'; ...
           'x_Model'; ...
           'y_Model'; ...
           'i_Scene'; ...
           'x_Scene'; ...
           'y_Scene'};
fprintf(fid, [repmat('%s,',1,numel(headers)-1), '%s\n'], headers{:});
for i=1:size(SR.ix, 1)
    if SR.ix(i) > 0
        fprintf(fid, [repmat('%f,',1,numel(headers)-1), '%f\n'], ...
            i, SR.mdp(i,1), SR.mdp(i,2), SR.ix(i), SR.scp(SR.ix(i),1), SR.scp(SR.ix(i),2));
    else
        fprintf(fid, [repmat('%f,',1,numel(headers)-1), '%f\n'], ...
            i, SR.mdp(i,1), SR.mdp(i,2), SR.ix(i), -1, -1);
    end
end
fclose(fid);

end

