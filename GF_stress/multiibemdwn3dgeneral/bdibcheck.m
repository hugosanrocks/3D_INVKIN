if isfield(para,'b_dib')%exist('b_dib','var')
    if ~isstruct(para.b_dib)
        para.b_dib=0;
    end
else
    para.b_dib=0;
end