if isfield(para.sortie,'UP')
    para.sortie.UPh = para.sortie.UP;
    para.sortie.USh = para.sortie.US;
    para.sortie.UIh = para.sortie.UI;
    para.sortie.UPt = 0;
    para.sortie.USt = 0;
    
    para.sortie     = rmfield(para.sortie,'UP');
    para.sortie     = rmfield(para.sortie,'US');
    para.sortie     = rmfield(para.sortie,'UI');
    fn0             = fieldnames(para.sortie);
end

if ~isfield(para.sortie,'szz')
    para.sortie.sxx =0;
    para.sortie.syy =0;
    para.sortie.szz =0;
    para.sortie.sxy =0;
    para.sortie.sxz =0;
    para.sortie.syz =0;
end