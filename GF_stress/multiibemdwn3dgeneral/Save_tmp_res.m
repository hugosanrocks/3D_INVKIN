function Save_tmp_res(para,uw,j) %#ok<INUSD>

name = para.nametmp;
save([name,'tmp',num2str(j)],'para','uw');
