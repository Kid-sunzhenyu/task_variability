function Func_write_func_gifti_32k( name_stem, alldata, OutPath, Lhdr, Rhdr )
%FUNC_WRITE_FUNC_GIFTI_32K Summary of this function goes here
%   Detailed explanation goes here

    load fsLR_32k_config.mat
    
    Lhdr.cdata = 0*Lhdr.cdata;
    Lhdr.cdata(Lvertlist) = alldata(1:29696);
    save(Lhdr, [OutPath '/' name_stem '_L.func.gii'])
    Rhdr.cdata = 0*Rhdr.cdata;
    Rhdr.cdata(Rvertlist) = alldata(29697:end)';
    save(Rhdr, [OutPath '/' name_stem '_R.func.gii'])

end

