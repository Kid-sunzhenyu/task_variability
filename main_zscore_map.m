clear,clc;

tasks = {'LANGUAGE', 'MOTOR', 'WM'};

for n = 1:3
    Lhdr = gifti(['InterSubject_Variability_' tasks{n} '_12mr_gsr_L.func.gii']);
    Rhdr = gifti(['InterSubject_Variability_' tasks{n} '_12mr_gsr_R.func.gii']);
    Ldata = Lhdr.cdata;
    Rdata = Rhdr.cdata;
    Lind = find(Ldata>0);
    Rind = find(Rdata>0);
    zdata = zscore([Ldata(Lind); Rdata(Rind)]);
    Lout = 0*Ldata; Rout = 0*Rdata;
    Lout(Lind) = zdata(1:Lind);
    Rout(Rind) = zdata(Lind+1:end);
    Lhdr.cdata = Lout;
    Rhdr.cdata = Rout;
    save(Lhdr, ['Zscore_InterSubject_Variability_' tasks{n} '_12mr_gsr_L.func.gii'])
    save(Rhdr, ['Zscore_InterSubject_Variability_' tasks{n} '_12mr_gsr_R.func.gii'])
end
