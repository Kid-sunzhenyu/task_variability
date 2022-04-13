clear,clc;

inpath = '/zfs/musc/david/HCP4variability/data/';
subs = textread('list_100Unrelated.txt','%s');

commonsubs = zeros(100,1);
for s = 1:100
    sub = subs{s}
    files = dir([inpath '/' sub '/']);
    if length(files) == 8
        commonsubs(s) = 1;
    end
end
sum(commonsubs)
saveas('list_common_subs.txt', commonsubs)