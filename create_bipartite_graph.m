%The raw_file contains all insects species and plants species

load raw_file.mat

C1= raw_file(:,1);
C2= raw_file(:,2);

[UniC1,ia1,iC1] = unique(C1);

[UniC2,ia2,iC2] = unique(C2);

A = zeros(length(UniC1),length(UniC2));

for i=1:length(iC1)
    H(iC1(i),iC2(i)) = 1
end
