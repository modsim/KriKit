nTest = 1e3;


tNChooseK = 0;
tVChooseK = 0;


for iTest=1:nTest
    t1 = tic;
    C1 = nchoosek(1:4,2);
    tNChooseK = toc(t1);
    
    t2 = tic;
    C2 = VChooseK(1:4,2);
    tVChooseK = toc(t2);
end
fprintf('t1 = %g\n',tNChooseK)
fprintf('t2 = %g\n',tVChooseK)
fprintf('t1/t2 = %g\n',tNChooseK/tVChooseK)