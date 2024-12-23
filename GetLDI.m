%% Get LDI


function LDI = GetLDI(Nresp)
    Nnew = sum(Nresp(1, :));
    Nsim = sum(Nresp(2, :));
    LDI = Nresp(2,2)/Nsim - Nresp(1,2)/Nnew;
end