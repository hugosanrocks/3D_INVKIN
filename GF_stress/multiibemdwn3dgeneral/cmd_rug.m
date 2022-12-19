med     =get(bouton.med,'value');
if med>1
    icont =get(bouton.cont,'value');
else
    icont = 1;
end

para.cont(med,icont).ruggeo=get(bouton.ruggeo,'value');