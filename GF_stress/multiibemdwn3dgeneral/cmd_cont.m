med     =get(bouton.med,'value');
set(bouton.xa     ,'string',para.cont(med,1).xa);
set(bouton.za     ,'string',para.cont(med,1).za);
set(bouton.a      ,'string',para.cont(med,1).a);
set(bouton.th     ,'string',para.cont(med,1).th);
if med>1
    icont   = get(bouton.cont,'value');
else
    icont   = 1;
end
set(bouton.contgeo,'value' ,para.cont(med,icont).geom);
set(bouton.base   ,'string',para.cont(med,icont).ba);
set(bouton.haut   ,'string',para.cont(med,icont).h);
set(bouton.ruggeo ,'value' ,para.cont(med,icont).ruggeo);
set(bouton.rugbase,'string',para.cont(med,icont).rba);
set(bouton.rughaut,'string',para.cont(med,icont).rh);

