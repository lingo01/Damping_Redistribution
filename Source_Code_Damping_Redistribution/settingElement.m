% 元件设置
function settingElement(filename)
    global type; 
    global Info Gen CDInv QDInv Load Cnct;
    global nNet;
    
    type = {};
    
    Info = readtable(['./data&figure/', filename, '.xlsx']);
    
    for rr = 1:nNet
        if isempty(find(rr == Info.Node))
            type(rr) = cellstr("Cnct");
            Cnct{rr} = 0;  
        else
            ii = find(rr == Info.Node);
            type(rr) = Info.Type(ii);
            if strcmp(Info.Type(ii), 'Gen')
                Gen{rr}.M = Info.Para1(ii);    Gen{rr}.D = Info.Para2(ii);   Gen{rr}.Td = Info.Para3(ii);  
                Gen{rr}.xd = Info.Para4(ii);   Gen{rr}.xdp = Info.Para5(ii);
                Gen{rr}.Pm = Info.Para6(ii);   Gen{rr}.Ef = Info.Para7(ii);
            elseif strcmp(Info.Type(ii), 'CDInv')
                CDInv{rr}.t1 = Info.Para1(ii);      CDInv{rr}.t2 = Info.Para2(ii);  
                CDInv{rr}.As = Info.Para3(ii);      CDInv{rr}.Vs = Info.Para4(ii); 
                CDInv{rr}.Ps = Info.Para5(ii);      CDInv{rr}.Qs = Info.Para6(ii); 
                CDInv{rr}.D1 = Info.Para7(ii);      CDInv{rr}.D2 = Info.Para8(ii);
            elseif strcmp(Info.Type(ii), 'QDInv')
                QDInv{rr}.t1 = Info.Para1(ii);      QDInv{rr}.t2 = Info.Para2(ii);  
                QDInv{rr}.As = Info.Para3(ii);      QDInv{rr}.Vs = Info.Para4(ii); 
                QDInv{rr}.Ps = Info.Para5(ii);      QDInv{rr}.Qs = Info.Para6(ii); 
                QDInv{rr}.D1 = Info.Para7(ii);      QDInv{rr}.D2 = Info.Para8(ii);
            elseif strcmp(Info.Type(ii), 'Load')
                Load{rr}.Ps = Info.Para1(ii);      Load{rr}.Qs = Info.Para2(ii);
            else
                error('Initial error: settingElement error');
            end
        end
    end
    
end

