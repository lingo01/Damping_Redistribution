% generate the reduced state dynamic matrices for transient stability analysis
function mat = funMatrixReduce(state0, deviceNumFlag)

    global nNet netB;
    global refPos LoadEquivTime LoadEquivDroop;
    global type CDInv Gen;

    global LoadEquivTimeTmp;

    if nargin == 1
        deviceNumFlag = 0;
    end

    if deviceNumFlag == 0
        refPos_tmp = 39;
    elseif deviceNumFlag == 1
        refPos_tmp = 54;
    end

    setV1 = [refPos_tmp];
    setV2 = [];
    setVc = [];
    for ii = 1:nNet
        if strcmp(type{ii}, "Gen") && ii ~= refPos_tmp
            setV1 = [setV1, ii];
        elseif strcmp(type{ii}, "Cnct")
            setVc = [setVc, ii];
        elseif ii ~= refPos_tmp
            setV2 = [setV2, ii];
        end
    end

    setV2 = [setV2, setVc];
    
    numV1 = length(setV1);      numV2 = length(setV2);
    setV  = [setV1, setV2];     numV  = length(setV);

    ang0 = state0(1:2:end);     ang0 = ang0(setV);
    vol0 = state0(2:2:end);     vol0 = vol0(setV);

    %% Power flow matrices
    netBRe = netB(setV,setV); % - netB(setV,setVc) /netB(setVc, setVc) * netB(setVc,setV);

    matA = zeros(numV, numV);
    matC1 = zeros(numV, numV);
    matC2 = zeros(numV, numV);
    matD = zeros(numV, numV);
    matE = zeros(numV, numV);
    
    for ii = 1:numV
        for jj = 1:numV
            if ii ~= jj
                matA(ii,jj) = -netBRe(ii,jj) * vol0(ii) * vol0(jj) * cos(ang0(ii)-ang0(jj));
            end
        end
        matA(ii,ii) = -sum(matA(ii,:));
    end
    
    for ii = 1:numV
        for jj = 1:numV
            if ii ~= jj
                matC1(ii,jj) = -netBRe(ii,jj) * cos(ang0(ii)-ang0(jj));
                matC2(ii,jj) = -netBRe(ii,jj) * vol0(ii) * cos(ang0(ii)-ang0(jj));
            else
                matC1(ii,ii) = -netBRe(ii,ii);
                matC2(ii,ii) = -2 * netBRe(ii,ii) * vol0(ii);
            end
        end
    end
    for ii = 1:numV
        for jj = 1:numV
            if ii ~= jj
                matC2(ii,ii) = matC2(ii,ii) - netBRe(ii,jj) * vol0(jj) * cos(ang0(ii)-ang0(jj));
            end
        end 
    end
    
    for ii = 1:numV
        for jj = 1:numV
            if ii ~= jj
                matD(ii,jj) = netBRe(ii,jj) * vol0(ii) * sin(ang0(ii)-ang0(jj));
                matE(ii,jj) = - netBRe(ii,jj) * vol0(ii) * vol0(jj) * sin(ang0(ii)-ang0(jj));
            end
        end
    end
    for ii = 1:numV
        matD(ii,ii) = -sum(matD(:,ii));
        matE(ii,ii) = -sum(matE(ii,:));
    end
    
    matA = matA(2:end, 2:end);
    matC1 = matC1;
    matC2 = matC2;
    matD = matD(2:end, :);
    matE = matE(:, 2:end);

    %% Device matrices
    matM = zeros(nNet, nNet);
    matd = zeros(nNet, nNet);
    matz = zeros(nNet, nNet);
    mate = zeros(nNet, nNet);
    
    matt = zeros(nNet, nNet);
    matr = zeros(nNet, nNet);
    matg = zeros(nNet, nNet);
    
    for ii = 1:nNet
        switch type{ii}
            case 'Gen'
                matM(ii,ii) = Gen{ii}.M;
                matd(ii,ii) = Gen{ii}.D;
                matz(ii,ii) = 0;
                mate(ii,ii) = 0;
    
                matt(ii,ii) = Gen{ii}.Td;
                matr(ii,ii) = Gen{ii}.xd - Gen{ii}.xdp;
                matg(ii,ii) = 0;
            case 'CDInv'
                matM(ii,ii) = 0;
                matd(ii,ii) = 0;
                matz(ii,ii) = CDInv{ii}.t1/CDInv{ii}.D1;
                mate(ii,ii) = 1/CDInv{ii}.D1;
    
                matt(ii,ii) = CDInv{ii}.t2;
                matr(ii,ii) = 0;
                matg(ii,ii) = CDInv{ii}.D2;
            case 'Load'
                matM(ii,ii) = 0;
                matd(ii,ii) = 0;
                matz(ii,ii) = LoadEquivTime/LoadEquivDroop;
                mate(ii,ii) = 1/LoadEquivDroop;
    
                matt(ii,ii) = LoadEquivTimeTmp;
                matr(ii,ii) = 0;
                matg(ii,ii) = LoadEquivDroop;
            case 'Cnct'
                matM(ii,ii) = 0;
                matd(ii,ii) = 0;
                matz(ii,ii) = LoadEquivTime/LoadEquivDroop;
                mate(ii,ii) = 1/LoadEquivDroop;
    
                matt(ii,ii) = LoadEquivTime;
                matr(ii,ii) = 0;
                matg(ii,ii) = LoadEquivDroop;
        end
    end
    
    matt = matt(setV, setV);
    matr = matr(setV, setV);
    matg = matg(setV, setV);

    setV1 = setV1(2:end);       numV1 = numV1-1;
    setV  = [setV1, setV2];     numV  = length(setV);

    matM = matM(setV1, setV1);
    matd = matd(setV1, setV1);
    matz = matz(setV2, setV2);
    mate = mate(setV2, setV2);

    %% Criteria matrices
    setV1 = 1:numV1;
    setV2 = (numV1+1):numV;

    % matTa
    matTa = [
      zeros(numV1,numV1),                eye(numV1),              zeros(numV1,numV2);
      -inv(matM)*matA(setV1,setV1),     -inv(matM)*matd,        -inv(matM)*matA(setV1,setV2);
      -inv(matz)*matA(setV2,setV1),     zeros(numV2,numV1),     -inv(matz)*(matA(setV2,setV2)+mate);
    ];
    
    % matTv
    matTv = -inv(matt) * (matr*matC1 + matg*matC2 + eye(numV+1));
    
    % matTav
    matTav = [zeros(numV1,numV+1); -inv(matM)*matD(setV1,:); -inv(matz)*matD(setV2,:)];
    
    % matTva
    matTmp = -inv(matt) * (matr*matD' + matg*matE);
    matTva = [matTmp(:,setV1), zeros(numV+1,numV1), matTmp(:,setV2)];   
    
    % matT
    matT = [matTa, matTav; matTva, matTv];
    matT_uncouple = [matTa, zeros(size(matTav)); zeros(size(matTva)), matTv];
    matT_couple   = matT - matT_uncouple;


    %% OUTPUT
    mat.A = matA;
    mat.C1 = matC1;
    mat.C2 = matC2;
    mat.D = matD;
    mat.E = matE;

    mat.M = matM;
    mat.d = matd;
    mat.k = matz;
    mat.t = matt;
    mat.r = matr;
    mat.g = matg;

    mat.Ta = matTa;
    mat.Tv = matTv;
    mat.Tav = matTav;
    mat.Tva = matTva;


    mat.T = matT;
    mat.T_uncouple = matT_uncouple;
    mat.T_couple = matT_couple;
end