function res = Place_FCs(n,Bidi_Ports)
    %% Post Processing
    res = [];
    for i=1:length(n)
        newmatrix = n{i};
        fc = create_fc(newmatrix,Bidi_Ports);
        res = [res fc];
    end
end

function fc = create_fc(x, Bidi_Ports)
    fc =cell(0);  
    base=cell(0);
    mat=cell(0);    %here all the available configurations of the specific power flow graph will be stored
    %bidirectional FC
    bidi_ind = find(x(5,:) ~= 0);
    uni_ind = find(x(5,:) == 0);
    pool = 2:(Bidi_Ports+1);
    original_list = 1:3;
    invpool = original_list(~ismember(original_list, pool));
    ind_set = [];
    I_exists=0;
    for i=1:size(x,2) % Check if a Type I SG exists in the PFGAM
        if x(4,i)==1
            I_exists=1;
        end
    end
    if I_exists
        ind_set=2;
        for i=1:length(bidi_ind)
            if x(4,bidi_ind(i))~=1 && Bidi_Ports==2
                for j=2:3
                    if x(j,i)==-1 %μαλλον θελει bidi_ind(i)13/05
                        ind_set=j;
                    end
                end
            end
        end
    else
        for i=1:length(pool)
            for j=1:length(invpool)
                ind_set = [ind_set; pool(i) invpool(j)]; %assign ports i, j as a possible pair for FC placement
            end
        end
        for i=1:length(pool)
            for j=1:length(pool)
                if i~=j
                    ind_set = [ind_set; pool(i) pool(j)];
                end
            end
        end
    end
    for ii=1:size(ind_set,1)%for each ind_set combination
        y = zeros(size(x,2)-length(bidi_ind)/2, size(x,2)+3);
        temp = x;
        for it=1:length(ind_set(1,:)) %for each element of current ind_set combination
            for ind=1:length(bidi_ind) %for each bidirectional SG
                y(it,bidi_ind(ind)) = x(ind_set(ii,it),bidi_ind(ind));
                y(it,size(x,2)+ind_set(ii,it)) = - x(ind_set(ii,it),bidi_ind(ind));
                temp(ind_set(ii,it),bidi_ind(ind))=0;
            end  
            %place a large value (10) where a bidirectional power flow to a
            %port is existent(+- 1 cannot be utilized)
            if (sum(y(it,1:size(x,2)))==0) && (sumabs(y(it,1:size(x,2)))~=0)
                y(it,size(x,2)+ind_set(ii,it)) = 10;
            end
        end
        base{ii}=temp;
        mat{ii}=y;
    end  
    %
    if isempty(uni_ind)
        uni_port=[];
        for i=1:size(x,2)    
            if x(4,i)~=1 && I_exists
                uni_ind=i;
                for k=1:3
                    if sum(x(k,:))~=0
                        uni_port=[uni_port k];
                    end
                end

            end
        end
        if I_exists
            for j=1:length(mat)
                mat{j}=[mat{j};zeros(1,(size(x,2)+3))];
            end
        end
    else
        uni_port=[1 2 3];
    end
    if isempty(mat)==1
        mat{1} = zeros(size(x,2)-length(bidi_ind)/2, size(x,2)+3);
        base{1}=x;
    end
    for j=1:length(uni_ind) %for each unidirectional subgraph (column of PFGAM)
        initial_mat_len=length(mat);     t=initial_mat_len+1;
        for i=1:length(uni_port) %for each port 
            if x(uni_port(i),uni_ind(j))~=0 && (x(4,uni_ind(j))~=1 || sumabs(x(1:uni_port(i),uni_ind(j)))==1)  %that is connected on the subgraph (non zero cell value)
                for k=1:initial_mat_len %work on the array already filled with the previous subgraphs (columns)
                    temp=mat{k};
                    if uni_port(i)==2 && length(uni_ind)==1 && I_exists && Bidi_Ports==1%not for TPS-EB
                        continue;
                    end
                    %if temp(length(bidi_ind)/2+j,size(x,2)+i)==0
                    if sumabs(temp(length(bidi_ind)/2+j,size(x,2)+uni_port(i)))==0
                        temp(length(bidi_ind)/2+j,uni_ind(j))=x(uni_port(i),uni_ind(j));      
                        temp(length(bidi_ind)/2+j,size(x,2)+uni_port(i))=-x(uni_port(i),uni_ind(j));
                        %tttt=base{t-initial_mat_len};
                        tttt=base{k};
                        tttt(uni_port(i),uni_ind(j))=0;
                        base{t}=tttt;
                        mat{t} = temp; t=t+1;
                    end
                end
            end
        end
        mat = mat(initial_mat_len+1:length(mat)); %remove matrices used to generate new matrices
        base = base(initial_mat_len+1:length(base));
    end
%     if length(base)==0
%         base{1}=x;
%     end
    for i=1:length(mat)
        temp=mat{1,i};
        temp2=base{i};
        mat{i}=[temp2(1:3,:) zeros(3,3);temp];
    end    
    fc=[fc mat];
end