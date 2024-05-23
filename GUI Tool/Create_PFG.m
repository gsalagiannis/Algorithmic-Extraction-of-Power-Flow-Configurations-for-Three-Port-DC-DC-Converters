%Create__PFG();

function res1 = Create__PFG(Bidi_Ports)
    global matrices enc_matrices N Nf SG h f n
    %% Solver Parameters
    N=1; 
    n=[];
    f=cell(0);
    %% Initialization
    tic
    numlist = [1]; %list of possible values of each cell

    [min_subgraph_num, max_subgraph_num] = find_num(Bidi_Ports);
    sg_list=min_subgraph_num:max_subgraph_num;

    %% Solution Set Creation and Rule Application
    for h=1:numel(sg_list)
        SG=sg_list(h); %current subgraph number
        matrices = cell(0); %Initialize an empty cell array
        current = zeros (N+2, SG); %Initialize an N+2 x SG zero matrix
        matrices{1} = current;
        m_i=2; %index for last position in matrices
        j=1;
        while (j<=SG) %for each column of the matrix
            ultra_len=length(matrices);
            for i=1:N+2 %for each row of the matrix        
                numlist = rule_1a(i, Bidi_Ports); %update the list of potential cell values based on port type
        %% -------------------------PFGAM Generation-------------------------------
                m_i=length(matrices)+1; %index for last position in matrices
                for l = 1:m_i-1 %do the below for all existent matrices. More matrices are created         
                    for k = 1:length(numlist)%in each cell put all the possible values
                        newmatrix = matrices{l};
                        newmatrix(i,j) = numlist(k);
                        matrices{m_i} = newmatrix;
                        m_i=m_i+1;                                     
                    end   
                end              
            end % column j is fully created here
        %--------------------------------------------------------------------------    
        % -------------------Optimization-Remove used matrices--------------------
            matrices = matrices(ultra_len+1:length(matrices)); %remove matrices used to generate new matrices
        % -------------------------------------------------------------------------
            if 1
                p=0;  enc_matrices = cell(0); e_p=1;
                while p<length(matrices)
                    clear_=0;     p=p+1;        newmatrix = matrices{p};

        %----------------------------Rule 2----------------------------------------
                    clear_ = rule_2(newmatrix,j);  p = erase(clear_, p); if clear_; continue, end       
                    if j>1
                        en_matrix = newmatrix(:,1:j) + 1;
                        enc_mm = encode_col(en_matrix);%+1 to get numbers 0,1,2
        %----------------------------Rule 3----------------------------------------                  
                        clear_ = rule_3(enc_mm); p = erase(clear_, p); if clear_; continue, end
        %----------------------------Rule 4----------------------------------------
                        clear_ = rule_4a(newmatrix,j,min_subgraph_num,Bidi_Ports);  p = erase(clear_, p); if clear_; continue, end
        %----------------------------Rule 4----------------------------------------                
                        clear_ = rule_4b(newmatrix,j);  p = erase(clear_, p); if clear_; continue, end  
        %----------------------------Rule 5----------------------------------------               
                        clear_ = rule_5(enc_mm,j); p = erase(clear_, p); if clear_; continue, end
                        e_p = add(clear_, e_p, enc_mm);                
                    end 
                 end
            end  
            j=j+1;
        end
        p=0;
        while p<length(matrices)
            clear_=0;
            p=p+1;newmatrix = matrices{p};
    %----------------------------Part B of Rule 1------------------------------
    %applied after creation because the whole row is needed
            clear_ = rule_1b(newmatrix, Bidi_Ports); p = erase(clear_, p); if clear_; continue, end
        %--------------------------------Rule 0------------------------------------
            clear_ = rule_0(newmatrix); p = erase(clear_, p); if clear_; continue, end
            type = find_type(newmatrix);
            pf = pf_find([newmatrix;type], Bidi_Ports);
            matrices{p}=[newmatrix;type;pf];
        end
        n=[n, matrices];
    end
    res1=n;
end

%--------------------------------------------------------------------------
function enc_m = encode_col(x)
    [rownum, colnum]=size(x);
    enc_m = zeros(1, colnum); %encoding by column
    for row = 1:rownum
        for column = 1:colnum
            enc_m(column) = enc_m(column)+x(row, column)*3^(row-1);
        end
    end
end
%--------------------------------------------------------------------------
function er = erase(rx, p)
    er=p;
    global matrices        
    if rx 
        matrices(p)=[]; 
        if p>=1
            er=p-1; %remove current matrix
        else
            er=1;
        end
    end
end
%--------------------------------------------------------------------------
function ad = add(rx, p, x)
    global enc_matrices
    ad = p;
    if rx==0
        enc_matrices{p} = x;
        ad=p+1;
    end
end
%--------------------------------------------------------------------------
function r0 = rule_0(x)
    global N
    r0 = 0;
    for i=1:N+2
        [s,~] = sumabs(x(i,:));
        if s==0
            r0 = 1; break;
        end
    end
end
%--------------------------------------------------------------------------
function r1a = rule_1a(row, Bidi_Ports)
    global N
    if row<=N
        r1a = [1];        
    elseif row==N+1
        r1a = [-1, 1];       
    elseif row==N+2 && Bidi_Ports==2
        r1a = [-1, 1];
    elseif row==N+2
        r1a = [-1];
    end
end
%--------------------------------------------------------------------------
function r1b = rule_1b(x, Bidi_Ports)
    global N
    r1b = 0;
    [s,~] = sumabs(x(N+1,:)); %Calculate the sum of cell absolute values of N+1 row.
    if s<=1 || s<=abs(sum(x(N+1,:))) %maybe first part of the condition is not required
       r1b = 1;
    end
    if Bidi_Ports==2
        [s,~] = sumabs(x(N+2,:)); %Calculate the sum of cell absolute values of N+1 row.
        if s<=1 || s<=abs(sum(x(N+2,:))) %maybe first part of the condition is not required
           r1b = 1;
        end 
    end
end
%--------------------------------------------------------------------------
function r2 = rule_2(x,col)
    r2 = 0;
    s_=sum(x(:,col));%sum of column 
    [sabs_,~] = sumabs(x(:,col)); %sum of abs elements of column 
    if ((abs(s_)==sabs_) )% application of Rule 2
        r2=1;
    end
end
%--------------------------------------------------------------------------
function r3 = rule_3(x)
    r3 = 0;   enc_m=[];
    for i=1:length(x)%Nf
        enc_m=[enc_m ; x(i)];
    end
    for h=1:height(enc_m)
        temp=enc_m(h,:);
        K=ismember(enc_m,temp,'rows'); %check if temp (encoding of one Xi-Xo pair) is existent in another Xi-Xo pair (has the same encoding)
        if nnz(K)>=2
            r3=1; break;
        end
    end
end
%--------------------------------------------------------------------------
function r4a = rule_4a(x,col,min_subgraph_num, Bidi_Ports)
    global h
    list = min_subgraph_num:-1:(min_subgraph_num-Bidi_Ports); %list of non-Type I subgraphs max required based on number of subgraphs
    r4a = 0;
    T=0; %number of non-Type I subgraphs
    for i=1:col %find T
        [s,~] = sumabs(x(:,i));
        if s>2
            T=T+1;
        end
    end
    if T>list(h)
        r4a = 1;
    end
end  
%--------------------------------------------------------------------------
function r4b = rule_4b(x, col)    
    r4b=0;
    tt=col-1; %clear identical power flow graphs
    for i=1:tt %for each column
        temp = repmat(x(:,i), 1, col);
        Q = (x(:,1:col) == temp); %all elements of x upt to column col (rest has not been updated yet)
        Q=Q.*x(:,1:col);%filter out zero element equality
        T=1;
        for o=1:col
            if i==o
                continue;
            end
            s=sumabs(Q(:,o));
            if s>1 && (sumabs(x(:,i))==2 || sumabs(x(:,o))==2)
                T=T+1;
            end
        end
        if T>1
            r4b=1;
        end
    end
end
%--------------------------------------------------------------------------
function r5 = rule_5(x, col)
    global enc_matrices
    r5 = 0;
    % Get the number of columns in x
    num_cols = size(x, 2);
    % Initialize an empty array for xi
    xi = []; 
    % Iterate through odd columns and append them to xi, similar for xo
    for col = 1:num_cols
        xi = [xi x(:, col)];
    end
    % Compute all permutations of the elements in the specified row
    perm_xi = perms(xi); 
    for i=1:length(enc_matrices)
        temp=enc_matrices{i};
        num_cols = size(temp, 2);
        tempi = []; tempo = [];
        for col = 1:num_cols
            tempi = [tempi temp(:, col)];
        end
        if ismember(tempi,perm_xi,'rows')
            r5=1; break;
        end    
    end   
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [minn, maxx] = find_num(Bidi_Ports)
    minn = 2;
%     if Bidi_Ports==2
%         maxx = minn+Bidi_Ports;
%     else 
%         maxx = 2;
%     end
    maxx = minn+Bidi_Ports;
end
%--------------------------------------------------------------------------
function type = find_type(x)
    type = zeros(1, size(x,2));
    for i=1:length(type)
        type(i)=sum(x(:,i));
    end
% Replace the numbers
    type(type == -1) = 3;%Type III
    type(type == 1) = 2; %Type II
    type(type == 0) = 1; %Type I
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function pf = pf_find(x, Bidi_Ports)
    ind = nchoosek(1:size(x,2),2); %create all possible combinations of powerflow subgraph pairs
    pf=zeros(1,size(x,2));
    for i=1:height(ind)
        if Bidi_Ports~=0
            temp= [];
            for it=1:Bidi_Ports
                temp = [temp; x(it+1,ind(i,1)) x(it+1,ind(i,2))];
            end
            for it=1:2-Bidi_Ports
                if sum(x(4,:))>size(x,2)
                    temp = [temp; -temp(1, 1) -temp(1, 2)];
                else
                    temp = [temp; 0 0];
                end    
            end
            prim_d_sum = trace(temp);
            sec_d_sum = trace(fliplr(temp));
            if prim_d_sum == -sec_d_sum && abs(prim_d_sum)>1
                pf(ind(i,1))=1; pf(ind(i,2))=1;
            end  
        end
    end
end
%--------------------------------------------------------------------------