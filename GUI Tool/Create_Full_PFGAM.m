function [res1, res2] = Create_Full_PFGAM(f)
    Colors=["blue" "red" "magenta"];
    res1 = cell(0); res2 = cell(0);
    for i=1:length(f)
            %Creating Adjacency Matrix
            %newmatrix = f{i};
            newmatrix = f{i};
            [row, column] = size(newmatrix);
            temp_matrix1 = zeros(3);
            temp_matrix1 = [temp_matrix1 newmatrix(1:3,1:column-3) zeros(3,row-3)];
            [row_indices, col_indices] = find(newmatrix(4:end, column-2:column).' == 10); 
            for p=1:length(row_indices)
                temp_matrix1(row_indices(p),column+col_indices(p)) = -10;
            end
            %temp_matrix2 = [-newmatrix(1:3,1:column-3).' zeros(column-3) newmatrix(4:end,1:column-3).'];
            temp_matrix2 = [zeros(column-3,3) zeros(column-3) zeros(column-3,row-3)];

            %temp_matrix3 = [newmatrix(4:end,:) zeros(row-3)];
            temp_matrix3 = [newmatrix(4:end,column-2:column) newmatrix(4:end,1:column-3) zeros(row-3)];
            A = [temp_matrix1;temp_matrix2;temp_matrix3];
            [row_indices, col_indices] = find(A == -1);
            for p=1:length(row_indices)
                A(row_indices(p),col_indices(p))=0;
                A(col_indices(p),row_indices(p))=1;
            end
            res1{i} = digraph(A);
            res2{i} = A;
            %nLabels = {'1','2','3','A', 'B'};%Ν=1,Νf=2
            P = {'1','2','3'};
            PF = {'A', 'B', 'C', 'D', 'E'};
            FC = {'a', 'b', 'c', 'd','e'};
            nLabels = [P PF(1:column-3) FC(1:row-3)];
            %plot(G,'Layout','force', 'NodeLabel',nLabels,'NodeFontSize',7,'ArrowSize',5)           
    end   
end

