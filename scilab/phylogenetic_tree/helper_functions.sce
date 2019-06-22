// BM4320 Genomic Signal Processing
//
// Operations
// 1. Construction of a WPGMA Tre1
// 1. Construction of a UPGMA Tree
// 1. Construction of a NJ Tree
//
// Objectives
// 1. Introduction to phylogenetic trees
// 2. Concept of the molecular clock
//
// Upeka Premaratne, ENTC, University of Moratuwa (upeka@uom.lk) 2019/05/22
// Free to use, distribute and modify for educational purposes with attribution

function result = remove_eols(text_in)
    // Remove EOLs of fasta file
    keys = find(text_in==10);
    if (isempty(keys)) then
        result = text_in;
    else
    text_out = text_in(1:(keys(1)-1));
    k_n = length(keys)-1;
    for k=1:k_n
        text_i = text_in((keys(k)+1):(keys(k+1)-1));
        text_out = [text_out, text_i];
    end
    result = [text_out,text_in((keys(k_n+1)+1):length(text_in))];
    end
endfunction

function comp=get_comp(base)
    if (base==65) then
        comp = 84;
    elseif (base==67) then
        comp = 71;
    elseif (base==71) then
        comp = 67;
    elseif (base==84) then
        comp = 65;
    end
endfunction

function gen_code=get_fasta_at(file_name,g_pos,g_end,strand)
   //Estimate the necessary overread to compensate for the EOL charactors of FASTA files
   g_len = g_end-g_pos;
   if g_len>0 then
          n_extra = floor(g_len/70);
          n_offset = floor(g_pos/70);
          file_details = fileinfo(file_name);
          file_len = file_details(1);
          fd = mopen(file_name,'rb');
          mseek(0,fd);
          header = mgetl(fd,1);
          g_start = length(header);
          mseek(g_start+g_pos+n_offset,fd);
          raw_code = mget(g_len+n_extra,'c',fd);
          mclose(fd);
          code_i = remove_eols(raw_code);
          if strand==1 then
              gen_code = code_i;
          else
              //get complementary strand
              len = length(code_i);
              code_c = [];
              for k=1:len
                  code_c = [code_c,get_comp(code_i(k))];
                  gen_code = code_c;
              end
          end
   else
       gen_code = [];
   end

endfunction

function y=e_dist(x1,x2)
    // Calculate the Euclidean distance
    y=sqrt(sum((x1-x2).^2));
endfunction

function [min_val,min_x,min_y]=find_min(dist_mat)
    // Find the minimum non-zero value of a 
    m=size(dist_mat,1)
    c_min=dist_mat(1,2)
    c_x = 1;
    c_y = 2;
    for k_x=1:m
        for k_y=(k_x+1):m
            if dist_mat(k_x,k_y)<c_min then
                c_min=dist_mat(k_x,k_y);
                c_x = k_x;
                c_y = k_y;
            end
        end
    end
    min_x=c_x;
    min_y=c_y;
    min_val=c_min;
endfunction

//function [new_mat,new_labels,new_min_val,new_min_lab]= c_join_wpgma(dist_mat,j_labels)
//    // Form two clusters using the WPGMA algorithm
//    m=size(dist_mat,1);
//    [min_val,min_x,min_y]=find_min(dist_mat); // Find the minimum
//        
//    m_row_keys = 1:m;
//    m_row_keys = m_row_keys((m_row_keys~=min_y)&(m_row_keys~=min_x));
//    m_row_x = dist_mat(min_x,m_row_keys);
//    m_row_y = dist_mat(min_y,m_row_keys);
//    m_res = dist_mat(m_row_keys,m_row_keys); // Remove the minimum nodes from the distance matrix
//    
//    g_label = sprintf('(%s,%s)',j_labels(min_x),j_labels(min_y)); // Create the new labels
//    n_labels = j_labels(m_row_keys);
//    
//    disp(g_label);
//    disp(n_labels);
//    
//    new_labels = [g_label; n_labels];
//    
//    n_row = [0,1/2*(m_row_x+m_row_y)]; // Calculate new distance according to WPGMA (half of the sum of the distance of each node to the two minimum values)
//    n_res = [zeros(m-2,1),triu(m_res)];
//    
//    n_mat = [n_row; n_res];
//   
//    new_mat=n_mat+n_mat';
//    
//    new_min_val = min_val/2; // Height of the parent node
//    new_min_lab = g_label;
//endfunction

function [new_mat,new_labels,new_min_val,new_min_lab,new_heights]= c_join_wpgma(dist_mat,j_labels,heights)
    // Form two clusters using the WPGMA algorithm
    m=size(dist_mat,1);
    [min_val,min_x,min_y]=find_min(dist_mat); // Find the minimum
        
    m_row_keys = 1:m;
    m_row_keys = m_row_keys((m_row_keys~=min_y)&(m_row_keys~=min_x));
    m_row_x = dist_mat(min_x,m_row_keys);
    m_row_y = dist_mat(min_y,m_row_keys);
    m_res = dist_mat(m_row_keys,m_row_keys); // Remove the minimum nodes from the distance matrix
    dis = min_val/2;
    g_heights_x = dis-heights(min_x);
    g_heights_y = dis-heights(min_y);
    //g_label = sprintf('(%s,%s)',j_labels(min_x),j_labels(min_y)); // Create the new labels
    g_label = sprintf('(%s:%.4f,%s:%.4f)',j_labels(min_x),g_heights_x,j_labels(min_y),g_heights_y);
    n_labels = j_labels(m_row_keys);
    n_heights = heights(m_row_keys);
    
    new_labels = [g_label; n_labels];
    new_heights =[dis; n_heights];
    n_row = [0,1/2*(m_row_x+m_row_y)]; // Calculate new distance according to WPGMA (half of the sum of the distance of each node to the two minimum values)
    n_res = [zeros(m-2,1),triu(m_res)];
    
    n_mat = [n_row; n_res];
   
    new_mat=n_mat+n_mat';
    
    new_min_val = min_val/2; // Height of the parent node
    //new_min_val = dis;
    new_min_lab = g_label;
    //new_min_lab = sprintf('%s:%.4f',g_label,new_min_val);
endfunction

//function print_wpgma_tree(dist_mat,j_labels)
//    // Iteratively construct the WPGMA tree and print it
//    m=size(dist_mat,1);
//    i_lab=j_labels;
//    i_mat=dist_mat;
//    for iter=1:(m-1)
//        [i_mat,i_lab,i_min_val,i_min_lab]=c_join_wpgma(i_mat,i_lab);
//        disp(i_mat);
//        disp(i_min_lab);
//        disp(i_min_val);
//    end
//endfunction

function print_wpgma_tree(dist_mat,j_labels)
    // Iteratively construct the WPGMA tree and print it
    m=size(dist_mat,1);
    i_lab=j_labels;
    i_mat=dist_mat;
    heights = zeros(dist_mat(1,:))';
    for iter=1:(m-1)
        [i_mat,i_lab,i_min_val,i_min_lab,heights]=c_join_wpgma(i_mat,i_lab,heights);
        //disp(i_mat);
        disp(i_min_lab);
        disp(i_min_val);
    end
endfunction

function [new_mat,new_labels,new_cardinality,new_min_val,new_min_lab]= c_join_upgma(dist_mat,j_labels,j_cardinality)
    // Form two clusters using the UPGMA algorithm
    // Start of with cardinality matrix of all ones
    m=size(dist_mat,1);
    [min_val,min_x,min_y]=find_min(dist_mat); // Find the minimum
        
    m_row_keys = 1:m;
    m_row_keys = m_row_keys((m_row_keys~=min_y)&(m_row_keys~=min_x));
    m_row_x = dist_mat(min_x,m_row_keys);
    m_row_y = dist_mat(min_y,m_row_keys);
    m_res = dist_mat(m_row_keys,m_row_keys); // Remove the minimum nodes from the distance matrix
    
    // Get cardinality of the two clusters
    cardinality_x = j_cardinality(min_x);
    cardinality_y = j_cardinality(min_y);
    
    cardinality_c = cardinality_x+cardinality_y; // Cardinality of cluster
    n_cardinality = j_cardinality(m_row_keys); // Remove cardinality of minimum nodes
    
    // Update labels
    g_label = sprintf('(%s,%s)',j_labels(min_x),j_labels(min_y));
    n_labels = j_labels(m_row_keys);
    
    new_labels = [g_label; n_labels]; // New labels
    new_cardinality = [cardinality_c; n_cardinality]; // New cardinality
    
    n_row = [0,1/cardinality_c*(m_row_x*cardinality_x+m_row_y*cardinality_y)]; // New distance vector based upon UPGMA cardinality scaling
    n_res = [zeros(m-2,1),triu(m_res)];
    
    n_mat = [n_row; n_res];
   
    new_mat=n_mat+n_mat';
    
    new_min_val = min_val/2; // Height of the parent node
    new_min_lab = g_label;
endfunction

function print_upgma_tree(dist_mat,j_labels)
    // Iteratively construct the UPGMA tree and print it
    m=size(dist_mat,1);
    i_cardinality=ones(m,1);
    i_lab=j_labels;
    i_mat=dist_mat;
    for iter=1:(m-1)
        [i_mat,i_lab,i_cardinality,i_min_val,i_min_lab]=c_join_upgma(i_mat,i_lab,i_cardinality);
        disp(i_mat);
        disp(i_cardinality');
        disp(i_min_lab);
        disp(i_min_val);
    end
endfunction

function [q_matrix,d_mat]=get_q_matrix(dist_mat)
    // Compute the Q matrix and new distance matrix
    m=size(dist_mat,1);
    q_mat=zeros(m,m);
    for k_x=1:m
        for k_y=1:m
            if k_x~=k_y then
                sum_d=0;
                for k_z=1:m
                    sum_d=sum_d-dist_mat(k_x,k_z)-dist_mat(k_y,k_z);
                end
                q_mat(k_x,k_y)=(m-2)*dist_mat(k_x,k_y)+sum_d;
            end
        end
    end
    [min_val,min_x,min_y]=find_min(q_mat); // Find minimum from Q matrix
    d_min=dist_mat(min_x,min_y); // Minimum distance
    sum_min=0;
    for k_z=1:m
        sum_min=sum_min+dist_mat(min_x,k_z)-dist_mat(min_y,k_z);
    end
    d_min_x = 1/2*d_min+sum_min/(2*(m-2)); // Minimum distance from parent node to x
    d_min_y = d_min-d_min_x; // Minimum distance from parent node to y
    q_matrix=q_mat;
    d_mat = [d_min_x d_min_y];
endfunction

function [new_mat,new_labels,new_d,new_min_lab]= c_join_nj(dist_mat,j_labels)
    // Cluster according to the Neighbor Joining algorithm
    m=size(dist_mat,1);
    [q_mat,d_mat] = get_q_matrix(dist_mat);
    [min_val,min_x,min_y]=find_min(q_mat);
        
    m_row_keys = 1:m;
    m_row_keys = m_row_keys((m_row_keys~=min_y)&(m_row_keys~=min_x));
    m_row_x = dist_mat(min_x,m_row_keys);
    m_row_y = dist_mat(min_y,m_row_keys);
    m_res = dist_mat(m_row_keys,m_row_keys);
    
    g_label = sprintf('(%s,%s)',j_labels(min_x),j_labels(min_y));
    n_labels = j_labels(m_row_keys);
    
    new_labels = [g_label; n_labels];
    
    dist_xy=dist_mat(min_x,min_y);
    
    n_row = [0,1/2*(m_row_x+m_row_y-dist_xy)]; // Compute new distances according to NJ formula
    n_res = [zeros(m-2,1),triu(m_res)];
    
    n_mat = [n_row; n_res];
   
    new_mat=n_mat+n_mat';
    
    new_d = d_mat;
    new_min_lab = g_label;
endfunction

function print_nj_tree(dist_mat,j_labels)
    // Iteratively construct the NJ tree and print it
    m=size(dist_mat,1);
    i_lab=j_labels;
    i_mat=dist_mat;
    for iter=1:(m-2)
        [i_mat,i_lab,i_min_val,i_min_lab]=c_join_nj(i_mat,i_lab);
        disp(i_min_lab);
        disp(i_min_val);
    end
    disp('Final Node');
    disp(i_mat(1,2));
    disp(i_lab(2));
endfunction

//////////////////////
//A=[0 1 16 14 15; 1 0 16 15 16; 16 16 0 16 16; 14 15 16 0 1; 15 16 16 1 0];
////A=[0 10 20 9 14; 0 0 19 8 16; 0 0 0 30 20; 0 0 0 0 12; 0 0 0 0 0];
//t_labels = ['A';'B';'C';'D';'E'];
//
//B = A;

//// Wikipedia example for UPGMA and WPGMA
//A=[0 17 21 31 23; 0 0 30 34 21; 0 0 0 28 39; 0 0 0 0 43; 0 0 0 0 0];
//B=A+A';

//disp('WPGMA Tree');
//print_wpgma_tree(B,t_labels);

//disp('WPGMA Tree');
//print_wpgma_tree(B,t_labels);

//// Wikipedia example for NJ
//A=[0 5 9 9 8; 0 0 10 10 9; 0 0 0 8 7; 0 0 0 0 3; 0 0 0 0 0];
//B=A+A';

//disp('NJ Tree;);
//print_nj_tree(B,t_labels);

//////////////////////
// Class LacZ example
//     1: AY746950.1     100      72      66      69      66      66
//     2: AY746951.1      72     100      67      67      65      67
//     3: AY746948.1      66      67     100      81      68      68
//     4: AY746947.1      69      67      81     100      68      69
//     5: AY746952.1      66      65      68      68     100      77
//     6: AY746949.1      66      67      68      69      77     100
//
//t_labels=['AY746950.1'; 'AY746951.1'; 'AY746948.1'; 'AY746947.1'; 'AY746952.1'; 'AY746949.1'];
////
//A=[100      72      66      69      66      66;
//72     100      67      67      65      67;
//66      67     100      81      68      68;
//69      67      81     100      68      69;
//66      65      68      68     100      77;
//66      67      68      69      77     100];

//t_lablels = ['EU872499.1(A)', 'KJ804868.1(B)', 'EU872498.1(C)','AF333829.1(D)','AF333821.1(E)', 'AF333827.1(F)', 'AF333822.1(G)','AY185543.1(H)'];
//
//A = [  0.000  0.000  0.003  0.003  0.007  0.007  0.007  0.011 0.010  0.016  0.008  0.007  0.015  0.011  0.013  0.013 
// 0.011  0.008  0.015  0.013  0.015  0.015  0.015  0.021 0.020  0.011  0.016  0.018  0.018  0.016  0.013  0.018 0.018  0.015  0.018  0.021  0.016  0.015  0.024  0.018 0.024  0.018  0.023  0.024  0.016  0.149  0.131  0;
//      0.000  0.000  0.003  0.003  0.006  0.006  0.006  0.011 0.010  0.016  0.008  0.006  0.015  0.011  0.013  0.013 0.011  0.008  0.014  0.013  0.014  0.014  0.014  0.021 0.019  0.011  0.016  0.018  0.018  0.016  0.013  0.018 0.018  0.014  0.018  0.021  0.016  0.014  0.024  0.018 0.024  0.018  0.023  0.024  0.016  0.149  0.130  0.197 0.197  0.195  0.195  0.197  0.203  0.204  0.525;
//       0.003  0.003  0.000  0.007  0.010  0.010  0.010  0.015 0.013  0.020  0.011  0.010  0.018  0.015  0.016  0.016 
//           0.015  0.011  0.018  0.016  0.018  0.018  0.018  0.024 
//           0.023  0.015  0.020  0.021  0.021  0.020  0.016  0.021 
//           0.021  0.018  0.021  0.024  0.020  0.018  0.028  0.021 
//           0.028  0.021  0.026  0.028  0.020  0.152  0.135  0.206 
//           0.206  0.204  0.204  0.206  0.207  0.206  0.519;
//        0.003  0.003  0.007  0.000  0.006  0.010  0.010  0.011 
//           0.010  0.013  0.008  0.006  0.014  0.011  0.013  0.013 
//           0.008  0.008  0.014  0.013  0.014  0.014  0.014  0.021 
//           0.019  0.011  0.016  0.018  0.018  0.016  0.013  0.018 
//           0.018  0.014  0.018  0.021  0.016  0.014  0.024  0.018 
//           0.024  0.018  0.022  0.024  0.016  0.147  0.131  0.199 
//           0.199  0.197  0.197  0.199  0.201  0.200  0.522;
//        0.007  0.006  0.010  0.006  0.000  0.013  0.010  0.014 
//           0.011  0.016  0.008  0.010  0.018  0.014  0.013  0.013 
//           0.011  0.011  0.014  0.016  0.018  0.018  0.018  0.024 
//           0.022  0.014  0.016  0.018  0.018  0.016  0.013  0.018 
//           0.018  0.018  0.018  0.021  0.016  0.014  0.024  0.018 
//           0.022  0.018  0.022  0.024  0.016  0.150  0.134  0.201 
//           0.201  0.199  0.199  0.201  0.203  0.202  0.523;
//        0.007  0.006  0.010  0.010  0.013  0.000  0.006  0.011 
//           0.013  0.016  0.008  0.006  0.014  0.011  0.011  0.013 
//           0.011  0.008  0.014  0.013  0.014  0.014  0.014  0.021 
//           0.019  0.011  0.016  0.018  0.018  0.016  0.013  0.018 
//           0.018  0.014  0.018  0.021  0.014  0.014  0.027  0.021 
//           0.024  0.018  0.019  0.021  0.013  0.148  0.128  0.203 
//           0.203  0.201  0.201  0.203  0.205  0.204  0.518;
//         0.007  0.006  0.010  0.010  0.010  0.006  0.000  0.014 
//           0.011  0.016  0.011  0.010  0.014  0.011  0.010  0.010 
//           0.008  0.008  0.011  0.013  0.018  0.018  0.018  0.024 
//           0.022  0.014  0.016  0.021  0.021  0.019  0.016  0.021 
//           0.021  0.018  0.021  0.024  0.019  0.018  0.030  0.024 
//           0.026  0.021  0.022  0.024  0.016  0.150  0.132  0.202 
//           0.202  0.200  0.200  0.202  0.205  0.204  0.520;
//         0.011  0.011  0.015  0.011  0.014  0.011  0.014  0.000 
//           0.011  0.016  0.010  0.008  0.016  0.013  0.014  0.014 
//           0.013  0.010  0.016  0.014  0.016  0.016  0.016  0.023 
//           0.021  0.013  0.018  0.019  0.019  0.018  0.014  0.019 
//           0.019  0.016  0.019  0.023  0.018  0.016  0.029  0.023 
//           0.026  0.019  0.023  0.026  0.018  0.147  0.133  0.206 
//           0.206  0.204  0.204  0.204  0.208  0.205  0.519            
//    ]
//    
//B = A;

//B = 100*ones(6,6)-A; // Convert similarity matrix into a distance matrix
//disp(B);

//disp('WPGMA Tree');
//print_wpgma_tree(B,t_labels);
//
//disp('UPGMA Tree');
//print_upgma_tree(B,t_labels);
//
//disp('NJ Tree');
//print_nj_tree(B,t_labels);
