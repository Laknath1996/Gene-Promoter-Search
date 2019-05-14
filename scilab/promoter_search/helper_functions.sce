// BM4321 Genomic Signal Processing
//
// Operations
// 1. Reading of a FASTA file
// 2. Extraction of coding and non-coding DNA
// 3. Reading of a GenBank Protein Table
// 4. Basic analysis of a coding and non-coding region
//
// Objectives
// 1. Familiarization with the coding regions of different organisms (Archaea, Bacteria, Eukaryota etc.)
// 2. Preliminary investigation of gene promoter regions
//
// Upeka Premaratne, ENTC, University of Moratuwa (upeka@uom.lk) 2016/12/05
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

function inc=base_inc(base,c_val)
    if (base==65) then
        inc=c_val+[1 0 0 0]';
    elseif (base==67) then
        inc=c_val+[0 1 0 0]';
    elseif (base==71) then
        inc=c_val+[0 0 1 0]';
    elseif (base==84) then
        inc=c_val+[0 0 0 1]';
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

function [validity] = verify(seq, que, thresh_len_down, y_len)
    // check whether the whole query is in the sequence
    idx = find(ascii(que)==ascii('W'));
    if length(idx) < y_len then
        validity = 0;
    elseif length(idx(idx > length(seq)-thresh_len_down)) ~= 0 then
        validity = 0;
    elseif length(seq(idx)(seq(idx)==ascii('C') | seq(idx)==ascii('G'))) ~= 0 then
        validity = 0;
    else
        validity = 1;
    end
endfunction

function dna_seq = get_dna_seq(n_key, fasta_file, gp, gn, thresh_up, thresh_down, st_genes)
    if n_key > st_genes then 
//        if gn(n_key-st_genes,2)-thresh_up >= 1 then      
        dna_seq = get_fasta_at(fasta_file,gn(n_key-st_genes,2)-thresh_down+1,gn(n_key-st_genes,2)+thresh_up+3,0);
        dna_seq = dna_seq(1:thresh_up+thresh_down);
        dna_seq = flipdim(dna_seq,2);
//        end
    else
//        if gp(n_key,1)-thresh_up >= 1 then 
        dna_seq = get_fasta_at(fasta_file,gp(n_key,1)-thresh_up,gp(n_key,1)+thresh_down+3,1);
        dna_seq = dna_seq(1:thresh_up+thresh_down);
//        end
    end
endfunction

function N = get_consecutive_Ws(ay, d_len, thresh_down)
    idx = find(ascii(ay) == ascii('W'))(1);
    N = 1;
    while (1)
        idx = idx + 1;
        if dna_seq(idx) == ascii('W') & idx <= d_len - thresh_down then //disregard Ws beyond d_len-thresh_dow
            N = N + 1;
        else 
            break;
        end
    end
endfunction

function result=compare_oligo(oligo_a,oligo_b)
    if (length(oligo_a)~=length(oligo_b)) then
        result = -1;
    else
        s = sum(abs(oligo_a-oligo_b))
        if (s==0) then
            result = 1;
        else
            result = 0;
        end
    end
endfunction

function [gene_array_p,gene_array_n,noncoding_array_p,noncoding_array_n]=get_protein_pos_array(filename)
    // Get the coding and non-coding DNA positions from a protein table
    fd = mopen(filename,'r');
    data = mgetl(fd,1);
    ga_p = [];
    ga_n = [];
    
    nca_p = [];
    nca_n = [];
    
    nc_prev_p = 0;
    nc_prev_n = 0;
    
    while (~meof(fd))
        if (data == []) then
            break;
        end  
        data = mgetl(fd,1);
        //disp(type(data));
        keys = strindex(data,ascii(9));
        p_data = strsplit(data,keys);
        pg_start = strtod(p_data(3));
        pg_stop = strtod(p_data(4));
        //disp(strcmp(p_data(5),'-'));
        //disp(strcmp(p_data(5),'+'));
        if (~isempty(p_data(5))) then
            if (strcmp(p_data(5),'-')==1) then
                ga_n = [ga_n; pg_start,pg_stop];
                nca_n = [nca_n; nc_prev_n, (pg_start-1)];
                nc_prev_n = pg_stop+1;
            else             
                ga_p = [ga_p; pg_start,pg_stop];
                nca_p = [nca_p; nc_prev_p, (pg_start-1)];
                nc_prev_p = pg_stop+1;
            end
        end
    end
    mclose(fd);
    gene_array_p=ga_p;
    noncoding_array_p=nca_p;
    gene_array_n=ga_n;
    noncoding_array_n=nca_n;
endfunction

function save_fasta(filename,header_line,gen_code)
    // Save results into a FASTA file
    fd = mopen(filename,'wc');
    mputl(header_line,fd);
    gen_len = length(gen_code);
    for key=1:gen_len
        mput(gen_code(key),'c');
        if pmodulo(key,70)==0 then
            mputl('',fd);
        end
    end
    mclose(fd)
endfunction

// BM4321 Genomic Signal Processing
//
// Standalone function file for alignment
//
// Operations
// 1. Global Search
// 2. Local Search
// 3. Semi-Global Search
//
// Objectives
// 1. Familiarization with basic sequence alignment algorithms
//
// Upeka Premaratne, ENTC, University of Moratuwa (upeka@uom.lk) 2016/12/05
// Free to use, distribute and modify for educational purposes with attribution

function g=gap_penalty(k)
    // Gap penalty function
    g_alpha = 0;
    g_beta = 2;
    g = -(g_alpha+g_beta*k);
endfunction

function result=prom_base(base_x,base_y,score_m,score_mm)
   // Match mismatch score
   // score_m - match score
   // score_mm - mismatch score
   if (base_x==base_y) then
       result = score_m;
   else
       result = score_mm;
   end
endfunction

function result=comp_base(base_x,base_y,score_m,score_mm)
    if (base_y == base_x) then
        result = score_m;
    else
        result = score_mm;
    end
endfunction

function matrix_result = scorematrix_global(seq_x,seq_y,score_m,score_mm,gap_penalty)
    // Function to perform global search
    len_x = length(seq_x);
    len_y = length(seq_y);
    basic_mat = zeros(len_x+1,len_y+1);
    
    //Initialize with gap penalties
    for k_x=1:(len_x+1)
        basic_mat(k_x,1)=gap_penalty(k_x-1);
    end
    
    for k_y=1:(len_y+1)
        basic_mat(1,k_y)=gap_penalty(k_y-1);
    end
    
    //Recurrance relation
    for k_x=2:(len_x+1)
        for k_y=2:(len_y+1)
            score_match = basic_mat(k_x-1,k_y-1)+comp_base(seq_x(k_x-1),seq_y(k_y-1),score_m,score_mm);
            score_gap_x = basic_mat(k_x,k_y-1)+gap_penalty(1);
            score_gap_y = basic_mat(k_x-1,k_y)+gap_penalty(1);
            basic_mat(k_x,k_y)=max([score_match score_gap_x score_gap_y]);
        end
    end
    matrix_result = basic_mat;
endfunction

function [align_x,align_y] = traceback_global(seq_x,seq_y,score_m,score_mm,gap_penalty)
    // Traceback for global alignment
    score_matrix = scorematrix_global(seq_x,seq_y,score_m,score_mm,gap_penalty);
    
    len_x = length(seq_x);
    len_y = length(seq_y);
    
    base_x = seq_x;
    base_y = seq_y;
    
    out_x = [];
    out_y = [];
    
    k_x = len_x;
    k_y = len_y;
    
    
    // Perform the traceback
    while (k_x>0&k_y>0) then
        k_gx = score_matrix(k_x,k_y+1);
        k_gy = score_matrix(k_x+1,k_y);
        k_diag = score_matrix(k_x,k_y);
        k_c = score_matrix(k_x+1,k_y+1);
        
        if k_c==k_diag+comp_base(seq_x(k_x),seq_y(k_y),score_m,score_mm) then
            out_x = [out_x,base_x(k_x)];
            out_y = [out_y,base_y(k_y)];
            k_x=k_x-1;
            k_y=k_y-1;
        elseif k_c==k_gx+gap_penalty(1) then
            out_x = [out_x,base_x(k_x)];
            out_y = [out_y,45];
            k_x=k_x-1;
        elseif k_c==k_gy+gap_penalty(1) then
            out_x = [out_x,45];
            out_y = [out_y,base_y(k_y)];
            k_y=k_y-1;
        end
    end
    
    // Write the output
    if (k_x>0&k_y==0) then
        pre_x = ascii(base_x(1:k_x));
        pre_y = ascii(45*ones(1,k_x));
        align_x = strcat([pre_x,strrev(ascii(out_x))]);
        align_y = strcat([pre_y,strrev(ascii(out_y))]);
    elseif (k_x==0&k_y>0) then
        pre_y = ascii(base_y(1:k_y));
        pre_x = ascii(45*ones(1,k_y));
        align_x = strcat([pre_x,strrev(ascii(out_x))]);
        align_y = strcat([pre_y,strrev(ascii(out_y))]);       
    else
        align_x = strrev(ascii(out_x));
        align_y = strrev(ascii(out_y));          
    end
endfunction

function matrix_result = scorematrix_local(seq_x,seq_y,score_m,score_mm,gap_penalty)
    //Function to perform local search
    len_x = length(seq_x);
    len_y = length(seq_y);
    basic_mat = zeros(len_x+1,len_y+1);
    
    //No need to initialize because gap penalty will always be negative
    //Recurrance relation
    for k_x=2:(len_x+1)
        for k_y=2:(len_y+1)
            score_match = basic_mat(k_x-1,k_y-1)+comp_base(seq_x(k_x-1),seq_y(k_y-1),score_m,score_mm);
            score_gap_x = basic_mat(k_x,k_y-1)+gap_penalty(1);
            score_gap_y = basic_mat(k_x-1,k_y)+gap_penalty(1);
            basic_mat(k_x,k_y)=max([score_match score_gap_x score_gap_y 0]); // Add the zero to global search
        end
    end
    matrix_result = basic_mat;
endfunction

function [align_x,align_y] = traceback_local(seq_x,seq_y,score_m,score_mm,gap_penalty)
    // Traceback for local alignment
    score_matrix = scorematrix_local(seq_x,seq_y,score_m,score_mm,gap_penalty);
    //disp(score_matrix');
    
    len_x = length(seq_x);
    len_y = length(seq_y);
    
    base_x = seq_x;
    base_y = seq_y;
    
    out_x = [];
    out_y = [];
    
    [m,n]=max(score_matrix);
    k_x = n(1,1)-1;
    k_y = n(1,2)-1;
    
    post_x = ascii(seq_x((k_x+1):len_x));
    post_y = ascii(seq_y((k_y+1):len_y));
    
    // Perform the traceback
    while (k_x>0&k_y>0) then
        k_gx = score_matrix(k_x,k_y+1);
        k_gy = score_matrix(k_x+1,k_y);
        k_diag = score_matrix(k_x,k_y);
        k_c = score_matrix(k_x+1,k_y+1);
        
        if k_c==k_diag+comp_base(seq_x(k_x),seq_y(k_y),score_m,score_mm) then
            out_x = [out_x,base_x(k_x)];
            out_y = [out_y,base_y(k_y)];
            k_x=k_x-1;
            k_y=k_y-1;
        elseif k_c==k_gx+gap_penalty(1) then
            out_x = [out_x,base_x(k_x)];
            out_y = [out_y,45];
            k_x=k_x-1;
        elseif k_c==k_gy+gap_penalty(1) then
            out_x = [out_x,45];
            out_y = [out_y,base_y(k_y)];
            k_y=k_y-1;
        elseif k_c==0;
            break;
        end
    end
    
    // Write the output
    if (k_x>0&k_y==0) then
        pre_x = ascii(base_x(1:k_x));
        pre_y = ascii(45*ones(1,k_x));
        align_x = strcat([pre_x,strrev(ascii(out_x)),post_x]);
        align_y = strcat([pre_y,strrev(ascii(out_y)),post_y]);
    elseif (k_x==0&k_y>0) then
        pre_y = ascii(base_y(1:k_y));
        pre_x = ascii(45*ones(1,k_y));
        align_x = strcat([pre_x,strrev(ascii(out_x)),post_x]);
        align_y = strcat([pre_y,strrev(ascii(out_y)),post_y]);       
    else
        align_x = strrev(ascii(out_x));
        align_y = strrev(ascii(out_y));          
    end     
endfunction

function matrix_result = scorematrix_prom(seq_x,seq_y,score_m,score_mm,gap_penalty)
    //Function to perform local search
    len_x = length(seq_x);
    len_y = length(seq_y);
    basic_mat = zeros(len_x+1,len_y+1);
    
    //No need to initialize because gap penalty will always be negative
    //Recurrance relation
    for k_x=2:(len_x+1)
        for k_y=2:(len_y+1)
            score_match = basic_mat(k_x-1,k_y-1)+prom_base(seq_x(k_x-1),seq_y(k_y-1),score_m,score_mm);
            score_gap_x = basic_mat(k_x,k_y-1)+gap_penalty(1);
            score_gap_y = basic_mat(k_x-1,k_y)+gap_penalty(1);
            basic_mat(k_x,k_y)=max([score_match score_gap_x score_gap_y 0]); // Add the zero to global search
        end
    end
    matrix_result = basic_mat;
endfunction

function [align_x,align_y] = traceback_prom(seq_x,seq_y,score_m,score_mm,gap_penalty)
    // Traceback for local alignment
    score_matrix = scorematrix_prom(seq_x,seq_y,score_m,score_mm,gap_penalty);
    disp(score_matrix');
    
    len_x = length(seq_x);
    len_y = length(seq_y);
    
    base_x = seq_x;
    base_y = seq_y;
    
    out_x = [];
    out_y = [];
    
    [m,n]=max(score_matrix);
    k_x = n(1,1)-1;
    k_y = n(1,2)-1;
    
    post_x = ascii(seq_x((k_x+1):len_x));
    post_y = ascii(seq_y((k_y+1):len_y));
    
    // Perform the traceback
    while (k_x>0&k_y>0) then
        k_gx = score_matrix(k_x,k_y+1);
        k_gy = score_matrix(k_x+1,k_y);
        k_diag = score_matrix(k_x,k_y);
        k_c = score_matrix(k_x+1,k_y+1);
        
        if k_c==k_diag+prom_base(seq_x(k_x),seq_y(k_y),score_m,score_mm) then
            out_x = [out_x,base_x(k_x)];
            out_y = [out_y,base_y(k_y)];
            k_x=k_x-1;
            k_y=k_y-1;
        elseif k_c==k_gx+gap_penalty(1) then
            out_x = [out_x,base_x(k_x)];
            out_y = [out_y,45];
            k_x=k_x-1;
        elseif k_c==k_gy+gap_penalty(1) then
            out_x = [out_x,45];
            out_y = [out_y,base_y(k_y)];
            k_y=k_y-1;
        elseif k_c==0;
            break;
        end
    end
    
    // Write the output
    if (k_x>0&k_y==0) then
        pre_x = ascii(base_x(1:k_x));
        pre_y = ascii(45*ones(1,k_x));
        align_x = strcat([pre_x,strrev(ascii(out_x)),post_x]);
        align_y = strcat([pre_y,strrev(ascii(out_y)),post_y]);
    elseif (k_x==0&k_y>0) then
        pre_y = ascii(base_y(1:k_y));
        pre_x = ascii(45*ones(1,k_y));
        align_x = strcat([pre_x,strrev(ascii(out_x)),post_x]);
        align_y = strcat([pre_y,strrev(ascii(out_y)),post_y]);       
    else
        align_x = strrev(ascii(out_x));
        align_y = strrev(ascii(out_y));          
    end     
endfunction

function matrix_result = scorematrix_semiglobal(seq_x,seq_y,score_m,score_mm,gap_penalty)
    //Function to perform local search
    len_x = length(seq_x);
    len_y = length(seq_y);
    basic_mat = zeros(len_x+1,len_y+1);
    
    //No need to initialize because gap penalty will always be negative
    //Recurrance relation
    for k_x=2:(len_x+1)
        for k_y=2:(len_y+1)
            score_match = basic_mat(k_x-1,k_y-1)+comp_base(seq_x(k_x-1),seq_y(k_y-1),score_m,score_mm);
            score_gap_x = basic_mat(k_x,k_y-1)+gap_penalty(1);
            score_gap_y = basic_mat(k_x-1,k_y)+gap_penalty(1);
            basic_mat(k_x,k_y)=max([score_match score_gap_x score_gap_y]); // Same as global search
        end
    end
    matrix_result = basic_mat;
endfunction

function [align_x,align_y] = traceback_semiglobal(seq_x,seq_y,score_m,score_mm,gap_penalty)
    // Traceback for semi-global alignment
    
    // Assuming x is longer than y to make code simpler
    
    score_matrix = scorematrix_semiglobal(seq_x,seq_y,score_m,score_mm,gap_penalty);
        
    len_x = length(seq_x);
    len_y = length(seq_y);
    
    base_x = seq_x;
    base_y = seq_y;
    
    out_x = [];
    out_y = [];
    
    k_x = len_x;
    [n_max,k_y] = max(score_matrix(len_x+1,:)); 
    k_y=k_y-1;
    
    post_x = ascii(seq_x((k_x+1):len_x));
    post_y = ascii(seq_y((k_y+1):len_y));
    
    // Perform the traceback
    while (k_x>0&k_y>0) then
        k_gx = score_matrix(k_x,k_y+1);
        k_gy = score_matrix(k_x+1,k_y);
        k_diag = score_matrix(k_x,k_y);
        k_c = score_matrix(k_x+1,k_y+1);
        
        if k_c==k_diag+comp_base(seq_x(k_x),seq_y(k_y),score_m,score_mm) then
            out_x = [out_x,base_x(k_x)];
            out_y = [out_y,base_y(k_y)];
            k_x=k_x-1;
            k_y=k_y-1;
        elseif k_c==k_gx+gap_penalty(1) then
            out_x = [out_x,base_x(k_x)];
            out_y = [out_y,45];
            k_x=k_x-1;
        elseif k_c==k_gy+gap_penalty(1) then
            out_x = [out_x,45];
            out_y = [out_y,base_y(k_y)];
            k_y=k_y-1;
        end
    end
    
    // Write the output
    if (k_x>0&k_y==0) then
        pre_x = ascii(base_x(1:k_x));
        pre_y = ascii(45*ones(1,k_x));
        align_x = strcat([pre_x,strrev(ascii(out_x)),post_x]);
        align_y = strcat([pre_y,strrev(ascii(out_y)),post_y]);
    elseif (k_x==0&k_y>0) then
        pre_y = ascii(base_y(1:k_y));
        pre_x = ascii(45*ones(1,k_y));
        align_x = strcat([pre_x,strrev(ascii(out_x)),post_x]);
        align_y = strcat([pre_y,strrev(ascii(out_y)),post_y]);       
    else
        align_x = strrev(ascii(out_x));
        align_y = strrev(ascii(out_y));          
    end     
endfunction

// BM4321 Genomic Signal Processing
//
// Operations
// 1. Position probability matrix
// 2. Statitical sequence alignment
// 3. Calculation of ppm entropy
//
// Objectives
// 1. Introduction to statistical sequence alignment
//
// Upeka Premaratne, ENTC, University of Moratuwa (upeka@uom.lk) 2019/04/03
// Free to use, distribute and modify for educational purposes with attribution

function bk=get_amino_key(amino)
// Keys of all amino acids and stop codons
    if (amino==ascii('A')) then
        bk=1; // L-Alanine (Ala/A)
    elseif (amino==ascii('C')) then
        bk=2; // L-Cysteine (Cys/C)
    elseif (amino==ascii('D')) then
        bk=3; // L-Aspartic acid (Asp/D)
    elseif (amino==ascii('E')) then
        bk=4; // L-Glutamic acid (Glu/E)
    elseif (amino==ascii('F')) then
        bk=5; // L-Phenylalanine (Phe/F)
    elseif (amino==ascii('G')) then
        bk=6; // Glycine (Gly/G)
    elseif (amino==ascii('H')) then
        bk=7; // L-Histidine (His/H)
    elseif (amino==ascii('I')) then
        bk=8; // L-Isoleucine (Ile/I)
    elseif (amino==ascii('K')) then
        bk=9; // L-Lysine (Lys/K)
    elseif (amino==ascii('L')) then
        bk=10; // L-Leucine (Leu/L)
    elseif (amino==ascii('M')) then
        bk=11; // L-Methionine (Met/M)
    elseif (amino==ascii('N')) then
        bk=12; // L-Asparagine (Asn/ N)
    elseif (amino==ascii('P')) then
        bk=13; // L-Proline (Pro/P)
    elseif (amino==ascii('Q')) then
        bk=14; // L-Glutamine (Gln/Q)
    elseif (amino==ascii('R')) then
        bk=15; // L-Arginine (Arg/R)
    elseif (amino==ascii('S')) then
        bk=16; // L-Serine (Ser/S)
    elseif (amino==ascii('T')) then
        bk=17; // L-Threonine (Thr/T)
    elseif (amino==ascii('V')) then
        bk=18; // L-Valine (Val/V)
    elseif (amino==ascii('W')) then
        bk=19; // L-Tryptophan (Trp/W)'
    elseif (amino==ascii('Y')) then
        bk=20; // L-Tyrosine (Tyr/Y)
    elseif (amino==ascii('+')|amino==ascii('U')|amino==ascii('O')) then
        bk=21; //Stop codons
    end
endfunction

function bk=get_base_key(base)
    // Base key as A=1, C=2, G=3, T=4 (alphabetic)
    if (base==65) then
        bk=1;
    elseif (base==67) then
        bk=2;
    elseif (base==71) then
        bk=3;
    elseif (base==84) then
        bk=4;
    // For multiple possibilities assign one for consistency
    // Evenly distributed as much as possible A=3,C=3, G=3, T=2
    elseif (base==ascii('R')) then
        bk=1; // Assign A (A or G)    
    elseif (base==ascii('Y')) then
        bk=2; // Assign C (C or T) 
    elseif (base==ascii('S')) then
        bk=3; // Assign G (C or G) 
    elseif (base==ascii('W')) then
        bk=4; // Assign T (A or T) 
    elseif (base==ascii('K')) then
        bk=3; // Assign G (G or T) 
    elseif (base==ascii('M')) then
        bk=1; // Assign A (A or C) 
    elseif (base==ascii('B')) then
        bk=2; // Assign C (C or G or T) 
    elseif (base==ascii('D')) then
        bk=3; // Assign G (A or G or T) 
    elseif (base==ascii('H')) then
        bk=4; // Assign T (A or C or T) 
    elseif (base==ascii('V')) then
        bk=1; // Assign A (A or C or G)
    elseif (base==ascii('N')) then
        bk=2; // Assign C (any base)
    end
endfunction

function bk=get_ct_key(base)
    // Base key according to the codon table T=1, C=2, A=3, G=4
    if (base==65) then
        bk=3;
    elseif (base==67) then
        bk=2;
    elseif (base==71) then
        bk=4;
    elseif (base==84) then
        bk=1;
    elseif (base==ascii('Y')|base==ascii('N')) then
        bk=2;
    elseif (base==ascii('R')|base==ascii('W')) then
        bk=3;
    // For multiple possibilities assign one for consistency
    // Evenly distributed as much as possible A=3,C=3, G=3, T=2
    elseif (base==ascii('R')) then
        bk=3; // Assign A (A or G)    
    elseif (base==ascii('Y')) then
        bk=2; // Assign C (C or T) 
    elseif (base==ascii('S')) then
        bk=4; // Assign G (C or G) 
    elseif (base==ascii('W')) then
        bk=1; // Assign T (A or T) 
    elseif (base==ascii('K')) then
        bk=4; // Assign G (G or T) 
    elseif (base==ascii('M')) then
        bk=3; // Assign A (A or C) 
    elseif (base==ascii('B')) then
        bk=2; // Assign C (C or G or T) 
    elseif (base==ascii('D')) then
        bk=3; // Assign G (A or G or T) 
    elseif (base==ascii('H')) then
        bk=1; // Assign T (A or C or T) 
    elseif (base==ascii('V')) then
        bk=3; // Assign A (A or C or G)
    elseif (base==ascii('N')) then
        bk=2; // Assign C (any base)
    end
endfunction

function base_hist=get_base_hist(seq)
    hist=ones(1,4)/4;
    l_seq = length(seq);
    for pos=1:l_seq
        key=get_base_key(seq(pos));
        hist(key)=hist(key)+1;
    end
    base_hist=hist;
endfunction

function ppm=get_ppm(motif)
    // Get the position probability matrix
    [m,n]=size(motif);
    hist_mat = [];
    for pos=1:n
       base_hist = get_base_hist(motif(1:m,pos));
       hist_mat = [hist_mat,base_hist'];
    end
    ppm=hist_mat/(m+1);
endfunction

function [w,su]=ppm_info(ppm,p0)
    [m,n]=size(ppm);
    p_mat = [];
    for k=1:m
        p_row = ones(1,n)/p0(k);
        p_mat = [p_mat; p_row];
    end
    w = log(ppm.*p_mat)/log(2);
    su = ppm.*w;
endfunction

function [s_mat,sm_pos,sm_seq]=get_motif_score(seq,ppm)
    l_seq = length(seq);
    [m,n]=size(ppm);
    score=[];
    max_score = 0;
    max_pos = 0;
    max_seq = [];
    for pos=1:(l_seq-n-1)
        sub_seq = seq(pos:(pos+n-1));
        ss_score=1;
        for pos_ss=1:n
            bk=get_base_key(sub_seq(pos_ss));
            ppm_s = ppm(bk,pos_ss);
            ss_score=ss_score*ppm_s;
        end
        score=[score,ss_score];
        if (ss_score>max_score) then
            max_score=ss_score;
            max_pos= pos;
            max_seq=sub_seq;
        end
    end
    s_mat = score;
    sm_pos = max_pos;
    sm_seq = max_seq;
endfunction


//    if n_key > size(gp,1) then 
//        if gn(n_key-size(gp,1),1)-thresh_up >= 1 then      
//            dna_seq = get_fasta_at(fasta_file,gn(n_key-size(gp,1),1)-thresh_up,gn(n_key-size(gp,1),1)+thresh_down,0);
//        end
//    else
//        if gp(n_key,1)-thresh_up >= 1 then 
//            dna_seq = get_fasta_at(fasta_file,gp(n_key,1)-thresh_up,gp(n_key,1)+thresh_down,1);
//        end
//    end

