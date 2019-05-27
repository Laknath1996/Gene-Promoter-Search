/////////////////////////////Execution starts here/////////////////////////////////

// instructions

// 1. sepcify the names of the protein table and the fasta file
// 2. make local_search = 1 if you need to perform local search to compute the ppm
// 3. make intact = 1 if you need to use intact query 
// 4. make intact = 0 if you need to use non-intact query
// 5. make local_search = 0 if you need to combine the ppms.

clear;
clc;

// import functions
exec('helper_functions.sce');

// define useful parameters
protein_table = 'ProteinTable10668_341867.txt';
fasta_file = 'Acetobacter tropicalis.fasta';
local_search = 0;
intact = 0;
y = ascii('WWWW'); // query
thresh_up = 50; //upstream bases
thresh_down = 3; // downstream bases

[gp,gn,ncp,ncn]= get_protein_pos_array(protein_table); //get the coding and non-coding arrays 

if intact ==1 then
    m = 1; // match score (intact : 1, non-intact : 3)
    mm = -1; // mismatch score (intact : -1, non-intact : -3)
else
    m = 3; // match score (intact : 1, non-intact : 3)
    mm = -3; // mismatch score (intact : -1, non-intact : -3)
end

num_st_genes = size(gp,1); // total number of sense strand genes
num_ast_genes = size(gn, 1); // total number of anti-sense strand genes
num_genes = num_st_genes + num_ast_genes; // total number of genes

disp('Number of Total Genes = ', num_genes)
disp('Number of Sense Strand Genes = ', num_st_genes)
disp('Number of Anti-Sense Strand Genes = ', num_ast_genes)

if local_search == 1 then
    consec_Ws = []; // record the number consective Ws
    // perform local search and get the number of consective Ws
    for n_key=1:num_genes
        dna_seq = get_dna_seq(n_key, fasta_file, gp, gn, thresh_up, thresh_down, num_st_genes); //get the dna sequence
        lg = dna_seq==ascii('A') | dna_seq==ascii('C') | dna_seq==ascii('G') | dna_seq==ascii('T'); // verify the sequence
        if length(dna_seq(~lg)) ~= 0 then
            continue;
        end
        d_len = length(dna_seq); // length of the sequence
        // disp(ascii(dna_seq));
        dna_seq(dna_seq==ascii('A') | dna_seq==ascii('T')) = ascii('W'); // W = A = T
        [ax,ay] = traceback_local(dna_seq,y,1,-1,gap_penalty); // local search
        // disp(ascii(dna_seq));
        // disp(ay);
        validity = verify(dna_seq, ay, thresh_down, length(y)); // check whether the query is valid
        if validity == 1 then
            consec_Ws = [consec_Ws, get_consecutive_Ws(ay, d_len, thresh_down)]; // record the consectives number of Ws

        end
    end
    
    N_max = max(consec_Ws); // maximum number of consecutive Ws
    
    acc_motif = []; // record the motifs
    key_s = []; // record the upstream position
    
    // get the motifs, the promoter percentage and upstream position
    for n_key=1:num_genes
        dna_seq = get_dna_seq(n_key, fasta_file, gp, gn, thresh_up, thresh_down, num_st_genes); //get the dna sequence
        lg = dna_seq==ascii('A') | dna_seq==ascii('C') | dna_seq==ascii('G') | dna_seq==ascii('T'); // verify the sequence
        if length(dna_seq(~lg)) ~= 0 then
            continue;
        end
        temp = dna_seq;
        disp(ascii(dna_seq));
        dna_seq(dna_seq==ascii('A') | dna_seq==ascii('T')) = ascii('W'); // W = A = T
        d_len = length(dna_seq);
        [ax,ay] = traceback_local(dna_seq,y,m,mm,gap_penalty);
        disp(ay);
        validity = verify(dna_seq, ay, thresh_down, length(y)); // check whether the query is valid
        if validity == 1 then
            pos = find(ascii(ay) == ascii('W'))(1);
            key_s = [key_s, thresh_up - pos + 1];
            if (pos + N_max - 1 <= d_len - thresh_down) then
                seq = temp(pos:pos+N_max-1);
                acc_motif = [acc_motif ; seq]; // record the motifs
            end
        end     
    end
    
    ppm = get_ppm(acc_motif); // position probability matrix
    if intact==1 then
        savematfile('ppm_intact','-v7.3', 'ppm');
    else
        savematfile('ppm_non_intact','-v7.3', 'ppm');
    end
    
else
    
    loadmatfile('ppm_intact.mat');
    temp = zeros(size(ppm, 1), size(ppm, 2));
    temp = temp + ppm;
    loadmatfile('ppm_non_intact.mat');
    ppm = (temp + ppm)/2;
    
end
    
// conduct the statistical alignment
stat_align_pos = []; // record the stat. alighnment upstream position
[w,su]=ppm_info(ppm,0.25);

//disp("Entropy of the PPM : ");
//disp(w);
//disp("Information content of the PPM : ");
//disp(su);

NI = sum(su,'r');
TI = sum(NI);

disp("Net information Content of each Position : ");
disp(NI);
disp("Total information Content : ");
disp(TI);

ppm = ppm(:,1:4); // only the first four bases have higher entropy values
N_max = size(ppm, 2);

imp = ['CCCC','GGGG']; // improbable promoter sequences
score_thresh = max(get_score(ascii(imp(1)),ppm),get_score(ascii(imp(2)),ppm)); // get the threshold for the statistical alignment

for n_key=1:num_genes
    dna_seq = get_dna_seq(n_key, fasta_file, gp, gn, thresh_up, thresh_down, num_st_genes); //get the dna sequence
    lg = dna_seq==ascii('A') | dna_seq==ascii('C') | dna_seq==ascii('G') | dna_seq==ascii('T'); // verify the sequence
    if length(dna_seq(~lg)) ~= 0 then
        continue;
    end
    d_len = length(dna_seq);
    [sm_score,sm_pos,sm_seq]=get_motif_score(dna_seq,ppm); // stat. alignment
    if (sm_pos + N_max - 1 <= d_len - thresh_down) then
        if sm_score > score_thresh  then
            // disp(ascii(sm_seq));
            stat_align_pos = [stat_align_pos, thresh_up - sm_pos + 1]; // record the stat. alighnment upstream positio
        end 
    end
end

// display results
if local_search==1 then
    disp(sprintf("Total number of genes : %d", num_genes));
    disp(sprintf("Number of Genes with Potentials Promoters : %d", length(key_s)));
    disp(sprintf("Percentage of Genes with Potential Promoters : %0.02f %%",length(key_s) / num_genes * 100 ));
    disp(sprintf("Upstream Position Distribution : %0.02f +/- %0.02f", mean(key_s), stdev(key_s)));
    disp(sprintf("Number of Consecutive Ws Distribution : %0.02f +/- %0.02f", mean(consec_Ws), stdev(consec_Ws)));
    disp(sprintf("Maximum Number of Consecutive Ws : %d", max(consec_Ws)));
end
disp("Position Probability Matrix");
disp(ppm);
disp(sprintf("Upstream Position Distribution for Stat. Align. : %0.02f +/- %0.02f", mean(stat_align_pos), stdev(stat_align_pos)));
disp(sprintf("Number of Genes with Potential Promoters (stat. alignment) : %d", length(stat_align_pos)));
disp(sprintf("Percentage of Genes with Potential Promoters (stat. alignment) : %0.02f %%",length(stat_align_pos) / num_genes * 100 ));

if local_search==1 then
    figure();
    histplot((1:1:thresh_up), key_s, style=2, normalization=%f); // upstream position histogram
    xlabel('upstream position','fontsize',5.5);
    ylabel('frequency','fontsize',5.5);
    if intact==1 then
        title("Upstream Position Histogram for Intact Query Local Search",'fontsize',5.5)
    else
        title("Upstream Position Histogram for Non-Intact Query Local Search",'fontsize',5.5)    
    end
    figure();
    histplot(5, consec_Ws, style=3, normalization=%f); // No. of consecutive Ws histogram
    xlabel('number of consecutive Ws', 'fontsize',5.5);
    ylabel('frequency', 'fontsize',5.5);
    if intact==1 then
        title("Consecutive Number of Ws Histogram for Intact Query Local Search", 'fontsize',5.5)
    else
        title("Consecutive Number of Ws Histogram for Non-Intact Query Local Search", 'fontsize',5.5)    
    end
    figure();
    histplot((1:1:thresh_up),stat_align_pos, style=5, normalization=%f); // upstream position histogram for stat. align.
    xlabel('upstream position', 'fontsize',5.5);
    ylabel('frequency', 'fontsize',5.5);
    if intact==1 then
        title("Upstream Position Histogram for Intact Statistical Alignment", 'fontsize',5.5)
    else
        title("Upstream Position Histogram for Non-Intact Statistical Alignment", 'fontsize',5.5)    
    end
else
    figure();
    histplot((1:1:thresh_up),stat_align_pos, style=5, normalization=%f);
    xlabel('upstream position', 'fontsize',5.5);
    ylabel('frequency', 'fontsize',5.5);
    title("Upstream Position Histogram for Combined Statistical Alignment", 'fontsize',5.5)  
end

