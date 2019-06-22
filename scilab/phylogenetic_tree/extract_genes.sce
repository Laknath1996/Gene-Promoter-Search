
clc;
clear all;
close();

exec('helper_functions.sce');

//proteinName = "diaminopimelate decarboxylase";
proteinName = 'DnaK';
numOfgenes = 22;

AssignmentFolder = '/Users/ashwin/Semester 7/Genomic Signal Processing/GSP_Codes/scilab/phylogenetic_tree'
ProteinFolderPath = '/Users/ashwin/Semester 7/Genomic Signal Processing/GSP_Codes/scilab/phylogenetic_tree/protein'
FastaFolderPath = '/Users/ashwin/Semester 7/Genomic Signal Processing/GSP_Codes/scilab/phylogenetic_tree/fasta'

//NZ_CP014692.1 Acetobacter aceti
//NC_013209.1/AP011121.1 Acetobacter pasteurianus
//NZ_CP023657.1/CP023657.1 Acetobacter pomorum
//NZ_CP022699.1/CP022699.1 Acetobacter tropicalis
//NZ_CP014687.1/CP014687.1 Acetobacter persici
//NZ_LN606600.1/LN606600.1 Acetobacter senegalensis
//NZ_CP011120.1/CP011120.1 Acetobacter oryzifermentans
//NZ_CP015164.1 Acetobacter ascendens
//NZ_CP015168.1 Acetobacter ascendens
//NZ_CP021524.1 Acetobacter ascendens
//NZ_CP022374.1/CP022374.1 Acetobacter oryzifermentans
//NZ_AP018515.1/AP018515.1 Acetobacter orientalis
//NZ_CP023189.1/CP023189.1 Acetobacter pomorum
//NC_017100.1/AP011128.1 Acetobacter pasteurianus
//NC_017121.1/AP011135.1 Acetobacter pasteurianus
//NC_017125.1/AP011142.1 Acetobacter pasteurianus
//NC_017146.1/AP011149.1 Acetobacter pasteurianus
//NZ_LN609302.1 Acetobacter ghanensis
//NC_017111.1/AP011156.1 Acetobacter pasteurianus
//NC_017150.1/AP011163.1 Acetobacter pasteurianus
//NC_017108.1/AP011170.1 Acetobacter pasteurianus
//NZ_AP014881.1/AP014881.1 Acetobacter pasteurianus

files =  findfiles(ProteinFolderPath,'*.txt');
fileNames = gsort(files,'g','i');

gene_pos=zeros(numOfgenes,3);
for i=1:numOfgenes
    fileName=ProteinFolderPath+"\"+fileNames(i);
    f=mopen(fileName,'r');
    disp(fileName);
    while(meof(f)==0)
       line=mgetl(f,1);
       indx=strindex(line,ascii(9));                            
       line=strsplit(line,indx);
       line=stripblanks(line,%t); //strip leading and trailing blanks including tabs
//           if (line(size(line)(1))==proteinName)||(line(size(line)(1))=='MULTISPECIES: '+proteinName) then
       if (grep(line(size(line)(1)),proteinName)==1)
           gPos(i,1)=strtod(line(3));
           gPos(i,2)=strtod(line(4));
           if (line(5)=='+') then               
               gPos(i,3)=1;
           else
               gPos(i,3)=0;
           end
           break
       end
    end
    mclose(f);
end

fileName = AssignmentFolder+"\"+proteinName+'_protein locations.txt';
f=mopen(fileName,'w');
for i=1:numOfgenes
    mputl((string(i)+ascii(9)+string(gPos(i,1))+ascii(9)+string(gPos(i,2))+ascii(9)+string(gPos(i,3))),f);
end
mclose(f);

fileName = AssignmentFolder+"\"+proteinName+'_protein locations.txt';
gPos=zeros(numOfgenes,3);
f=mopen(fileName,'r');
for i=1:numOfgenes
    line=mgetl(f,1);
    indx=strindex(line,ascii(9));                            
    line=strsplit(line,indx);
    line=stripblanks(line,%t);
    gPos(i,1)=strtod(line(2));
    gPos(i,2)=strtod(line(3));
    gPos(i,3)=strtod(line(4));
end   
mclose(f)

f_files =  findfiles(FastaFolderPath,'*.fasta');
f_fileNames = gsort(f_files,'g','i');

geneSet={};
organisms = {};

//NZ_CP014692.1
//NC_013209.1/AP011121.1
//NZ_CP023657.1/CP023657.1
//NZ_CP022699.1/CP022699.1
//NZ_CP014687.1/CP014687.1
//NZ_LN606600.1/LN606600.1
//NZ_CP011120.1/CP011120.1
//NZ_CP015164.1
//NZ_CP015168.1
//NZ_CP021524.1
//NZ_CP022374.1/CP022374.1
//NZ_AP018515.1/AP018515.1
//NZ_CP023189.1/CP023189.1
//NC_017100.1/AP011128.1
//NC_017121.1/AP011135.1
//NC_017125.1/AP011142.1
//NC_017146.1/AP011149.1
//NZ_LN609302.1
//NC_017111.1/AP011156.1
//NC_017150.1/AP011163.1
//NC_017108.1/AP011170.1
//NZ_AP014881.1/AP014881.1

gene_lengths = [];
for i=1:numOfgenes
    fileName=FastaFolderPath+"\"+f_fileNames(i);
    organisms{i} = part(f_fileNames(i), [1:$-6]);
    f=mopen(fileName,'r');
    disp(fileName);
    gene=get_fasta_at(fileName,gPos(i,1),gPos(i,2),gPos(i,3));
    disp(length(gene));
    disp(ascii(gene));
    gene_lengths = [gene_lengths length(gene)];
    geneSet{i} = gene;
end

idx = [21, 22, 19, 17, 16, 15, 14, 20, 2, 11, 13, 7, 18, 6, 1, 8, 9, 10, 3, 4, 5, 12];
lengths = gene_lengths(idx);
disp(lengths);

f=mopen(AssignmentFolder+"\"+'set_of_extracted_genes_dnak.txt','w');

for i=1:numOfgenes
//    header='>'+ascii(65-1+i);
    organism = organisms{i}
    header='>'+organism;
    //header = ">" + ;
    mputl(header,f);
    mputl(ascii(geneSet{i}),f);
end
mclose(f);

