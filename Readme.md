🚀 Como Funciona na Prática (Linha de Comando e R)1.Alinhamento com HISAT2:Terminal / Linux.Mapeie suas sequências de RNA (reads) contra o genoma de referência previamente indexado:

Bash# Mapeamento de reads paired-end e conversão direta para BAM ordenado
hisat2 -x /path/to/genome_index \
       -1 amostra1_R1.fastq.gz \
       -2 amostra1_R2.fastq.gz \
       -S amostra1.sam

# Converter SAM para BAM e ordenar (necessário para economizar espaço e memória)
samtools sort -@ 4 -o amostra1.bam amostra1.sam
2.Quantificação com featureCounts:Terminal / Subread package.Converta o alinhamento físico em uma matriz de contagens absolutas por gene usando o arquivo de anotação GTF:BashfeatureCounts -p -T 4 \
              -a anotacao.gtf \
              -o contagens_matriz.txt \
              amostra1.bam amostra2.bam amostra3.bam
3.Análise Diferencial com DESeq2:Linguagem R.No R, importe a matriz do featureCounts para estimar a variação biológica e identificar a expressão diferencial:Rlibrary(DESeq2)

# Importar matriz e metadata
cts <- as.matrix(read.table("contagens_matriz.txt", header=TRUE, row.names=1))
colData <- read.csv("metadados.csv", row.names=1) # contém os grupos (ex: Controle vs Tratado)

# Criar objeto DESeq
dds <- DESeqDataSetFromMatrix(countData = cts, colData = colData, design = ~ condicao)

# Executar a pipeline estatística
dds <- DESeq(dds)
res <- results(dds)

# Filtrar genes estatisticamente significativos (p-ajustado < 0.05)
genes_sig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)


Dica: O DESeq2 espera receber contagens brutas (raw counts), e não valores normalizados como RPKM/TPM. O próprio pacote calcula internamente os fatores de normalização para corrigir variações na profundidade do sequenciamento.
