# download Proteobacterial genomes (1839 at the time)

datasets summary genome taxon Proteobacteria --assembly-source RefSeq --reference --as-json-lines --exclude-atypical --assembly-level complete

datasets download genome taxon Proteobacteria --assembly-source RefSeq --reference --dehydrated --exclude-atypical --include genome

# predict and translate coding sequences

for i in *.fasta; do prodigal -a faa/${i%%.fasta}.faa -i $i -q; done

# create tree and look for KO-targets

GToTree -A fasta.files -H Proteobacteria -G 0 -K ../../flag.kegg -n 7 -j 7 -M 7 -c 0.33 -o proto_out

# look again with hmmsearch

parallel -j 3 hmmsearch --tblout proto_out/KO_search_results/hmmsearch_results/{.}.tblout -T 0 --domT 0 proto_out/KO_search_results/FAP.hmm {} ::: *.faa

