## Downloading Satpathy data from ENA database

## definition of the urls
srr_ids <- readLines('../data/satpathy2016/SRR_Acc_List.txt')
base_url <- 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/' 
dir1 <- substr(srr_ids, 1, 6)
last_char <- sapply(srr_ids, function(s) substr(s, nchar(s), nchar(s)))
urls <- sapply(seq_along(srr_ids), 
              function(i) paste0(base_url,
                                 dir1[i],s
                                 '/00',
                                 last_char[i],
                                 '/', srr_ids[i]))

params <- 'wget -r --tries=50 -P ../data/satpathy2016/fastq/ '
comms <-  paste0(params, urls, '/')
comms

writeLines(comms, '../scripts/download_satpathy_ena.sh')
