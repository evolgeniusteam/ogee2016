## data download
ftp://ftp.informatics.jax.org/pub/reports/index.html#pheno

## definition of essential genes
Genes associated with MP terms (see file: MPheno_OBO.ontology) whose names contain any of the following strings are essential:
* lethality
* infertility
* premature death

## get a list of essential and non-essential genes
```perl
perl get_mouse_essential_genes.pl -o mouse_essential_genes.txt
```

## prepare the essential genes for loading into database
```shell
more mouse_essential_genes.txt | awk 'BEGIN{FS="\t";OFS="\t"} {print 349, $1, $2, 0, 10090, 0, 0, 0, -1}' > mouse.ess
```

## update database / load genes into database
```mysql
update data_summary set dateadded = now() where datasetid = 349;

delete from essential_genes where datasetid = 349;
select count(*) from essential_genes where datasetid = 349;
load data local infile "/Volumes/data/WebServers/ogee2016/perlscripts/mouse.ess" INTO table essential_genes fields terminated by "\t";
```
