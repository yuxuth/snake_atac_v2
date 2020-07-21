# snakme-ChIPseq 


snakemake -p -j 99 --cluster-config cluster.json --cluster "sbatch -J {cluster.job} --mem={cluster.mem} -N 1 -n {threads} -o {cluster.out} -e {cluster.err} " &> log &


adapt from Ming Tang's repo [pyflow-CHIPseq](https://github.com/crazyhottommy/pyflow-ChIPseq).


