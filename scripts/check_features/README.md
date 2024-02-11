Example:

```shell
cd ~

wfmash /lizardfs/guarracino/wfmash-paper/assemblies/scerevisiae/scerevisiae7.fa.gz -t 48 --lower-triangular > scerevisiae7.paf

bash /lizardfs/guarracino/wfmash-paper/scripts/check_features/feature_checker.sh scerevisiae7.paf /lizardfs/guarracino/wfmash-paper/busco/scerevisiae/S288C/saccharomycetes_odb10/S288C.busco-genes.bed /scratch

ls scerevisiae7.*.tsv -lh
    -rw-r--r-- 1 guarracino guarracino 1.3M Feb 11 10:27 scerevisiae7.report.features.tsv
    -rw-r--r-- 1 guarracino guarracino  15K Feb 11 10:27 scerevisiae7.report.genomes.tsv

head scerevisiae7.report.features.tsv -n 5 | column -t
    feature.name  query             query.feature.start  query.feature.end  query.strand  target            target.feature.start  target.feature.end  aligned.bases  not.aligned.in.query.bp  not.aligned.in.target.bp
    12947at4891   DBVPG6765#1#chrI  147566               148823             +             DBVPG6044#1#chrI  141399                142647              1248           9                        0
    13507at4891   DBVPG6765#1#chrI  157130               158267             +             DBVPG6044#1#chrI  156929                158066              1137           0                        0
    15784at4891   DBVPG6765#1#chrI  61639                62743              +             DBVPG6044#1#chrI  60782                 61886               1104           0                        0
    17975at4891   DBVPG6765#1#chrI  95485                96889              +             DBVPG6044#1#chrI  94618                 96022               1404           0                        0


head scerevisiae7.report.genomes.tsv -n 5 | column -t
    query               num.features.in.query  target              num.feature.in.target  num.features.in.common  aligned.bases  not.aligned.in.query.bp  not.aligned.in.target.bp
    DBVPG6765#1#chrI    30                     DBVPG6044#1#chrI    30                     30                      51090          762                      75
    DBVPG6765#1#chrII   146                    DBVPG6044#1#chrII   146                    145                     224399         519                      4115
    DBVPG6765#1#chrIII  46                     DBVPG6044#1#chrIII  46                     46                      72119          354                      1146
    DBVPG6765#1#chrIV   284                    DBVPG6044#1#chrIV   284                    284                     431509         17461                    3017
```
