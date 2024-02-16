Example:

```shell
DIR_WF_PAPER=/lizardfs/guarracino/wfmash-paper/

cd ~

wfmash $DIR_WF_PAPER/assemblies/scerevisiae/scerevisiae7.fa.gz -t 48 --lower-triangular > $DIR_WF_PAPER/alignment/scerevisiae/scerevisiae7.paf

export PATH="/home/guarracino/tools/feature_level_report/target/release:$PATH"

bash $DIR_WF_PAPER/scripts/check_features/check_features.sh $DIR_WF_PAPER/alignment/scerevisiae/scerevisiae7.paf $DIR_BASE/alignment/scerevisiae/scerevisiae.busco-genes.single.bed 50 ./scerevisiae7.m50 48 /scratch

ls scerevisiae7.m50.*.tsv -lh
    -rw-r--r-- 1 guarracino guarracino 1.3M Feb 14 15:45 scerevisiae7.m50.report.features.tsv
    -rw-r--r-- 1 guarracino guarracino  15K Feb 14 15:45 scerevisiae7.m50.report.genomes.tsv


head scerevisiae7.m50.report.features.tsv -n 5 | column -t
    feature.name  query             query.feature.start  query.feature.end  query.strand  target            target.feature.start  target.feature.end  aligned.bp  not.aligned.in.query.bp  not.aligned.in.target.bp
    12947at4891   DBVPG6765#1#chrI  147566               148823             +             DBVPG6044#1#chrI  141399                142647              1248        9                        0
    13507at4891   DBVPG6765#1#chrI  157130               158267             +             DBVPG6044#1#chrI  156929                158066              1137        0                        0
    15784at4891   DBVPG6765#1#chrI  61639                62743              +             DBVPG6044#1#chrI  60782                 61886               1104        0                        0
    17975at4891   DBVPG6765#1#chrI  95485                96889              +             DBVPG6044#1#chrI  94618                 96022               1404        0                        0

head scerevisiae7.m50.report.genomes.tsv -n 5 | column -t
    query               num.features.in.query  target              num.feature.in.target  num.features.in.common  aligned.bases  not.aligned.in.query.bp  not.aligned.in.target.bp
    DBVPG6765#1#chrI    30                     DBVPG6044#1#chrI    30                     30                      51090          747                      72
    DBVPG6765#1#chrII   146                    DBVPG6044#1#chrII   146                    145                     224399         498                      4065
    DBVPG6765#1#chrIII  46                     DBVPG6044#1#chrIII  46                     46                      72119          348                      1143
    DBVPG6765#1#chrIV   284                    DBVPG6044#1#chrIV   284                    284                     431509         17309                    2966
```
