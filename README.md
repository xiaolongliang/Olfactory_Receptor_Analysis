### Fst statistics of olfactory receptor genes and bitter taste receptor genes  
#### 1. calculated the wondows Fst of gene and chosen randomly Non-coding region that located at least 200 kb away from known gene  
```python
python Generate.Prir.Fst.py > fst.sh
sh fst.sh
```
#### 2.  Pickup the mean Fst of gene and Non-coding region  
```shell
sh Datadeal.sh
```
#### 3. T test to compare the Fst value of gene and Non-coding region  
```R
ls File/*.gene |while read line;do id=${line#*/};id=${id%%.*};Rscript T-test.R $id >> T-test.res;done
```
#### 4. plot
```R
Rscript T-test.res.plot.R
```

