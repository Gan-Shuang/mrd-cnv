# mrd-cnv
Ultra low WGS analyze CNV infer tumor purity in human plasma pipline
## Introduction  
base on ichorCNA ,with built up environment.  
## Method  
> build environment
```
cd mrd-cnv
singularity build ichorcna.sif docker://ganshuang0925/ichorcna_load
```
```
singularity exec -B /mnt:/mnt ichorcna.sif python3 run_ichorcna.py -i [PATH/bam] -o [/PATH_to_result_dir/]
```
## Results
ex:test_example
![WCUCADLYF0708_genomeWide](https://user-images.githubusercontent.com/50703435/197149192-23ec41e0-5e4e-402d-a834-db0f184fa9f4.png)
![chr1](https://user-images.githubusercontent.com/50703435/197149202-fc75002a-69c1-4ebd-bccd-d175072dda5d.png)
