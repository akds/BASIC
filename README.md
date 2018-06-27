## BASIC: BCR (and now TCR) assembly from single cells


### Pre-requisites and installation
* BASIC was developed and tested using Python 3.4.4 and 2.7.10
* BASIC requires Bowtie2 to run.
* The `db/` folder must remain with the BASIC.py file

### Basic usage example

```python
python BASIC.py -b <path_to_Bowtie2> -PE_1 <R1.fastq.gz> -PE_2 <R2.fastq.gz> -g <species> -i <receptor>
```

**Parameters**
+ **species**: "human" or "mouse"

    Note: other species are possible by adding the appropriate bowtie2 indices and following the existing `db/` directory structure
+ **receptor**: "BCR" or "TCR"

For more advanced usage or other use cases e.g. single end reads:
```python
python BASIC.py -h
```

### More information
http://ttic.uchicago.edu/~aakhan/BASIC/

### Version
1.3.1 (2018/06/27)

### Contact
Please contact: Aly Azeem Khan <aakhan@ttic.edu> for any questions or comments.
More recent updates have been contributed by: Derek Croote <dcroote@stanford.edu>

### License
Software provided to academic users under MIT License
