

```python
%run setup.ipynb
```


<style type="text/css">
.container {
    width: 100%;
}
#maintoolbar {
    display: none;
}
#header-container {
    display: none;
}
div#notebook {
    padding-top: 0;
}
</style>



```python
region_vgsc = SeqFeature('2L', 2358158, 2431617)
region_vgsc
```




    <__main__.SeqFeature at 0x7f72c4cfa668>




```python
region_vgsc.region_str
```




    '2L:2358158-2431617'




```python
!bcftools --help
```

    
    Program: bcftools (Tools for variant calling and manipulating VCFs and BCFs)
    Version: 1.3.1 (using htslib 1.3.1)
    
    Usage:   bcftools [--version|--version-only] [--help] <command> <argument>
    
    Commands:
    
     -- Indexing
        index        index VCF/BCF files
    
     -- VCF/BCF manipulation
        annotate     annotate and edit VCF/BCF files
        concat       concatenate VCF/BCF files from the same set of samples
        convert      convert VCF/BCF files to different formats and back
        isec         intersections of VCF/BCF files
        merge        merge VCF/BCF files files from non-overlapping sample sets
        norm         left-align and normalize indels
        plugin       user-defined plugins
        query        transform VCF/BCF into user-defined formats
        reheader     modify VCF/BCF header, change sample names
        view         VCF/BCF conversion, view, subset and filter VCF/BCF files
    
     -- VCF/BCF analysis
        call         SNP/indel calling
        consensus    create consensus sequence by applying VCF variants
        cnv          HMM CNV calling
        filter       filter VCF/BCF files using fixed thresholds
        gtcheck      check sample concordance, detect sample swaps and contamination
        roh          identify runs of autozygosity (HMM)
        stats        produce VCF/BCF stats
    
     Most commands accept VCF, bgzipped VCF, and BCF with the file type detected
     automatically even when streaming from a pipe. Indexed VCF and BCF will work
     in all situations. Un-indexed VCF and BCF and streams will work in most but
     not all situations.
    



```python
input_filename = '../ngs.sanger.ac.uk/production/ag1000g/phase1/AR3/variation/main/vcf/ag1000g.phase1.ar3.%s.vcf.gz' % region_vgsc.seqid
output_filename = '../data/ag1000g.phase1.ar3.%s.%s.%s.vcf' % tuple(region_vgsc)
output_filename
```




    '../data/ag1000g.phase1.ar3.2L.2358158.2431617.vcf'




```python
!bcftools view
```

    
    About:   VCF/BCF conversion, view, subset and filter VCF/BCF files.
    Usage:   bcftools view [options] <in.vcf.gz> [region1 [...]]
    
    Output options:
        -G,   --drop-genotypes              drop individual genotype information (after subsetting if -s option set)
        -h/H, --header-only/--no-header     print the header only/suppress the header in VCF output
        -l,   --compression-level [0-9]     compression level: 0 uncompressed, 1 best speed, 9 best compression [-1]
              --no-version                  do not append version and command line to the header
        -o,   --output-file <file>          output file name [stdout]
        -O,   --output-type <b|u|z|v>       b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]
        -r, --regions <region>              restrict to comma-separated list of regions
        -R, --regions-file <file>           restrict to regions listed in a file
        -t, --targets [^]<region>           similar to -r but streams rather than index-jumps. Exclude regions with "^" prefix
        -T, --targets-file [^]<file>        similar to -R but streams rather than index-jumps. Exclude regions with "^" prefix
            --threads <int>                 number of extra output compression threads [0]
    
    Subset options:
        -a, --trim-alt-alleles        trim alternate alleles not seen in the subset
        -I, --no-update               do not (re)calculate INFO fields for the subset (currently INFO/AC and INFO/AN)
        -s, --samples [^]<list>       comma separated list of samples to include (or exclude with "^" prefix)
        -S, --samples-file [^]<file>  file of samples to include (or exclude with "^" prefix)
            --force-samples           only warn about unknown subset samples
    
    Filter options:
        -c/C, --min-ac/--max-ac <int>[:<type>]      minimum/maximum count for non-reference (nref), 1st alternate (alt1), least frequent
                                                       (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]
        -f,   --apply-filters <list>                require at least one of the listed FILTER strings (e.g. "PASS,.")
        -g,   --genotype [^]<hom|het|miss>          require one or more hom/het/missing genotype or, if prefixed with "^", exclude sites with hom/het/missing genotypes
        -i/e, --include/--exclude <expr>            select/exclude sites for which the expression is true (see man page for details)
        -k/n, --known/--novel                       select known/novel sites only (ID is not/is '.')
        -m/M, --min-alleles/--max-alleles <int>     minimum/maximum number of alleles listed in REF and ALT (e.g. -m2 -M2 for biallelic sites)
        -p/P, --phased/--exclude-phased             select/exclude sites where all samples are phased
        -q/Q, --min-af/--max-af <float>[:<type>]    minimum/maximum frequency for non-reference (nref), 1st alternate (alt1), least frequent
                                                       (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]
        -u/U, --uncalled/--exclude-uncalled         select/exclude sites without a called genotype
        -v/V, --types/--exclude-types <list>        select/exclude comma-separated list of variant types: snps,indels,mnps,other [null]
        -x/X, --private/--exclude-private           select/exclude sites where the non-reference alleles are exclusive (private) to the subset samples
    



```python
!bcftools view --output-file {output_filename} --output-type v --regions {region_vgsc.region_str} {input_filename}
```


```python
!ls -lh ../data
```

    total 6.1M
    -rw-r--r-- 1 aliman kwiat-cluster-users 142M Feb 27 12:34 ag1000g.phase1.ar3.2L.2358158.2431617.vcf
    -rw-rw-r-- 1 aliman                2834 2.6M Feb 23 09:42 Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.4.gff3.gz
    -rw-rw-r-- 1 aliman                2834  47K Feb 23 09:42 davies_vgsc_model_20170125.gff3
    -rw-rw-r-- 1 aliman                2834 1.5M Feb 27 11:55 demo.npy
    -rw-r--r-- 1 aliman kwiat-cluster-users  19K Feb 27 11:45 domestica_gambiae_map.txt
    -rw-r--r-- 1 aliman kwiat-cluster-users 4.3K Feb 27 11:45 domestica_gambiae_PROT_MEGA.fas
    -rw-r--r-- 1 aliman kwiat-cluster-users 4.2K Feb 27 11:45 domestica_gambiae_PROT_MEGA.fas.flat
    -rw-r--r-- 1 aliman kwiat-cluster-users   66 Feb 27 11:45 domestica_gambiae_PROT_MEGA.fas.gdx
    -rw-rw-r-- 1 aliman                2834  11M Feb 27 12:29 tbl_variants_phase1.pkl
    -rw-rw-r-- 1 aliman                2834 4.3M Feb 27 12:29 tbl_variants_phase1.txt



```python
!bgzip -f {output_filename}
```


```python
!tabix -p vcf -f {output_filename}.gz
```


```python
!ls -lh ../data
```

    total 34M
    -rw-r--r-- 1 aliman kwiat-cluster-users  28M Feb 27 12:34 ag1000g.phase1.ar3.2L.2358158.2431617.vcf.gz
    -rw-r--r-- 1 aliman kwiat-cluster-users  187 Feb 27 12:34 ag1000g.phase1.ar3.2L.2358158.2431617.vcf.gz.tbi
    -rw-rw-r-- 1 aliman                2834 2.6M Feb 23 09:42 Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.4.gff3.gz
    -rw-rw-r-- 1 aliman                2834  47K Feb 23 09:42 davies_vgsc_model_20170125.gff3
    -rw-rw-r-- 1 aliman                2834 1.5M Feb 27 11:55 demo.npy
    -rw-r--r-- 1 aliman kwiat-cluster-users  19K Feb 27 11:45 domestica_gambiae_map.txt
    -rw-r--r-- 1 aliman kwiat-cluster-users 4.3K Feb 27 11:45 domestica_gambiae_PROT_MEGA.fas
    -rw-r--r-- 1 aliman kwiat-cluster-users 4.2K Feb 27 11:45 domestica_gambiae_PROT_MEGA.fas.flat
    -rw-r--r-- 1 aliman kwiat-cluster-users   66 Feb 27 11:45 domestica_gambiae_PROT_MEGA.fas.gdx
    -rw-rw-r-- 1 aliman                2834  11M Feb 27 12:29 tbl_variants_phase1.pkl
    -rw-rw-r-- 1 aliman                2834 4.3M Feb 27 12:29 tbl_variants_phase1.txt



```python

```
