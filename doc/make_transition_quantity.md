# How to make Transition-quantity

## Requirements
* alignment database: https://cdn.andgo.io/database/processed/files.txt
* MTRAP version >= 2.0.0

## Make Transition-quantity
1. Store alignment database.
2. Edit configures of `[project]/utils/make_tq.sh`.
3. Run `[project]/utils/make_tq.sh`.

## Notes

### Database structure

```
./[database]/ref/*.msf             reference alignment (MSF format)
./[database]/ref/*.xml             reference alignment (XML format, used for BaliScore)
./[database]/ref/*.ref_fasta       reference alignemnt (Captal Letter means Homologous region, used for QScore)
./[database]/inputdata/*.fasta     input sequences (All capital letters, no gap)

./[database]/*.dnd junks w.r.t. ClustalW
./[database]/*.tfa junks w.r.t. BAliBASE
```

Notes:
```
./homstrad/ref  All captal letter (No homologous region annotations)
./prefab4/ref   All Captal letter (No homologous region annotations)
./SABmark_*/    Unmaintained
```

## References
* original database: https://cdn.andgo.io/database/original/files.txt
* test training infomation: https://cdn.andgo.io/database/original/files.txt
