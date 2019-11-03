# TODO

1. Create all possible HLA structs at compile time.
  ```
     if you do it in build.rs, you don't need any serialization
     the data would not be stored in json or protobuf format but embedded into your tool
     in your build.rs, you would have code that
     1. reads the file and parses it (should be rather simple if it's tab-separated values)
     2. creates a bunch of rust code like const DATA : MyStruct { foo: 1234, bar: 5678, … }
     3. writes that to a file like generated/data.rs
  ```

# DATA

- NetMHCpan 4.0 allele list: http://www.cbs.dtu.dk/services/NetMHCpan/MHC_allele_names.txt
- Obtain table of ligands:
  - A: `curl -s "https://www.ebi.ac.uk/cgi-bin/ipd/kir/retrieve_ligands.cgi?A*" | grep -i -e '</\?table\|</\?td\|</\?tr\|</\?th' | sed 's/^[\ \t]*//g' | tr -d '\n' | sed 's/<\/tr[^>]*>/\n/Ig'  | sed 's/<\/\?\(table\|tr\)[^>]*>//Ig' | sed 's/^<t[dh][^>]*>\|<\/\?t[dh][^>]*>$//Ig' | sed 's/<\/t[dh][^>]*><t[dh][^>]*>/\t/Ig' | gzip > b_ligand_groups.tsv.gz`
  - B: `curl -s "https://www.ebi.ac.uk/cgi-bin/ipd/kir/retrieve_ligands.cgi?B*" | grep -i -e '</\?table\|</\?td\|</\?tr\|</\?th' | sed 's/^[\ \t]*//g' | tr -d '\n' | sed 's/<\/tr[^>]*>/\n/Ig'  | sed 's/<\/\?\(table\|tr\)[^>]*>//Ig' | sed 's/^<t[dh][^>]*>\|<\/\?t[dh][^>]*>$//Ig' | sed 's/<\/t[dh][^>]*><t[dh][^>]*>/\t/Ig' | gzip > b_ligand_groups.tsv.gz`
  - C: `curl -s "https://www.ebi.ac.uk/cgi-bin/ipd/kir/retrieve_ligands.cgi?C*" | grep -i -e '</\?table\|</\?td\|</\?tr\|</\?th' | sed 's/^[\ \t]*//g' | tr -d '\n' | sed 's/<\/tr[^>]*>/\n/Ig'  | sed 's/<\/\?\(table\|tr\)[^>]*>//Ig' | sed 's/^<t[dh][^>]*>\|<\/\?t[dh][^>]*>$//Ig' | sed 's/<\/t[dh][^>]*><t[dh][^>]*>/\t/Ig' | gzip > b_ligand_groups.tsv.gz`

