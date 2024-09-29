## clinvar-ancestry

**A little command-line tool to process clinvar.vcf together with AncestryDNA.txt.**

The tool filters clinvar.vcf on GENEINFO and then matches up all AncestryDNA.txt entries which refer to that rsID. For each match it outputs a combined line to stdout. It also parses the INFO part of the clinvar.vcf line.

Matching lines are at the end written to the file "filtered-DNA.json" in JSON format, which can be read with other tools, such as python pandas.

An initial commit with a command-line utility that takes a few arguments and filters the rsID list in AncestryDNA.txt and outputs matching data from clinvar.vcf.

### Examples:

Output all rsIDs where GENEINFO begins with MTHFR:

`clinvar-ancestry clinvar.vcf AncestryDNA.txt --gene "MTHFR.*"`

Output rsIDs where the clinvar GENEINFO begins with CBS or MTHFR and at least one allele in the AncestryDNA.txt entry is different from the REF allele:

`clinvar-ancestry clinvar.vcf AncestryDNA.txt --gene "(CBS|MTHFR).*" --discrepancies`

### Options:

`--gene <regex>` : Output only rsIDs related to entries where GENEINFO matches this regex.

`--discrepancies` : Output only heterozygous or homozygous rsIDs.

`--unmatched-rsid` : Output also entries where there is no matching rsID in the clinvar.vcf.
