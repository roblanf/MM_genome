
Setting up folders for genome_assembly branch readme: (Names to be added/changed as I determine the key sections further)
```bash
mkdir -p 01_qc/{scripts,config}
mkdir -p 02_filtering/{scripts,config}
mkdir -p 03_assembly/{scripts,config}
mkdir -p 04_binning/{scripts,config}
mkdir -p 06_comparison/{scripts,config}
```
Adding readme stub for each stage
```bash
for d in 01_qc 02_filtering 03_assembly 04_binning 06_comparison; do
	echo "# ${d}" > "${d}/README.md"
	echo "" >> "${d}/README.md"
	echo "## Method" >> "${d}/README.md"
	echo "" >> "${d}/README.md"
	echo "## Status" >> "${d}/README.md"
done
```

