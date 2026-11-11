# Organelle

This README is purposed to store additional details required to conduct the `Assembly of organelle genomes` step from the main <a href="../../README.md">`README.md`</a>.

## Download the organelles database for oatk
```
db_dir="04_organelle/db"
base_url="https://raw.githubusercontent.com/c-zhou/OatkDB/main/v20230921"

files=(
    "magnoliopsida_mito.fam"
    "magnoliopsida_mito.fam.h3f"
    "magnoliopsida_mito.fam.h3i"
    "magnoliopsida_mito.fam.h3m"
    "magnoliopsida_mito.fam.h3p"
    "magnoliopsida_pltd.fam"
    "magnoliopsida_pltd.fam.h3f"
    "magnoliopsida_pltd.fam.h3i"
    "magnoliopsida_pltd.fam.h3m"
    "magnoliopsida_pltd.fam.h3p"
)

# download individual files
for file in "${files[@]}"; do
    wget -q "${base_url}/${file}" -O "${db_dir}/${file}"
    echo "  Downloaded: ${file}"
done
```