# synthetic_data_fix

Add `GQ`, `DP` and `AD` to synthetic `GT`-only VCF files. 


## Installation

### Direct Download
```bash
curl -sSL https://raw.githubusercontent.com/hmartiniano/synthetic_data_fix/main/fix_synthetic_data.py -o fix_synthetic_data.py
```

### Pip Install
```bash
pip install git+https://github.com/hmartiniano/synthetic_data_fix.git
```

### Local Clone & uv
```bash
git clone git@github.com:hmartiniano/synthetic_data_fix.git
cd synthetic_data_fix
uv sync
```

## Usage

### 1. Streaming VCF input and output to bcftools (recommended to get bgzipped output VCF). All samples get the same values.
```bash
bcftools view synthetic.vcf.gz | python fix_synthetic_data.py | bgzip -@ 4 -c > synthetic_fixed.vcf.gz
bcftools index -t synthetic_fixed.vcf.gz
```

### 2. Same as above, but the outut VCF has 0 to 5% simulated low-quality genotypes
```bash
bcftools view synthetic.vcf.gz | python fix_synthetic_data.py \
    --simulate-low-quality --max-low-quality-fraction 0.05 --seed 42 | bgzip -@ 4 -c > synthetic_fixed.vcf.gz
bcftools index -t synthetic_fixed.vcf.gz
```

