# After-Loter: LAI-Results-Analysis-and-Visualizatioon-
Loter Segment Analyzer 是一个用于分析Loter（Local Ancestry Inference Tool）输出结果的Python工具。Loter是一种用于局部祖源推断的算法，其输出为一个二维NumPy数组，表示每个单倍型在每个SNP位点上的祖源赋值。然而，该数组仅使用SNP索引，缺乏物理位置信息。本工具通过将SNP索引映射回VCF文件中的物理位置（染色体和位置），从而能够提取具有实际基因组位置的单倍型片段。此外，工具还提供了片段统计、祖源比例计算以及可视化功能。

# 项目描述

## Loter Segment Analyzer

Loter Segment Analyzer 是一个用于分析Loter（Local Ancestry Inference Tool）局部祖源推断结果的Python工具包。该工具解决了Loter输出结果中仅包含SNP索引而缺乏物理位置信息的问题，通过将SNP索引映射回VCF文件的物理位置，实现了对单倍型祖源片段的精确识别、提取和分析。

### 主要功能：
- **SNP索引到物理位置的映射**：从VCF文件或预提取的TSV文件中提取SNP的染色体和位置信息
- **局部祖源片段识别**：根据Loter的祖源赋值结果识别连续的单倍型片段
- **跨染色体片段处理**：智能检测和处理跨染色体的片段，确保每个片段仅位于单一染色体
- **统计分析**：计算各祖源群体的全局比例、片段长度分布和染色体分布
- **结果导出**：生成详细的表格报告和可视化图表
- **可视化支持**：提供片段长度分布和染色体分布的可视化脚本

### 应用场景：
- 群体遗传学研究中的局部祖源推断分析
- 混合群体（如维吾尔族、塔吉克族、羌族）的祖源成分精细解析
- 比较不同参考群体对目标群体祖源推断的影响
- 研究自然选择和基因流在基因组中的分布模式

### 技术特点：
- 支持标准VCF文件格式输入
- 处理大规模基因组数据的高效算法
- 可自定义最小片段长度阈值
- 完整的统计报告和可视化输出
- 适用于多种祖源推断研究设计

---

# README.md

# Loter Segment Analyzer

A Python toolkit for analyzing local ancestry inference results from Loter (Local Ancestry Inference Tool). This tool bridges the gap between SNP indices in Loter output and physical genomic positions, enabling precise identification, extraction, and analysis of haplotype ancestry segments.

## Features

- **SNP Index to Physical Position Mapping**: Extract chromosome and position information from VCF files or pre-extracted TSV files
- **Local Ancestry Segment Identification**: Identify continuous haplotype segments based on Loter ancestry assignments
- **Cross-Chromosome Segment Handling**: Intelligent detection and processing of segments spanning multiple chromosomes
- **Statistical Analysis**: Calculate global ancestry proportions, segment length distributions, and chromosomal distributions
- **Comprehensive Output**: Generate detailed tabular reports and visualization plots
- **Visualization Support**: Built-in scripts for visualizing segment length distributions and chromosomal distributions

## Requirements

### Software Dependencies
- Python 3.6+
- bcftools (for extracting SNP positions from VCF files)

### Python Packages
- numpy >= 1.18.0
- pandas >= 1.0.0
- matplotlib >= 3.1.0
- seaborn >= 0.10.0

## Installation

### 1. Clone the Repository
```bash
git clone https://github.com/yourusername/loter-segment-analyzer.git
cd loter-segment-analyzer
```

### 2. Install Python Dependencies
```bash
pip install numpy pandas matplotlib seaborn
```

### 3. Install bcftools (for VCF processing)
```bash
# Ubuntu/Debian
sudo apt-get install bcftools

# CentOS/RHEL
sudo yum install bcftools

# macOS
brew install bcftools
```

### 4. Verify Installation
```bash
python loter_segment_analyzer.py --help
bcftools --version
```

## Quick Start

### Step 1: Extract SNP Positions (One-time Operation)
```bash
bcftools query -f '%CHROM\t%POS\n' \
  /path/to/your/vcf/file.vcf.gz \
  > snp_positions.tsv
```

### Step 2: Run Loter Segment Analysis
```bash
python loter_segment_analyzer.py \
  --loter-npy /path/to/loter_output.npy \
  --group-names "Han,Dai,Pumi,Wa" \
  --snp-map-tsv snp_positions.tsv \
  --output-prefix jino_analysis
```

### Step 3: Visualize Results (Optional)
```bash
python plot_segment_lengths.py jino_analysis_segments.tsv
```

## Detailed Usage

### Main Script: `loter_segment_analyzer.py`

#### Basic Command Structure:
```bash
python loter_segment_analyzer.py \
  --loter-npy <path_to_loter_npy> \
  --group-names <comma_separated_group_names> \
  [--vcf-file <path_to_vcf> | --snp-map-tsv <path_to_tsv>] \
  [options]
```

#### Required Arguments:
| Argument | Description | Example |
|----------|-------------|---------|
| `--loter-npy` | Path to Loter output .npy file | `loter_output.npy` |
| `--group-names` | Comma-separated reference population names | `"Han,Dai,Pumi,Wa"` |

#### Optional Arguments:
| Argument | Description | Default |
|----------|-------------|---------|
| `--vcf-file` | VCF file for SNP position extraction | None |
| `--snp-map-tsv` | Pre-extracted SNP positions TSV file | None |
| `--output-prefix` | Prefix for output files | `loter_analysis` |
| `--min-snps` | Minimum SNPs per segment (filter short segments) | 3 |
| `--handle-cross-chr` | Split segments that span multiple chromosomes | False |
| `--verbose` | Display detailed progress information | False |

#### Example Commands:

**Example 1: Basic Analysis with Pre-extracted SNP Positions**
```bash
python loter_segment_analyzer.py \
  --loter-npy /data/loter_results/Jino_vs_4refs.npy \
  --group-names "Han,Dai,Pumi,Wa" \
  --snp-map-tsv /data/snp_positions.tsv \
  --output-prefix Jino_vs_HDPW
```

**Example 2: Direct VCF Processing with Advanced Options**
```bash
python loter_segment_analyzer.py \
  --loter-npy /data/loter_results/Jino_vs_4refs.npy \
  --group-names "Han,Dai,Pumi,Wa,Miao,Yao" \
  --vcf-file /data/reference/9pops.subset.vcf.gz \
  --handle-cross-chr \
  --min-snps 5 \
  --verbose \
  --output-prefix Jino_vs_6refs
```

### Visualization Script: `plot_segment_lengths.py`

```bash
python plot_segment_lengths.py <segments_tsv_file>
```

**Example:**
```bash
python plot_segment_lengths.py Jino_vs_HDPW_segments.tsv
```

## Output Files

After successful execution, the following files will be generated:

### 1. **Ancestry Proportions**
- **File:** `<prefix>_ancestry_proportions.tsv`
- **Format:** Tab-separated values
- **Columns:** `ANCESTRY`, `PROPORTION`
- **Description:** Global ancestry proportions for each reference population

### 2. **Segment Information**
- **File:** `<prefix>_segments.tsv`
- **Format:** Tab-separated values
- **Columns:**
  - `HAPLOTYPE_ID`: Haplotype index (0-based)
  - `ANCESTRY`: Assigned ancestry population
  - `START_SNP`: Starting SNP index
  - `END_SNP`: Ending SNP index
  - `LENGTH_SNP`: Segment length in SNPs
  - `CHR_START`: Starting chromosome
  - `CHR_END`: Ending chromosome
  - `POS_START`: Starting physical position (bp)
  - `POS_END`: Ending physical position (bp)
  - `LENGTH_BP`: Segment length in base pairs
  - `CHR_STATUS`: Chromosome status (chromosome name or 'cross_chromosome')

### 3. **Segment Statistics**
- **File:** `<prefix>_segment_statistics.tsv`
- **Format:** Tab-separated values
- **Columns:**
  - `ANCESTRY`: Ancestry population
  - `N_SEGMENTS`: Number of segments
  - `MEAN_LENGTH_SNP`: Mean segment length (SNPs)
  - `MEDIAN_LENGTH_SNP`: Median segment length (SNPs)
  - `MEAN_LENGTH_BP`: Mean segment length (bp)
  - `MEDIAN_LENGTH_BP`: Median segment length (bp)
  - `TOTAL_LENGTH_BP`: Total segment length (bp)
  - `MAX_LENGTH_SNP`: Maximum segment length (SNPs)
  - `MAX_LENGTH_BP`: Maximum segment length (bp)
  - `MIN_LENGTH_SNP`: Minimum segment length (SNPs)
  - `MIN_LENGTH_BP`: Minimum segment length (bp)

### 4. **Summary Report**
- **File:** `<prefix>_summary_report.txt`
- **Format:** Human-readable text
- **Contents:**
  - Global ancestry proportions
  - Segment statistics summary
  - Chromosomal distribution of segments
  - List of generated files

### 5. **Visualization Output (from plot_segment_lengths.py)**
- **File:** `segment_length_distribution.png`
- **Contents:** Four-panel figure showing:
  1. Haplotype segment length distribution (SNPs) - KDE plot
  2. Haplotype segment length distribution (bp) - KDE plot
  3. Segment length comparison (SNPs) - Box plot
  4. Segment length comparison (bp) - Box plot

- **File:** `chromosome_distribution.png`
- **Contents:** Stacked bar chart showing distribution of ancestry segments across chromosomes

## Directory Structure

A typical analysis project might have the following structure:

```
loter_analysis_project/
│
├── data/
│   ├── vcf_files/
│   │   └── 9pops.subset.vcf.gz          # Input VCF file
│   ├── loter_results/
│   │   └── Jino_vs_4refs.npy            # Loter output file
│   └── intermediate/
│       └── snp_positions.tsv            # Extracted SNP positions
│
├── scripts/
│   ├── loter_segment_analyzer.py        # Main analysis script
│   └── plot_segment_lengths.py          # Visualization script
│
├── results/
│   ├── Jino_vs_HDPW_ancestry_proportions.tsv
│   ├── Jino_vs_HDPW_segments.tsv
│   ├── Jino_vs_HDPW_segment_statistics.tsv
│   ├── Jino_vs_HDPW_summary_report.txt
│   └── visualization/
│       ├── segment_length_distribution.png
│       └── chromosome_distribution.png
│
├── README.md                            # This file
├── requirements.txt                     # Python dependencies
└── run_analysis.sh                      # Example analysis script
```

## Example Workflow

### Case Study: Jino Population Ancestry Analysis

#### Step 1: Prepare SNP Position File
```bash
# Extract SNP positions from the reference VCF
bcftools query -f '%CHROM\t%POS\n' \
  /home/litianxing/100My_Jino/110.Relate/data/chrvcf/9pops.subset.vcf.gz \
  > snp_positions.tsv

# Verify the extraction
wc -l snp_positions.tsv
head -5 snp_positions.tsv
```

#### Step 2: Analyze Loter Results
```bash
# Run analysis with four reference populations
python loter_segment_analyzer.py \
  --loter-npy loter_output_Jino_vs_4refs.npy \
  --group-names "Han,Dai,Pumi,Wa" \
  --snp-map-tsv snp_positions.tsv \
  --handle-cross-chr \
  --min-snps 5 \
  --output-prefix Jino_vs_HDPW
```

#### Step 3: Examine Results
```bash
# Check the summary report
cat Jino_vs_HDPW_summary_report.txt

# Look at segment statistics
head Jino_vs_HDPW_segment_statistics.tsv

# Count segments by ancestry
cut -f2 Jino_vs_HDPW_segments.tsv | sort | uniq -c
```

#### Step 4: Create Visualizations
```bash
# Generate distribution plots
python plot_segment_lengths.py Jino_vs_HDPW_segments.tsv

# Open the generated images
open segment_length_distribution.png
open chromosome_distribution.png
```

## Advanced Usage

### Batch Processing Multiple Loter Results

Create a batch processing script `batch_analyze.sh`:

```bash
#!/bin/bash
# batch_analyze.sh

SNP_MAP="snp_positions.tsv"

# Array of analyses to run
declare -A analyses=(
    ["Jino_vs_HDPW"]="Han,Dai,Pumi,Wa"
    ["Jino_vs_HMY"]="Han,Miao,Yao"
    ["Jino_vs_all"]="Han,Dai,Pumi,Wa,Miao,Yao,Tibetan"
)

for prefix in "${!analyses[@]}"; do
    echo "Processing: $prefix"
    echo "Reference groups: ${analyses[$prefix]}"
    
    python loter_segment_analyzer.py \
        --loter-npy "loter_output_${prefix}.npy" \
        --group-names "${analyses[$prefix]}" \
        --snp-map-tsv "$SNP_MAP" \
        --handle-cross-chr \
        --min-snps 5 \
        --output-prefix "results/${prefix}"
    
    echo "Completed: $prefix"
    echo "----------------------------------------"
done
```

### Customizing Visualization

Modify `plot_segment_lengths.py` to customize plots:

```python
# Example: Change color palette
sns.set_palette("husl")  # Different color scheme

# Example: Adjust figure size
plt.figure(figsize=(16, 10))

# Example: Add custom title
ax.set_title("Jino Population Ancestry Segments", fontsize=16, weight='bold')
```

## Troubleshooting

### Common Issues and Solutions

**Issue 1: bcftools not found**
```
Error: bcftools command not found
```
**Solution:** Install bcftools using package manager (apt, yum, brew) or compile from source.

**Issue 2: SNP count mismatch**
```
Warning: SNP count mismatch: Loter=500000, SNP_map=600000
```
**Solution:** Ensure you're using the same VCF file that was used for Loter analysis. Check if VCF was filtered before Loter analysis.

**Issue 3: Memory error with large datasets**
```
MemoryError: Unable to allocate array with shape...
```
**Solution:** Process chromosomes separately or increase system memory. Consider using memory-efficient data types.

**Issue 4: Cross-chromosome segments**
```
Found 15 cross-chromosome segments
```
**Solution:** Use `--handle-cross-chr` flag to automatically split these segments.

### Log File Interpretation

The script provides verbose output when using `--verbose` flag:

```
[INFO] 从VCF文件提取SNP位置信息: 9pops.subset.vcf.gz
[INFO] 成功加载 1234567 个SNP位置
[INFO] 染色体分布:
chr1    234567
chr2    210987
...
[INFO] 处理单倍型 100/2000...
```

## Citation

If you use this tool in your research, please cite:

```
Loter Segment Analyzer v1.0 - A tool for mapping Loter SNP indices to physical positions
Author: Jino Ancestry Inference Project
GitHub: https://github.com/yourusername/loter-segment-analyzer
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

For questions, issues, or feature requests, please:
1. Check the [Issues](https://github.com/yourusername/loter-segment-analyzer/issues) page
2. Create a new issue if needed

---

**Note:** This tool is specifically designed for use with Loter output. Ensure your Loter results are generated using the same VCF file provided to this tool for accurate SNP index mapping.
