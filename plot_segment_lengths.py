#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
可视化单倍型片段长度分布
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys


def plot_length_distribution(segment_file):
    """绘制片段长度分布图"""

    # 加载数据
    df = pd.read_csv(segment_file, sep='\t')

    # 设置绘图样式
    sns.set_style("whitegrid")
    sns.set_context("talk")

    # 创建图形
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    axes = axes.flatten()

    # 1. 各祖源片段长度分布（SNP数量）
    ax1 = axes[0]
    for ancestry in df['ANCESTRY'].unique():
        anc_data = df[df['ANCESTRY'] == ancestry]['LENGTH_SNP']
        # 使用KDE绘制
        sns.kdeplot(data=anc_data, label=ancestry, ax=ax1, fill=True, alpha=0.3)

    ax1.set_xlabel('Segment Length (SNPs)')
    ax1.set_ylabel('Density')
    ax1.set_title('Haplotype Segment Length Distribution (SNPs)')
    ax1.legend()
    ax1.set_xlim(0, df['LENGTH_SNP'].quantile(0.95))  # 去除极端值

    # 2. 各祖源片段长度分布（物理长度）
    ax2 = axes[1]
    for ancestry in df['ANCESTRY'].unique():
        anc_data = df[df['ANCESTRY'] == ancestry]['LENGTH_BP']
        # 过滤极端值
        q95 = anc_data.quantile(0.95)
        anc_data_filtered = anc_data[anc_data <= q95]
        sns.kdeplot(data=anc_data_filtered, label=ancestry, ax=ax2, fill=True, alpha=0.3)

    ax2.set_xlabel('Segment Length (bp)')
    ax2.set_ylabel('Density')
    ax2.set_title('Haplotype Segment Length Distribution (bp)')
    ax2.legend()

    # 3. 箱线图比较（SNP长度）
    ax3 = axes[2]
    sns.boxplot(data=df, x='ANCESTRY', y='LENGTH_SNP', ax=ax3)
    ax3.set_ylabel('Segment Length (SNPs)')
    ax3.set_title('Segment Length Comparison (SNPs)')
    ax3.tick_params(axis='x', rotation=45)

    # 4. 箱线图比较（物理长度）
    ax4 = axes[3]
    # 过滤极端值
    df_filtered = df.copy()
    for ancestry in df['ANCESTRY'].unique():
        anc_mask = df['ANCESTRY'] == ancestry
        q95 = df.loc[anc_mask, 'LENGTH_BP'].quantile(0.95)
        df_filtered.loc[anc_mask & (df['LENGTH_BP'] > q95), 'LENGTH_BP'] = q95

    sns.boxplot(data=df_filtered, x='ANCESTRY', y='LENGTH_BP', ax=ax4)
    ax4.set_ylabel('Segment Length (bp)')
    ax4.set_title('Segment Length Comparison (bp)')
    ax4.tick_params(axis='x', rotation=45)

    plt.tight_layout()
    plt.savefig('segment_length_distribution.png', dpi=300, bbox_inches='tight')
    plt.show()

    print(f"图形已保存为 'segment_length_distribution.png'")


def plot_chromosome_distribution(segment_file):
    """绘制染色体分布图"""

    df = pd.read_csv(segment_file, sep='\t')

    # 按染色体和祖源统计
    chr_anc_counts = df.groupby(['CHR_START', 'ANCESTRY']).size().unstack(fill_value=0)

    # 绘制堆叠条形图
    plt.figure(figsize=(12, 6))
    chr_anc_counts.plot(kind='bar', stacked=True)
    plt.xlabel('Chromosome')
    plt.ylabel('Number of Segments')
    plt.title('Distribution of Ancestry Segments Across Chromosomes')
    plt.legend(title='Ancestry')
    plt.tight_layout()
    plt.savefig('chromosome_distribution.png', dpi=300, bbox_inches='tight')
    plt.show()


def main():
    if len(sys.argv) < 2:
        print("用法: python plot_segment_lengths.py <segment_file.tsv>")
        sys.exit(1)

    segment_file = sys.argv[1]

    print(f"加载数据: {segment_file}")
    df = pd.read_csv(segment_file, sep='\t')
    print(f"数据形状: {df.shape}")
    print("\n数据预览:")
    print(df.head())

    # 绘制图形
    plot_length_distribution(segment_file)
    plot_chromosome_distribution(segment_file)


if __name__ == "__main__":
    main()