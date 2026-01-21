#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Loter单倍型片段分析工具
将Loter SNP索引映射回VCF物理位置，提取祖源片段信息
作者: 基诺族祖源推断项目
"""

import numpy as np
import pandas as pd
import argparse
import sys
from pathlib import Path
import warnings

warnings.filterwarnings('ignore')


class LoterSegmentAnalyzer:
    """Loter结果分析器，支持SNP索引到物理位置的映射"""

    def __init__(self, verbose=True):
        self.verbose = verbose
        self.snp_map = None  # 存储SNP索引到(CHR, POS)的映射
        self.group_names = None

    def load_snp_map_from_vcf(self, vcf_file: str) -> None:
        """
        从VCF文件中提取SNP位置信息

        参数:
        -----------
        vcf_file : str
            VCF文件路径（支持.gz压缩）
        """
        print(f"[INFO] 从VCF文件提取SNP位置信息: {vcf_file}")

        # 使用bcftools提取CHR和POS信息
        import subprocess

        # 构建bcftools命令
        cmd = f"bcftools query -f '%CHROM\\t%POS\\n' {vcf_file}"

        try:
            # 执行命令并读取输出
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"[ERROR] bcftools命令执行失败: {result.stderr}")
                sys.exit(1)

            # 解析输出
            lines = result.stdout.strip().split('\n')
            snp_data = []

            for i, line in enumerate(lines):
                if not line:
                    continue
                chrom, pos = line.split('\t')
                snp_data.append({
                    'SNP_INDEX': i,
                    'CHR': chrom,
                    'POS': int(pos)
                })

            self.snp_map = pd.DataFrame(snp_data)

            if self.verbose:
                print(f"[INFO] 成功加载 {len(self.snp_map)} 个SNP位置")
                print(f"[INFO] 染色体分布:")
                print(self.snp_map['CHR'].value_counts().to_string())

        except FileNotFoundError:
            print("[ERROR] bcftools未安装，请安装bcftools或使用备用方法")
            sys.exit(1)

    def load_snp_map_from_tsv(self, tsv_file: str) -> None:
        """
        从TSV文件加载SNP映射（如果已经提前提取）

        参数:
        -----------
        tsv_file : str
            TSV文件，格式：CHROM\\tPOS
        """
        print(f"[INFO] 从TSV文件加载SNP映射: {tsv_file}")

        self.snp_map = pd.read_csv(
            tsv_file,
            sep='\t',
            names=['CHR', 'POS']
        )
        self.snp_map['SNP_INDEX'] = range(len(self.snp_map))

        if self.verbose:
            print(f"[INFO] 成功加载 {len(self.snp_map)} 个SNP位置")
            print(f"[INFO] 前5个SNP位置:")
            print(self.snp_map.head())

    def load_loter_result(self, npy_file: str) -> np.ndarray:
        """
        加载Loter输出结果

        参数:
        -----------
        npy_file : str
            Loter输出的.npy文件路径

        返回:
        --------
        np.ndarray : Loter矩阵，形状为 (n_haps, n_snps)
        """
        print(f"[INFO] 加载Loter结果: {npy_file}")

        res = np.load(npy_file)

        if res.ndim != 2:
            raise ValueError(f"Loter矩阵应为二维数组，实际维度: {res.ndim}")

        n_haps, n_snps = res.shape
        print(f"[INFO] 单倍型数量: {n_haps}, SNP数量: {n_snps}")

        # 检查SNP数量是否匹配
        if self.snp_map is not None and n_snps != len(self.snp_map):
            print(f"[WARNING] SNP数量不匹配: Loter={n_snps}, SNP映射={len(self.snp_map)}")
            # 如果VCF中的SNP更多，我们可以截取前n_snps个
            if len(self.snp_map) > n_snps:
                print(f"[INFO] 截取前{n_snps}个SNP")
                self.snp_map = self.snp_map.iloc[:n_snps].copy()

        return res

    def calculate_ancestry_proportions(self, loter_mat: np.ndarray) -> pd.Series:
        """
        计算全局祖源比例

        参数:
        -----------
        loter_mat : np.ndarray
            Loter矩阵

        返回:
        --------
        pd.Series : 每个祖源的比例
        """
        total_sites = loter_mat.size

        proportions = {}
        for idx, group in enumerate(self.group_names):
            prop = np.sum(loter_mat == idx) / total_sites
            proportions[group] = prop

        prop_series = pd.Series(proportions)
        prop_series['total_sites'] = total_sites
        prop_series['n_haplotypes'] = loter_mat.shape[0]

        if self.verbose:
            print("\n[INFO] 全局祖源比例:")
            for group, prop in proportions.items():
                print(f"  {group}: {prop:.4f} ({prop * 100:.2f}%)")

        return prop_series

    def extract_segments_with_positions(self, loter_mat: np.ndarray,
                                        min_snps: int = 3) -> pd.DataFrame:
        """
        提取单倍型片段并添加物理位置信息

        参数:
        -----------
        loter_mat : np.ndarray
            Loter矩阵
        min_snps : int
            最小SNP数量，用于过滤短片段

        返回:
        --------
        pd.DataFrame : 片段信息数据框
        """
        if self.snp_map is None:
            raise ValueError("请先加载SNP位置映射")

        n_haps, n_snps = loter_mat.shape
        records = []

        print(f"[INFO] 开始提取片段信息 (min_snps={min_snps})...")

        for hap_id in range(n_haps):
            if self.verbose and hap_id % 100 == 0:
                print(f"[INFO] 处理单倍型 {hap_id}/{n_haps}...")

            hap = loter_mat[hap_id, :]

            # # 找到祖源变化的位置
            # change_points = np.where(hap[:-1] != hap[1:])[0] + 1
            # start_idx = 0

            # 找到祖源变化的位置（祖源变化 或 染色体边界）
            change_points = []

            for i in range(n_snps - 1):
                # 祖源发生变化
                if hap[i] != hap[i + 1]:
                    change_points.append(i + 1)
                    continue

                # 物理位置回退（新染色体）
                pos_i = self.snp_map.loc[i, 'POS']
                pos_j = self.snp_map.loc[i + 1, 'POS']
                chr_i = self.snp_map.loc[i, 'CHR']
                chr_j = self.snp_map.loc[i + 1, 'CHR']

                if (chr_i != chr_j) or (pos_j < pos_i):
                    change_points.append(i + 1)

            change_points = np.array(change_points)
            start_idx = 0


            for end_idx in change_points:
                ancestry_idx = int(hap[start_idx])
                segment_length_snps = end_idx - start_idx

                # 过滤短片段
                if segment_length_snps >= min_snps:
                    records.append(self._create_segment_record(
                        hap_id, ancestry_idx, start_idx, end_idx - 1
                    ))

                start_idx = end_idx

            # 处理最后一个片段
            ancestry_idx = int(hap[start_idx])
            segment_length_snps = n_snps - start_idx

            if segment_length_snps >= min_snps:
                records.append(self._create_segment_record(
                    hap_id, ancestry_idx, start_idx, n_snps - 1
                ))

        # 转换为DataFrame
        segment_df = pd.DataFrame(records)

        if self.verbose:
            print(f"[INFO] 提取了 {len(segment_df)} 个片段")
            print(f"[INFO] 各祖源片段数量:")
            print(segment_df['ANCESTRY'].value_counts().to_string())

        return segment_df

    def _create_segment_record(self, hap_id: int, ancestry_idx: int,
                               start_idx: int, end_idx: int) -> dict:
        """创建单个片段的记录"""
        # 获取染色体和位置信息
        chr_start = self.snp_map.loc[start_idx, 'CHR']
        pos_start = self.snp_map.loc[start_idx, 'POS']
        chr_end = self.snp_map.loc[end_idx, 'CHR']
        pos_end = self.snp_map.loc[end_idx, 'POS']

        # 检查是否跨染色体
        if chr_start != chr_end:
            # 跨染色体的情况 - 需要分割，这里先标记
            chr_status = 'cross_chromosome'
        else:
            chr_status = chr_start

        return {
            'HAPLOTYPE_ID': hap_id,
            'ANCESTRY': self.group_names[ancestry_idx],
            'START_SNP': start_idx,
            'END_SNP': end_idx,
            'LENGTH_SNP': end_idx - start_idx + 1,
            'CHR_START': chr_start,
            'CHR_END': chr_end,
            'POS_START': pos_start,
            'POS_END': pos_end,
            'LENGTH_BP': pos_end - pos_start,
            'CHR_STATUS': chr_status
        }

    def handle_cross_chromosome_segments(self, segment_df: pd.DataFrame) -> pd.DataFrame:
        """
        处理跨染色体的片段（分割成多个片段）

        参数:
        -----------
        segment_df : pd.DataFrame
            包含跨染色体片段的原始数据框

        返回:
        --------
        pd.DataFrame : 处理后的数据框
        """
        cross_mask = segment_df['CHR_STATUS'] == 'cross_chromosome'
        cross_segments = segment_df[cross_mask].copy()

        if len(cross_segments) == 0:
            print("[INFO] 没有跨染色体的片段")
            return segment_df

        print(f"[INFO] 发现 {len(cross_segments)} 个跨染色体片段，正在处理...")

        new_segments = []

        for _, seg in cross_segments.iterrows():
            start_idx = seg['START_SNP']
            end_idx = seg['END_SNP']
            ancestry = seg['ANCESTRY']
            hap_id = seg['HAPLOTYPE_ID']

            # 找到染色体边界
            current_chr = self.snp_map.loc[start_idx, 'CHR']
            segment_start_idx = start_idx

            for idx in range(start_idx + 1, end_idx + 1):
                chr_at_idx = self.snp_map.loc[idx, 'CHR']

                if chr_at_idx != current_chr:
                    # 染色体边界，保存当前片段
                    new_segments.append(self._create_segment_record(
                        hap_id, self.group_names.index(ancestry),
                        segment_start_idx, idx - 1
                    ))

                    # 开始新片段
                    current_chr = chr_at_idx
                    segment_start_idx = idx

            # 保存最后一个片段
            new_segments.append(self._create_segment_record(
                hap_id, self.group_names.index(ancestry),
                segment_start_idx, end_idx
            ))

        # 创建新的DataFrame（不包含跨染色体片段）
        non_cross_df = segment_df[~cross_mask].copy()
        new_segments_df = pd.DataFrame(new_segments)

        # 合并
        result_df = pd.concat([non_cross_df, new_segments_df], ignore_index=True)

        print(f"[INFO] 跨染色体片段处理完成，片段数量从 {len(segment_df)} 变为 {len(result_df)}")

        return result_df

    def calculate_segment_statistics(self, segment_df: pd.DataFrame) -> pd.DataFrame:
        """
        计算片段的统计信息

        参数:
        -----------
        segment_df : pd.DataFrame
            片段数据框

        返回:
        --------
        pd.DataFrame : 统计信息
        """
        stats = []

        for ancestry in segment_df['ANCESTRY'].unique():
            anc_df = segment_df[segment_df['ANCESTRY'] == ancestry]

            stats.append({
                'ANCESTRY': ancestry,
                'N_SEGMENTS': len(anc_df),
                'MEAN_LENGTH_SNP': anc_df['LENGTH_SNP'].mean(),
                'MEDIAN_LENGTH_SNP': anc_df['LENGTH_SNP'].median(),
                'MEAN_LENGTH_BP': anc_df['LENGTH_BP'].mean(),
                'MEDIAN_LENGTH_BP': anc_df['LENGTH_BP'].median(),
                'TOTAL_LENGTH_BP': anc_df['LENGTH_BP'].sum(),
                'MAX_LENGTH_SNP': anc_df['LENGTH_SNP'].max(),
                'MAX_LENGTH_BP': anc_df['LENGTH_BP'].max(),
                'MIN_LENGTH_SNP': anc_df['LENGTH_SNP'].min(),
                'MIN_LENGTH_BP': anc_df['LENGTH_BP'].min(),
            })

        stats_df = pd.DataFrame(stats)

        if self.verbose:
            print("\n[INFO] 片段统计信息:")
            print(stats_df.to_string())

        return stats_df

    def save_results(self, proportions: pd.Series,
                     segment_df: pd.DataFrame,
                     stats_df: pd.DataFrame,
                     output_prefix: str):
        """
        保存分析结果

        参数:
        -----------
        proportions : pd.Series
            祖源比例
        segment_df : pd.DataFrame
            片段数据框
        stats_df : pd.DataFrame
            统计信息
        output_prefix : str
            输出文件前缀
        """
        # 保存祖源比例
        prop_df = pd.DataFrame({
            'ANCESTRY': proportions.index,
            'PROPORTION': proportions.values
        })
        prop_file = f"{output_prefix}_ancestry_proportions.tsv"
        prop_df.to_csv(prop_file, sep='\t', index=False)
        print(f"[INFO] 祖源比例保存到: {prop_file}")

        # 保存片段信息
        segment_file = f"{output_prefix}_segments.tsv"
        segment_df.to_csv(segment_file, sep='\t', index=False)
        print(f"[INFO] 片段信息保存到: {segment_file}")

        # 保存统计信息
        stats_file = f"{output_prefix}_segment_statistics.tsv"
        stats_df.to_csv(stats_file, sep='\t', index=False)
        print(f"[INFO] 统计信息保存到: {stats_file}")

        # 保存摘要报告
        self._save_summary_report(proportions, segment_df, stats_df, output_prefix)

    def _save_summary_report(self, proportions: pd.Series,
                             segment_df: pd.DataFrame,
                             stats_df: pd.DataFrame,
                             output_prefix: str):
        """保存摘要报告"""
        report_file = f"{output_prefix}_summary_report.txt"

        with open(report_file, 'w') as f:
            f.write("=" * 60 + "\n")
            f.write("LOTER单倍型片段分析报告\n")
            f.write("=" * 60 + "\n\n")

            f.write("1. 全局祖源比例\n")
            f.write("-" * 40 + "\n")
            for idx, (anc, prop) in enumerate(proportions.items()):
                if anc not in ['total_sites', 'n_haplotypes']:
                    f.write(f"{anc:10s}: {prop:.4f} ({prop * 100:.2f}%)\n")

            f.write(f"\n总位点数: {proportions.get('total_sites', 'N/A')}\n")
            f.write(f"单倍型数: {proportions.get('n_haplotypes', 'N/A')}\n")

            f.write("\n2. 片段统计摘要\n")
            f.write("-" * 40 + "\n")
            f.write(stats_df.to_string() + "\n")

            f.write("\n3. 各染色体片段分布\n")
            f.write("-" * 40 + "\n")
            chr_counts = segment_df.groupby(['CHR_START', 'ANCESTRY']).size().unstack(fill_value=0)
            f.write(chr_counts.to_string() + "\n")

            f.write("\n4. 文件列表\n")
            f.write("-" * 40 + "\n")
            f.write(f"祖源比例: {output_prefix}_ancestry_proportions.tsv\n")
            f.write(f"片段信息: {output_prefix}_segments.tsv\n")
            f.write(f"统计信息: {output_prefix}_segment_statistics.tsv\n")
            f.write(f"摘要报告: {output_prefix}_summary_report.txt\n")

        print(f"[INFO] 摘要报告保存到: {report_file}")


def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='Loter单倍型片段分析工具 - 将SNP索引映射到物理位置'
    )

    parser.add_argument('--loter-npy', required=True,
                        help='Loter输出的.npy文件路径')

    parser.add_argument('--group-names', required=True,
                        help='祖源群体名称，用逗号分隔，如: Han,Dai,Pumi,Wa')

    parser.add_argument('--vcf-file',
                        help='VCF文件路径（用于提取SNP位置）')

    parser.add_argument('--snp-map-tsv',
                        help='SNP映射TSV文件（如果已提取）')

    parser.add_argument('--output-prefix', default='loter_analysis',
                        help='输出文件前缀（默认: loter_analysis）')

    parser.add_argument('--min-snps', type=int, default=3,
                        help='最小片段SNP数量（默认: 3）')

    parser.add_argument('--handle-cross-chr', action='store_true',
                        help='处理跨染色体的片段')

    parser.add_argument('--verbose', action='store_true',
                        help='显示详细输出信息')

    return parser.parse_args()


def main():
    """主函数"""
    args = parse_arguments()

    # 初始化分析器
    analyzer = LoterSegmentAnalyzer(verbose=args.verbose)

    # 解析祖源群体名称
    analyzer.group_names = args.group_names.split(',')
    print(f"[INFO] 祖源群体: {analyzer.group_names}")

    # 加载SNP位置映射
    if args.snp_map_tsv:
        analyzer.load_snp_map_from_tsv(args.snp_map_tsv)
    elif args.vcf_file:
        analyzer.load_snp_map_from_vcf(args.vcf_file)
    else:
        print("[ERROR] 请提供--vcf-file或--snp-map-tsv参数")
        sys.exit(1)

    # 加载Loter结果
    loter_mat = analyzer.load_loter_result(args.loter_npy)

    # 计算全局祖源比例
    proportions = analyzer.calculate_ancestry_proportions(loter_mat)

    # 提取片段信息（带物理位置）
    segment_df = analyzer.extract_segments_with_positions(
        loter_mat, min_snps=args.min_snps
    )

    # 处理跨染色体片段（如果指定）
    if args.handle_cross_chr:
        segment_df = analyzer.handle_cross_chromosome_segments(segment_df)

    # 计算统计信息
    stats_df = analyzer.calculate_segment_statistics(segment_df)

    # 保存结果
    analyzer.save_results(proportions, segment_df, stats_df, args.output_prefix)

    print("\n" + "=" * 60)
    print("分析完成！")
    print("=" * 60)


if __name__ == "__main__":

    main()
