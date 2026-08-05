#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""

Loter单倍型片段分析工具 - 优化版本

将Loter SNP索引映射回VCF物理位置，提取祖源片段信息

作者: your-username

"""



import numpy as np

import pandas as pd

import argparse

import sys

from pathlib import Path

import warnings

import time

from typing import Dict, List, Tuple, Optional



warnings.filterwarnings('ignore')





class LoterSegmentAnalyzer:

    """Loter结果分析器，支持SNP索引到物理位置的映射 - 优化版本"""



    def __init__(self, verbose=True):

        self.verbose = verbose

        self.snp_map = None  # 存储SNP索引到(CHR, POS)的映射

        self.chr_array = None  # 预处理的染色体数组

        self.pos_array = None  # 预处理的位置数组

        self.group_names = None



    def load_snp_map_from_vcf(self, vcf_file: str) -> None:

        """

        从VCF文件中提取SNP位置信息 - 优化版本



        参数:

        -----------

        vcf_file : str

            VCF文件路径（支持.gz压缩）

        """

        print(f"[INFO] 从VCF文件提取SNP位置信息: {vcf_file}")

        start_time = time.time()



        import subprocess



        # 构建bcftools命令

        cmd = f"bcftools query -f '%CHROM\\t%POS\\n' {vcf_file}"



        try:

            # 执行命令并读取输出

            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

            if result.returncode != 0:

                print(f"[ERROR] bcftools命令执行失败: {result.stderr}")

                sys.exit(1)



            # 使用pandas快速读取

            from io import StringIO

            self.snp_map = pd.read_csv(

                StringIO(result.stdout),

                sep='\t',

                names=['CHR', 'POS'],

                dtype={'CHR': str, 'POS': np.int64}

            )

            self.snp_map['SNP_INDEX'] = np.arange(len(self.snp_map))



        except FileNotFoundError:

            print("[ERROR] bcftools未安装，请安装bcftools或使用备用方法")

            sys.exit(1)



        # 预处理数组

        self._preprocess_snp_arrays()

        

        elapsed = time.time() - start_time

        if self.verbose:

            print(f"[INFO] 成功加载 {len(self.snp_map)} 个SNP位置 (耗时: {elapsed:.2f}s)")

            print(f"[INFO] 染色体分布:")

            print(self.snp_map['CHR'].value_counts().to_string())



    def load_snp_map_from_tsv(self, tsv_file: str) -> None:

        """

        从TSV文件加载SNP映射（如果已经提前提取）



        参数:

        -----------

        tsv_file : str

            TSV文件，格式：CHROM\\tPOS

        """

        print(f"[INFO] 从TSV文件加载SNP映射: {tsv_file}")

        start_time = time.time()



        self.snp_map = pd.read_csv(

            tsv_file,

            sep='\t',

            names=['CHR', 'POS'],

            dtype={'CHR': str, 'POS': np.int64}

        )

        self.snp_map['SNP_INDEX'] = np.arange(len(self.snp_map))

        

        # 预处理数组

        self._preprocess_snp_arrays()

        

        elapsed = time.time() - start_time

        if self.verbose:

            print(f"[INFO] 成功加载 {len(self.snp_map)} 个SNP位置 (耗时: {elapsed:.2f}s)")

            print(f"[INFO] 前5个SNP位置:")

            print(self.snp_map.head())



    def _preprocess_snp_arrays(self) -> None:

        """预处理SNP数组以提高性能"""

        # 转换为numpy数组以提高访问速度

        self.chr_array = self.snp_map['CHR'].values

        self.pos_array = self.snp_map['POS'].values

        

        # 预计算染色体边界（为了更快的跨染色体检测）

        self.chr_boundaries = np.where(self.chr_array[1:] != self.chr_array[:-1])[0] + 1

        self.chr_boundaries = np.concatenate([[0], self.chr_boundaries, [len(self.chr_array)]])



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

        start_time = time.time()



        res = np.load(npy_file, mmap_mode='r')



        if res.ndim != 2:

            raise ValueError(f"Loter矩阵应为二维数组，实际维度: {res.ndim}")



        n_haps, n_snps = res.shape

        

        elapsed = time.time() - start_time

        print(f"[INFO] 单倍型数量: {n_haps}, SNP数量: {n_snps} (耗时: {elapsed:.2f}s)")



        # 检查SNP数量是否匹配

        if self.snp_map is not None and n_snps != len(self.snp_map):

            print(f"[WARNING] SNP数量不匹配: Loter={n_snps}, SNP映射={len(self.snp_map)}")

            # 如果VCF中的SNP更多，我们可以截取前n_snps个

            if len(self.snp_map) > n_snps:

                print(f"[INFO] 截取前{n_snps}个SNP")

                self.snp_map = self.snp_map.iloc[:n_snps].copy()

                self._preprocess_snp_arrays()



        # 如果内存充足，转换为内存数组以提高速度

        if res.flags['C_CONTIGUOUS']:

            return np.asarray(res)

        else:

            return res.copy()



    def calculate_ancestry_proportions(self, loter_mat: np.ndarray) -> pd.Series:

        """

        计算全局祖源比例 - 优化版本



        参数:

        -----------

        loter_mat : np.ndarray

            Loter矩阵



        返回:

        --------

        pd.Series : 每个祖源的比例

        """

        print(f"[INFO] 计算祖源比例...")

        start_time = time.time()



        total_sites = loter_mat.size

        proportions = {}



        for idx, group in enumerate(self.group_names):

            # 使用向量化计算，比循环快得多

            prop = np.sum(loter_mat == idx) / total_sites

            proportions[group] = prop



        prop_series = pd.Series(proportions)

        prop_series['total_sites'] = total_sites

        prop_series['n_haplotypes'] = loter_mat.shape[0]



        elapsed = time.time() - start_time

        if self.verbose:

            print(f"[INFO] 全局祖源比例 (耗时: {elapsed:.2f}s):")

            for group, prop in proportions.items():

                print(f"  {group}: {prop:.4f} ({prop * 100:.2f}%)")



        return prop_series



    def extract_segments_with_positions(self, loter_mat: np.ndarray,

                                        min_snps: int = 3) -> pd.DataFrame:

        """

        提取单倍型片段并添加物理位置信息 - 优化版本



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

        

        print(f"[INFO] 开始提取片段信息 (min_snps={min_snps})...")

        start_time = time.time()



        records = []

        

        # 预处理染色体边界检测

        # 创建染色体变化掩码（比逐点比较快得多）

        chr_change_mask = np.zeros(n_snps - 1, dtype=bool)

        chr_change_mask[self.chr_boundaries[1:-1] - 1] = True

        

        # 位置回退掩码（检测位置是否减少）

        pos_decrease_mask = self.pos_array[1:] <= self.pos_array[:-1]



        for hap_id in range(n_haps):

            if self.verbose and hap_id % 500 == 0 and hap_id > 0:

                elapsed = time.time() - start_time

                rate = hap_id / elapsed if elapsed > 0 else 0

                print(f"[INFO] 处理单倍型 {hap_id}/{n_haps}... (速度: {rate:.1f} haps/s)")



            hap = loter_mat[hap_id, :]

            

            # 使用向量化操作找到所有变化点

            # 1. 祖源变化点

            ancestry_changes = np.where(hap[:-1] != hap[1:])[0] + 1

            

            # 2. 染色体变化点或位置回退点

            chrom_or_pos_changes = np.where(chr_change_mask | pos_decrease_mask)[0] + 1

            

            # 合并所有变化点

            all_changes = np.unique(np.concatenate([ancestry_changes, chrom_or_pos_changes]))

            all_changes.sort()

            

            # 如果没有任何变化，整个单倍型就是一个片段

            if len(all_changes) == 0:

                if n_snps >= min_snps:

                    records.append(self._create_segment_record_fast(

                        hap_id, hap[0], 0, n_snps - 1

                    ))

                continue

            

            # 处理第一个片段（从0到第一个变化点）

            start_idx = 0

            end_idx = all_changes[0] - 1

            if end_idx - start_idx + 1 >= min_snps:

                records.append(self._create_segment_record_fast(

                    hap_id, hap[start_idx], start_idx, end_idx

                ))

            

            # 处理中间的片段

            for i in range(len(all_changes) - 1):

                start_idx = all_changes[i]

                end_idx = all_changes[i + 1] - 1

                segment_length = end_idx - start_idx + 1

                

                if segment_length >= min_snps:

                    records.append(self._create_segment_record_fast(

                        hap_id, hap[start_idx], start_idx, end_idx

                    ))

            

            # 处理最后一个片段

            start_idx = all_changes[-1]

            end_idx = n_snps - 1

            if end_idx - start_idx + 1 >= min_snps:

                records.append(self._create_segment_record_fast(

                    hap_id, hap[start_idx], start_idx, end_idx

                ))



        # 转换为DataFrame

        segment_df = pd.DataFrame(records)



        elapsed = time.time() - start_time

        if self.verbose:

            print(f"[INFO] 提取了 {len(segment_df)} 个片段 (耗时: {elapsed:.2f}s)")

            print(f"[INFO] 各祖源片段数量:")

            print(segment_df['ANCESTRY'].value_counts().to_string())



        return segment_df



    def _create_segment_record_fast(self, hap_id: int, ancestry_idx: int,

                                    start_idx: int, end_idx: int) -> dict:

        """创建单个片段的记录 - 快速版本"""

        # 直接从预处理的数组中获取信息，避免loc操作

        chr_start = self.chr_array[start_idx]

        pos_start = self.pos_array[start_idx]

        chr_end = self.chr_array[end_idx]

        pos_end = self.pos_array[end_idx]



        # 验证片段是否在同一染色体上

        if chr_start != chr_end:

            # 这不应该发生，因为我们在提取时已经处理了跨染色体

            print(f"[WARNING] 跨染色体片段: hap={hap_id}, {chr_start}:{pos_start}-{chr_end}:{pos_end}")

            # 可以在这里进一步处理或记录日志



        length_bp = pos_end - pos_start

        

        # 确保长度非负

        if length_bp < 0:

            print(f"[ERROR] 负长度片段: hap={hap_id}, {chr_start}:{pos_start}-{pos_end}")

            length_bp = 0  # 设置为0或考虑其他处理方式



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

            'LENGTH_BP': length_bp,

        }



    def calculate_segment_statistics(self, segment_df: pd.DataFrame) -> pd.DataFrame:

        """

        计算片段的统计信息 - 优化版本



        参数:

        -----------

        segment_df : pd.DataFrame

            片段数据框



        返回:

        --------

        pd.DataFrame : 统计信息

        """

        print(f"[INFO] 计算片段统计信息...")

        start_time = time.time()



        stats = []

        

        # 使用groupby和agg进行向量化计算

        grouped = segment_df.groupby('ANCESTRY')

        

        for ancestry, anc_df in grouped:

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



        elapsed = time.time() - start_time

        if self.verbose:

            print(f"\n[INFO] 片段统计信息 (耗时: {elapsed:.2f}s):")

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

        print(f"[INFO] 保存结果...")

        start_time = time.time()



        # 保存祖源比例

        prop_df = pd.DataFrame({

            'ANCESTRY': proportions.index,

            'PROPORTION': proportions.values

        })

        prop_file = f"{output_prefix}_ancestry_proportions.tsv"

        prop_df.to_csv(prop_file, sep='\t', index=False, float_format='%.6f')

        print(f"[INFO] 祖源比例保存到: {prop_file}")



        # 保存片段信息

        segment_file = f"{output_prefix}_segments.tsv"

        # 只保存关键列，减少文件大小

        segment_df[['HAPLOTYPE_ID', 'ANCESTRY', 'CHR_START', 'POS_START', 'POS_END', 'LENGTH_BP', 'LENGTH_SNP']].to_csv(

            segment_file, sep='\t', index=False

        )

        print(f"[INFO] 片段信息保存到: {segment_file}")



        # 保存完整片段信息（如果需要）

        full_segment_file = f"{output_prefix}_segments_full.tsv"

        segment_df.to_csv(full_segment_file, sep='\t', index=False)

        print(f"[INFO] 完整片段信息保存到: {full_segment_file}")



        # 保存统计信息

        stats_file = f"{output_prefix}_segment_statistics.tsv"

        stats_df.to_csv(stats_file, sep='\t', index=False, float_format='%.2f')

        print(f"[INFO] 统计信息保存到: {stats_file}")



        # 保存摘要报告

        self._save_summary_report(proportions, segment_df, stats_df, output_prefix)



        elapsed = time.time() - start_time

        print(f"[INFO] 结果保存完成 (耗时: {elapsed:.2f}s)")



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

            f.write(f"完整片段信息: {output_prefix}_segments_full.tsv\n")

            f.write(f"统计信息: {output_prefix}_segment_statistics.tsv\n")

            f.write(f"摘要报告: {output_prefix}_summary_report.txt\n")



        print(f"[INFO] 摘要报告保存到: {report_file}")





def parse_arguments():

    """解析命令行参数"""

    parser = argparse.ArgumentParser(

        description='Loter单倍型片段分析工具 - 将SNP索引映射到物理位置（优化版本）'

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



    parser.add_argument('--verbose', action='store_true',

                        help='显示详细输出信息')



    parser.add_argument('--chunk-size', type=int, default=100,

                        help='批量处理单倍型的大小（默认: 100）')



    return parser.parse_args()





def main():

    """主函数"""

    args = parse_arguments()

    

    total_start_time = time.time()



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



    # 计算统计信息

    stats_df = analyzer.calculate_segment_statistics(segment_df)



    # 保存结果

    analyzer.save_results(proportions, segment_df, stats_df, args.output_prefix)



    total_elapsed = time.time() - total_start_time

    print("\n" + "=" * 60)

    print(f"分析完成！总耗时: {total_elapsed:.2f} 秒")

    print("=" * 60)





if __name__ == "__main__":

    main()



