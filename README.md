# After-Loter---LAI-Results-Analysis-and-Visualizatioon-
Loter Segment Analyzer 是一个用于分析Loter（Local Ancestry Inference Tool）输出结果的Python工具。Loter是一种用于局部祖源推断的算法，其输出为一个二维NumPy数组，表示每个单倍型在每个SNP位点上的祖源赋值。然而，该数组仅使用SNP索引，缺乏物理位置信息。本工具通过将SNP索引映射回VCF文件中的物理位置（染色体和位置），从而能够提取具有实际基因组位置的单倍型片段。此外，工具还提供了片段统计、祖源比例计算以及可视化功能。
