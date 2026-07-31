#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import argparse
import sys
from pathlib import Path

def get_category(presence_count, n_species):
    """标准 pan-genome 分类"""
    if n_species <= 0:
        raise ValueError("n_species must be > 0")
    
    perc = (presence_count / n_species) * 100
    if presence_count == n_species:
        return "core"
    elif perc > 90:
        return "softcore"
    elif perc >= 10:
        return "shell"
    else:
        return "cloud"


def main():
    parser = argparse.ArgumentParser(description="Orthogroups pan-genome classification (core/softcore/shell/cloud)")
    parser.add_argument("-c", "--count", required=True, help="Orthogroups.GeneCount.tsv")
    parser.add_argument("-t", "--tsv", required=True, help="Orthogroups.tsv")
    parser.add_argument("-o", "--output", required=True, help="Output prefix")
    args = parser.parse_args()

    count_file = Path(args.count)
    og_file = Path(args.tsv)
    prefix = args.output

    # 文件检查
    for f in [count_file, og_file]:
        if not f.exists():
            print(f"ERROR: File not found: {f}", file=sys.stderr)
            sys.exit(1)

    print("Loading files...")
    count_df = pd.read_csv(count_file, sep="\t", low_memory=False)
    og_df = pd.read_csv(og_file, sep="\t", low_memory=False)

    # 获取物种列
    species_cols = [col for col in count_df.columns if col not in ["Orthogroup", "Total"]]
    n_species = len(species_cols)
    print(f"Detected {n_species} genomes: {species_cols}")

    # 数值转换 + 分类
    count_df[species_cols] = count_df[species_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    presence_matrix = count_df[species_cols] > 0
    count_df["PresenceCount"] = presence_matrix.sum(axis=1)
    count_df["Category"] = count_df["PresenceCount"].apply(lambda x: get_category(x, n_species))

    # 输出1：OG 分类结果
    count_df[["Orthogroup", "PresenceCount", "Category"]].to_csv(
        f"{prefix}.orthogroup_category.tsv", sep="\t", index=False
    )

    # 统计各类别数量
    print("\nOrthogroup category counts:")
    for cat in ["core", "softcore", "shell", "cloud"]:
        num = (count_df["Category"] == cat).sum()
        print(f"  {cat:10}: {num:6d}  ({num/n_species*100:5.2f}%)")

    # ====================== 基因映射 ======================
    print("Mapping genes to categories for each genome...")

    # 构建 Orthogroup -> Category 的字典（加速查找）
    og_to_cat = dict(zip(count_df["Orthogroup"], count_df["Category"]))

    summary_list = []
    detail_list = []

    for genome in species_cols:
        cat_data = {"core": [], "softcore": [], "shell": [], "cloud": []}

        for _, row in og_df[["Orthogroup", genome]].iterrows():
            og_id = row["Orthogroup"]
            genes_str = str(row[genome]).strip() if pd.notna(row[genome]) else ""
            
            if not genes_str or genes_str == "NaN":
                continue
                
            gene_list = [g.strip() for g in genes_str.split(",") if g.strip()]
            category = og_to_cat.get(og_id)
            
            if category:
                cat_data[category].extend(gene_list)

        # 汇总
        for cat in ["core", "softcore", "shell", "cloud"]:
            genes = cat_data[cat]
            gene_count = len(genes)
            summary_list.append([genome, cat, gene_count])
            
            # 只在 detail 中保存基因列表（可选：如果基因太多可注释掉）
            gene_str = ",".join(genes) if genes else ""
            detail_list.append([genome, cat, gene_count, gene_str])

    # 输出2：统计表
    pd.DataFrame(summary_list, columns=["Genome", "Category", "GeneCount"]).to_csv(
        f"{prefix}.summary.tsv", sep="\t", index=False
    )

    # 输出3：详细基因列表
    pd.DataFrame(detail_list, columns=["Genome", "Category", "GeneCount", "GeneIDs"]).to_csv(
        f"{prefix}.detail.tsv", sep="\t", index=False
    )

    print("\n✅ Analysis complete!")
    print(f"Output files:")
    print(f"   1. {prefix}.orthogroup_category.tsv")
    print(f"   2. {prefix}.summary.tsv")
    print(f"   3. {prefix}.detail.tsv")


if __name__ == "__main__":
    main()
