#!/usr/bin/env python3

import json
import argparse
import sys
from typing import List, Dict, Union, Optional
from dataclasses import dataclass, field, asdict
from pathlib import Path

@dataclass
class ProteinModification:
    ptmType: str
    ptmPosition: int

@dataclass
class Protein:
    id: Union[str, List[str]]
    sequence: str
    modifications: Optional[List[ProteinModification]] = None
    unpairedMsa: Optional[str] = None
    pairedMsa: Optional[str] = None
    templates: Optional[List[Dict]] = None

@dataclass
class RNA:
    id: Union[str, List[str]]
    sequence: str
    modifications: Optional[List[Dict]] = None
    unpairedMsa: Optional[str] = None
    pairedMsa: Optional[str] = None

@dataclass
class DNA:
    id: Union[str, List[str]]
    sequence: str
    modifications: Optional[List[Dict]] = None

@dataclass
class Ligand:
    id: str
    ccd: Optional[str] = None
    smiles: Optional[str] = None

@dataclass
class AlphaFold3Input:
    name: str
    modelSeeds: List[int]
    sequences: List[Dict]
    bondedAtomPairs: Optional[List[Dict]] = None
    userCCD: Optional[str] = None
    dialect: str = "alphafold3"
    version: int = 1

def read_fasta(fasta_path: str) -> tuple[str, str]:
    """读取FASTA文件，返回序列名和序列。"""
    with open(fasta_path) as f:
        lines = f.readlines()
    
    name = lines[0].strip().lstrip(">")
    sequence = "".join(line.strip() for line in lines[1:])
    return name, sequence

def read_a3m(a3m_path: str) -> str:
    """读取A3M文件，返回MSA内容。"""
    with open(a3m_path) as f:
        return f.read()

def generate_protein_input():
    """生成蛋白质输入。"""
    print("\n=== 蛋白质输入 ===")
    
    # 询问用户是输入序列还是从文件读取
    use_file = input("是否从FASTA文件读取序列？(y/n): ").lower() == 'y'
    
    if use_file:
        fasta_path = input("请输入FASTA文件路径: ")
        name, sequence = read_fasta(fasta_path)
    else:
        sequence = input("请输入氨基酸序列: ")
    
    # 链ID
    chain_id = input("请输入链ID (多个ID用逗号分隔，如 A 或 A,B,C): ")
    chain_ids = [id.strip() for id in chain_id.split(",")]
    id = chain_ids[0] if len(chain_ids) == 1 else chain_ids
    
    # MSA
    use_msa = input("是否提供MSA文件？(y/n): ").lower() == 'y'
    unpairedMsa = None
    if use_msa:
        msa_path = input("请输入MSA文件路径 (A3M格式): ")
        unpairedMsa = read_a3m(msa_path)
    
    # PTM修饰
    use_ptm = input("是否添加PTM修饰？(y/n): ").lower() == 'y'
    modifications = []
    if use_ptm:
        while True:
            ptm_type = input("请输入PTM类型 (如 HY3，输入空行结束): ")
            if not ptm_type:
                break
            ptm_pos = int(input("请输入修饰位点 (1-based): "))
            modifications.append(ProteinModification(ptm_type, ptm_pos))
    
    protein = Protein(
        id=id,
        sequence=sequence,
        modifications=modifications if modifications else None,
        unpairedMsa=unpairedMsa
    )
    
    return {"protein": asdict(protein)}

def generate_rna_input():
    """生成RNA输入。"""
    print("\n=== RNA输入 ===")
    
    # 询问用户是输入序列还是从文件读取
    use_file = input("是否从FASTA文件读取序列？(y/n): ").lower() == 'y'
    
    if use_file:
        fasta_path = input("请输入FASTA文件路径: ")
        name, sequence = read_fasta(fasta_path)
    else:
        sequence = input("请输入RNA序列: ")
    
    # 链ID
    chain_id = input("请输入链ID (多个ID用逗号分隔，如 R 或 R,S,T): ")
    chain_ids = [id.strip() for id in chain_id.split(",")]
    id = chain_ids[0] if len(chain_ids) == 1 else chain_ids
    
    # MSA
    use_msa = input("是否提供MSA文件？(y/n): ").lower() == 'y'
    unpairedMsa = None
    if use_msa:
        msa_path = input("请输入MSA文件路径 (A3M格式): ")
        unpairedMsa = read_a3m(msa_path)
    
    rna = RNA(
        id=id,
        sequence=sequence,
        unpairedMsa=unpairedMsa
    )
    
    return {"rna": asdict(rna)}

def generate_dna_input():
    """生成DNA输入。"""
    print("\n=== DNA输入 ===")
    
    # 询问用户是输入序列还是从文件读取
    use_file = input("是否从FASTA文件读取序列？(y/n): ").lower() == 'y'
    
    if use_file:
        fasta_path = input("请输入FASTA文件路径: ")
        name, sequence = read_fasta(fasta_path)
    else:
        sequence = input("请输入DNA序列: ")
    
    # 链ID
    chain_id = input("请输入链ID (多个ID用逗号分隔，如 D 或 D,E,F): ")
    chain_ids = [id.strip() for id in chain_id.split(",")]
    id = chain_ids[0] if len(chain_ids) == 1 else chain_ids
    
    dna = DNA(
        id=id,
        sequence=sequence
    )
    
    return {"dna": asdict(dna)}

def generate_ligand_input():
    """生成配体输入。"""
    print("\n=== 配体输入 ===")
    
    ligand_id = input("请输入配体ID: ")
    
    # 询问用户使用CCD还是SMILES
    use_ccd = input("使用CCD代码还是SMILES？(ccd/smiles): ").lower()
    
    if use_ccd == "ccd":
        ccd = input("请输入CCD代码: ")
        ligand = Ligand(id=ligand_id, ccd=ccd)
    else:
        smiles = input("请输入SMILES字符串: ")
        ligand = Ligand(id=ligand_id, smiles=smiles)
    
    return {"ligand": asdict(ligand)}

def main():
    print("=== AlphaFold 3 输入生成器 ===")
    
    # 工作名称
    name = input("请输入工作名称: ")
    
    # 随机种子
    seeds_input = input("请输入随机种子 (多个种子用逗号分隔，如 1,2,3): ")
    modelSeeds = [int(seed.strip()) for seed in seeds_input.split(",")]
    
    # 序列
    sequences = []
    while True:
        print("\n请选择要添加的分子类型：")
        print("1. 蛋白质")
        print("2. RNA")
        print("3. DNA")
        print("4. 配体")
        print("0. 完成")
        
        choice = input("请选择 (0-4): ")
        
        if choice == "0":
            break
        elif choice == "1":
            sequences.append(generate_protein_input())
        elif choice == "2":
            sequences.append(generate_rna_input())
        elif choice == "3":
            sequences.append(generate_dna_input())
        elif choice == "4":
            sequences.append(generate_ligand_input())
    
    # 生成最终的输入对象
    af3_input = AlphaFold3Input(
        name=name,
        modelSeeds=modelSeeds,
        sequences=sequences
    )
    
    # 保存到文件
    output_path = input("\n请输入输出文件路径 (默认: input.json): ").strip() or "input.json"
    with open(output_path, 'w') as f:
        json.dump(asdict(af3_input), f, indent=2)
    
    print(f"\n输入文件已保存到: {output_path}")

if __name__ == "__main__":
    main()
