# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Functions for getting MSA and calculating alignment features."""

from collections.abc import MutableMapping, Sequence
import dataclasses
import string
from typing import Self, Union, Optional

from absl import logging
from alphafold3.constants import mmcif_names
from alphafold3.data import msa_config
from alphafold3.data import msa_features
from alphafold3.data import msa_identifiers
from alphafold3.data.tools import jackhmmer
from alphafold3.data.tools import mmseqs2
from alphafold3.data.tools import nhmmer
import numpy as np


@dataclasses.dataclass(frozen=True, slots=True)
class Msa:
  """Class representing an MSA."""
  sequences: list[str]
  deletion_matrix: list[list[int]]
  descriptions: list[str]

  @classmethod
  def from_a3m(
      cls,
      query_sequence: str,
      chain_poly_type: str,
      a3m: str,
      max_depth: Optional[int] = None,
      deduplicate: bool = False,
  ) -> Self:
    """Creates an MSA from an a3m string."""
    sequences = []
    deletion_matrix = []
    descriptions = []
    
    for line in a3m.splitlines():
      if not line:
        continue
      if line.startswith('>'):
        descriptions.append(line[1:].strip())
      else:
        sequence = ''
        deletion_row = []
        for char in line:
          if char.isupper():
            sequence += char
            deletion_row.append(0)
          elif char.islower():
            sequence += char.upper()
            deletion_row.append(1)
        sequences.append(sequence)
        deletion_matrix.append(deletion_row)
    
    if deduplicate:
      unique_sequences = {}
      for i, sequence in enumerate(sequences):
        if sequence not in unique_sequences:
          unique_sequences[sequence] = i
      
      sequences = [sequences[i] for i in unique_sequences.values()]
      deletion_matrix = [deletion_matrix[i] for i in unique_sequences.values()]
      descriptions = [descriptions[i] for i in unique_sequences.values()]
    
    if max_depth is not None:
      sequences = sequences[:max_depth]
      deletion_matrix = deletion_matrix[:max_depth]
      descriptions = descriptions[:max_depth]
    
    return cls(sequences=sequences, deletion_matrix=deletion_matrix, descriptions=descriptions)

  def to_a3m(self) -> str:
    """Converts the MSA to an a3m string."""
    a3m_lines = []
    for sequence, deletion_row, description in zip(
        self.sequences, self.deletion_matrix, self.descriptions
    ):
      a3m_lines.append(f">{description}")
      a3m_sequence = ""
      for aa, deletion in zip(sequence, deletion_row):
        if deletion:
          a3m_sequence += aa.lower()
        else:
          a3m_sequence += aa
      a3m_lines.append(a3m_sequence)
    return "\n".join(a3m_lines)

  @classmethod
  def from_empty(cls, query_sequence: str, chain_poly_type: str) -> Self:
    """Creates an empty MSA containing just the query sequence."""
    return cls(
        sequences=[query_sequence],
        deletion_matrix=[[0] * len(query_sequence)],
        descriptions=["Original query"],
    )

  def featurize(self) -> MutableMapping[str, np.ndarray]:
    """Featurizes the MSA and returns a map of feature names to features."""
    msa_array, deletion_matrix = msa_features.extract_msa_features(
        msa_sequences=self.sequences,
        chain_poly_type=mmcif_names.PROTEIN_CHAIN,
    )
    species_ids = msa_features.extract_species_ids(self.descriptions)
    
    return {
        'msa': msa_array,
        'deletion_matrix_int': deletion_matrix,
        'msa_species_identifiers': np.array(species_ids, dtype=object),
        'num_alignments': np.array(len(self.sequences), dtype=np.int32),
    }


def merge_msas(msas: Sequence[Msa], deduplicate: bool = False) -> Msa:
  """Merges multiple MSAs into a single MSA."""
  all_sequences = []
  all_deletion_matrix = []
  all_descriptions = []
  
  for msa in msas:
    all_sequences.extend(msa.sequences)
    all_deletion_matrix.extend(msa.deletion_matrix)
    all_descriptions.extend(msa.descriptions)
  
  if deduplicate:
    unique_sequences = {}
    for i, sequence in enumerate(all_sequences):
      if sequence not in unique_sequences:
        unique_sequences[sequence] = i
    
    all_sequences = [all_sequences[i] for i in unique_sequences.values()]
    all_deletion_matrix = [all_deletion_matrix[i] for i in unique_sequences.values()]
    all_descriptions = [all_descriptions[i] for i in unique_sequences.values()]
  
  return Msa(
      sequences=all_sequences,
      deletion_matrix=all_deletion_matrix,
      descriptions=all_descriptions,
  )


def get_msa_tool(config: Union[msa_config.MMseqs2Config, msa_config.JackhmmerConfig, msa_config.NhmmerConfig]) -> Union[mmseqs2.MMseqs2, jackhmmer.Jackhmmer, nhmmer.Nhmmer]:
  """Returns the appropriate MSA tool based on the config type."""
  if isinstance(config, msa_config.MMseqs2Config):
    return mmseqs2.MMseqs2(
        binary_path=config.binary_path,
        database_path=config.database_path,
        n_cpu=config.n_cpu,
        e_value=config.e_value,
        max_sequences=config.max_sequences,
        sensitivity=config.sensitivity,
    )
  elif isinstance(config, msa_config.JackhmmerConfig):
    return jackhmmer.Jackhmmer(
        binary_path=config.binary_path,
        database_path=config.database_path,
        n_cpu=config.n_cpu,
        e_value=config.e_value,
        max_sequences=config.max_sequences,
    )
  elif isinstance(config, msa_config.NhmmerConfig):
    return nhmmer.Nhmmer(
        binary_path=config.binary_path,
        database_path=config.database_config.path,
        n_cpu=config.n_cpu,
        e_value=config.e_value,
        max_sequences=config.max_sequences,
        alphabet=config.alphabet,
    )
  else:
    raise ValueError(f"Unsupported MSA tool config type: {type(config)}")


def get_msa(
    target_sequence: str,
    run_config: msa_config.RunConfig,
    deduplicate: bool = False,
) -> Msa:
  """Computes the MSA for a given query sequence."""
  return Msa.from_a3m(
      query_sequence=target_sequence,
      chain_poly_type=run_config.chain_poly_type,
      a3m=get_msa_tool(run_config.config).query(target_sequence).a3m,
      max_depth=run_config.crop_size,
      deduplicate=deduplicate,
  )


def get_protein_msa_config(data_pipeline_config, database_path: str) -> msa_config.RunConfig:
  """Creates the appropriate MSA configuration based on the selected tool."""
  if data_pipeline_config.msa_tool == "mmseqs2":
    config = msa_config.MMseqs2Config(
        binary_path=data_pipeline_config.mmseqs2_binary_path,
        database_path=database_path,
        n_cpu=data_pipeline_config.mmseqs2_n_cpu,
        e_value=0.0001,
        max_sequences=10000,
        sensitivity=7.5,
        gpu_devices=data_pipeline_config.mmseqs2_gpu_devices,  # Pass GPU devices parameter
    )
  else:  # jackhmmer
    config = msa_config.JackhmmerConfig(
        binary_path=data_pipeline_config.jackhmmer_binary_path,
        database_path=database_path,
        n_cpu=data_pipeline_config.jackhmmer_n_cpu,
        e_value=0.0001,
        max_sequences=10000,
    )
  
  return msa_config.RunConfig(
      config=config,
      chain_poly_type=mmcif_names.PROTEIN_CHAIN,
      crop_size=None,
  )
