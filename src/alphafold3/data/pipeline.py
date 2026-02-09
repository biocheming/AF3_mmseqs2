# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Functions for running the MSA and template tools for the AlphaFold model."""

from concurrent import futures
import dataclasses
import datetime
import functools
import logging
import sys
import time
from typing import List, Dict, Any

from alphafold3.common import folding_input
from alphafold3.constants import mmcif_names
from alphafold3.data import msa
from alphafold3.data import msa_config
from alphafold3.data import structure_stores
from alphafold3.data import templates
from alphafold3.data.tools import mmseqs2


def _get_parallel_msa_results(sequence: str,
                             msa_configs: List[msa_config.RunConfig]) -> Dict[str, msa.Msa]:
  """Get MSA results in parallel using MMseqs2's parallel method if available."""
  try:
    # Try to use MMseqs2's parallel method for better GPU/CPU coordination
    mmseqs2_tool = None
    db_configs = []

    for config in msa_configs:
      if isinstance(config.config, msa_config.MMseqs2Config):
        if mmseqs2_tool is None:
          # Create MMseqs2 tool instance to access parallel method
          mmseqs2_tool = msa.get_msa_tool(config.config).__class__(
            binary_path=config.config.binary_path,
            database_path=config.config.database_path,
            n_cpu=config.config.n_cpu,
            e_value=config.config.e_value,
            max_sequences=config.config.max_sequences,
            sensitivity=config.config.sensitivity,
            gpu_devices=getattr(config.config, 'gpu_devices', None),
            use_gpu=getattr(config.config, 'use_gpu', False)
          )

        db_name = config.config.database_path.split('/')[-1].split('_')[0]
        db_configs.append({
          'binary_path': config.config.binary_path,
          'database_path': config.config.database_path,
          'n_cpu': config.config.n_cpu,
          'e_value': config.config.e_value,
          'max_sequences': config.config.max_sequences,
          'sensitivity': config.config.sensitivity,
          'gpu_devices': getattr(config.config, 'gpu_devices', None),
          'use_gpu': getattr(config.config, 'use_gpu', False)
        })

    if mmseqs2_tool and hasattr(mmseqs2_tool, 'query_multiple_databases_parallel'):
      logging.info("🚀 Using parallel MSA generation across databases")
      parallel_results = mmseqs2_tool.query_multiple_databases_parallel(sequence, db_configs)

      # Convert to Msa objects
      msa_results = {}
      for db_name, result in parallel_results.items():
        msa_results[db_name] = msa.Msa.from_a3m(
          query_sequence=sequence,
          chain_poly_type=msa_configs[0].chain_poly_type,  # Use first config's poly type
          a3m=result.a3m,
          max_depth=msa_configs[0].crop_size
        )
      return msa_results

  except Exception as e:
    logging.warning(f"Parallel MSA failed, falling back to sequential: {e}")

  # Fallback to sequential processing
  logging.info("⚡ Using sequential MSA generation")
  msa_results = {}
  for config in msa_configs:
    db_name = config.config.database_path.split('/')[-1].split('_')[0]
    msa_results[db_name] = msa.get_msa(sequence, config)

  return msa_results


# Cache to avoid re-running the MSA tools for the same sequence in homomers.
@functools.cache
def _get_protein_msa_and_templates(
    sequence: str,
    uniref90_msa_config: msa_config.RunConfig,
    mgnify_msa_config: msa_config.RunConfig,
    small_bfd_msa_config: msa_config.RunConfig,
    uniprot_msa_config: msa_config.RunConfig,
    templates_config: msa_config.TemplatesConfig,
    pdb_database_path: str,
) -> tuple[msa.Msa, msa.Msa, list[templates.Hit]]:
  """Gets the MSA and templates for a protein sequence."""
  # Try parallel MSA generation first (GPT Scheme 3)
  msa_configs = [uniref90_msa_config, mgnify_msa_config, small_bfd_msa_config, uniprot_msa_config]
  parallel_results = _get_parallel_msa_results(sequence, msa_configs)

  # Extract results - parallel MSA should have all results
  required_dbs = ['uniref90', 'mgnify', 'bfd', 'uniprot']
  missing_dbs = [db for db in required_dbs if db not in parallel_results]

  if missing_dbs:
    logging.warning(f"Parallel MSA missing results for: {missing_dbs}, falling back to sequential")
    raise ValueError(f"Missing MSA results for databases: {missing_dbs}")

  uniref90_msa = parallel_results['uniref90']
  mgnify_msa = parallel_results['mgnify']
  small_bfd_msa = parallel_results['bfd']
  uniprot_msa = parallel_results['uniprot']

  # Now run template search using the uniref90 MSA
  # If uniref90_msa is empty, create a minimal MSA with just the query sequence
  if not uniref90_msa.sequences:
    logging.warning('UniRef90 MSA is empty, using query sequence only for template search')
    template_search_msa = msa.Msa(
        sequences=[sequence],
        deletion_matrix=[[0] * len(sequence)],
        descriptions=['query'],
    )
  else:
    template_search_msa = uniref90_msa

  templates_obj = templates.Templates.from_seq_and_a3m(
      query_sequence=sequence,
      msa_a3m=template_search_msa.to_a3m(),
      max_template_date=templates_config.filter_config.max_template_date,
      database_path=templates_config.template_tool_config.database_path,
      hmmsearch_config=templates_config.template_tool_config.hmmsearch_config,
      max_a3m_query_sequences=None,  # Use all sequences
      structure_store=structure_stores.StructureStore(pdb_database_path),
  )

  # Filter templates according to the config
  templates_obj = templates_obj.filter(
      max_subsequence_ratio=templates_config.filter_config.max_subsequence_ratio,
      min_align_ratio=templates_config.filter_config.min_align_ratio,
      min_hit_length=templates_config.filter_config.min_hit_length,
      deduplicate_sequences=templates_config.filter_config.deduplicate_sequences,
      max_hits=templates_config.filter_config.max_hits,
  )

  # Get both hits and their structures
  template_hits = templates_obj.get_hits_with_structures()

  # Create unpaired MSA: merge uniref90, small_bfd, mgnify (with deduplication)
  # This matches AlphaFold3's logic: unpaired_msa from multiple sources with deduplication
  unpaired_msa = msa.merge_msas(
      msas=[uniref90_msa, small_bfd_msa, mgnify_msa],
      deduplicate=True,
  )

  # Create paired MSA: just uniprot (without deduplication)
  # This matches AlphaFold3's logic: paired_msa from single source without deduplication
  paired_msa = msa.merge_msas(
      msas=[uniprot_msa],
      deduplicate=False,
  )

  # Ensure unpaired_msa is not empty and contains query sequence
  if not unpaired_msa.sequences:
    logging.warning('Unpaired MSA is empty, using query sequence only')
    unpaired_msa = msa.Msa(
        sequences=[sequence],
        deletion_matrix=[[0] * len(sequence)],
        descriptions=['query'],
    )
  elif unpaired_msa.sequences[0] != sequence:
    logging.warning('Query sequence not found in unpaired MSA, adding it as first sequence')
    unpaired_msa = msa.Msa(
        sequences=[sequence] + unpaired_msa.sequences,
        deletion_matrix=[[0] * len(sequence)] + unpaired_msa.deletion_matrix,
        descriptions=['query'] + unpaired_msa.descriptions,
    )

  # Ensure paired_msa is not empty and contains query sequence
  if not paired_msa.sequences:
    logging.warning('Paired MSA is empty, using query sequence only')
    paired_msa = msa.Msa(
        sequences=[sequence],
        deletion_matrix=[[0] * len(sequence)],
        descriptions=['query'],
    )
  elif paired_msa.sequences[0] != sequence:
    logging.warning('Query sequence not found in paired MSA, adding it as first sequence')
    paired_msa = msa.Msa(
        sequences=[sequence] + paired_msa.sequences,
        deletion_matrix=[[0] * len(sequence)] + paired_msa.deletion_matrix,
        descriptions=['query'] + paired_msa.descriptions,
    )

  # Final debug: ensure both MSAs are valid before returning
  logging.warning(f'Final MSA check - unpaired: {len(unpaired_msa.sequences)} seqs, paired: {len(paired_msa.sequences)} seqs')
  if len(unpaired_msa.sequences) == 0:
    logging.error('CRITICAL: Unpaired MSA is still empty!')
  if len(paired_msa.sequences) == 0:
    logging.error('CRITICAL: Paired MSA is still empty!')

  return unpaired_msa, paired_msa, template_hits


# Cache to avoid re-running the Nhmmer for the same sequence in homomers.
@functools.cache
def _get_rna_msa(
    sequence: str,
    nt_rna_msa_config: msa_config.NhmmerConfig,
    rfam_msa_config: msa_config.NhmmerConfig,
    rnacentral_msa_config: msa_config.NhmmerConfig,
) -> msa.Msa:
  """Processes a single RNA chain."""
  logging.info('Getting RNA MSAs for sequence %s', sequence)
  rna_msa_start_time = time.time()
  # Run various MSA tools in parallel. Use a ThreadPoolExecutor because
  # they're not blocked by the GIL, as they're sub-shelled out.
  with futures.ThreadPoolExecutor() as executor:
    nt_rna_msa_future = executor.submit(
        msa.get_msa,
        target_sequence=sequence,
        run_config=nt_rna_msa_config,
        chain_poly_type=mmcif_names.RNA_CHAIN,
    )
    rfam_msa_future = executor.submit(
        msa.get_msa,
        target_sequence=sequence,
        run_config=rfam_msa_config,
        chain_poly_type=mmcif_names.RNA_CHAIN,
    )
    rnacentral_msa_future = executor.submit(
        msa.get_msa,
        target_sequence=sequence,
        run_config=rnacentral_msa_config,
        chain_poly_type=mmcif_names.RNA_CHAIN,
    )
  nt_rna_msa = nt_rna_msa_future.result()
  rfam_msa = rfam_msa_future.result()
  rnacentral_msa = rnacentral_msa_future.result()
  logging.info(
      'Getting RNA MSAs took %.2f seconds for sequence %s',
      time.time() - rna_msa_start_time,
      sequence,
  )

  return msa.Msa.from_multiple_msas(
      msas=[rfam_msa, rnacentral_msa, nt_rna_msa],
      deduplicate=True,
  )


@dataclasses.dataclass(frozen=True, kw_only=True, slots=True)
class DataPipelineConfig:
  """Configuration for the data pipeline.

  Attributes:
    jackhmmer_binary_path: Path to the jackhmmer binary.
    nhmmer_binary_path: Path to the nhmmer binary.
    hmmalign_binary_path: Path to the hmmalign binary.
    hmmsearch_binary_path: Path to the hmmsearch binary.
    hmmbuild_binary_path: Path to the hmmbuild binary.
    mmseqs2_binary_path: Path to the MMseqs2 binary.
    small_bfd_database_path: Path to the small BFD database.
    mgnify_database_path: Path to the MGnify database.
    uniprot_cluster_annot_database_path: Path to the UniProt cluster annotation database.
    uniref90_database_path: Path to the UniRef90 database.
    ntrna_database_path: Path to the NTRNA database.
    rfam_database_path: Path to the Rfam database.
    rna_central_database_path: Path to the RNAcentral database.
    seqres_database_path: Path to the PDB seqres database.
    pdb_database_path: Path to the PDB database.
    jackhmmer_n_cpu: Number of CPUs to use for jackhmmer.
    nhmmer_n_cpu: Number of CPUs to use for nhmmer.
    mmseqs2_n_cpu: Number of CPUs to use for MMseqs2.
    msa_tool: Tool to use for MSA generation (mmseqs2 or jackhmmer).
    mmseqs2_gpu_devices: List of GPU devices to use for MMseqs2 searches.
  """

  jackhmmer_binary_path: str
  nhmmer_binary_path: str
  hmmalign_binary_path: str
  hmmsearch_binary_path: str
  hmmbuild_binary_path: str
  mmseqs2_binary_path: str
  small_bfd_database_path: str
  mgnify_database_path: str
  uniprot_cluster_annot_database_path: str
  uniref90_database_path: str
  ntrna_database_path: str
  rfam_database_path: str
  rna_central_database_path: str
  seqres_database_path: str
  pdb_database_path: str
  jackhmmer_n_cpu: int = 8
  nhmmer_n_cpu: int = 8
  mmseqs2_n_cpu: int = 8
  msa_tool: str = "mmseqs2"
  mmseqs2_gpu_devices: tuple[str, ...] | None = None
  mmseqs2_use_gpu: bool = False

class DataPipeline:
  """Runs the alignment tools and assembles the input features."""

  def __init__(self, data_pipeline_config: DataPipelineConfig):
    """Initializes the data pipeline with default configurations."""
    self._uniref90_msa_config = msa.get_protein_msa_config(
        data_pipeline_config,
        data_pipeline_config.uniref90_database_path,
    )
    self._mgnify_msa_config = msa.get_protein_msa_config(
        data_pipeline_config,
        data_pipeline_config.mgnify_database_path,
    )
    self._small_bfd_msa_config = msa.get_protein_msa_config(
        data_pipeline_config,
        data_pipeline_config.small_bfd_database_path,
    )
    self._uniprot_msa_config = msa.get_protein_msa_config(
        data_pipeline_config,
        data_pipeline_config.uniprot_cluster_annot_database_path,
    )
    self._nt_rna_msa_config = msa_config.RunConfig(
        config=msa_config.NhmmerConfig(
            binary_path=data_pipeline_config.nhmmer_binary_path,
            hmmalign_binary_path=data_pipeline_config.hmmalign_binary_path,
            hmmbuild_binary_path=data_pipeline_config.hmmbuild_binary_path,
            database_config=msa_config.DatabaseConfig(
                name='nt_rna',
                path=data_pipeline_config.ntrna_database_path,
            ),
            n_cpu=data_pipeline_config.nhmmer_n_cpu,
            e_value=1e-3,
            alphabet='rna',
            max_sequences=10000,
        ),
        chain_poly_type=mmcif_names.RNA_CHAIN,
        crop_size=None,
    )
    self._rfam_msa_config = msa_config.RunConfig(
        config=msa_config.NhmmerConfig(
            binary_path=data_pipeline_config.nhmmer_binary_path,
            hmmalign_binary_path=data_pipeline_config.hmmalign_binary_path,
            hmmbuild_binary_path=data_pipeline_config.hmmbuild_binary_path,
            database_config=msa_config.DatabaseConfig(
                name='rfam_rna',
                path=data_pipeline_config.rfam_database_path,
            ),
            n_cpu=data_pipeline_config.nhmmer_n_cpu,
            e_value=1e-3,
            alphabet='rna',
            max_sequences=10000,
        ),
        chain_poly_type=mmcif_names.RNA_CHAIN,
        crop_size=None,
    )
    self._rnacentral_msa_config = msa_config.RunConfig(
        config=msa_config.NhmmerConfig(
            binary_path=data_pipeline_config.nhmmer_binary_path,
            hmmalign_binary_path=data_pipeline_config.hmmalign_binary_path,
            hmmbuild_binary_path=data_pipeline_config.hmmbuild_binary_path,
            database_config=msa_config.DatabaseConfig(
                name='rna_central_rna',
                path=data_pipeline_config.rna_central_database_path,
            ),
            n_cpu=data_pipeline_config.nhmmer_n_cpu,
            e_value=1e-3,
            alphabet='rna',
            max_sequences=10000,
        ),
        chain_poly_type=mmcif_names.RNA_CHAIN,
        crop_size=None,
    )

    self._templates_config = msa_config.TemplatesConfig(
        template_tool_config=msa_config.TemplateToolConfig(
            database_path=data_pipeline_config.seqres_database_path,
            chain_poly_type=mmcif_names.PROTEIN_CHAIN,
            hmmsearch_config=msa_config.HmmsearchConfig(
                hmmsearch_binary_path=data_pipeline_config.hmmsearch_binary_path,
                hmmbuild_binary_path=data_pipeline_config.hmmbuild_binary_path,
                filter_f1=0.1,
                filter_f2=0.1,
                filter_f3=0.1,
                e_value=100,
                inc_e=100,
                dom_e=100,
                incdom_e=100,
                alphabet='amino',
            ),
        ),
        filter_config=msa_config.TemplateFilterConfig(
            max_subsequence_ratio=0.95,
            min_align_ratio=0.1,
            min_hit_length=10,
            deduplicate_sequences=True,
            max_hits=4,
            # By default, use the date from AF3 paper.
            max_template_date=datetime.date(2021, 9, 30),
        ),
    )
    self._pdb_database_path = data_pipeline_config.pdb_database_path

  def process_protein_chain(
      self, chain: folding_input.ProteinChain
  ) -> folding_input.ProteinChain:
    """Processes a single protein chain."""
    if chain.unpaired_msa or chain.paired_msa or chain.templates:
      if (
          chain.unpaired_msa is None
          or chain.paired_msa is None
          or chain.templates is None
      ):
        raise ValueError(
            f'Protein chain {chain.id} has unpaired MSA, paired MSA, or'
            ' templates set only partially. If you want to run the pipeline'
            ' with custom MSA/templates, you need to set all of them. You can'
            ' set MSA to empty string and templates to empty list to signify'
            ' that they should not be used and searched for.'
        )
      logging.info(
          'Skipping MSA and template search for protein chain %s because it '
          'already has MSAs and templates.',
          chain.id,
      )
      return chain

    unpaired_msa, paired_msa, template_hits = _get_protein_msa_and_templates(
        sequence=chain.sequence,
        uniref90_msa_config=self._uniref90_msa_config,
        mgnify_msa_config=self._mgnify_msa_config,
        small_bfd_msa_config=self._small_bfd_msa_config,
        uniprot_msa_config=self._uniprot_msa_config,
        templates_config=self._templates_config,
        pdb_database_path=self._pdb_database_path,
    )

    unpaired_a3m = unpaired_msa.to_a3m()
    paired_a3m = paired_msa.to_a3m()

    # Create a new ProteinChain instance with updated MSA and templates
    new_chain = folding_input.ProteinChain(
        id=chain.id,
        sequence=chain.sequence,
        ptms=chain.ptms,
        unpaired_msa=unpaired_a3m,
        paired_msa=paired_a3m,
        templates=[
            folding_input.Template(
                mmcif=struc.to_mmcif(),
                query_to_template_map=hit.query_to_hit_mapping,
            )
            for hit, struc in template_hits
        ],
    )


    return new_chain

  def process_rna_chain(
      self, chain: folding_input.RnaChain
  ) -> folding_input.RnaChain:
    """Processes a single RNA chain."""
    if chain.unpaired_msa:
      # Don't run MSA tools if the chain already has an MSA.
      logging.info(
          'Skipping MSA search for RNA chain %s because it already has MSA.',
          chain.id,
      )
      return chain

    rna_msa = _get_rna_msa(
        sequence=chain.sequence,
        nt_rna_msa_config=self._nt_rna_msa_config,
        rfam_msa_config=self._rfam_msa_config,
        rnacentral_msa_config=self._rnacentral_msa_config,
    )
    return dataclasses.replace(chain, unpaired_msa=rna_msa.to_a3m())

  def process(self, fold_input: folding_input.Input) -> folding_input.Input:
    """Runs MSA and template tools and returns a new Input with the results."""
    processed_chains = []
    for chain in fold_input.chains:
      print(f'Processing chain {chain.id}')
      process_chain_start_time = time.time()
      match chain:
        case folding_input.ProteinChain():
          processed_chains.append(self.process_protein_chain(chain))
        case folding_input.RnaChain():
          processed_chains.append(self.process_rna_chain(chain))
        case _:
          processed_chains.append(chain)
      print(
          f'Processing chain {chain.id} took'
          f' {time.time() - process_chain_start_time:.2f} seconds',
      )

    return dataclasses.replace(fold_input, chains=processed_chains)
