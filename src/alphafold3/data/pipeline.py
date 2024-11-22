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
import time

from alphafold3.common import folding_input
from alphafold3.constants import mmcif_names
from alphafold3.data import msa
from alphafold3.data import msa_config
from alphafold3.data import structure_stores
from alphafold3.data import templates
from alphafold3.data.tools import mmseqs2


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
  # We first run the MSA tools in parallel.
  with futures.ThreadPoolExecutor() as executor:
    uniref90_msa_future = executor.submit(
        msa.get_msa, sequence, uniref90_msa_config
    )
    mgnify_msa_future = executor.submit(
        msa.get_msa, sequence, mgnify_msa_config
    )
    small_bfd_msa_future = executor.submit(
        msa.get_msa, sequence, small_bfd_msa_config
    )
    uniprot_msa_future = executor.submit(
        msa.get_msa, sequence, uniprot_msa_config
    )

  # Get MSA results first since we need uniref90_msa for template search
  uniref90_msa = uniref90_msa_future.result()
  mgnify_msa = mgnify_msa_future.result()
  small_bfd_msa = small_bfd_msa_future.result()
  uniprot_msa = uniprot_msa_future.result()

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

  # Then we merge the MSAs.
  unpaired_msa = msa.merge_msas(
      msas=[uniref90_msa, mgnify_msa, small_bfd_msa, uniprot_msa],
      deduplicate=True,
  )

  # If the merged MSA is empty, create a minimal MSA with just the query sequence
  if not unpaired_msa.sequences:
    logging.warning('All MSAs are empty, using query sequence only')
    unpaired_msa = msa.Msa(
        sequences=[sequence],
        deletion_matrix=[[0] * len(sequence)],
        descriptions=['query'],
    )
  # If the query sequence is not the first sequence, add it
  elif unpaired_msa.sequences[0] != sequence:
    logging.warning('Query sequence not found in MSA, adding it as first sequence')
    unpaired_msa = msa.Msa(
        sequences=[sequence] + unpaired_msa.sequences,
        deletion_matrix=[[0] * len(sequence)] + unpaired_msa.deletion_matrix,
        descriptions=['query'] + unpaired_msa.descriptions,
    )

  return unpaired_msa, None, template_hits


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

    unpaired_msa, _, template_hits = _get_protein_msa_and_templates(
        sequence=chain.sequence,
        uniref90_msa_config=self._uniref90_msa_config,
        mgnify_msa_config=self._mgnify_msa_config,
        small_bfd_msa_config=self._small_bfd_msa_config,
        uniprot_msa_config=self._uniprot_msa_config,
        templates_config=self._templates_config,
        pdb_database_path=self._pdb_database_path,
    )

    return dataclasses.replace(
        chain,
        unpaired_msa=unpaired_msa.to_a3m(),
        paired_msa="",
        templates=[
            folding_input.Template(
                mmcif=struc.to_mmcif(),
                query_to_template_map=hit.query_to_hit_mapping,
            )
            for hit, struc in template_hits
        ],
    )

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
