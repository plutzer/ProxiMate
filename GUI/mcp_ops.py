"""The operations an agent can reach through ``call_tool``.

Each function is a thin, fully typed wrapper over ``backend`` or ``cytoscape_ctl``
that passes actor ``'mcp'``, takes every setting as an explicit argument, and returns
plain JSON.  Read the module docstring of ``mcp_registry`` for what the modes mean.
"""

import math
import os

import numpy as np

import backend
import cytoscape_ctl as ctl
import cytoscape_net
import dataset_store as store
import provenance
from mcp_registry import register
from setup_datasets import CORUM_FILENAME, ORGANISMS

ACTOR = 'mcp'


def _json(value):
    """A JSON-safe copy: NaN and numpy scalars become None and Python numbers."""
    if isinstance(value, dict):
        return {str(k): _json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return None if math.isnan(value) else float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def _records(frame, limit=None):
    n = int(len(frame))
    shown = frame if limit is None else frame.head(limit)
    return {'n': n, 'returned': int(len(shown)), 'columns': list(map(str, frame.columns)),
            'rows': _json(shown.to_dict('records'))}


# --- read ------------------------------------------------------------------------------------

@register('read', "The datasets in the session with their input type, quant type and scored state.",
          tags=('list', 'session'))
def list_datasets() -> dict:
    """Every row of the "Datasets in this Session" table."""
    return {'datasets': _json(store.table().to_dict('records'))}


@register('read', "One dataset: its row, which result files exist, its baits, annotation settings.",
          tags=('dataset', 'baits', 'info'))
def get_dataset_info(dataset: str) -> dict:
    """The session row, presence of ED.csv / merged.csv / annotated_scores.csv /
    Feature_enrichment.csv / run.json, the bait names (from annotated_scores.csv), the
    organism and Human Cell Map setting it was annotated with, and whether a job is
    running on it."""
    return _json(backend.dataset_info(dataset))


@register('read', "Server state: output directory, datasets, running jobs, drawn Cytoscape network, recent MCP activity.",
          tags=('status', 'health', 'jobs'))
def server_status() -> dict:
    import mcp_registry
    snap = ctl.snapshot()
    return _json({'out_dir': backend.OUT_DIR, 'datasets': store.names(),
                  'scored': store.scored_names(), 'jobs': backend.running_jobs(),
                  'cytoscape': {'dataset': snap['dataset'], 'n_nodes': snap['n_nodes'],
                                'n_edges': snap['n_edges'], 'thresholds': snap['thresholds'],
                                'busy': snap['busy']},
                  'activity': mcp_registry.activity(10)})


# --- sandbox ---------------------------------------------------------------------------------

@register('sandbox', "QC metrics at a threshold set: median network size, BioGRID enrichment, mean degree.",
          tags=('thresholds', 'qc', 'metrics', 'known', 'degree'))
def threshold_metrics(dataset: str, thresholds: dict, bait: str = None) -> dict:
    """The Data Thresholding tab's metrics for ``thresholds`` = {SaintScore, BFDR, WD,
    WDFDR}: interactions before and after filtering, known (In.BioGRID) before and
    after, median passing preys per bait, enrichment of known interactions, and the
    mean number of BioGRID partners each passing prey has among the other passing
    preys.  ``bait`` restricts everything to one bait.  Call it at several threshold
    sets to compare them; the GUI's sliders never move."""
    return _json(backend.threshold_metrics(dataset, thresholds, bait=bait))


@register('sandbox', "Protein feature enrichment (GO, UniProt features) of each bait's passing preys.",
          tags=('feature', 'enrichment', 'GO', 'domains', 'hypergeometric'))
def feature_analysis(dataset: str, thresholds: dict, feature_types: list = None,
                     top_n: int = 200) -> dict:
    """Runs the Protein Feature Analysis tab's test at ``thresholds`` and returns the
    table (Bait, Feature, Feature_type, k, n, K, M, p_value, enrichment, adj_p) sorted
    by adjusted p-value, at most ``top_n`` rows.  ``feature_types`` is a subset of
    GO_CC, GO_BP, GO_MF, Motifs, Regions, Repeats, Compositions, Domains (default all).
    Nothing is written: Feature_enrichment.csv stays as the GUI last wrote it."""
    result = backend.feature_analysis(dataset, thresholds, feature_types)
    if len(result):
        result = result.sort_values(['adj_p', 'p_value'], kind='stable')
    out = _records(result, top_n)
    out.update(dataset=dataset, thresholds=backend.validate_thresholds(thresholds))
    return out


@register('sandbox', "Compare two baits: volcano data and the Venn gene lists at separate thresholds.",
          tags=('compare', 'volcano', 'venn', 'baits', 'fold change'))
def compare_networks(dataset: str, bait_a: str, bait_b: str, thresholds_a: dict,
                     thresholds_b: dict, top_n: int = 500) -> dict:
    """The Network Comparison tab's data: ``genes`` holds the sorted gene lists passing
    thresholds in only A, only B and both; ``volcano`` holds every prey seen under
    either bait with status (shared / a_only / b_only), mean intensities, log2 fold
    change, adjusted p-value and the -log10 BFDR for one-sided preys, at most
    ``top_n`` rows ordered by adjusted p-value."""
    out = backend.compare_networks(dataset, bait_a, bait_b, thresholds_a, thresholds_b)
    volcano = out.pop('volcano')
    if len(volcano) and 'pval_adj' in volcano.columns:
        volcano = volcano.sort_values('pval_adj', kind='stable', na_position='last')
    out['volcano'] = _records(volcano, top_n)
    return _json(out)


# --- dataset ----------------------------------------------------------------------------------

@register('dataset', "Parse raw quantification files into a new dataset (MaxQuant, DIA-NN, Pioneer, FragPipe, MSstats or SAINT).",
          tags=('parse', 'import', 'upload', 'maxquant', 'diann', 'saint', 'experimental design'))
def parse_dataset(dataset: str, input_format: str, files: dict, quant_type: str = 'Intensity') -> dict:
    """Creates ``<out_dir>/<dataset>/`` with the SAINT and CompPASS inputs and adds the
    dataset to the session.  ``files`` maps file keys to paths readable by the server:
    MaxQuant {pg, ed}; DIA-NN and Pioneer {matrix, ed}; FragPipe {fp, ed}; MSstats
    {msstats, ed}; SAINT {bait, prey, interaction}.  ``ed`` is the experimental design
    CSV (Experiment Name, Type, Bait, Replicate, Bait ID).  ``quant_type`` is
    Intensity, LFQ or Spectral Counts; DIA-NN, Pioneer and MSstats always use
    Intensity.  The name must be new and use letters, digits and underscores only."""
    return _json(backend.run_parse(dataset, input_format, files, quant_type, actor=ACTOR))


@register('dataset', "Score a parsed dataset with SAINTexpress and CompPASS, then annotate it.",
          tags=('score', 'saint', 'comppass', 'annotate', 'imputation', 'wdfdr'))
def score_dataset(dataset: str, imputation: int, wdfdr_iterations: int, organism: str,
                  exclude_hcm: bool, pi_method: str = None, pi_bait: str = None,
                  seed: int = None) -> dict:
    """Runs score.py and annotator.py on the dataset, as the Scoring card does.
    ``imputation``: 0 none, 1 prey-specific AFT, 2 refactored AFT, 3 one-component
    AFT.  ``wdfdr_iterations``: permutations for the WD FDR (0 skips it).
    ``organism``: human, mouse or yeast.  ``exclude_hcm`` removes Human Cell Map
    evidence from BioGRID (human only).  ``pi_method`` (weighted_average or
    single_bait, with ``pi_bait``) applies to imputation 2.  ``seed`` fixes the
    CompPASS permutations.  Blocks until both stages finish (minutes); refused while
    the GUI or another call is working on the dataset."""
    return _json(backend.run_score(dataset, imputation, wdfdr_iterations, organism, exclude_hcm,
                                   pi_method=pi_method, pi_bait=pi_bait, seed=seed, actor=ACTOR))


@register('dataset', "Replace the whole session with the datasets in a session zip. Destructive; needs confirm=true.",
          tags=('session', 'restore', 'load', 'zip'))
def load_session(zip_path: str, confirm: bool = False) -> dict:
    """Removes every dataset in the output directory and unpacks the archive (one
    written by Download Session) in its place.  Refused unless ``confirm`` is true and
    no job is running.  Tell the user before calling this."""
    if not confirm:
        raise ValueError("load_session replaces every dataset in the session; pass confirm=true "
                         "after the user has agreed")
    table = backend.load_session(zip_path, actor=ACTOR)
    return {'datasets': _json(table.to_dict('records'))}


# --- cytoscape --------------------------------------------------------------------------------

def _network():
    snap = ctl.snapshot()
    return {'dataset': snap['dataset'], 'title': snap['title'], 'n_nodes': snap['n_nodes'],
            'n_edges': snap['n_edges'], 'n_hidden': snap['n_hidden'],
            'thresholds': snap['thresholds'], 'style': snap['style'], 'busy': snap['busy']}


@register('read', "Whether Cytoscape answers, and what ProxiMate network is drawn at which thresholds.",
          tags=('cytoscape', 'connect', 'health', 'status'))
def cytoscape_status(probe: bool = True) -> dict:
    """``probe`` contacts CyREST (a few seconds when Cytoscape is down); the drawn
    network's dataset, counts, thresholds and style come from the controller state."""
    out = {'network': _network(), 'activity': ctl.snapshot()['log'][-10:]}
    if probe:
        out['cytoscape'] = ctl.health()
    return _json(out)


@register('cytoscape', "Draw a dataset's thresholded network in Cytoscape, replacing the previous one.",
          tags=('cytoscape', 'send', 'draw', 'network', 'layout', 'biogrid', 'corum'))
def cytoscape_send(dataset: str, thresholds: dict, baits: list = None, prey_prey: bool = True,
                   label_policy: str = 'all', layout: str = 'force-directed',
                   width_source: str = 'abundance', literature_weighted: bool = False,
                   biogrid_scope: str = 'all', corum: bool = False, corum_min_members: int = 3,
                   corum_min_fraction: float = 0.5) -> dict:
    """As the Cytoscape tab's Send button, with every option explicit.  ``baits``
    limits the drawing (default all).  ``prey_prey`` adds BioGRID prey-prey edges;
    ``biogrid_scope`` all or multivalidated; ``literature_weighted`` thickens them by
    publications.  ``corum`` adds complex edges (human datasets only) under the two
    criteria.  ``label_policy`` all, baits or none.  ``width_source`` abundance,
    SaintScore, WD, FoldChange or uniform.  ``layout`` is a Cytoscape layout name
    (force-directed, grid, circular, ...).  The thresholds and options are recorded in
    the dataset's run.json; the GUI's own sliders keep the user's values."""
    thresholds = backend.validate_thresholds(thresholds)
    if width_source not in cytoscape_net.EDGE_WIDTH_SOURCES:
        raise ValueError(f"width_source must be one of {cytoscape_net.EDGE_WIDTH_SOURCES}")
    if label_policy not in cytoscape_net.LABEL_POLICIES:
        raise ValueError(f"label_policy must be one of {cytoscape_net.LABEL_POLICIES}")
    if biogrid_scope not in cytoscape_net.BIOGRID_SCOPES:
        raise ValueError(f"biogrid_scope must be one of {cytoscape_net.BIOGRID_SCOPES}")
    dataset_path = backend.dataset_dir(dataset)
    scores_path = backend.results_path(dataset)
    if not os.path.isfile(scores_path):
        raise FileNotFoundError(f"dataset {dataset!r} has no annotated_scores.csv; score it first")
    settings = provenance.annotation_settings(dataset_path)
    organism = settings['organism']
    biogrid_path = provenance.biogrid_summary_path(organism, exclude_hcm=settings['exclude_hcm'])
    corum_path = None
    if corum:
        if not ORGANISMS[organism]['has_corum']:
            raise ValueError(f"CORUM covers human complexes only; this dataset is {organism}")
        corum_path = os.path.join(provenance.DEFAULT_DATASETS_DIR, CORUM_FILENAME)
    params = {'thresholds': thresholds, 'baits': baits, 'prey_prey': prey_prey,
              'label_policy': label_policy, 'layout': layout, 'width_source': width_source,
              'literature_weighted': literature_weighted, 'biogrid_scope': biogrid_scope,
              'corum': corum, 'corum_min_members': corum_min_members,
              'corum_min_fraction': corum_min_fraction, 'actor': ACTOR}
    with provenance.stage(dataset_path, 'cytoscape', entrypoint='mcp.cytoscape_send',
                          params=params) as record:
        record.add_input(scores_path, role='scores')
        ctl.draw(dataset, scores_path, thresholds, baits=list(baits or []), prey_prey=prey_prey,
                 biogrid_path=biogrid_path, label_policy=label_policy, layout=layout,
                 width_source=width_source, literature_weighted=literature_weighted,
                 biogrid_scope=biogrid_scope, corum_path=corum_path,
                 corum_min_members=int(corum_min_members),
                 corum_min_fraction=float(corum_min_fraction), actor=ACTOR)
        network = _network()
        record.metric('n_nodes', network['n_nodes'])
        record.metric('n_edges', network['n_edges'])
    return _json({'network': network})


@register('cytoscape', "Hide edges failing tighter thresholds on the drawn network without redrawing.",
          tags=('cytoscape', 'thresholds', 'rethreshold', 'hide'))
def cytoscape_apply_thresholds(thresholds: dict) -> dict:
    """Only thresholds at least as tight as the drawn ones apply in place; looser
    ones need cytoscape_send.  Returns the number of hidden edges."""
    hidden = ctl.apply_thresholds(backend.validate_thresholds(thresholds), actor=ACTOR)
    return {'hidden': int(hidden), 'network': _json(_network())}


@register('cytoscape', "Restyle the drawn network: edge width source, BioGRID scope, publication weighting.",
          tags=('cytoscape', 'style', 'edge width', 'biogrid'))
def cytoscape_restyle(width_source: str = 'abundance', literature_weighted: bool = False,
                      biogrid_scope: str = 'all') -> dict:
    """A style update in place; nothing moves.  Returns how many edges changed."""
    changed = ctl.restyle_edges(width_source, literature_weighted, biogrid_scope, actor=ACTOR)
    return {'changed': int(changed)}


@register('read', "The nodes selected in Cytoscape with the scores of their edges.",
          tags=('cytoscape', 'selection', 'read'))
def cytoscape_read_selection() -> dict:
    """Fails when nothing is selected.  The selection is the user's: ask before acting on it."""
    nodes, edges = ctl.read_selection()
    return _json({'nodes': nodes.to_dict('records'), 'edges': edges.to_dict('records')})


@register('cytoscape', "Select nodes by accession or gene symbol, replacing or adding to the selection.",
          tags=('cytoscape', 'select', 'nodes'))
def cytoscape_select_nodes(ids: list, add: bool = False) -> dict:
    """``ids`` are UniProt accessions or gene symbols as drawn; unknown names fail."""
    return {'selected': ctl.select_nodes(ids, add=add, actor=ACTOR), 'added': bool(add)}


@register('cytoscape', "Select the nodes a relation names for a seed: interactors, singletons, partners, cocomplex.",
          tags=('cytoscape', 'select', 'relation', 'interactors', 'partners'))
def cytoscape_select_related(seed: str, relation: str, add: bool = False, min_saint: float = None,
                             max_bfdr: float = None, min_abundance: float = None,
                             min_publications: int = None) -> dict:
    """``interactors`` and ``singletons`` need a bait seed and take the SAINT, BFDR and
    abundance cuts; ``partners`` are BioGRID or complex neighbors of any node
    (``min_publications`` applies); ``cocomplex`` needs the CORUM layer drawn."""
    ids = ctl.select_related(seed, relation, add=add, actor=ACTOR, min_saint=min_saint,
                             max_bfdr=max_bfdr, min_abundance=min_abundance,
                             min_publications=min_publications)
    return {'selected': list(ids), 'added': bool(add)}


@register('cytoscape', "With one bait selected, select it with the preys whose only neighbor it is.",
          tags=('cytoscape', 'select', 'loners', 'bait'))
def cytoscape_select_loners() -> dict:
    return {'selected': list(ctl.select_loners(actor=ACTOR))}


@register('cytoscape', "With one bait selected, select it with its own preys and the two-bait preys nearer to it.",
          tags=('cytoscape', 'select', 'satellites', 'bait'))
def cytoscape_select_satellites() -> dict:
    return {'selected': list(ctl.select_satellites(actor=ACTOR))}


@register('cytoscape', "Hide or show edges against the selection: hide_selected, show_selected, hide_unselected, show_all.",
          tags=('cytoscape', 'edges', 'hide', 'show', 'visibility'))
def cytoscape_set_edge_visibility(action: str) -> dict:
    """A column update; nothing moves and nothing is rebuilt."""
    return {'changed': int(ctl.set_edge_visibility(action, actor=ACTOR))}


@register('read', "Node positions in Cytoscape, all drawn nodes or the named ones.",
          tags=('cytoscape', 'positions', 'layout', 'read'))
def cytoscape_get_positions(ids: list = None) -> dict:
    return {'positions': _json(ctl.get_positions(ids))}


@register('cytoscape', "Move nodes to absolute positions {id or symbol: [x, y]}.",
          tags=('cytoscape', 'move', 'positions', 'layout'))
def cytoscape_move_nodes(positions: dict) -> dict:
    """Plain position writes, so the nodes stay draggable.  A hand layout in progress
    is affected: say what will move before calling this."""
    return {'moved': int(ctl.move_nodes(positions, actor=ACTOR))}


@register('cytoscape', "Leiden over the selected nodes, colored by community and re-packed inside their box.",
          tags=('cytoscape', 'cluster', 'leiden', 'community', 'repack'))
def cytoscape_cluster_selection(resolution: float = 1.0, seed: int = 17,
                                literature_weight: float = 1.0) -> dict:
    """Only the selected nodes move; the community number lands in the node table."""
    return _json(ctl.cluster_selection(resolution, seed, literature_weight, actor=ACTOR))


@register('cytoscape', "Record the positions the user dragged nodes to into the controller state.",
          tags=('cytoscape', 'sync', 'positions'))
def cytoscape_sync_positions() -> dict:
    return {'recorded': int(ctl.sync_positions(actor=ACTOR))}


@register('cytoscape', "Export the drawn network as a PNG under the dataset's cytoscape folder.",
          tags=('cytoscape', 'export', 'png', 'image'))
def cytoscape_export_image(height: int = 2000) -> dict:
    """Returns the image path; recorded in the dataset's run.json with the drawn thresholds."""
    snap = ctl.snapshot()
    if snap['dataset'] is None:
        raise RuntimeError("no ProxiMate network is in Cytoscape; send one first")
    dataset_path = backend.dataset_dir(snap['dataset'])
    with provenance.stage(dataset_path, 'cytoscape', entrypoint='mcp.cytoscape_export_image',
                          params={'height': height, 'actor': ACTOR}) as record:
        path = ctl.export_image(dataset_path, height=int(height), actor=ACTOR)
        record.add_output(path, role='image')
        record.extra(thresholds=snap['thresholds'])
    return {'path': path}


@register('cytoscape', "Release any lock a Cytoscape view picked up, restoring pan and zoom.",
          tags=('cytoscape', 'unlock', 'view'))
def cytoscape_unlock() -> dict:
    return {'released': list(ctl.unlock(actor=ACTOR))}
