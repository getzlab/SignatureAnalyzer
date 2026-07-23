from importlib import import_module

__version__ = '0.0.9'

__all__ = [
    'ardnmf',
    'postprocess_msigs',
    'consensus_cluster',
    'run_maf',
    'run_spectra',
    'run_matrix',
    'supervised_ardnmf',
    'pl',
    'pw',
]


def __getattr__(name):
    try:
        if name == 'pl':
            return import_module('.plotting', __name__)
        if name == 'pw':
            return import_module('.pathways', __name__)
        if name == 'ardnmf':
            return getattr(import_module('.bnmf', __name__), 'ardnmf')
        if name == 'postprocess_msigs':
            return getattr(import_module('.utils', __name__), 'postprocess_msigs')
        if name == 'consensus_cluster':
            return getattr(import_module('.consensus', __name__), 'consensus_cluster')
        if name == 'run_maf':
            return getattr(import_module('.signatureanalyzer', __name__), 'run_maf')
        if name == 'run_spectra':
            return getattr(import_module('.signatureanalyzer', __name__), 'run_spectra')
        if name == 'run_matrix':
            return getattr(import_module('.signatureanalyzer', __name__), 'run_matrix')
        if name == 'supervised_ardnmf':
            return getattr(import_module('.supervised_bnmf', __name__), 'supervised_ardnmf')
    except (ModuleNotFoundError, ImportError):
        raise AttributeError(f'module {__name__!r} has no attribute {name!r}') from None
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')
