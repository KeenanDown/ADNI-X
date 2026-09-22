from dataclasses import replace
from pathlib import Path
import json

import numpy as np
import pytest
from scipy.io import savemat

from adni_x.config import Config, Cohort
from adni_x.inputs import acquisition_metadata, entities, fingerprints, inventory
from adni_x.preparation import prepare


def latest(config):
    return sorted(Path(config.inputs.tensors).glob('*ses-M012*.mat'))


def set_tr(tensor, value):
    sidecar = Path(acquisition_metadata(tensor)['bold_json'])
    sidecar.write_text(json.dumps({'RepetitionTime': value}))
    return sidecar


def test_tr_boundaries_and_outlier_before_visit_selection(config):
    tensors = latest(config)[:5]
    values = [2.875, 3.150, np.nextafter(2.875, 0), np.nextafter(3.150, 4), 6.015]
    for tensor, tr in zip(tensors, values):
        set_tr(tensor, tr)
    # The manuscript's outlier passes the length rule but fails acquisition TR.
    series = next(Path(config.inputs.series).glob(f'*{entities(tensors[-1])["sub"]}*M012*'))
    savemat(series, {'activity': np.ones((4, 197))})
    inv = inventory(config)
    for tensor, month, reason, tr in zip(tensors, [12, 12, 0, 0, 0],
            ['eligible', 'eligible', 'below_tr_window', 'above_tr_window', 'above_tr_window'], values):
        rows = inv.scans.loc[inv.scans.subject == entities(tensor)['sub']]
        assert rows.loc[rows.selected, 'month'].item() == month
        scan = rows.loc[rows.month == 12].iloc[0]
        assert scan.reason == reason and scan.TR_seconds == tr
        assert scan.tr_status == 'ok'
    assert inv.summary()['repetition_time_rule']['inclusive'] is True
    for visit in ('baseline', 'M012'):
        audit = inventory(replace(config, cohort=replace(config.cohort, visit=visit))).scans
        assert audit.loc[audit.selected, 'TR_seconds'].between(2.875, 3.150).all()
        assert not audit.loc[audit.reason == 'above_tr_window', 'eligible'].any()


@pytest.mark.parametrize('value', [None, True, '3.0', 0, -3, float('nan'), float('inf'), [], {}])
def test_invalid_tr_excluded(config, value):
    tensor = latest(config)[0]
    set_tr(tensor, value)
    audit = inventory(config).scans
    row = audit.loc[audit.tensor_file == str(tensor)].iloc[0]
    assert row.reason == 'invalid_repetition_time'
    assert not row.eligible


@pytest.mark.parametrize('case, expected', [
    ('missing_tensor', 'missing_tensor_json'), ('broken_tensor', 'unreadable_tensor_json'),
    ('list_tensor', 'invalid_tensor_json'), ('no_source', 'missing_or_invalid_source_bold'),
    ('missing_bold', 'missing_bold_json'), ('broken_bold', 'unreadable_bold_json'),
    ('no_tr', 'missing_repetition_time'), ('wrong_subject', 'conflicting_source_bold_identity'),
    ('wrong_run', 'conflicting_source_bold_identity'), ('extension', 'invalid_source_bold_extension'),
])
def test_metadata_failures_have_separate_audit_reasons(config, case, expected):
    tensor = latest(config)[0]
    tensor_json = tensor.with_suffix('.json')
    bold_json = Path(acquisition_metadata(tensor)['bold_json'])
    if case == 'missing_tensor':
        tensor_json.unlink()
    elif case == 'broken_tensor':
        tensor_json.write_text('{')
    elif case == 'list_tensor':
        tensor_json.write_text('[]')
    elif case == 'no_source':
        tensor_json.write_text('{}')
    elif case == 'missing_bold':
        bold_json.unlink()
    elif case == 'broken_bold':
        bold_json.write_text('{')
    elif case == 'no_tr':
        bold_json.write_text('{}')
    else:
        source = json.loads(tensor_json.read_text())['source_bold']
        if case == 'wrong_subject':
            source = source.replace(entities(tensor)['sub'], 'ADNI999S9999')
        elif case == 'wrong_run':
            renamed = tensor.with_name(tensor.name.replace('task-rest', 'task-rest_run-1'))
            tensor.rename(renamed)
            tensor_json.rename(renamed.with_suffix('.json'))
            tensor, tensor_json = renamed, renamed.with_suffix('.json')
            source = source.replace('task-rest', 'task-rest_run-2')
        else:
            source = source.replace('.nii.gz', '.mat')
        tensor_json.write_text(json.dumps({'source_bold': source}))
    inv = inventory(config)
    row = inv.scans.loc[inv.scans.tensor_file == str(tensor)].iloc[0]
    assert row.reason == expected and row.tr_status == expected
    assert not row.eligible
    # Unresolved metadata must not abort provenance or unrelated eligible scans.
    assert inv.scans.selected.sum() == 24
    fingerprints(config)


@pytest.mark.parametrize('absolute, extension', [(False, '.nii'), (False, '.nii.gz'), (True, '.nii.gz')])
def test_explicit_source_paths_without_loading_images(config, absolute, extension):
    tensor = latest(config)[0]
    record = acquisition_metadata(tensor)
    source = Path(record['source_bold']).with_name(Path(record['bold_json']).stem + extension)
    value = str(source) if absolute else source.name
    tensor.with_suffix('.json').write_text(json.dumps({'source_bold': value}))
    assert not source.exists()  # Metadata-only review; no raw image read needed.
    assert acquisition_metadata(tensor)['TR_seconds'] == 3.0


def test_tr_metadata_and_settings_invalidate_preparation(config):
    tensor = latest(config)[0]
    first = prepare(config)
    bold_json = set_tr(tensor, 6.015)
    second = prepare(config)
    assert second.path != first.path
    assert str(bold_json) in second.manifest['sources']
    assert str(tensor.with_suffix('.json')) in second.manifest['sources']
    for name in ('visits', 'selected', 'subjects'):
        assert str(tensor) not in second.table(name).tensor_file.values
        assert 'TR_seconds' in second.table(name)
    wider = prepare(replace(config, cohort=replace(config.cohort, tr_max_seconds=6.1)))
    assert wider.path != second.path
    assert str(tensor) in wider.table('selected').tensor_file.values
    bold_json.unlink()
    absent = prepare(config)
    assert absent.path != second.path
    set_tr(tensor, 3.0)
    assert prepare(config).path == first.path
    # Changing the explicit link to an equivalent sidecar is also a dependency.
    other = bold_json.with_name(bold_json.stem + '_desc-copy.json')
    other.write_text(bold_json.read_text())
    tensor.with_suffix('.json').write_text(json.dumps({'source_bold': other.with_suffix('.nii').name}))
    assert prepare(config).path != first.path


def test_tensor_symlink_uses_original_sidecar_and_relative_source(config, tmp_path):
    tensor = latest(config)[0]
    pilot = tmp_path / 'pilot'
    pilot.mkdir()
    link = pilot / tensor.name
    link.symlink_to(tensor)
    assert acquisition_metadata(link) == acquisition_metadata(tensor)


def test_unreadable_metadata_is_audited_and_fingerprinted(config, monkeypatch):
    import builtins
    tensor = latest(config)[0]
    bold_json = Path(acquisition_metadata(tensor)['bold_json'])
    original = Path.open
    def restricted(path, *args, **kwargs):
        if path == bold_json:
            raise PermissionError(13, 'synthetic unreadable metadata')
        return original(path, *args, **kwargs)
    monkeypatch.setattr(Path, 'open', restricted)
    original_builtin = builtins.open
    def restricted_builtin(path, *args, **kwargs):
        if Path(path) == bold_json:
            raise PermissionError(13, 'synthetic unreadable metadata')
        return original_builtin(path, *args, **kwargs)
    monkeypatch.setattr(builtins, 'open', restricted_builtin)
    audit = inventory(config).scans
    assert audit.loc[audit.tensor_file == str(tensor), 'reason'].item() == 'unreadable_bold_json'
    assert fingerprints(config)[str(bold_json)]['unavailable'] == 'PermissionError'


@pytest.mark.parametrize('lo, hi', [(0, 3), (4, 3), (float('nan'), 3), (2, float('inf')), (True, 3), ('2', 3), (None, 3)])
def test_invalid_tr_configuration(lo, hi):
    with pytest.raises(ValueError, match='tr_min_seconds'):
        Config(cohort=Cohort(tr_min_seconds=lo, tr_max_seconds=hi))
