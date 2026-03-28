from __future__ import annotations

import benchmarks.fasta_backend_probe as probe


def test_probe_evaluate_without_pyfastx_reports_contract_gap(monkeypatch) -> None:
    monkeypatch.setattr(probe, "pyfastx", None)

    result = probe.evaluate(
        records=32,
        seq_len=80,
        iterations=1,
        random_access_queries=8,
    )

    assert result.parity["pyfaidx_matches_sg"] is True
    assert result.pyfastx is None
    assert result.contracts["sg_supports_file_handles"] is True
    assert result.contracts["pyfastx_supports_file_handles"] is False
    assert result.recommendation == "reject_default_pyfastx_not_installed"
