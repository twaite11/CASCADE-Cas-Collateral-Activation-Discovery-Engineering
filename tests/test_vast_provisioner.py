"""Provisioner tests with the vastai CLI mocked via a fake binary on PATH."""
import json
import os
import stat
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from dashboard_backend.vast.provisioner import VastCliError, VastProvisioner
from dashboard_backend.vast.runner import RunCommand, SshEndpoint, rsync_cli_invocation


def _write_fake_vastai_stdout(dir_: Path, stdout_payload: str, *, exit_code: int = 0, to_stderr: bool = False) -> Path:
    """Create a fake ``vastai`` executable that prints ``stdout_payload``.

    To avoid shell-quoting headaches (especially Windows cmd.exe with JSON),
    write the payload to a sidecar file and have the stub ``type``/``cat`` it
    back to stdout (or stderr). Cross-platform and quote-safe.
    """
    payload_file = dir_ / "payload.txt"
    payload_file.write_text(stdout_payload, encoding="utf-8")

    if os.name == "nt":
        binpath = dir_ / "vastai.cmd"
        redirect = "1>&2" if to_stderr else ""
        # type prints file content. Preserves all JSON characters verbatim.
        binpath.write_text(
            f'@echo off\r\ntype "{payload_file}" {redirect}\r\nexit /b {exit_code}\r\n',
            encoding="ascii",
        )
    else:
        binpath = dir_ / "vastai"
        redirect = ">&2" if to_stderr else ""
        binpath.write_text(
            f'#!/usr/bin/env bash\ncat "{payload_file}" {redirect}\nexit {exit_code}\n',
            encoding="utf-8",
        )
        st = os.stat(binpath)
        os.chmod(binpath, st.st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return binpath


@pytest.mark.asyncio
async def test_search_offers_parses_raw(tmp_path):
    payload = json.dumps(
        [
            {
                "id": 12345,
                "gpu_name": "A100_SXM4",
                "num_gpus": 1,
                "gpu_ram": 81920,
                "cpu_cores": 16,
                "cpu_ram": 131072,
                "disk_space": 200.0,
                "dph_total": 1.42,
                "inet_down": 800,
                "inet_up": 400,
                "datacenter": "us-east",
                "reliability2": 0.987,
            }
        ]
    )
    binpath = _write_fake_vastai_stdout(tmp_path, payload)

    prov = VastProvisioner(bin_path=str(binpath))
    offers = await prov.search_offers(limit=1)
    assert len(offers) == 1
    assert offers[0].id == 12345
    assert offers[0].gpu_name == "A100_SXM4"
    assert offers[0].gpu_ram_gb == pytest.approx(80.0, rel=0.01)
    assert offers[0].dph_total == 1.42


@pytest.mark.asyncio
async def test_create_instance_returns_id(tmp_path):
    payload = json.dumps({"success": True, "new_contract": 99887766})
    binpath = _write_fake_vastai_stdout(tmp_path, payload)

    prov = VastProvisioner(bin_path=str(binpath))
    inst_id = await prov.create_instance(
        offer_id=12345,
        image="ghcr.io/x/cascade-orchestrator:test",
        onstart_cmd="echo hi",
        env_vars={"CASCADE_RUN_ID": "abc"},
        label="run-abc",
    )
    assert inst_id == 99887766


@pytest.mark.asyncio
async def test_cli_error_propagates(tmp_path):
    binpath = _write_fake_vastai_stdout(tmp_path, "kaboom\n", exit_code=2, to_stderr=True)
    prov = VastProvisioner(bin_path=str(binpath))
    with pytest.raises(VastCliError) as excinfo:
        await prov.search_offers()
    assert "kaboom" in str(excinfo.value)
    assert excinfo.value.returncode == 2


def test_run_command_onstart_cmd_shape():
    cmd = RunCommand(
        run_id="abc123",
        baseline_ids=["B1", "B2"],
        crrna_lookup_ids=["L1", "L2"],
        max_generations=4,
        variants_per_gen=3,
        workers=2,
    ).onstart_cmd()
    assert "conda activate cascade" in cmd
    assert "--run-id abc123" in cmd
    assert "--baseline-ids B1,B2" in cmd
    assert "--crrna-lookup-ids L1,L2" in cmd
    assert "--max-generations 4" in cmd
    assert "--variants-per-gen 3" in cmd
    assert "--workers 2" in cmd
    assert "tee" in cmd
    assert "run-abc123.log" in cmd


def test_rsync_cli_invocation_includes_endpoint():
    ep = SshEndpoint(host="1.2.3.4", port=22022, user="root")
    cmd = rsync_cli_invocation("./local", "root@1.2.3.4:/remote", ep)
    assert "-azP" in cmd
    # -e arg carries the ssh invocation with our port + strict-host flags
    idx = cmd.index("-e")
    ssh_str = cmd[idx + 1]
    assert "-p 22022" in ssh_str
    assert "StrictHostKeyChecking=no" in ssh_str
