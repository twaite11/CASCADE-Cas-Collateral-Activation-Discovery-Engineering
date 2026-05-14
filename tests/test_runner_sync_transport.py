import shutil
from pathlib import Path

import pytest

from dashboard_backend.vast import runner as runner_mod


def test_sync_exe_prefers_rsync(monkeypatch):
    def fake(which: str) -> str | None:
        if which == "rsync":
            return "/usr/bin/rsync"
        if which == "scp":
            return "/usr/bin/scp"
        return None

    monkeypatch.setattr(shutil, "which", fake)
    assert runner_mod._sync_executable() == ("/usr/bin/rsync", "rsync")


def test_sync_exe_falls_back_to_scp(monkeypatch):
    def fake(which: str) -> str | None:
        if which == "rsync":
            return None
        if which == "scp":
            return r"C:\Windows\System32\OpenSSH\scp.exe"
        return None

    monkeypatch.setattr(shutil, "which", fake)
    assert runner_mod._sync_executable()[1] == "scp"


def test_sync_exe_neither_installed(monkeypatch):
    monkeypatch.setattr(shutil, "which", lambda _n: None)
    with pytest.raises(RuntimeError, match="Neither"):
        runner_mod._sync_executable()


def test_cascade_sync_bin_full_file_path(monkeypatch, tmp_path):
    shim = tmp_path / "mysync.bat"
    shim.write_text("@exit /b\n", encoding="ascii")
    monkeypatch.delenv("CASCADE_SYNC_BIN", raising=False)
    monkeypatch.setenv("CASCADE_SYNC_BIN", str(shim))
    monkeypatch.setattr(shutil, "which", lambda _n: None)
    exe, kind = runner_mod._sync_executable()
    assert exe == str(shim.resolve())
    assert kind == "scp"
