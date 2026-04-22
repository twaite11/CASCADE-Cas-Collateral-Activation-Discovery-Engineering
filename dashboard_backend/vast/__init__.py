"""Vast.ai provisioner + asyncssh runner + runs catalog.

The controller uses these modules to turn a `POST /api/runs` request into a
freshly provisioned A100 VPS, stream stdout back to the dashboard, and rsync
artifacts home on completion. No external services required (all state lives
in a single SQLite `runs` table).
"""
from .runs_store import Run, RunsStore, RunStatus  # noqa: F401
from .provisioner import (  # noqa: F401
    VastProvisioner,
    VastOffer,
    VastInstance,
    VastCliError,
)
