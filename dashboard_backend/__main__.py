"""Entrypoint for ``python -m dashboard_backend``.

Why this exists
---------------
On Windows + Python 3.14, uvicorn's default event-loop policy ends up being
``WindowsSelectorEventLoop``, which raises ``NotImplementedError`` on
``asyncio.create_subprocess_exec``.  Our ``/api/vast/*`` and ``/api/runs/*``
routes shell out to the ``vastai`` CLI via subprocess_exec, so the entire
Vast.ai integration is dead under that policy.

uvicorn creates its event loop BEFORE importing the FastAPI app module, so a
policy assignment inside ``dashboard_backend.main`` runs too late.  The only
reliable fix is a tiny shim that:

  1. Sets ``WindowsProactorEventLoopPolicy`` (supports subprocess) BEFORE
     uvicorn is imported.
  2. Calls ``uvicorn.run`` programmatically with the same defaults the
     dashboard launcher script uses.

Usage
-----
Replace ``python -m uvicorn dashboard_backend.main:app ...`` with::

    python -m dashboard_backend [--host HOST] [--port PORT] [--no-reload]

Same CLI shape as ``uvicorn`` for the flags we accept; everything else uses
sane localhost-dev defaults.
"""
from __future__ import annotations

import argparse
import asyncio
import sys


def _install_windows_proactor_policy() -> None:
    if sys.platform != "win32":
        return
    try:
        asyncio.set_event_loop_policy(asyncio.WindowsProactorEventLoopPolicy())
    except Exception as exc:  # noqa: BLE001 -- not fatal on weird embeddings
        print(f"warning: could not install WindowsProactorEventLoopPolicy: {exc}", file=sys.stderr)


def main() -> int:
    _install_windows_proactor_policy()

    parser = argparse.ArgumentParser(
        prog="python -m dashboard_backend",
        description="CASCADE dashboard backend launcher (subprocess-capable on Windows).",
    )
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument(
        "--no-reload",
        dest="reload",
        action="store_false",
        default=True,
        help="Disable hot reload (default: on for dev).",
    )
    parser.add_argument("--log-level", default="info", choices=["critical", "error", "warning", "info", "debug", "trace"])
    parser.add_argument("--workers", type=int, default=1, help="Uvicorn worker count (reload requires 1).")
    args = parser.parse_args()

    import uvicorn  # imported AFTER the policy is set
    uvicorn.run(
        "dashboard_backend.main:app",
        host=args.host,
        port=args.port,
        reload=args.reload,
        log_level=args.log_level,
        workers=args.workers if not args.reload else 1,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
